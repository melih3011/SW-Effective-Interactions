/*
 * sw_interface.c -- planar liquid-vapor interface of a binary square-well
 * mixture in classical DFT, and the effective one-body potential the
 * interface exerts on a dilute solute.
 *
 *   M. Gül and R. Roth, "Effective one-body interactions due to the presence
 *   of a liquid-vapor interface", Physica A 674, 130686 (2025).
 *   arXiv:2608.09375
 *
 * Model
 *   Component 1 (solvent) and component 2 (dilute solute, "nano particles")
 *   are hard spheres with square-well attractions. Cross interactions follow
 *   from the Lorentz-Berthelot rules eps_12 = sqrt(eps_1 eps_2),
 *   lambda_12 sigma_12 = (lambda_1 sigma_1 + lambda_2 sigma_2)/2, which are
 *   built into the functional.
 *
 * Functional
 *   F = F_id + F_hs + F_sw, with F_hs the White-Bear FMT functional and F_sw
 *   the RPA of the SW attraction written in FMT form (Eq. 17 of the paper):
 *       beta Phi_sw = n0 n3 + n1 n2 - n1v . n2v ,
 *   with SW weights of radius lambda_i R_i that carry sqrt(beta eps_i).
 *
 * What the program does
 *   1. Solves the binary bulk coexistence P^l = P^v, mu_i^l = mu_i^v for
 *      (rho_1^l, rho_1^v, rho_2^v) at fixed rho_2^l = rho_{2,0}; with -s it
 *      instead finds the "symmetric" eps_2 for which rho_2^l = rho_2^v.
 *   2. Minimises the grand potential for the free planar interface (vapor at
 *      z<0, liquid at z>0) by Picard iteration and computes the surface tension.
 *   3. Computes the effective one-body potential in two ways:
 *        mixture route  beta V(z) = -ln(rho_2(z)/rho_{2,0})
 *        Widom route    beta V(z) = c2^(1)(+inf) - c2^(1)(z), evaluated with
 *                       the solvent profile and rho_2 -> 0.
 *
 * Units
 *   Internally lengths are measured in R1 = sigma_1/2 and energies in kT.
 *   All output is in units of sigma_1: z/sigma_1, rho*sigma_1^3,
 *   beta gamma sigma_1^2.
 */
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "white_bear.h"

#define SIGMA3 8.0 /* (sigma_1/R1)^3: converts rho R1^3 to rho sigma_1^3 */

/* ------------------------------------------------------------------------ */
/* parameters                                                                */
/* ------------------------------------------------------------------------ */

typedef struct {
  double q;       /* size ratio sigma_2/sigma_1                          */
  double eps1;    /* beta eps_1                                           */
  double lam1;    /* lambda_1                                             */
  double eps2;    /* beta eps_2                                           */
  double lam2;    /* lambda_2 (<0: equal well width, 1+(lam1-1)/q)        */
  double eta2;    /* solute packing fraction in the liquid, rho_{2,0}     */
  int grid;       /* grid points per R1 = sigma_1/2                       */
  double alpha;   /* Picard mixing parameter                              */
  long max_iter;  /* maximum number of Picard iterations                  */
  double tol;     /* convergence: sum_i (Delta rho_i)^2 / 2 < tol         */
  int symmetric;  /* solve for eps_2 with rho_2^l = rho_2^v               */
  const char *outdir;
} params_t;

typedef struct {
  double R, lam, seps, Rsw; /* radius, SW range, sqrt(beta eps), lam*R     */
  int nR, nsw;              /* kernel half-widths in grid points           */
  double *w3, *w2, *w2v;    /* hard-sphere weights (length 2 nR + 1)       */
  double *s3, *s2, *s2v;    /* square-well weights (length 2 nsw + 1)      */
  double rho_l, rho_v;      /* coexisting bulk densities                   */
  double mu_l, mu_v;        /* bulk excess chemical potentials             */
  double *rho, *dep;        /* profile and dep = -c1 - beta mu_ex^l        */
  double *d2, *dv2;         /* contracted HS derivatives                   */
  double *ds2, *ds2v;       /* contracted SW derivatives                   */
} species_t;

static species_t sp[2];
static double dz;          /* grid spacing in units of R1                   */
static int N, Rb, Rmax;    /* grid size, interface position, widest kernel  */
static int win;            /* Picard update window: Rb - win ... Rb + win   */
static double p_bulk;      /* coexistence pressure                          */

/* weighted densities and the species-independent derivative dPhi/dn3 */
static double *n3, *n2, *n1, *n0, *n2v, *n1v;
static double *nsw3, *nsw2, *nsw1, *nsw0, *nsw2v, *nsw1v;
static double *d3;

/* ------------------------------------------------------------------------ */
/* helpers                                                                   */
/* ------------------------------------------------------------------------ */

static void *xcalloc(size_t n, size_t size)
{
  void *p = calloc(n, size);
  if (!p) {
    fprintf(stderr, "out of memory\n");
    exit(EXIT_FAILURE);
  }
  return p;
}

/* Integral of f over [i0, i1] with the extended closed formula, Eq. (4.1.14)
 * of Numerical Recipes (needs i1 - i0 > 6). */
static double integrate(const double *f, int i0, int i1)
{
  double sum = 3. / 8. * (f[i0] + f[i1]) + 7. / 6. * (f[i0 + 1] + f[i1 - 1]) + 23. / 24 * (f[i0 + 2] + f[i1 - 2]);
  for (int i = i0 + 3; i <= i1 - 3; i++) sum += f[i];
  return sum * dz;
}

enum { W3, W2, W2V };

/* Planar weight function of a sphere of radius a on the grid, including the
 * quadrature weights of the extended closed formula at both ends:
 *   w3(z) = pi (a^2 - z^2),  w2(z) = 2 pi a,  w2v(z) = 2 pi z,  |z| <= a.   */
static double *make_kernel(int n, int type, double a, double pref)
{
  double *k = xcalloc(2 * n + 1, sizeof(double));

  for (int i = 0; i <= 2 * n; i++) {
    double z = (i - n) * dz;
    switch (type) {
    case W3: k[i] = pref * PI * dz * (a * a - z * z); break;
    case W2: k[i] = pref * PI_2 * dz * a; break;
    case W2V: k[i] = pref * PI_2 * dz * z; break;
    }
  }
  k[0] *= 3. / 8.;      k[2 * n] *= 3. / 8.;
  k[1] *= 7. / 6.;      k[2 * n - 1] *= 7. / 6.;
  k[2] *= 23. / 24.;    k[2 * n - 2] *= 23. / 24.;

  return k;
}

/* (f * k)(i) = sum_j f_j k(z_i - z_j) */
static inline double conv(const double *f, int i, const double *k, int n)
{
  double s = 0.0;
  for (int j = i - n; j <= i + n; ++j) s += f[j] * k[n + i - j];
  return s;
}


/* ------------------------------------------------------------------------ */
/* bulk thermodynamics                                                       */
/* ------------------------------------------------------------------------ */

/* Excess chemical potentials and pressure of the bulk mixture (r1, r2). */
static void bulk(double r1, double r2, double *mu1, double *mu2, double *p)
{
  const species_t *a = &sp[0], *b = &sp[1];
  double R1 = a->R, R2 = b->R, L1 = a->Rsw, L2 = b->Rsw;

  double x0 = r1 + r2;
  double x1 = r1 * R1 + r2 * R2;
  double x2 = PI_4 * (r1 * R1 * R1 + r2 * R2 * R2);
  double x3 = PI_4 * r1 / 3. * pow(R1, 3) + PI_4 * r2 / 3. * pow(R2, 3);

  double s0 = -(a->seps * r1 + b->seps * r2);
  double s1 = -(a->seps * r1 * L1 + b->seps * r2 * L2);
  double s2 = PI_4 * (a->seps * r1 * L1 * L1 + b->seps * r2 * L2 * L2);
  double s3 = a->seps * PI_4 / 3. * r1 * pow(L1, 3) + b->seps * PI_4 / 3. * r2 * pow(L2, 3);

  double musw1 = sw_mu_ex(s0, s1, s2, s3, a->seps, L1);
  double musw2 = sw_mu_ex(s0, s1, s2, s3, b->seps, L2);
  double psw = -(s0 * s3 + s1 * s2) + r1 * musw1 + r2 * musw2;

  *mu1 = wb_mu_ex(x0, x1, x2, x3, R1) + musw1;
  *mu2 = wb_mu_ex(x0, x1, x2, x3, R2) + musw2;
  *p = wb_pressure(x0, x1, x2, x3) + psw;
}

static void set_eps2(double eps2)
{
  sp[1].seps = sqrt(eps2);
}

/* Coexistence conditions P^l = P^v, mu_i^l = mu_i^v.
 * Unknowns x = (ln rho_1^l, ln rho_1^v[, y]) with y = ln rho_2^v, or
 * y = beta eps_2 in symmetric mode (rho_2^v = rho_2^l). With n = 2 the solute
 * is absent (pure solvent). */
static struct {
  int n, symmetric;
  double r2l;
} cx;

static void coex_residual(const double *x, double *F)
{
  double r1l = exp(x[0]), r1v = exp(x[1]), r2l = cx.n == 2 ? 0.0 : cx.r2l, r2v = r2l;
  double m1l, m2l, pl, m1v, m2v, pv;

  if (cx.n == 3) {
    if (cx.symmetric)
      set_eps2(x[2]);
    else
      r2v = exp(x[2]);
  }

  bulk(r1l, r2l, &m1l, &m2l, &pl);
  bulk(r1v, r2v, &m1v, &m2v, &pv);

  F[0] = (pl - pv) / pl;
  F[1] = m1l + log(r1l) - m1v - log(r1v);
  if (cx.n == 3) F[2] = m2l + log(r2l) - m2v - log(r2v);
}

static int solve_linear(int n, double A[3][3], double *b)
{
  for (int c = 0; c < n; c++) {
    int piv = c;
    for (int r = c + 1; r < n; r++)
      if (fabs(A[r][c]) > fabs(A[piv][c])) piv = r;
    if (fabs(A[piv][c]) < 1e-300) return -1;
    for (int k = 0; k < n; k++) { double t = A[c][k]; A[c][k] = A[piv][k]; A[piv][k] = t; }
    double t = b[c]; b[c] = b[piv]; b[piv] = t;
    for (int r = c + 1; r < n; r++) {
      double f = A[r][c] / A[c][c];
      for (int k = c; k < n; k++) A[r][k] -= f * A[c][k];
      b[r] -= f * b[c];
    }
  }
  for (int r = n - 1; r >= 0; r--) {
    for (int k = r + 1; k < n; k++) b[r] -= A[r][k] * b[k];
    b[r] /= A[r][r];
  }
  return 0;
}

static double norm(int n, const double *F)
{
  double s = 0;
  for (int k = 0; k < n; k++) s += fabs(F[k]);
  return s;
}

/* Damped Newton iteration with a finite-difference Jacobian; returns the
 * final residual. Steps are limited to 0.2 in each (logarithmic) variable
 * and halved until the residual is finite. */
static double newton(double *x)
{
  int n = cx.n;
  double F[3], Fp[3], Fm[3], J[3][3], xt[3];

  for (int it = 0; it < 100; it++) {
    coex_residual(x, F);
    for (int k = 0; k < n; k++) {
      double h = 1e-7 * fmax(1.0, fabs(x[k])), xs = x[k];
      x[k] = xs + h; coex_residual(x, Fp);
      x[k] = xs - h; coex_residual(x, Fm);
      x[k] = xs;
      for (int r = 0; r < n; r++) J[r][k] = (Fp[r] - Fm[r]) / (2 * h);
    }
    double dx[3];
    for (int k = 0; k < n; k++) dx[k] = -F[k];
    if (solve_linear(n, J, dx)) break;

    double big = 0, fac = 1.0;
    for (int k = 0; k < n; k++) big = fmax(big, fabs(dx[k]));
    if (big > 0.2) fac = 0.2 / big;
    for (int tries = 0; tries < 30; tries++, fac *= 0.5) {
      for (int k = 0; k < n; k++) xt[k] = x[k] + fac * dx[k];
      coex_residual(xt, Fp);
      if (isfinite(norm(n, Fp))) break;
    }
    for (int k = 0; k < n; k++) x[k] = xt[k];
    if (fac * big < 1e-14) break;
  }
  coex_residual(x, F);
  return norm(n, F);
}

/* Pressure and total chemical potential of the pure solvent at packing
 * fraction eta. */
static void pure_state(double eta, double *p, double *mu)
{
  double r = 3. / PI_4 * eta / pow(sp[0].R, 3), m2;
  bulk(r, 0.0, mu, &m2, p);
  *mu += log(r);
}

static double pure_root(double p, double a, double b)
{
  /* P(eta) - p changes sign once on [a, b] */
  double pa, pm, mu;
  pure_state(a, &pa, &mu);
  for (int k = 0; k < 200; k++) {
    double m = 0.5 * (a + b);
    pure_state(m, &pm, &mu);
    if ((pm - p > 0) == (pa - p > 0)) { a = m; pa = pm; } else b = m;
  }
  return 0.5 * (a + b);
}

/* Liquid-vapor coexistence of the pure solvent: locate the van der Waals
 * loop of P(eta), then bisect on the pressure until mu^l = mu^v. Writes
 * x = (ln rho^l, ln rho^v); returns nonzero if there is no loop. */
static int pure_coexistence(double *x)
{
  const int M = 20000;
  const double etamax = 0.7;
  double p, mu, pprev = 0, ea = -1, eb = -1;

  for (int k = 1; k <= M; k++) {
    double e = etamax * k / M;
    pure_state(e, &p, &mu);
    if (k > 1 && ea < 0 && p < pprev) ea = etamax * (k - 1) / M; /* local maximum */
    if (ea > 0 && eb < 0 && p > pprev) eb = etamax * (k - 1) / M; /* local minimum */
    pprev = p;
  }
  if (ea < 0 || eb < 0) return -1;

  double pmax, pmin, ev = 0, el = 0;
  pure_state(ea, &pmax, &mu);
  pure_state(eb, &pmin, &mu);
  double plo = fmax(pmin, 1e-12 * pmax), phi = pmax;

  for (int k = 0; k < 200; k++) {
    double pm = 0.5 * (plo + phi), ml, mv, pp;
    ev = pure_root(pm, 1e-12, ea);
    el = pure_root(pm, eb, etamax);
    pure_state(ev, &pp, &mv);
    pure_state(el, &pp, &ml);
    if (ml > mv) plo = pm; else phi = pm; /* d(mu^l - mu^v)/dp = 1/rho^l - 1/rho^v < 0 */
  }
  x[0] = log(3. / PI_4 * el / pow(sp[0].R, 3));
  x[1] = log(3. / PI_4 * ev / pow(sp[0].R, 3));
  return 0;
}

/* Binary liquid-vapor coexistence at fixed rho_2^l = rho_{2,0}: first the
 * pure solvent, then the dilute solute is added (Eq. 24 of the paper). */
static void coexistence(params_t *P)
{
  double x[3] = {log(3. / PI_4 * 0.35), log(3. / PI_4 * 0.01), 0.0};
  double m1l, m2l, pl, m1v, m2v, pv;

  cx.r2l = 3. / PI_4 * P->eta2 / pow(sp[1].R, 3);
  cx.symmetric = P->symmetric;

  /* 1. pure solvent */
  cx.n = 2;
  if (pure_coexistence(x)) {
    fprintf(stderr, "no liquid-vapor coexistence for beta*eps1 = %g, lambda1 = %g\n", P->eps1, P->lam1);
    exit(EXIT_FAILURE);
  }
  double res = newton(x);
  if (!(res < 1e-10) || x[0] - x[1] < 0.1) {
    fprintf(stderr, "no liquid-vapor coexistence found for the pure solvent (residual %.3e)\n", res);
    exit(EXIT_FAILURE);
  }

  /* 2. initial guess for the solute from the dilute limit rho_2 -> 0 */
  bulk(exp(x[0]), 0.0, &m1l, &m2l, &pl);
  bulk(exp(x[1]), 0.0, &m1v, &m2v, &pv);
  if (P->symmetric) {
    /* For rho_2 -> 0, mu_ex,2 = mu_hs,2 + sqrt(eps2) g_sw is linear in
     * sqrt(eps2): evaluate it at eps2 = 0 and 1 and solve mu^l = mu^v. */
    set_eps2(0.0);
    bulk(exp(x[0]), 0.0, &m1l, &m2l, &pl);
    bulk(exp(x[1]), 0.0, &m1v, &m2v, &pv);
    double hl = m2l, hv = m2v;
    set_eps2(1.0);
    bulk(exp(x[0]), 0.0, &m1l, &m2l, &pl);
    bulk(exp(x[1]), 0.0, &m1v, &m2v, &pv);
    double s = -(hl - hv) / ((m2l - hl) - (m2v - hv));
    x[2] = s * s;
  } else {
    x[2] = log(cx.r2l) + m2l - m2v;
  }

  /* 3. full binary coexistence */
  cx.n = 3;
  res = newton(x);
  if (!(res < 1e-10)) {
    fprintf(stderr, "binary coexistence solver failed (residual %.3e)\n", res);
    exit(EXIT_FAILURE);
  }

  sp[0].rho_l = exp(x[0]);
  sp[0].rho_v = exp(x[1]);
  sp[1].rho_l = cx.r2l;
  if (P->symmetric) {
    P->eps2 = x[2];
    sp[1].rho_v = cx.r2l;
  } else {
    sp[1].rho_v = exp(x[2]);
  }
  set_eps2(P->eps2);

  bulk(sp[0].rho_l, sp[1].rho_l, &sp[0].mu_l, &sp[1].mu_l, &pl);
  bulk(sp[0].rho_v, sp[1].rho_v, &sp[0].mu_v, &sp[1].mu_v, &pv);
  p_bulk = 0.5 * (pl + pv);
}

/* ------------------------------------------------------------------------ */
/* DFT                                                                       */
/* ------------------------------------------------------------------------ */

static void init(params_t *P)
{
  dz = 1.0 / P->grid;

  sp[0].R = 1.0;
  sp[0].nR = P->grid;
  sp[1].R = P->q;
  sp[1].nR = (int)lround(P->q * P->grid);
  if (fabs(sp[1].nR * dz - P->q) > 1e-12)
    fprintf(stderr, "warning: sigma2/sigma1 = %g is not a multiple of 1/%d, using %g\n", P->q, P->grid,
            sp[1].nR * dz);

  sp[0].lam = P->lam1;
  sp[1].lam = P->lam2 > 0 ? P->lam2 : 1. + (P->lam1 - 1.) / P->q;
  sp[0].seps = sqrt(P->eps1);
  set_eps2(P->eps2);

  for (int s = 0; s < 2; s++) {
    sp[s].R = sp[s].nR * dz;
    sp[s].Rsw = sp[s].R * sp[s].lam;
    sp[s].nsw = (int)lround(sp[s].nR * sp[s].lam);
  }

  coexistence(P);

  Rmax = 0;
  for (int s = 0; s < 2; s++) {
    if (sp[s].nR > Rmax) Rmax = sp[s].nR;
    if (sp[s].nsw > Rmax) Rmax = sp[s].nsw;
  }
  /* Picard updates are done in Rb +- win; outside, the densities are kept at
   * their bulk values and act as reservoirs. Two further kernel widths are
   * needed on either side to evaluate the weighted densities. */
  win = 20 * Rmax;
  Rb = win + 2 * Rmax + 1;
  N = 2 * Rb + 1;

  double **arrs[] = {&n3, &n2, &n1, &n0, &n2v, &n1v, &nsw3, &nsw2, &nsw1, &nsw0, &nsw2v, &nsw1v, &d3};
  for (size_t k = 0; k < sizeof arrs / sizeof *arrs; k++) *arrs[k] = xcalloc(N, sizeof(double));

  for (int s = 0; s < 2; s++) {
    species_t *c = &sp[s];
    c->w3 = make_kernel(c->nR, W3, c->R, 1.0);
    c->w2 = make_kernel(c->nR, W2, c->R, 1.0);
    c->w2v = make_kernel(c->nR, W2V, c->R, 1.0);
    c->s3 = make_kernel(c->nsw, W3, c->Rsw, c->seps);
    c->s2 = make_kernel(c->nsw, W2, c->Rsw, c->seps);
    c->s2v = make_kernel(c->nsw, W2V, c->Rsw, c->seps);

    c->rho = xcalloc(N, sizeof(double));
    c->dep = xcalloc(N, sizeof(double));
    c->d2 = xcalloc(N, sizeof(double));
    c->dv2 = xcalloc(N, sizeof(double));
    c->ds2 = xcalloc(N, sizeof(double));
    c->ds2v = xcalloc(N, sizeof(double));

    /* step profile: vapor for z < 0, liquid for z >= 0 */
    for (int i = 0; i < N; i++) c->rho[i] = i < Rb ? c->rho_v : c->rho_l;
  }
}

/* Weighted densities on [Rb - win - Rmax, Rb + win + Rmax]. */
static void weighted_densities(void)
{
  const species_t *a = &sp[0], *b = &sp[1];
  const double La = a->Rsw, Lb = b->Rsw;

#pragma omp parallel for schedule(static)
  for (int i = Rb - win - Rmax; i <= Rb + win + Rmax; ++i) {
    double n3a = conv(a->rho, i, a->w3, a->nR), n3b = conv(b->rho, i, b->w3, b->nR);
    double n2a = conv(a->rho, i, a->w2, a->nR), n2b = conv(b->rho, i, b->w2, b->nR);
    double nva = conv(a->rho, i, a->w2v, a->nR), nvb = conv(b->rho, i, b->w2v, b->nR);
    double n1a = n2a / (PI_4 * a->R), n1b = n2b / (PI_4 * b->R);
    double n0a = n1a / a->R, n0b = n1b / b->R;

    n3[i] = n3a + n3b;
    n2[i] = n2a + n2b;
    n1[i] = n1a + n1b;
    n0[i] = n0a + n0b;
    n2v[i] = nva + nvb;
    n1v[i] = nva / (PI_4 * a->R) + nvb / (PI_4 * b->R);

    /* SW weights: w1 = -w2/(4 pi L), w0 = w1/L, w1v = -w2v/(4 pi L) */
    double s3a = conv(a->rho, i, a->s3, a->nsw), s3b = conv(b->rho, i, b->s3, b->nsw);
    double s2a = conv(a->rho, i, a->s2, a->nsw), s2b = conv(b->rho, i, b->s2, b->nsw);
    double sva = conv(a->rho, i, a->s2v, a->nsw), svb = conv(b->rho, i, b->s2v, b->nsw);
    double s1a = -s2a / (PI_4 * La), s1b = -s2b / (PI_4 * Lb);

    nsw3[i] = s3a + s3b;
    nsw2[i] = s2a + s2b;
    nsw1[i] = s1a + s1b;
    nsw0[i] = s1a / La + s1b / Lb;
    nsw2v[i] = sva + svb;
    nsw1v[i] = -sva / (PI_4 * La) + -svb / (PI_4 * Lb);
  }
}

/* dep_s(z) = -c1_s(z) - beta mu_ex,s^l on the Picard window. */
static void derivative(void)
{
  const int lo = Rb - win - Rmax, hi = Rb + win + Rmax;

#pragma omp parallel for schedule(static)
  for (int i = lo; i <= hi; ++i) {
    d3[i] = wb_d3(n0[i], n1[i], n2[i], n3[i], n1v[i], n2v[i]);
    for (int s = 0; s < 2; s++) {
      species_t *c = &sp[s];
      double L = c->Rsw;
      c->d2[i] = wb_d2(n0[i], n1[i], n2[i], n3[i], n1v[i], n2v[i], c->R);
      c->dv2[i] = wb_dv2(n0[i], n1[i], n2[i], n3[i], n1v[i], n2v[i], c->R);
      /* dPhi_sw/dn_alpha contracted with the SW weights of species s */
      c->ds2[i] = nsw1[i] - nsw2[i] / (PI_4 * L) - nsw3[i] / (PI_4 * L * L);
      c->ds2v[i] = nsw1v[i] - nsw2v[i] / (PI_4 * L);
    }
  }

#pragma omp parallel for schedule(static)
  for (int i = Rb - win; i < Rb + win; ++i) {
    for (int s = 0; s < 2; s++) {
      species_t *c = &sp[s];
      c->dep[i] = conv(d3, i, c->w3, c->nR) + conv(c->d2, i, c->w2, c->nR) - conv(c->dv2, i, c->w2v, c->nR)
                  + conv(nsw0, i, c->s3, c->nsw) + conv(c->ds2, i, c->s2, c->nsw) + conv(c->ds2v, i, c->s2v, c->nsw)
                  - c->mu_l;
    }
  }
}

/* Picard iteration rho <- (1-alpha) rho + alpha rho^l exp(-dep). */
static long minimize(const params_t *P)
{
  long it;
  double res = 1.0;

  for (it = 0; it < P->max_iter && res > P->tol; it++) {
    weighted_densities();
    derivative();

    res = 0.0;
#pragma omp parallel for schedule(static) reduction(+ : res)
    for (int i = Rb - win; i < Rb + win; i++) {
      double h1 = sp[0].rho[i], h2 = sp[1].rho[i];
      sp[0].rho[i] = (1. - P->alpha) * h1 + P->alpha * sp[0].rho_l * exp(-sp[0].dep[i]);
      sp[1].rho[i] = (1. - P->alpha) * h2 + P->alpha * sp[1].rho_l * exp(-sp[1].dep[i]);
      res += 0.5 * (h1 - sp[0].rho[i]) * (h1 - sp[0].rho[i]) + 0.5 * (h2 - sp[1].rho[i]) * (h2 - sp[1].rho[i]);
    }

    if (it % 500 == 0 && isatty(fileno(stderr))) {
      fprintf(stderr, "\r  iteration %7ld   residual %.3e", it, res);
      fflush(stderr);
    }
  }
  if (isatty(fileno(stderr))) fprintf(stderr, "\r  iteration %7ld   residual %.3e\n", it, res);
  if (res > P->tol) fprintf(stderr, "warning: not converged after %ld iterations\n", it);

  return it;
}

/* Surface tension beta gamma = (Omega + P V)/A in units of R1^-2. */
static double surface_tension(void)
{
  const int lo = Rb - 18 * Rmax, hi = Rb + 18 * Rmax;
  double *w = xcalloc(N, sizeof(double));

  weighted_densities();
  for (int i = lo; i <= hi; i++) {
    double f = wb_phi(n0[i], n1[i], n2[i], n3[i], n1v[i], n2v[i]) + p_bulk;
    f += nsw0[i] * nsw3[i] + nsw2[i] * nsw1[i] - nsw2v[i] * nsw1v[i];
    for (int s = 0; s < 2; s++) {
      const species_t *c = &sp[s];
      double rb = i < Rb ? c->rho_v : c->rho_l, mu = i < Rb ? c->mu_v : c->mu_l;
      f += c->rho[i] * (log(c->rho[i] / rb) - 1. - mu);
    }
    w[i] = f;
  }
  double g = integrate(w, lo, hi);
  free(w);
  return g;
}

/* ------------------------------------------------------------------------ */
/* output                                                                    */
/* ------------------------------------------------------------------------ */

static void write_header(FILE *fp, const params_t *P, double gamma, long iters)
{
  double s3 = SIGMA3;
  fprintf(fp, "# binary square-well mixture, planar free liquid-vapor interface (White Bear + FMT-like SW RPA)\n");
  fprintf(fp, "# Gul & Roth, Physica A 674, 130686 (2025)\n");
  fprintf(fp, "# sigma2/sigma1 = %.6f  beta*eps1 = %.6f  lambda1 = %.6f  beta*eps2 = %.10f  lambda2 = %.6f\n", P->q,
          P->eps1, sp[0].lam, P->eps2, sp[1].lam);
  fprintf(fp, "# rho1^l s1^3 = %.10f  rho1^v s1^3 = %.10f\n", s3 * sp[0].rho_l, s3 * sp[0].rho_v);
  fprintf(fp, "# rho2^l s1^3 = %.10e  rho2^v s1^3 = %.10e  (rho_2,0 = rho2^l)\n", s3 * sp[1].rho_l,
          s3 * sp[1].rho_v);
  fprintf(fp, "# beta*P s1^3 = %.10e\n", s3 * p_bulk);
  fprintf(fp, "# beta*gamma*s1^2 = %.6f\n", 4. * gamma);
  fprintf(fp, "# grid = %d points per sigma1/2, alpha = %g, iterations = %ld\n", P->grid, P->alpha, iters);
}

static FILE *open_out(const params_t *P, const char *kind)
{
  char name[1024];
  snprintf(name, sizeof name, "%s/%s_q%.2f_eps%.4f.dat", P->outdir, kind, P->q, P->eps2);
  FILE *fp = fopen(name, "w");
  if (!fp) {
    fprintf(stderr, "cannot open %s: %s\n", name, strerror(errno));
    exit(EXIT_FAILURE);
  }
  printf("  wrote %s\n", name);
  return fp;
}

static void usage(const char *prog)
{
  fprintf(stderr,
          "usage: %s [options]\n"
          "  -q ratio   size ratio sigma2/sigma1                (default 2.0)\n"
          "  -e eps2    beta*eps2 of the solute                 (default 1.3918)\n"
          "  -s         find the symmetric eps2 (rho2^l = rho2^v), -e is the initial guess\n"
          "  -l lam2    SW range of the solute     (default 1 + (lambda1 - 1)/q)\n"
          "  -E eps1    beta*eps1 of the solvent                (default 1.0)\n"
          "  -L lam1    SW range of the solvent                 (default 1.5)\n"
          "  -x eta2    solute packing fraction in the liquid   (default 1e-6)\n"
          "  -g grid    grid points per sigma1/2                (default 60)\n"
          "  -a alpha   Picard mixing parameter                 (default 0.01)\n"
          "  -n iter    maximum number of iterations            (default 75000)\n"
          "  -t tol     convergence tolerance                   (default 1e-14)\n"
          "  -o dir     output directory                        (default .)\n",
          prog);
}

int main(int argc, char **argv)
{
  params_t P = {.q = 2.0, .eps1 = 1.0, .lam1 = 1.5, .eps2 = 1.3918, .lam2 = -1, .eta2 = 1e-6, .grid = 60,
                .alpha = 0.01, .max_iter = 75000, .tol = 1e-14, .symmetric = 0, .outdir = "."};
  int opt;

  while ((opt = getopt(argc, argv, "q:e:sl:E:L:x:g:a:n:t:o:h")) != -1) {
    switch (opt) {
    case 'q': P.q = atof(optarg); break;
    case 'e': P.eps2 = atof(optarg); break;
    case 's': P.symmetric = 1; break;
    case 'l': P.lam2 = atof(optarg); break;
    case 'E': P.eps1 = atof(optarg); break;
    case 'L': P.lam1 = atof(optarg); break;
    case 'x': P.eta2 = atof(optarg); break;
    case 'g': P.grid = atoi(optarg); break;
    case 'a': P.alpha = atof(optarg); break;
    case 'n': P.max_iter = atol(optarg); break;
    case 't': P.tol = atof(optarg); break;
    case 'o': P.outdir = optarg; break;
    default: usage(argv[0]); return opt == 'h' ? EXIT_SUCCESS : EXIT_FAILURE;
    }
  }
  if (P.q <= 0 || P.grid < 4 || P.alpha <= 0 || P.alpha > 1) {
    usage(argv[0]);
    return EXIT_FAILURE;
  }
  mkdir(P.outdir, 0755);

  init(&P);

  double s3 = SIGMA3;
  printf("sigma2/sigma1 = %g, beta*eps2 = %.10f%s, lambda2 = %.6f\n", P.q, P.eps2,
         P.symmetric ? " (symmetric)" : "", sp[1].lam);
  printf("coexistence: rho1 s1^3 = %.6f (l) / %.6f (v), rho2^v/rho2^l = %.6f\n", s3 * sp[0].rho_l,
         s3 * sp[0].rho_v, sp[1].rho_v / sp[1].rho_l);

  long iters = minimize(&P);
  double gamma = surface_tension();
  printf("beta*gamma*sigma1^2 = %.6f\n", 4. * gamma);

  /* density profiles */
  FILE *fp = open_out(&P, "profile");
  write_header(fp, &P, gamma, iters);
  fprintf(fp, "# z/s1  rho1*s1^3  rho2*s1^3  eta(z)=n3\n");
  for (int i = Rb - 10 * Rmax; i < Rb + 10 * Rmax; i++)
    fprintf(fp, "%10.5f  %.10e %.10e %.10e\n", (i - Rb) * dz / 2., s3 * sp[0].rho[i], s3 * sp[1].rho[i], n3[i]);
  fclose(fp);

  /* effective one-body potential: mixture route, then Widom route (rho2 -> 0) */
  double *vmix = xcalloc(N, sizeof(double));
  for (int i = Rb - win; i < Rb + win; i++) vmix[i] = -log(sp[1].rho[i] / sp[1].rho_l);
  for (int i = 0; i < N; i++) sp[1].rho[i] = 0.0;
  weighted_densities();
  derivative();

  fp = open_out(&P, "veff");
  write_header(fp, &P, gamma, iters);
  fprintf(fp, "# z/s1  betaV_eff(mixture route)  betaV_eff(Widom route)\n");
  for (int i = Rb - 15 * Rmax; i <= Rb + 15 * Rmax; i++)
    fprintf(fp, "%10.5f  %.10e %.10e\n", (i - Rb) * dz / 2., vmix[i], sp[1].dep[i]);
  fclose(fp);
  free(vmix);

  return EXIT_SUCCESS;
}
