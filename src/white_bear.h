/*
 * white_bear.h -- White-Bear FMT and FMT-like square-well (RPA) free-energy
 * densities, their partial derivatives and bulk thermodynamics.
 *
 * Part of the code accompanying
 *   M. Gül and R. Roth, "Effective one-body interactions due to the presence
 *   of a liquid-vapor interface", Physica A 674, 130686 (2025).
 */
#ifndef WHITE_BEAR_H
#define WHITE_BEAR_H

#define PI   3.14159265358979323846
#define PI_2 6.2831853071795864769  /* 2 pi */
#define PI_4 12.5663706143591729539 /* 4 pi */

/* ---- hard spheres: White-Bear (mark I) functional ---------------------- */

/* Free-energy density Phi_hs(n_alpha). */
double wb_phi(double n0, double n1, double n2, double n3, double n1v, double n2v);

/* Derivatives of Phi_hs contracted with the scalar weights of a species with
 * radius hR, i.e. the coefficient that multiplies w2 after using
 * w1 = w2/(4 pi R), w0 = w2/(4 pi R^2):
 *   wb_d2  = dPhi/dn2  + dPhi/dn1 /(4 pi R) + dPhi/dn0 /(4 pi R^2)
 *   wb_d3  = dPhi/dn3
 *   wb_dv2 = dPhi/dn2v + dPhi/dn1v/(4 pi R)                              */
double wb_d2(double n0, double n1, double n2, double n3, double n1v, double n2v, double hR);
double wb_d3(double n0, double n1, double n2, double n3, double n1v, double n2v);
double wb_dv2(double n0, double n1, double n2, double n3, double n1v, double n2v, double hR);

/* Bulk (BMCSL) pressure and excess chemical potential of a species with
 * radius hR, in terms of the bulk scaled-particle variables n0..n3.        */
double wb_pressure(double n0, double n1, double n2, double n3);
double wb_mu_ex(double n0, double n1, double n2, double n3, double hR);

/* ---- square well: RPA in FMT form, Phi_sw = n0 n3 + n1 n2 - n1v.n2v ---- */

/* Bulk SW excess chemical potential of a species with sqrt(beta eps) = seps
 * and SW radius hR = lambda R.                                              */
double sw_mu_ex(double n0, double n1, double n2, double n3, double seps, double hR);

#endif
