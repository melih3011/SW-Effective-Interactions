/*
 * white_bear.c -- White-Bear FMT free-energy density, its derivatives and the
 * corresponding bulk (BMCSL) thermodynamics, plus the bulk chemical potential
 * of the FMT-like square-well functional.
 *
 * Roth, Evans, Lang & Kahl, J. Phys.: Condens. Matter 14, 12063 (2002);
 * Yu & Wu, J. Chem. Phys. 117, 10156 (2002).
 */
#include <math.h>

#include "white_bear.h"

/* Below this packing fraction the closed forms suffer from cancellation and
 * Taylor series are used instead. */
#define ETA_SMALL 1.e-3

double wb_pressure(double n0, double n1, double n2, double n3)
{
  double hc = 1. / (1. - n3);

  return hc * (n0 + hc * (n1 * n2 + hc * (n2 * n2 * n2 * (1. - n3 / 3.)) / (12. * PI)));
}

double wb_mu_ex(double n0, double n1, double n2, double n3, double hR)
{
  double hc = 1. - n3;

  return -log(hc) + hR * n2 / hc
         + PI_4 * hR * hR * (n1 / hc + n2 * n2 * (n3 + hc * hc * log(hc)) / (12. * PI * hc * hc * n3 * n3))
         + PI_4 * hR * hR * hR / 3.0
               * (n0 / hc + n1 * n2 / (hc * hc)
                  - n2 * n2 * n2 * (n3 * (2. - 5. * n3 + n3 * n3) + 2 * hc * hc * hc * log(hc))
                        / (36. * PI * hc * hc * hc * n3 * n3 * n3));
}

double sw_mu_ex(double n0, double n1, double n2, double n3, double seps, double hR)
{
  return -n3 * seps + n0 * seps * PI_4 / 3. * hR * hR * hR - n2 * seps * hR + n1 * seps * PI_4 * hR * hR;
}

double wb_phi(double n0, double n1, double n2, double n3, double n1v, double n2v)
{
  double n13 = 1. - n3, hc = n13 * n13;
  double phi1 = -n0 * log(n13);
  double phi2 = (n1 * n2 - n2v * n1v) / n13;
  double phi3 = 0.0;

  if (n3 != 0) phi3 = n2 * (n2 * n2 - 3 * n2v * n2v) * (n3 + hc * log(n13)) / (36.0 * PI * n3 * n3 * hc);

  return phi1 + phi2 + phi3;
}

double wb_d2(double n0, double n1, double n2, double n3, double n1v, double n2v, double hR)
{
  double d21, d22, d23, ln3 = 1. - n3, hc = ln3 * n3;
  (void)n0;
  (void)n1v;

  if (n3 < ETA_SMALL)
    d21 = n1 / ln3 + (n2 * n2 - n2v * n2v) / PI * (1. / 8. + 2. * n3 / 9. + 5. * n3 * n3 / 16. + 2. * n3 * n3 * n3 / 5.);
  else
    d21 = n1 / ln3 + (n2 * n2 - n2v * n2v) * (n3 + ln3 * ln3 * log(ln3)) / (12. * PI * hc * hc);

  d22 = n2 / ln3;
  d23 = -log(ln3);

  return d21 + d22 / (PI_4 * hR) + d23 / (PI_4 * hR * hR);
}

double wb_d3(double n0, double n1, double n2, double n3, double n1v, double n2v)
{
  double d13, d23, ln3 = 1. - n3, hc = ln3 * n3;

  d13 = n0 / ln3 + (n1 * n2 - n1v * n2v) / (ln3 * ln3);

  if (n3 < ETA_SMALL)
    d23 = n2 * (n2 * n2 - 3. * n2v * n2v) / PI * (2. / 27. + 5. * n3 / 48. + 2. * n3 * n3 / 15. + 35. * n3 * n3 * n3 / 54.);
  else
    d23 = -n2 * (n2 * n2 - 3. * n2v * n2v) * (n3 * (2. - 5. * n3 + n3 * n3) + 2. * ln3 * ln3 * ln3 * log(ln3))
          / (36. * PI * hc * hc * hc);

  return d13 + d23;
}

double wb_dv2(double n0, double n1, double n2, double n3, double n1v, double n2v, double hR)
{
  double dv21, dv22, ln3 = 1. - n3, hc = ln3 * n3;
  (void)n0;
  (void)n1;

  if (n3 < ETA_SMALL)
    dv21 = n1v / ln3 + 2. * n2 * n2v / PI * (1. / 8. + 2. * n3 / 9. + 5. * n3 * n3 / 16. + 2. * n3 * n3 * n3 / 5.);
  else
    dv21 = n1v / ln3 + n2 * n2v * (n3 + ln3 * ln3 * log(ln3)) / (6. * PI * hc * hc);

  dv22 = n2v / ln3;

  return -(dv21 + dv22 / (PI_4 * hR));
}
