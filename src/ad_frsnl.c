
#include <math.h>
#include <float.h>
#include <stdio.h>
#include "ad_frsnl.h"


static double Fresnel (double n_i, double n_t, double mu_i);

static double R1 (double ni, double nt);



double
Cos_Critical_Angle (double ni, double nt)
{
  double x;

  if (nt >= ni)
    return 0.0;
  else
    {
      x = nt / ni;
      x = sqrt (1.0 - x * x);
      return x;
    }
}




double
Cos_Snell (double n_i, double mu_i, double n_t)
{
  double temp;

  if (mu_i == 1.0)
    return 1.0;

  if (n_i == n_t)
    return mu_i;

  temp = n_i / n_t;
  temp = 1.0 - temp * temp * (1.0 - mu_i * mu_i);
  if (temp < 0)
    return 0.0;
  else
    return (sqrt (temp));
}




static double
Fresnel (double n_i, double n_t, double mu_i)
{
  double mu_t, ratio, temp, temp1;

  if (n_i == n_t)
    return 0.0;

  if (mu_i == 1.0)
    {
      temp = (n_i - n_t) / (n_i + n_t);
      return (temp * temp);
    }

  if (mu_i == 0.0)
    return 1.0;

  mu_t = Cos_Snell (n_i, mu_i, n_t);
  if (mu_t == 0.0)
    return 1.0;
  ratio = n_i / n_t;
  temp = ratio * mu_t;
  temp1 = (mu_i - temp) / (mu_i + temp);
  temp = ratio * mu_i;
  temp = (mu_t - temp) / (mu_t + temp);
  return ((temp1 * temp1 + temp * temp) / 2);
}




double
Glass (double n_i, double n_g, double n_t, double mu_i)
{
  double r1, r2, mu_g, temp;

  if (n_i == n_g)
    return (Fresnel (n_g, n_t, mu_i));

  r1 = Fresnel (n_i, n_g, mu_i);
  if (r1 >= 1.0 || n_g == n_t)
    return r1;

  mu_g = Cos_Snell (n_i, mu_i, n_g);
  r2 = Fresnel (n_g, n_t, mu_g);
  temp = r1 * r2;
  temp = (r1 + r2 - 2 * temp) / (1 - temp);
  return temp;
}




void
Absorbing_Glass_RT (double n_i, double n_g, double n_t, double mu_i, double b,
		    double *r, double *t)
{
  double r1, r2, mu_g, expo, denom;
  *t = 0;

  *r = Fresnel (n_i, n_g, mu_i);
  if (*r >= 1.0 || b == HUGE_VAL || mu_i == 0.0)
    return;

  mu_g = Cos_Snell (n_i, mu_i, n_g);
  r1 = *r;
  r2 = Fresnel (n_g, n_t, mu_g);

  if (b == 0.0)
    {
      *r = (r1 + r2 - 2.0 * r1 * r2) / (1 - r1 * r2);
      *t = 1.0 - (*r);
    }
  else
    {
      expo = -b / mu_g;
      if (2 * expo <= DBL_MIN_10_EXP * 2.3025851)
	return;
      expo = exp (expo);

      denom = 1.0 - r1 * r2 * expo * expo;
      *r = (r1 + (1.0 - 2.0 * r1) * r2 * expo * expo) / denom;
      *t = (1.0 - r1) * (1.0 - r2) * expo / denom;
    }
}





static double
R1 (double ni, double nt)
{
  double m, m2, m4, mm1, mp1, r, temp;

  if (ni == nt)
    return 0.0;

  if (ni < nt)
    m = nt / ni;
  else
    m = ni / nt;

  m2 = m * m;
  m4 = m2 * m2;
  mm1 = m - 1;
  mp1 = m + 1;
  temp = (m2 - 1) / (m2 + 1);

  r = 0.5 + mm1 * (3 * m + 1) / 6 / mp1 / mp1;
  r += m2 * temp * temp / (m2 + 1) * log (mm1 / mp1);
  r -= 2 * m * m2 * (m2 + 2 * m - 1) / (m2 + 1) / (m4 - 1);
  r += 8 * m4 * (m4 + 1) / (m2 + 1) / (m4 - 1) / (m4 - 1) * log (m);

  if (ni < nt)
    return r;
  else
    return (1 - (1 - r) / m2);
}




void
Sp_mu_RT (double n_top, double n_slab, double n_bottom,
	  double tau_top, double tau_slab, double tau_bottom, double mu,
	  double *r, double *t)
{
  double r_top, r_bottom, t_top, t_bottom, mu_slab, beer, denom, temp,
    mu_in_slab;
  *r = 0;
  *t = 0;
  Absorbing_Glass_RT (1.0, n_top, n_slab, mu, tau_top, &r_top, &t_top);

  mu_in_slab = Cos_Snell (1.0, mu, n_slab);
  Absorbing_Glass_RT (n_slab, n_bottom, 1.0, mu_in_slab, tau_bottom,
		      &r_bottom, &t_bottom);


  mu_slab = Cos_Snell (1.0, mu, n_slab);

  if (mu_slab == 0)
    beer = 0.0;
  else if (tau_slab == HUGE_VAL)
    beer = 0.0;
  else
    {
      temp = -tau_slab / mu_slab;
      if (2 * temp <= DBL_MIN_10_EXP * 2.3025851)
	beer = 0.0;
      else
	beer = exp (temp);
    }



  if (beer == 0.0)
    {
      *r = r_top;
    }
  else
    {
      temp = t_top * beer;
      denom = 1 - r_top * r_bottom * beer * beer;
      *r = r_top + r_bottom * temp * temp / denom;
      *t = t_bottom * temp / denom;
    }


}




void
Sp_mu_RT_Flip (int flip, double n_top, double n_slab, double n_bottom,
	       double tau_top, double tau_slab, double tau_bottom, double mu,
	       double *r, double *t)
{
  Sp_mu_RT (n_top, n_slab, n_bottom, tau_top, tau_slab, tau_bottom, mu, r, t);
  if (flip && n_top != n_bottom && tau_top != tau_bottom)
    {
      double correct_r = *r;
      Sp_mu_RT (n_bottom, n_slab, n_top, tau_bottom, tau_slab, tau_top, mu, r,
		t);
      *r = correct_r;
    }
}




double
Diffuse_Glass_R (double nair, double nslide, double nslab)
{
  double rairglass, rglasstissue, rtemp;

  rairglass = R1 (nair, nslide);
  rglasstissue = R1 (nslide, nslab);
  rtemp = rairglass * rglasstissue;
  if (rtemp >= 1)
    return 1.0;
  else
    return ((rairglass + rglasstissue - 2 * rtemp) / (1 - rtemp));
}
