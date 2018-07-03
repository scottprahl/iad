
#include <math.h>
#include <float.h>
#include <stdio.h>
#include "nr_util.h"
#include "ad_globl.h"
#include "ad_frsnl.h"
#include "ad_bound.h"
#include "iad_type.h"
#include "iad_calc.h"
#include "iad_pub.h"
#include "iad_util.h"
unsigned long g_util_debugging = 0;

#define BIG_A_VALUE 999999.0
#define SMALL_A_VALUE 0.000001 \





double
What_Is_B (struct AD_slab_type slab, double Tc)
{
  double r1, r2, t1, t2, mu_in_slab;


  Absorbing_Glass_RT (1.0, slab.n_top_slide, slab.n_slab,
		      slab.cos_angle, slab.b_top_slide, &r1, &t1);

  mu_in_slab = Cos_Snell (1.0, slab.cos_angle, slab.n_slab);

  Absorbing_Glass_RT (slab.n_slab, slab.n_bottom_slide, 1.0,
		      mu_in_slab, slab.b_bottom_slide, &r2, &t2);




  if (Tc <= 0)
    return (HUGE_VAL);

  if (Tc >= t1 * t2 / (1 - r1 * r2))
    return (0.001);



  if (r1 == 0 || r2 == 0)
    return (-slab.cos_angle * log (Tc / t1 / t2));





  {
    double B;

    B = t1 * t2;
    return (-slab.cos_angle *
	    log (2 * Tc / (B + sqrt (B * B + 4 * Tc * Tc * r1 * r2))));
  }


}





void
Estimate_RT (struct measure_type m, struct invert_type r, double *rt,
	     double *tt, double *rd, double *rc, double *td, double *tc)
{


  Calculate_Minimum_MR (m, r, rc, tc);




  if (m.fraction_of_rc_in_mr)
    {
      *rt = m.m_r;
      *rd = *rt - m.fraction_of_rc_in_mr * (*rc);
      if (Debug (DEBUG_SEARCH))
	{
	  fprintf (stderr, "        rt = %.5f\n", *rt);
	  fprintf (stderr, "    est rd = %.5f\n", *rd);
	}
      if (*rd < 0)
	{
	  *rd = 0;
	  *rc = *rt;
	}
    }
  else
    {
      *rd = m.m_r;
      *rt = *rd + *rc;
    }



  if (m.num_measures == 1)
    {
      *tt = 0.0;
      *td = 0.0;
    }
  else if (m.fraction_of_tc_in_mt)
    {
      *tt = m.m_t;
      *td = *tt - *tc;
      if (*td < 0)
	{
	  *tc = *tt;
	  *td = 0;
	}
    }
  else
    {
      *td = m.m_t;
      *tt = *td + *tc;
    }


}




double
a2acalc (double a)
{
  if (a <= 0)
    return -BIG_A_VALUE;

  if (a >= 1)
    return BIG_A_VALUE;
  return ((2 * a - 1) / a / (1 - a));
}




double
acalc2a (double acalc)
{
  if (acalc == BIG_A_VALUE)
    return 1.0;
  else if (acalc == -BIG_A_VALUE)
    return 0.0;
  else if (fabs (acalc) < SMALL_A_VALUE)
    return 0.5;
  else
    return ((-2 + acalc + sqrt (acalc * acalc + 4)) / (2 * acalc));
}




double
g2gcalc (double g)
{
  if (g <= -1)
    return (-HUGE_VAL);

  if (g >= 1)
    return (HUGE_VAL);

  return (g / (1 - fabs (g)));
}




double
gcalc2g (double gcalc)
{
  if (gcalc == -HUGE_VAL)
    return -1.0;
  if (gcalc == HUGE_VAL)
    return 1.0;
  return (gcalc / (1 + fabs (gcalc)));
}




double
b2bcalc (double b)
{
  if (b == HUGE_VAL)
    return HUGE_VAL;
  if (b <= 0)
    return 0.0;
  return (log (b));
}




double
bcalc2b (double bcalc)
{
  if (bcalc == HUGE_VAL)
    return HUGE_VAL;
  if (bcalc > 2.3 * DBL_MAX_10_EXP)
    return HUGE_VAL;
  return (exp (bcalc));
}




void
twoprime (double a, double b, double g, double *ap, double *bp)
{
  if (a == 1 && g == 1)
    *ap = 0.0;
  else
    *ap = (1 - g) * a / (1 - a * g);

  if (b == HUGE_VAL)
    *bp = HUGE_VAL;
  else
    *bp = (1 - a * g) * b;
}




void
twounprime (double ap, double bp, double g, double *a, double *b)
{
  *a = ap / (1 - g + ap * g);
  if (bp == HUGE_VAL)
    *b = HUGE_VAL;
  else
    *b = (1 + ap * g / (1 - g)) * bp;
}




void
abgg2ab (double a1, double b1, double g1, double g2, double *a2, double *b2)
{
  double a, b;

  twoprime (a1, b1, g1, &a, &b);
  twounprime (a, b, g2, a2, b2);
}




void
abgb2ag (double a1, double b1, double b2, double *a2, double *g2)
{
  if (b1 == 0 || b2 == 0)
    {
      *a2 = a1;
      *g2 = 0;
    }

  if (b2 < b1)
    b2 = b1;

  if (a1 == 0)
    *a2 = 0.0;
  else
    {
      if (a1 == 1)
	*a2 = 1.0;
      else
	{
	  if (b1 == 0 || b2 == HUGE_VAL)
	    *a2 = a1;
	  else
	    *a2 = 1 + b1 / b2 * (a1 - 1);
	}
    }
  if (*a2 == 0 || b2 == 0 || b2 == HUGE_VAL)
    *g2 = 0.5;
  else
    *g2 = (1 - b1 / b2) / (*a2);
}




void
quick_guess (struct measure_type m, struct invert_type r, double *a,
	     double *b, double *g)
{
  double UR1, UT1, rd, td, tc, rc, bprime, aprime, alpha, beta, logr;

  Estimate_RT (m, r, &UR1, &UT1, &rd, &rc, &td, &tc);

  if (UT1 == 1)
    aprime = 1.0;
  else if (rd / (1 - UT1) >= 0.1)
    {
      double tmp = (1 - rd - UT1) / (1 - UT1);
      aprime = 1 - 4.0 / 9.0 * tmp * tmp;
    }
  else if (rd < 0.05 && UT1 < 0.4)
    aprime = 1 - (1 - 10 * rd) * (1 - 10 * rd);
  else if (rd < 0.1 && UT1 < 0.4)
    aprime = 0.5 + (rd - 0.05) * 4;
  else
    {
      double tmp = (1 - 4 * rd - UT1) / (1 - UT1);
      aprime = 1 - tmp * tmp;
    }



  switch (m.num_measures)
    {
    case 1:

      *g = r.default_g;
      *a = aprime / (1 - *g + aprime * (*g));
      *b = HUGE_VAL;


      break;
    case 2:


      if (rd < 0.01)
	{
	  bprime = What_Is_B (r.slab, UT1);
	  fprintf (stderr, "low rd<0.01! ut1=%f aprime=%f bprime=%f\n", UT1,
		   aprime, bprime);
	}
      else if (UT1 <= 0)
	bprime = HUGE_VAL;
      else if (UT1 > 0.1)
	bprime = 2 * exp (5 * (rd - UT1) * log (2.0));
      else
	{
	  alpha = 1 / log (0.05 / 1.0);
	  beta = log (1.0) / log (0.05 / 1.0);
	  logr = log (UR1);
	  bprime = log (UT1) - beta * log (0.05) + beta * logr;
	  bprime /= alpha * log (0.05) - alpha * logr - 1;
	}



      *g = r.default_g;
      *a = aprime / (1 - *g + aprime ** g);
      *b = bprime / (1 - *a ** g);


      break;
    case 3:


      switch (r.search)
	{
	case FIND_A:


	  *g = r.default_g;
	  *a = aprime / (1 - *g + aprime ** g);
	  *b = What_Is_B (r.slab, m.m_u);


	  break;
	case FIND_B:


	  *g = r.default_g;
	  *a = 0.0;
	  *b = What_Is_B (r.slab, m.m_u);


	  break;
	case FIND_AB:


	  *g = r.default_g;

	  if (*g == 1)
	    *a = 0.0;
	  else
	    *a = aprime / (1 - *g + aprime ** g);


	  if (rd < 0.01)
	    {
	      bprime = What_Is_B (r.slab, UT1);
	      fprintf (stderr, "low rd<0.01! ut1=%f aprime=%f bprime=%f\n",
		       UT1, aprime, bprime);
	    }
	  else if (UT1 <= 0)
	    bprime = HUGE_VAL;
	  else if (UT1 > 0.1)
	    bprime = 2 * exp (5 * (rd - UT1) * log (2.0));
	  else
	    {
	      alpha = 1 / log (0.05 / 1.0);
	      beta = log (1.0) / log (0.05 / 1.0);
	      logr = log (UR1);
	      bprime = log (UT1) - beta * log (0.05) + beta * logr;
	      bprime /= alpha * log (0.05) - alpha * logr - 1;
	    }


	  if (bprime == HUGE_VAL || *a ** g == 1)
	    *b = HUGE_VAL;
	  else
	    *b = bprime / (1 - *a ** g);


	  break;
	case FIND_AG:

	  *b = What_Is_B (r.slab, m.m_u);
	  if (*b == HUGE_VAL || *b == 0)
	    {
	      *a = aprime;
	      *g = r.default_g;
	    }
	  else
	    {

	      if (rd < 0.01)
		{
		  bprime = What_Is_B (r.slab, UT1);
		  fprintf (stderr,
			   "low rd<0.01! ut1=%f aprime=%f bprime=%f\n", UT1,
			   aprime, bprime);
		}
	      else if (UT1 <= 0)
		bprime = HUGE_VAL;
	      else if (UT1 > 0.1)
		bprime = 2 * exp (5 * (rd - UT1) * log (2.0));
	      else
		{
		  alpha = 1 / log (0.05 / 1.0);
		  beta = log (1.0) / log (0.05 / 1.0);
		  logr = log (UR1);
		  bprime = log (UT1) - beta * log (0.05) + beta * logr;
		  bprime /= alpha * log (0.05) - alpha * logr - 1;
		}


	      *a = 1 + bprime * (aprime - 1) / (*b);
	      if (*a < 0.1)
		*g = 0.0;
	      else
		*g = (1 - bprime / (*b)) / (*a);
	    }


	  break;
	}


      break;
    }


  if (*a < 0)
    *a = 0.0;
  if (*g < 0)
    *g = 0.0;
  else if (*g >= 1)
    *g = 0.5;


}




void
Set_Debugging (unsigned long debug_level)
{
  g_util_debugging = debug_level;
}




int
Debug (unsigned long mask)
{
  if (g_util_debugging & mask)
    return 1;
  else
    return 0;
}




void
Print_Invert_Type (struct invert_type r)
{
  fprintf (stderr, "\n");
  fprintf (stderr, "default  a=%10.5f   b=%10.5f    g=%10.5f\n",
	   r.default_a, r.default_b, r.default_g);
  fprintf (stderr, "slab     a=%10.5f   b=%10.5f    g=%10.5f\n",
	   r.slab.a, r.slab.b, r.slab.g);
  fprintf (stderr, "n      top=%10.5f mid=%10.5f  bot=%10.5f\n",
	   r.slab.n_top_slide, r.slab.n_slab, r.slab.n_bottom_slide);
  fprintf (stderr, "thick  top=%10.5f cos=%10.5f  bot=%10.5f\n",
	   r.slab.b_top_slide, r.slab.cos_angle, r.slab.b_bottom_slide);
  fprintf (stderr, "search = %d quadrature points = %d\n", r.search,
	   r.method.quad_pts);
}




void
Print_Measure_Type (struct measure_type m)
{
  fprintf (stderr, "\n");
  fprintf (stderr, "#                        Beam diameter = %7.1f mm\n",
	   m.d_beam);
  fprintf (stderr, "#                     Sample thickness = %7.1f mm\n",
	   m.slab_thickness);
  fprintf (stderr, "#                  Top slide thickness = %7.1f mm\n",
	   m.slab_top_slide_thickness);
  fprintf (stderr, "#               Bottom slide thickness = %7.1f mm\n",
	   m.slab_bottom_slide_thickness);
  fprintf (stderr, "#           Sample index of refraction = %7.3f\n",
	   m.slab_index);
  fprintf (stderr, "#        Top slide index of refraction = %7.3f\n",
	   m.slab_top_slide_index);
  fprintf (stderr, "#     Bottom slide index of refraction = %7.3f\n",
	   m.slab_bottom_slide_index);
  fprintf (stderr, "#    Fraction unscattered light in M_R = %7.1f %%\n",
	   m.fraction_of_rc_in_mr * 100);
  fprintf (stderr, "#    Fraction unscattered light in M_T = %7.1f %%\n",
	   m.fraction_of_tc_in_mt * 100);
  fprintf (stderr, "# \n");
  fprintf (stderr, "# Reflection sphere\n");
  fprintf (stderr, "#                      sphere diameter = %7.1f mm\n",
	   m.d_sphere_r);
  fprintf (stderr, "#                 sample port diameter = %7.1f mm\n",
	   2 * m.d_sphere_r * sqrt (m.as_r));
  fprintf (stderr, "#               entrance port diameter = %7.1f mm\n",
	   2 * m.d_sphere_r * sqrt (m.ae_r));
  fprintf (stderr, "#               detector port diameter = %7.1f mm\n",
	   2 * m.d_sphere_r * sqrt (m.ad_r));
  fprintf (stderr, "#                     wall reflectance = %7.1f %%\n",
	   m.rw_r * 100);
  fprintf (stderr, "#                 standard reflectance = %7.1f %%\n",
	   m.rstd_r * 100);
  fprintf (stderr, "#                 detector reflectance = %7.1f %%\n",
	   m.rd_r * 100);
  fprintf (stderr, "#                              spheres = %7d\n",
	   m.num_spheres);
  fprintf (stderr, "#                             measures = %7d\n",
	   m.num_measures);
  fprintf (stderr, "#                               method = %7d\n",
	   m.method);
  fprintf (stderr, "area_r as=%10.5f  ad=%10.5f    ae=%10.5f  aw=%10.5f\n",
	   m.as_r, m.ad_r, m.ae_r, m.aw_r);
  fprintf (stderr, "refls  rd=%10.5f  rw=%10.5f  rstd=%10.5f   f=%10.5f\n",
	   m.rd_r, m.rw_r, m.rstd_r, m.f_r);
  fprintf (stderr, "area_t as=%10.5f  ad=%10.5f    ae=%10.5f  aw=%10.5f\n",
	   m.as_t, m.ad_t, m.ae_t, m.aw_t);
  fprintf (stderr, "refls  rd=%10.5f  rw=%10.5f  rstd=%10.5f   f=%10.5f\n",
	   m.rd_t, m.rw_t, m.rstd_t, m.f_t);
  fprintf (stderr, "lost  ur1=%10.5f ut1=%10.5f   uru=%10.5f  utu=%10.5f\n",
	   m.ur1_lost, m.ut1_lost, m.utu_lost, m.utu_lost);
}
