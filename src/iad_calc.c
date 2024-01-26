

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "nr_util.h"
#include "nr_zbrent.h"
#include "ad_globl.h"
#include "ad_frsnl.h"
#include "ad_prime.h"
#include "iad_type.h"
#include "iad_util.h"
#include "iad_calc.h"

#define ABIT 1e-6
#define A_COLUMN 1
#define B_COLUMN 2
#define G_COLUMN 3
#define URU_COLUMN 4
#define UTU_COLUMN 5
#define UR1_COLUMN 6
#define UT1_COLUMN 7
#define REFLECTION_SPHERE 1
#define TRANSMISSION_SPHERE 0
#define GRID_SIZE 101
#define T_TRUST_FACTOR 2

static int CALCULATING_GRID = 1;
static struct measure_type MM;
static struct invert_type RR;
static struct measure_type MGRID;
static struct invert_type RGRID;
static double **The_Grid = NULL;
static double GG_a;
static double GG_b;
static double GG_g;
static double GG_bs;
static double GG_ba;
static boolean_type The_Grid_Initialized = FALSE;
static boolean_type The_Grid_Search = -1;



void
Set_Calc_State (struct measure_type m, struct invert_type r)
{
  memcpy (&MM, &m, sizeof (struct measure_type));
  memcpy (&RR, &r, sizeof (struct invert_type));
  if (Debug (DEBUG_ITERATIONS) && !CALCULATING_GRID)
    {
      fprintf (stderr, "UR1 loss=%g, UT1 loss=%g  ", m.ur1_lost, m.ut1_lost);
      fprintf (stderr, "URU loss=%g, UTU loss=%g\n", m.uru_lost, m.utu_lost);
    }
}




void
Get_Calc_State (struct measure_type *m, struct invert_type *r)
{
  memcpy (m, &MM, sizeof (struct measure_type));
  memcpy (r, &RR, sizeof (struct invert_type));
}




boolean_type
Same_Calc_State (struct measure_type m, struct invert_type r)
{
  if (The_Grid == NULL)
    return FALSE;
  if (!The_Grid_Initialized)
    return FALSE;

  if (r.search != RR.search)
    return FALSE;
  if (r.method.quad_pts != RR.method.quad_pts)
    return FALSE;
  if (r.slab.a != RR.slab.a)
    return FALSE;
  if (r.slab.b != RR.slab.b)
    return FALSE;
  if (r.slab.g != RR.slab.g)
    return FALSE;
  if (r.slab.phase_function != RR.slab.phase_function)
    return FALSE;
  if (r.slab.n_slab != RR.slab.n_slab)
    return FALSE;
  if (r.slab.n_top_slide != RR.slab.n_top_slide)
    return FALSE;
  if (r.slab.n_bottom_slide != RR.slab.n_bottom_slide)
    return FALSE;
  if (r.slab.b_top_slide != RR.slab.b_top_slide)
    return FALSE;
  if (r.slab.b_bottom_slide != RR.slab.b_bottom_slide)
    return FALSE;
  if (r.slab.cos_angle != RR.slab.cos_angle)
    return FALSE;
  if ((m.num_measures == 3) && (m.m_u != MGRID.m_u))
    return (FALSE);
  return TRUE;
}




void Fill_AB_Grid (struct measure_type m, struct invert_type r);

void Fill_AG_Grid (struct measure_type m, struct invert_type r);



void
RT_Flip (int flip, int n, struct AD_slab_type *slab, double *UR1, double *UT1,
	 double *URU, double *UTU)
{
  double swap, correct_UR1, correct_URU;
  RT (n, slab, UR1, UT1, URU, UTU);
  if (flip)
    {
      correct_UR1 = *UR1;
      correct_URU = *URU;
      swap = slab->n_top_slide;
      slab->n_top_slide = slab->n_bottom_slide;
      slab->n_bottom_slide = swap;
      swap = slab->b_top_slide;
      slab->b_top_slide = slab->b_bottom_slide;
      slab->b_bottom_slide = swap;
      RT (n, slab, UR1, UT1, URU, UTU);
      swap = slab->n_top_slide;
      slab->n_top_slide = slab->n_bottom_slide;
      slab->n_bottom_slide = swap;
      swap = slab->b_top_slide;
      slab->b_top_slide = slab->b_bottom_slide;
      slab->b_bottom_slide = swap;
      *UR1 = correct_UR1;
      *URU = correct_URU;
    }
}




void
Allocate_Grid (search_type s)
{
  The_Grid = dmatrix (0, GRID_SIZE * GRID_SIZE, 1, 7);
  if (The_Grid == NULL)
    AD_error ("unable to allocate the grid matrix");
  The_Grid_Initialized = FALSE;
}




boolean_type
Valid_Grid (struct measure_type m, search_type s)
{

  if (The_Grid == NULL)
    {
      if (Debug (DEBUG_GRID))
	fprintf (stderr, "GRID: Fill because NULL\n");
      return (FALSE);
    }
  if (!The_Grid_Initialized)
    {
      if (Debug (DEBUG_GRID))
	fprintf (stderr, "GRID: Fill because not initialized\n");
      return (FALSE);
    }


  if (The_Grid_Search != s)
    {
      if (Debug (DEBUG_GRID))
	fprintf (stderr, "GRID: Fill because search type changed\n");
      return (FALSE);
    }



  if ((m.num_measures == 3) && (m.m_u != MGRID.m_u))
    {
      if (Debug (DEBUG_GRID))
	fprintf (stderr, "GRID: Fill because unscattered light changed\n");
      return (FALSE);
    }


  if (m.slab_index != MGRID.slab_index)
    {
      if (Debug (DEBUG_GRID))
	fprintf (stderr, "GRID: Fill slab refractive index changed\n");
      return (FALSE);
    }
  if (m.slab_cos_angle != MGRID.slab_cos_angle)
    {
      if (Debug (DEBUG_GRID))
	fprintf (stderr, "GRID: Fill incident light changed\n");
      return (FALSE);
    }

  if (m.slab_top_slide_index != MGRID.slab_top_slide_index)
    {
      if (Debug (DEBUG_GRID))
	fprintf (stderr, "GRID: Fill top slide refractive index changed\n");
      return (FALSE);
    }

  if (m.slab_bottom_slide_index != MGRID.slab_bottom_slide_index)
    {
      if (Debug (DEBUG_GRID))
	fprintf (stderr,
		 "GRID: Fill bottom slide refractive index changed\n");
      return (FALSE);
    }




  return (TRUE);
}



static void
fill_grid_entry (int i, int j)
{
  double ur1, ut1, uru, utu;

  if (RR.slab.b <= 1e-6)
    RR.slab.b = 1e-6;
  if (Debug (DEBUG_EVERY_CALC))
    {
      if (!CALCULATING_GRID)
	fprintf (stderr, "a=%8.5f b=%10.5f g=%8.5f ", RR.slab.a, RR.slab.b,
		 RR.slab.g);
      else
	{
	  if (j == 0)
	    fprintf (stderr, ".");
	  if (i + 1 == GRID_SIZE && j == 0)
	    fprintf (stderr, "\n");
	}
    }

  RT_Flip (MM.flip_sample, RR.method.quad_pts, &RR.slab, &ur1, &ut1, &uru,
	   &utu);

  if (Debug (DEBUG_EVERY_CALC) && !CALCULATING_GRID)
    fprintf (stderr, "ur1=%8.5f ut1=%8.5f\n", ur1, ut1);

  The_Grid[GRID_SIZE * i + j][A_COLUMN] = RR.slab.a;
  The_Grid[GRID_SIZE * i + j][B_COLUMN] = RR.slab.b;
  The_Grid[GRID_SIZE * i + j][G_COLUMN] = RR.slab.g;
  The_Grid[GRID_SIZE * i + j][UR1_COLUMN] = ur1;
  The_Grid[GRID_SIZE * i + j][UT1_COLUMN] = ut1;
  The_Grid[GRID_SIZE * i + j][URU_COLUMN] = uru;
  The_Grid[GRID_SIZE * i + j][UTU_COLUMN] = utu;

  if (Debug (DEBUG_GRID_CALC))
    {
      fprintf (stderr, "+ %2d %2d ", i, j);
      fprintf (stderr, "%10.5f %10.5f %10.5f |", RR.slab.a, RR.slab.b,
	       RR.slab.g);
      fprintf (stderr, "%10.5f %10.5f |", MM.m_r, uru);
      fprintf (stderr, "%10.5f %10.5f \n", MM.m_t, utu);
    }
}




void
Fill_Grid (struct measure_type m, struct invert_type r, int force_new)
{
  if (force_new || !Same_Calc_State (m, r))
    {
      switch (r.search)
	{
	case FIND_AB:
	  if (Debug (DEBUG_SEARCH))
	    fprintf (stderr, "filling AB Grid\n");
	  Fill_AB_Grid (m, r);
	  break;
	case FIND_AG:
	  if (Debug (DEBUG_SEARCH))
	    fprintf (stderr, "filling AG Grid\n");
	  Fill_AG_Grid (m, r);
	  break;
	case FIND_BG:
	  if (Debug (DEBUG_SEARCH))
	    fprintf (stderr, "filling BG Grid\n");
	  Fill_BG_Grid (m, r);
	  break;
	case FIND_BaG:
	  if (Debug (DEBUG_SEARCH))
	    fprintf (stderr, "filling BaG Grid\n");
	  Fill_BaG_Grid (m, r);
	  break;
	case FIND_BsG:
	  if (Debug (DEBUG_SEARCH))
	    fprintf (stderr, "filling BsG Grid\n");
	  Fill_BsG_Grid (m, r);
	  break;
	default:
	  AD_error ("Attempt to fill grid for unusual search case.");
	}
    }

  Get_Calc_State (&MGRID, &RGRID);
}




void
Near_Grid_Points (double r, double t, search_type s, int *i_min, int *j_min)
{
  int i, j;
  double fval;
  double smallest = 10.0;
  struct measure_type old_mm;
  struct invert_type old_rr;

  Get_Calc_State (&old_mm, &old_rr);

  *i_min = 0;
  *j_min = 0;
  for (i = 0; i < GRID_SIZE; i++)
    {
      for (j = 0; j < GRID_SIZE; j++)
	{

	  CALCULATING_GRID = 1;
	  fval = Calculate_Grid_Distance (i, j);
	  CALCULATING_GRID = 0;

	  if (fval < smallest)
	    {
	      *i_min = i;
	      *j_min = j;
	      smallest = fval;
	    }
	}
    }

  Set_Calc_State (old_mm, old_rr);
}




void
Fill_AB_Grid (struct measure_type m, struct invert_type r)
{
  int i, j;
  double a;
  double min_b = -8;
  double max_b = +8;

  if (Debug (Debug (DEBUG_GRID)))
    fprintf (stderr, "Filling AB grid\n");

  if (The_Grid == NULL)
    Allocate_Grid (r.search);

  GG_a = 0.0;
  GG_b = 0.0;
  GG_g = 0.0;
  GG_bs = 0.0;
  GG_ba = 0.0;




  Set_Calc_State (m, r);

  GG_g = RR.slab.g;
  for (i = 0; i < GRID_SIZE; i++)
    {
      double x = (double) i / (GRID_SIZE - 1.0);

      RR.slab.b = exp (min_b + (max_b - min_b) * x);
      for (j = 0; j < GRID_SIZE; j++)
	{

	  a = (double) j / (GRID_SIZE - 1.0);
	  if (a < 0.25)
	    RR.slab.a = 1.0 - a * a;
	  else if (a > 0.75)
	    RR.slab.a = (1.0 - a) * (1.0 - a);
	  else
	    RR.slab.a = 1 - a;


	  a = (double) j / (GRID_SIZE - 1.0);
	  RR.slab.a = (1.0 - a * a) * (1.0 - a) + (1.0 - a) * (1.0 - a) * a;



	  fill_grid_entry (i, j);
	}
    }

  The_Grid_Initialized = TRUE;
  The_Grid_Search = FIND_AB;
}




void
Fill_AG_Grid (struct measure_type m, struct invert_type r)
{
  int i, j;
  double a;

  if (Debug (Debug (DEBUG_GRID)))
    fprintf (stderr, "Filling AG grid\n");

  if (The_Grid == NULL)
    Allocate_Grid (r.search);

  GG_a = 0.0;
  GG_b = 0.0;
  GG_g = 0.0;
  GG_bs = 0.0;
  GG_ba = 0.0;




  Set_Calc_State (m, r);
  GG_b = r.slab.b;
  for (i = 0; i < GRID_SIZE; i++)
    {
      RR.slab.g = 0.9999 * (2.0 * i / (GRID_SIZE - 1.0) - 1.0);
      for (j = 0; j < GRID_SIZE; j++)
	{


	  a = (double) j / (GRID_SIZE - 1.0);
	  if (a < 0.25)
	    RR.slab.a = 1.0 - a * a;
	  else if (a > 0.75)
	    RR.slab.a = (1.0 - a) * (1.0 - a);
	  else
	    RR.slab.a = 1 - a;


	  a = (double) j / (GRID_SIZE - 1.0);
	  RR.slab.a = (1.0 - a * a) * (1.0 - a) + (1.0 - a) * (1.0 - a) * a;


	  fill_grid_entry (i, j);
	}
    }

  The_Grid_Initialized = TRUE;
  The_Grid_Search = FIND_AG;
}




void
Fill_BG_Grid (struct measure_type m, struct invert_type r)
{
  int i, j;

  if (The_Grid == NULL)
    Allocate_Grid (r.search);

  GG_a = 0.0;
  GG_b = 0.0;
  GG_g = 0.0;
  GG_bs = 0.0;
  GG_ba = 0.0;




  if (Debug (Debug (DEBUG_GRID)))
    fprintf (stderr, "Filling BG grid\n");

  Set_Calc_State (m, r);
  RR.slab.b = 1.0 / 32.0;
  RR.slab.a = RR.default_a;
  GG_a = RR.slab.a;

  for (i = 0; i < GRID_SIZE; i++)
    {
      RR.slab.b *= 2;
      for (j = 0; j < GRID_SIZE; j++)
	{
	  RR.slab.g = 0.9999 * (2.0 * j / (GRID_SIZE - 1.0) - 1.0);
	  fill_grid_entry (i, j);
	}
    }

  The_Grid_Initialized = TRUE;
  The_Grid_Search = FIND_BG;
}




void
Fill_BaG_Grid (struct measure_type m, struct invert_type r)
{
  int i, j;
  double bs, ba;

  if (The_Grid == NULL)
    Allocate_Grid (r.search);

  GG_a = 0.0;
  GG_b = 0.0;
  GG_g = 0.0;
  GG_bs = 0.0;
  GG_ba = 0.0;




  if (Debug (Debug (DEBUG_GRID)))
    fprintf (stderr, "Filling BaG grid\n");

  Set_Calc_State (m, r);
  ba = 1.0 / 32.0;
  bs = RR.default_bs;
  GG_bs = bs;
  for (i = 0; i < GRID_SIZE; i++)
    {
      ba *= 2;
      ba = exp ((double) i / (GRID_SIZE - 1.0) * log (1024.0)) / 16.0;
      RR.slab.b = ba + bs;
      if (RR.slab.b > 0)
	RR.slab.a = bs / RR.slab.b;
      else
	RR.slab.a = 0;
      for (j = 0; j < GRID_SIZE; j++)
	{
	  RR.slab.g = 0.9999 * (2.0 * j / (GRID_SIZE - 1.0) - 1.0);
	  fill_grid_entry (i, j);
	}
    }

  The_Grid_Initialized = TRUE;
  The_Grid_Search = FIND_BaG;
}




void
Fill_BsG_Grid (struct measure_type m, struct invert_type r)
{
  int i, j;
  double bs, ba;

  if (The_Grid == NULL)
    Allocate_Grid (r.search);

  GG_a = 0.0;
  GG_b = 0.0;
  GG_g = 0.0;
  GG_bs = 0.0;
  GG_ba = 0.0;




  Set_Calc_State (m, r);
  bs = 1.0 / 32.0;
  ba = RR.default_ba;
  GG_ba = ba;
  for (i = 0; i < GRID_SIZE; i++)
    {
      bs *= 2;
      RR.slab.b = ba + bs;
      if (RR.slab.b > 0)
	RR.slab.a = bs / RR.slab.b;
      else
	RR.slab.a = 0;
      for (j = 0; j < GRID_SIZE; j++)
	{
	  RR.slab.g = 0.9999 * (2.0 * j / (GRID_SIZE - 1.0) - 1.0);
	  fill_grid_entry (i, j);
	}
    }

  The_Grid_Initialized = TRUE;
  The_Grid_Search = FIND_BsG;
}




void
Grid_ABG (int i, int j, guess_type *guess)
{
  if (0 <= i && i < GRID_SIZE && 0 <= j && j < GRID_SIZE)
    {
      guess->a = The_Grid[GRID_SIZE * i + j][A_COLUMN];
      guess->b = The_Grid[GRID_SIZE * i + j][B_COLUMN];
      guess->g = The_Grid[GRID_SIZE * i + j][G_COLUMN];
      guess->distance = Calculate_Grid_Distance (i, j);
    }
  else
    {
      guess->a = 0.5;
      guess->b = 0.5;
      guess->g = 0.5;
      guess->distance = 999;
    }
}




double
Gain (int sphere, struct measure_type m, double URU)
{
  double G, tmp;

  if (sphere == REFLECTION_SPHERE)
    {
      tmp =
	m.rw_r * (m.aw_r + (1 - m.ae_r) * (m.ad_r * m.rd_r + m.as_r * URU));
      if (tmp == 1.0)
	G = 1;
      else
	G = 1.0 + tmp / (1.0 - tmp);

    }
  else
    {
      tmp =
	m.rw_t * (m.aw_t + (1 - m.ae_t) * (m.ad_t * m.rd_t + m.as_t * URU));
      if (tmp == 1.0)
	G = 1;
      else
	G = 1.0 + tmp / (1.0 - tmp);
    }

  return G;
}




double
Gain_11 (struct measure_type m, double URU, double tdiffuse)
{
  double G, GP, G11;

  G = Gain (REFLECTION_SPHERE, m, URU);
  GP = Gain (TRANSMISSION_SPHERE, m, URU);

  G11 =
    G / (1 -
	 m.as_r * m.as_t * m.aw_r * m.aw_t * (1 - m.ae_r) * (1 -
							     m.ae_t) * G *
	 GP * tdiffuse * tdiffuse);

  return G11;
}




double
Gain_22 (struct measure_type m, double URU, double tdiffuse)
{
  double G, GP, G22;

  G = Gain (REFLECTION_SPHERE, m, URU);
  GP = Gain (TRANSMISSION_SPHERE, m, URU);

  G22 =
    GP / (1 -
	  m.as_r * m.as_t * m.aw_r * m.aw_t * (1 - m.ae_r) * (1 -
							      m.ae_t) * G *
	  GP * tdiffuse * tdiffuse);

  return G22;
}




double
Two_Sphere_R (struct measure_type m,
	      double UR1, double URU, double UT1, double UTU)
{
  double x, GP;
  GP = Gain (TRANSMISSION_SPHERE, m, URU);

  x = m.ad_r * (1 - m.ae_r) * m.rw_r * Gain_11 (m, URU, UTU);
  x *=
    (1 - m.f_r) * UR1 + m.rw_r * m.f_r + (1 - m.f_r) * m.as_t * (1 -
								 m.ae_t) *
    m.rw_t * UT1 * UTU * GP;
  return x;
}




double
Two_Sphere_T (struct measure_type m,
	      double UR1, double URU, double UT1, double UTU)
{
  double x, G;
  G = Gain (REFLECTION_SPHERE, m, URU);
  x = m.ad_t * (1 - m.ae_t) * m.rw_t * Gain_22 (m, URU, UTU);
  x *=
    (1 - m.f_r) * UT1 + (1 -
			 m.ae_r) * m.rw_r * m.as_r * UTU * (m.f_r * m.rw_r +
							    (1 -
							     m.f_r) * UR1) *
    G;
  return x;
}





void
Calculate_Distance_With_Corrections (double UR1, double UT1,
				     double Rc, double Tc,
				     double URU, double UTU,
				     double *M_R, double *M_T, double *dev)
{
  double R_direct, T_direct, R_diffuse, T_diffuse;

  R_diffuse = URU - MM.uru_lost;
  T_diffuse = UTU - MM.utu_lost;

  R_direct = UR1 - MM.ur1_lost - (1.0 - MM.fraction_of_rc_in_mr) * Rc;
  T_direct = UT1 - MM.ut1_lost - (1.0 - MM.fraction_of_tc_in_mt) * Tc;

  switch (MM.num_spheres)
    {
    case 0:

      *M_R = R_direct;
      *M_T = T_direct;


      break;

    case 1:
    case -2:
      if (MM.method == COMPARISON)

	{
	  *M_R =
	    (1 - MM.f_r) * R_direct / ((1 - MM.f_r) +
				       MM.f_r * MM.rw_r / MM.rstd_r);
	  *M_T = T_direct;
	}


      else
	{
	  double P_std, P_d, P_0;
	  double G, G_0, G_std, GP_std, GP;

	  G_0 = Gain (REFLECTION_SPHERE, MM, 0.0);
	  G = Gain (REFLECTION_SPHERE, MM, R_diffuse);
	  G_std = Gain (REFLECTION_SPHERE, MM, MM.rstd_r);

	  P_d = G * (R_direct * (1 - MM.f_r) + MM.f_r * MM.rw_r);
	  P_std = G_std * (MM.rstd_r * (1 - MM.f_r) + MM.f_r * MM.rw_r);
	  P_0 = G_0 * (MM.f_r * MM.rw_r);
	  *M_R = MM.rstd_r * (P_d - P_0) / (P_std - P_0);

	  GP = Gain (TRANSMISSION_SPHERE, MM, R_diffuse);
	  GP_std = Gain (TRANSMISSION_SPHERE, MM, 0.0);
	  *M_T = T_direct * GP / GP_std;
	}


      break;

    case 2:

      {
	double R_0, T_0;
	R_0 = Two_Sphere_R (MM, 0, 0, 0, 0);
	T_0 = Two_Sphere_T (MM, 0, 0, 0, 0);

	*M_R =
	  MM.rstd_r *
	  (Two_Sphere_R (MM, R_direct, R_diffuse, T_direct, T_diffuse) -
	   R_0) / (Two_Sphere_R (MM, MM.rstd_r, MM.rstd_r, 0, 0) - R_0);
	*M_T =
	  (Two_Sphere_T (MM, R_direct, R_diffuse, T_direct, T_diffuse) -
	   T_0) / (Two_Sphere_T (MM, 0, 0, 1, 1) - T_0);
      }


      break;
    }



  if (RR.search == FIND_A || RR.search == FIND_G || RR.search == FIND_B ||
      RR.search == FIND_Bs || RR.search == FIND_Ba)
    {


      if (MM.m_t > 0)
	{
	  if (RR.metric == RELATIVE)
	    *dev = fabs (MM.m_t - *M_T) / (MM.m_t + ABIT);
	  else
	    *dev = fabs (MM.m_t - *M_T);
	}
      else
	{
	  if (RR.metric == RELATIVE)
	    *dev = fabs (MM.m_r - *M_R) / (MM.m_r + ABIT);
	  else
	    *dev = fabs (MM.m_r - *M_R);
	}


    }
  else
    {


      if (RR.metric == RELATIVE)
	{
	  *dev = 0;
	  if (MM.m_t > ABIT)
	    *dev = T_TRUST_FACTOR * fabs (MM.m_t - *M_T) / (MM.m_t + ABIT);
	  if (RR.default_a != 0)
	    *dev += fabs (MM.m_r - *M_R) / (MM.m_r + ABIT);
	}
      else
	{
	  *dev = T_TRUST_FACTOR * fabs (MM.m_t - *M_T);
	  if (RR.default_a != 0)
	    *dev += fabs (MM.m_r - *M_R);
	}



    }



  if ((Debug (DEBUG_ITERATIONS) && !CALCULATING_GRID) ||
      (Debug (DEBUG_GRID_CALC) && CALCULATING_GRID))
    {
      static int once = 0;

      if (once == 0)
	{
	  fprintf (stderr, "%10s %10s %10s |%10s %10s |%10s %10s |%10s\n",
		   "a", "b", "g", "m_r", "fit", "m_t", "fit", "delta");
	  once = 1;
	}

      fprintf (stderr, "%10.5f %10.5f %10.5f |", RR.slab.a, RR.slab.b,
	       RR.slab.g);
      fprintf (stderr, "%10.5f %10.5f |", MM.m_r, *M_R);
      fprintf (stderr, "%10.5f %10.5f |", MM.m_t, *M_T);
      fprintf (stderr, "%10.5f \n", *dev);
    }


}




double
Calculate_Grid_Distance (int i, int j)
{
  double ur1, ut1, uru, utu, Rc, Tc, b, dev, LR, LT;

  if (Debug (DEBUG_GRID_CALC))
    fprintf (stderr, "g %2d %2d ", i, j);

  b = The_Grid[GRID_SIZE * i + j][B_COLUMN];
  ur1 = The_Grid[GRID_SIZE * i + j][UR1_COLUMN];
  ut1 = The_Grid[GRID_SIZE * i + j][UT1_COLUMN];
  uru = The_Grid[GRID_SIZE * i + j][URU_COLUMN];
  utu = The_Grid[GRID_SIZE * i + j][UTU_COLUMN];
  RR.slab.a = The_Grid[GRID_SIZE * i + j][A_COLUMN];
  RR.slab.b = The_Grid[GRID_SIZE * i + j][B_COLUMN];
  RR.slab.g = The_Grid[GRID_SIZE * i + j][G_COLUMN];

  Sp_mu_RT_Flip (MM.flip_sample,
		 RR.slab.n_top_slide, RR.slab.n_slab, RR.slab.n_bottom_slide,
		 RR.slab.b_top_slide, b, RR.slab.b_bottom_slide,
		 RR.slab.cos_angle, &Rc, &Tc);

  CALCULATING_GRID = 1;
  Calculate_Distance_With_Corrections (ur1, ut1, Rc, Tc, uru, utu, &LR, &LT,
				       &dev);
  CALCULATING_GRID = 0;

  return dev;
}




void
Calculate_Distance (double *M_R, double *M_T, double *deviation)
{
  double Rc, Tc, ur1, ut1, uru, utu;

  if (RR.slab.b <= 1e-6)
    RR.slab.b = 1e-6;
  if (Debug (DEBUG_EVERY_CALC))
    fprintf (stderr, "a=%8.5f b=%10.5f g=%8.5f ", RR.slab.a, RR.slab.b,
	     RR.slab.g);

  RT_Flip (MM.flip_sample, RR.method.quad_pts, &RR.slab, &ur1, &ut1, &uru,
	   &utu);

  if (Debug (DEBUG_EVERY_CALC))
    fprintf (stderr, "ur1=%8.5f ut1=%8.5f (not M_R and M_T!)\n", ur1, ut1);

  Sp_mu_RT_Flip (MM.flip_sample,
		 RR.slab.n_top_slide, RR.slab.n_slab, RR.slab.n_bottom_slide,
		 RR.slab.b_top_slide, RR.slab.b, RR.slab.b_bottom_slide,
		 RR.slab.cos_angle, &Rc, &Tc);

  if ((!CALCULATING_GRID && Debug (DEBUG_ITERATIONS)) ||
      (CALCULATING_GRID && Debug (DEBUG_GRID_CALC)))
    fprintf (stderr, "        ");

  Calculate_Distance_With_Corrections (ur1, ut1, Rc, Tc, uru, utu, M_R, M_T,
				       deviation);
}




void
abg_distance (double a, double b, double g, guess_type *guess)
{
  double m_r, m_t, distance;
  struct measure_type old_mm;
  struct invert_type old_rr;

  Get_Calc_State (&old_mm, &old_rr);

  RR.slab.a = a;
  RR.slab.b = b;
  RR.slab.g = g;

  Calculate_Distance (&m_r, &m_t, &distance);

  Set_Calc_State (old_mm, old_rr);

  guess->a = a;
  guess->b = b;
  guess->g = g;
  guess->distance = distance;
}





double
Find_AG_fn (double x[])
{
  double m_r, m_t, deviation;
  RR.slab.a = acalc2a (x[1]);
  RR.slab.g = gcalc2g (x[2]);
  Calculate_Distance (&m_r, &m_t, &deviation);
  return deviation;
}




double
Find_AB_fn (double x[])
{
  double m_r, m_t, deviation;
  RR.slab.a = acalc2a (x[1]);
  RR.slab.b = bcalc2b (x[2]);
  Calculate_Distance (&m_r, &m_t, &deviation);
  return deviation;
}




double
Find_Ba_fn (double x)
{
  double m_r, m_t, deviation, ba, bs;

  bs = RR.slab.b;
  ba = bcalc2b (x);
  RR.slab.b = ba + bs;
  RR.slab.a = bs / (ba + bs);

  Calculate_Distance (&m_r, &m_t, &deviation);

  RR.slab.b = bs;
  return deviation;
}




double
Find_Bs_fn (double x)
{
  double m_r, m_t, deviation, ba, bs;

  ba = RR.slab.b;
  bs = bcalc2b (x);
  RR.slab.b = ba + bs;
  RR.slab.a = bs / (ba + bs);

  Calculate_Distance (&m_r, &m_t, &deviation);

  RR.slab.b = ba;
  return deviation;
}




double
Find_A_fn (double x)
{
  double m_r, m_t, deviation;
  RR.slab.a = acalc2a (x);
  Calculate_Distance (&m_r, &m_t, &deviation);
  return deviation;
}




double
Find_B_fn (double x)
{
  double m_r, m_t, deviation;
  RR.slab.b = bcalc2b (x);
  Calculate_Distance (&m_r, &m_t, &deviation);
  return deviation;
}




double
Find_G_fn (double x)
{
  double m_r, m_t, deviation;
  RR.slab.g = gcalc2g (x);
  Calculate_Distance (&m_r, &m_t, &deviation);
  return deviation;
}




double
Find_BG_fn (double x[])
{
  double m_r, m_t, deviation;
  RR.slab.b = bcalc2b (x[1]);
  RR.slab.g = gcalc2g (x[2]);
  RR.slab.a = RR.default_a;
  Calculate_Distance (&m_r, &m_t, &deviation);
  return deviation;
}




double
Find_BaG_fn (double x[])
{
  double m_r, m_t, deviation;

  RR.slab.b = bcalc2b (x[1]) + RR.default_bs;
  if (RR.slab.b <= 0)
    RR.slab.a = 0;
  else
    RR.slab.a = RR.default_bs / RR.slab.b;

  RR.slab.g = gcalc2g (x[2]);

  Calculate_Distance (&m_r, &m_t, &deviation);
  return deviation;
}




double
Find_BsG_fn (double x[])
{
  double m_r, m_t, deviation;

  RR.slab.b = bcalc2b (x[1]) + RR.default_ba;
  if (RR.slab.b <= 0)
    RR.slab.a = 0;
  else
    RR.slab.a = 1.0 - RR.default_ba / RR.slab.b;

  RR.slab.g = gcalc2g (x[2]);
  Calculate_Distance (&m_r, &m_t, &deviation);
  return deviation;
}




double
maxloss (double f)
{
  struct measure_type m_old;
  struct invert_type r_old;
  double m_r, m_t, deviation;

  Get_Calc_State (&m_old, &r_old);

  RR.slab.a = 1.0;
  MM.ur1_lost *= f;
  MM.ut1_lost *= f;

  Calculate_Distance (&m_r, &m_t, &deviation);

  Set_Calc_State (m_old, r_old);
  deviation = ((MM.m_r + MM.m_t) - (m_r + m_t));

  return deviation;
}




void
Max_Light_Loss (struct measure_type m, struct invert_type r,
		double *ur1_loss, double *ut1_loss)
{
  struct measure_type m_old;
  struct invert_type r_old;

  *ur1_loss = m.ur1_lost;
  *ut1_loss = m.ut1_lost;

  if (Debug (DEBUG_LOST_LIGHT))
    fprintf (stderr, "\nlost before ur1=%7.5f, ut1=%7.5f\n", *ur1_loss,
	     *ut1_loss);

  Get_Calc_State (&m_old, &r_old);

  Set_Calc_State (m, r);

  if (maxloss (1.0) * maxloss (0.0) < 0)
    {
      double frac;
      frac = zbrent (maxloss, 0.00, 1.0, 0.001);

      *ur1_loss = m.ur1_lost * frac;
      *ut1_loss = m.ut1_lost * frac;
    }

  Set_Calc_State (m_old, r_old);
  if (Debug (DEBUG_LOST_LIGHT))
    fprintf (stderr, "lost after  ur1=%7.5f, ut1=%7.5f\n", *ur1_loss,
	     *ut1_loss);
}
