
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ad_globl.h"
#include "nr_util.h"
#include "nr_mnbrk.h"
#include "nr_brent.h"
#include "nr_amoeb.h"
#include "iad_type.h"
#include "iad_util.h"
#include "iad_calc.h"
#include "iad_find.h"
#define NUMBER_OF_GUESSES 10

guess_type guess[NUMBER_OF_GUESSES];

int
compare_guesses (const void *p1, const void *p2)
{
  guess_type *g1 = (guess_type *) p1;
  guess_type *g2 = (guess_type *) p2;

  if (g1->distance < g2->distance)
    return -1;
  else if (g1->distance == g2->distance)
    return 0;
  else
    return 1;
}



void
U_Find_Ba (struct measure_type m, struct invert_type *r)
{
  double ax, bx, cx, fa, fb, fc, ba;

  if (Debug (DEBUG_SEARCH))
    {
      fprintf (stderr, "In U_Find_Bs");
      fprintf (stderr, " (mu=%6.4f)", r->slab.cos_angle);
      if (r->default_bs != UNINITIALIZED)
	fprintf (stderr, "  default_bs = %8.5f", r->default_bs);
      if (r->default_g != UNINITIALIZED)
	fprintf (stderr, "  default_g = %8.5f", r->default_g);
      fprintf (stderr, "\n");
    }

  r->slab.a = 0;
  r->slab.g = (r->default_g == UNINITIALIZED) ? 0 : r->default_g;
  r->slab.b = (r->default_bs == UNINITIALIZED) ? HUGE_VAL : r->default_bs;

  Set_Calc_State (m, *r);

  ax = b2bcalc (0.1);
  bx = b2bcalc (1.0);
  mnbrak (&ax, &bx, &cx, &fa, &fb, &fc, Find_Ba_fn);
  r->final_distance = brent (ax, bx, cx, Find_Ba_fn, r->tolerance, &ba);


  r->slab.a = (r->slab.b) / (bcalc2b (ba) + r->slab.b);
  r->slab.b = bcalc2b (ba) + r->slab.b;
  Set_Calc_State (m, *r);


  r->a = r->slab.a;
  r->b = r->slab.b;
  r->g = r->slab.g;
  r->found = (r->tolerance <= r->final_distance);


}




void
U_Find_Bs (struct measure_type m, struct invert_type *r)
{
  double ax, bx, cx, fa, fb, fc, bs;

  if (Debug (DEBUG_SEARCH))
    {
      fprintf (stderr, "In U_Find_Bs");
      fprintf (stderr, " (mu=%6.4f)", r->slab.cos_angle);
      if (r->default_ba != UNINITIALIZED)
	fprintf (stderr, "  default_ba = %8.5f", r->default_ba);
      if (r->default_g != UNINITIALIZED)
	fprintf (stderr, "  default_g = %8.5f", r->default_g);
      fprintf (stderr, "\n");
    }

  r->slab.a = 0;
  r->slab.g = (r->default_g == UNINITIALIZED) ? 0 : r->default_g;
  r->slab.b = (r->default_ba == UNINITIALIZED) ? HUGE_VAL : r->default_ba;

  Set_Calc_State (m, *r);

  ax = b2bcalc (0.1);
  bx = b2bcalc (1.0);
  mnbrak (&ax, &bx, &cx, &fa, &fb, &fc, Find_Bs_fn);
  r->final_distance = brent (ax, bx, cx, Find_Bs_fn, r->tolerance, &bs);


  r->slab.a = bcalc2b (bs) / (bcalc2b (bs) + r->slab.b);
  r->slab.b = bcalc2b (bs) + r->slab.b;
  Set_Calc_State (m, *r);


  r->a = r->slab.a;
  r->b = r->slab.b;
  r->g = r->slab.g;
  r->found = (r->tolerance <= r->final_distance);


}




void
U_Find_A (struct measure_type m, struct invert_type *r)
{
  double Rt, Tt, Rd, Rc, Td, Tc;

  if (Debug (DEBUG_SEARCH))
    {
      fprintf (stderr, "In U_Find_A");
      fprintf (stderr, " (mu=%6.4f)", r->slab.cos_angle);
      if (r->default_b != UNINITIALIZED)
	fprintf (stderr, "  default_b = %8.5f", r->default_b);
      if (r->default_g != UNINITIALIZED)
	fprintf (stderr, "  default_g = %8.5f", r->default_g);
      fprintf (stderr, "\n");
    }

  Estimate_RT (m, *r, &Rt, &Tt, &Rd, &Rc, &Td, &Tc);

  r->slab.g = (r->default_g == UNINITIALIZED) ? 0 : r->default_g;
  r->slab.b = (r->default_b == UNINITIALIZED) ? HUGE_VAL : r->default_b;
  r->slab.a = 0.0;
  r->final_distance = 0.0;
  Set_Calc_State (m, *r);

  if (Rt > 0.99999)
    r->final_distance = Find_A_fn (a2acalc (1.0));
  else
    {
      double x, ax, bx, cx, fa, fb, fc;

      ax = a2acalc (0.3);
      bx = a2acalc (0.5);

      mnbrak (&ax, &bx, &cx, &fa, &fb, &fc, Find_A_fn);
      r->final_distance = brent (ax, bx, cx, Find_A_fn, r->tolerance, &x);
      r->slab.a = acalc2a (x);
    }


  r->a = r->slab.a;
  r->b = r->slab.b;
  r->g = r->slab.g;
  r->found = (r->tolerance <= r->final_distance);


}




void
U_Find_B (struct measure_type m, struct invert_type *r)
{
  double Rt, Tt, Rd, Rc, Td, Tc;

  if (Debug (DEBUG_SEARCH))
    {
      fprintf (stderr, "In U_Find_B");
      fprintf (stderr, " (mu=%6.4f)", r->slab.cos_angle);
      if (r->default_a != UNINITIALIZED)
	fprintf (stderr, "  default_a = %8.5f", r->default_a);
      if (r->default_g != UNINITIALIZED)
	fprintf (stderr, "  default_g = %8.5f", r->default_g);
      fprintf (stderr, "\n");
    }

  Estimate_RT (m, *r, &Rt, &Tt, &Rd, &Rc, &Td, &Tc);

  r->slab.g = (r->default_g == UNINITIALIZED) ? 0 : r->default_g;
  r->slab.a = (r->default_a == UNINITIALIZED) ? 0 : r->default_a;
  r->slab.b = 0.5;
  r->final_distance = 0.0;
  Set_Calc_State (m, *r);


  {
    double x, ax, bx, cx, fa, fb, fc;

    ax = b2bcalc (0.1);
    bx = b2bcalc (10);

    mnbrak (&ax, &bx, &cx, &fa, &fb, &fc, Find_B_fn);
    r->final_distance = brent (ax, bx, cx, Find_B_fn, r->tolerance, &x);
    r->slab.b = bcalc2b (x);
    Set_Calc_State (m, *r);
  }





  r->a = r->slab.a;
  r->b = r->slab.b;
  r->g = r->slab.g;
  r->found = (r->tolerance <= r->final_distance);



  if (Debug (DEBUG_SEARCH))
    {
      fprintf (stderr, "In U_Find_B final (a,b,g) = ");
      fprintf (stderr, "(%8.5f,%8.5f,%8.5f)\n", r->a, r->b, r->g);
    }
}




void
U_Find_G (struct measure_type m, struct invert_type *r)
{
  double Rt, Tt, Rd, Rc, Td, Tc;

  if (Debug (DEBUG_SEARCH))
    {
      fprintf (stderr, "In U_Find_G");
      fprintf (stderr, " (mu=%6.4f)", r->slab.cos_angle);
      if (r->default_a != UNINITIALIZED)
	fprintf (stderr, "  default_a = %8.5f", r->default_a);
      if (r->default_b != UNINITIALIZED)
	fprintf (stderr, "  default_b = %8.5f", r->default_b);
      fprintf (stderr, "\n");
    }

  Estimate_RT (m, *r, &Rt, &Tt, &Rd, &Rc, &Td, &Tc);

  r->slab.a = (r->default_a == UNINITIALIZED) ? 0.5 : r->default_a;
  r->slab.b = (r->default_b == UNINITIALIZED) ? HUGE_VAL : r->default_b;
  r->slab.g = 0.0;
  r->final_distance = 0.0;
  Set_Calc_State (m, *r);

  if (Rd > 0.0)
    {
      double x, ax, bx, cx, fa, fb, fc;

      ax = g2gcalc (-0.99);
      bx = g2gcalc (0.99);

      mnbrak (&ax, &bx, &cx, &fa, &fb, &fc, Find_G_fn);
      r->final_distance = brent (ax, bx, cx, Find_G_fn, r->tolerance, &x);
      r->slab.g = gcalc2g (x);
      Set_Calc_State (m, *r);
    }


  r->a = r->slab.a;
  r->b = r->slab.b;
  r->g = r->slab.g;
  r->found = (r->tolerance <= r->final_distance);


}





void
U_Find_AG (struct measure_type m, struct invert_type *r)
{

  int i, i_best, j_best;
  double *x, *y, **p;

  x = dvector (1, 2);
  y = dvector (1, 3);
  p = dmatrix (1, 3, 1, 2);



  if (Debug (DEBUG_SEARCH))
    {
      fprintf (stderr, "In U_Find_AG");
      fprintf (stderr, " (mu=%6.4f)", r->slab.cos_angle);
      if (r->default_b != UNINITIALIZED)
	fprintf (stderr, "  default_b = %8.5f", r->default_b);
      fprintf (stderr, "\n");
    }

  if (m.num_measures == 3)
    r->slab.b = What_Is_B (r->slab, m.m_u);
  else if (r->default_b == UNINITIALIZED)
    r->slab.b = 1;
  else
    r->slab.b = r->default_b;

  Set_Calc_State (m, *r);

  {

    size_t count = NUMBER_OF_GUESSES;

    abg_distance (r->slab.a, r->slab.b, r->slab.g, &(guess[0]));

    if (!Valid_Grid (m, r->search))
      Fill_Grid (m, *r);


    Near_Grid_Points (m.m_r, m.m_t, r->search, &i_best, &j_best);
    Grid_ABG (i_best, j_best, &(guess[1]));
    Grid_ABG (i_best + 1, j_best, &(guess[2]));
    Grid_ABG (i_best - 1, j_best, &(guess[3]));
    Grid_ABG (i_best, j_best + 1, &(guess[4]));
    Grid_ABG (i_best, j_best - 1, &(guess[5]));
    Grid_ABG (i_best + 1, j_best + 1, &(guess[6]));
    Grid_ABG (i_best - 1, j_best - 1, &(guess[7]));
    Grid_ABG (i_best + 1, j_best - 1, &(guess[8]));
    Grid_ABG (i_best - 1, j_best + 1, &(guess[9]));

    qsort ((void *) guess, count, sizeof (guess_type), compare_guesses);

    if (Debug (DEBUG_BEST_GUESS))
      {
	int k;
	fprintf (stderr, "after\n");
	for (k = 0; k <= 6; k++)
	  {
	    fprintf (stderr, "%3d  ", k);
	    fprintf (stderr, "%10.5f ", guess[k].a);
	    fprintf (stderr, "%10.5f ", guess[k].b);
	    fprintf (stderr, "%10.5f ", guess[k].g);
	    fprintf (stderr, "%10.5f\n", guess[k].distance);
	  }
      }
  }



  {
    int k, kk;

    p[1][1] = a2acalc (guess[0].a);
    p[1][2] = g2gcalc (guess[0].g);

    for (k = 1; k < 7; k++)
      {
	if (guess[0].a != guess[k].a)
	  break;
      }

    p[2][1] = a2acalc (guess[k].a);
    p[2][2] = g2gcalc (guess[k].g);

    for (kk = 1; kk < 7; kk++)
      {
	if (guess[0].g != guess[kk].g && guess[k].g != guess[kk].g)
	  break;
      }
    p[3][1] = a2acalc (guess[kk].a);
    p[3][2] = g2gcalc (guess[kk].g);

    if (Debug (DEBUG_BEST_GUESS))
      {
	fprintf (stderr, "guess 1");
	fprintf (stderr, "%10.5f ", guess[0].a);
	fprintf (stderr, "%10.5f ", guess[0].b);
	fprintf (stderr, "%10.5f ", guess[0].g);
	fprintf (stderr, "%10.5f\n", guess[0].distance);
	fprintf (stderr, "guess 2");
	fprintf (stderr, "%10.5f ", guess[k].a);
	fprintf (stderr, "%10.5f ", guess[k].b);
	fprintf (stderr, "%10.5f ", guess[k].g);
	fprintf (stderr, "%10.5f\n", guess[k].distance);
	fprintf (stderr, "guess 3");
	fprintf (stderr, "%10.5f ", guess[kk].a);
	fprintf (stderr, "%10.5f ", guess[kk].b);
	fprintf (stderr, "%10.5f ", guess[kk].g);
	fprintf (stderr, "%10.5f\n", guess[kk].distance);
      }
  }




  for (i = 1; i <= 3; i++)
    {
      x[1] = p[i][1];
      x[2] = p[i][2];
      y[i] = Find_AG_fn (x);
    }



  amoeba (p, y, 2, r->tolerance, Find_AG_fn, &r->iterations);

  r->final_distance = 10;
  for (i = 1; i <= 3; i++)
    {
      if (y[i] < r->final_distance)
	{
	  r->slab.a = acalc2a (p[i][1]);
	  r->slab.g = gcalc2g (p[i][2]);
	  r->final_distance = y[i];
	}
    }




  free_dvector (x, 1, 2);
  free_dvector (y, 1, 3);
  free_dmatrix (p, 1, 3, 1, 2);




  r->a = r->slab.a;
  r->b = r->slab.b;
  r->g = r->slab.g;
  r->found = (r->tolerance <= r->final_distance);


}




void
U_Find_AB (struct measure_type m, struct invert_type *r)
{

  int i, i_best, j_best;
  double *x, *y, **p;

  x = dvector (1, 2);
  y = dvector (1, 3);
  p = dmatrix (1, 3, 1, 2);



  if (Debug (DEBUG_SEARCH))
    {
      fprintf (stderr, "In U_Find_AB");
      fprintf (stderr, " (mu=%6.4f)", r->slab.cos_angle);
      if (r->default_g != UNINITIALIZED)
	fprintf (stderr, "  default_g = %8.5f", r->default_g);
      fprintf (stderr, "\n");
    }

  r->slab.g = (r->default_g == UNINITIALIZED) ? 0 : r->default_g;
  Set_Calc_State (m, *r);


  {

    size_t count = NUMBER_OF_GUESSES;

    abg_distance (r->slab.a, r->slab.b, r->slab.g, &(guess[0]));

    if (!Valid_Grid (m, r->search))
      Fill_Grid (m, *r);


    Near_Grid_Points (m.m_r, m.m_t, r->search, &i_best, &j_best);
    Grid_ABG (i_best, j_best, &(guess[1]));
    Grid_ABG (i_best + 1, j_best, &(guess[2]));
    Grid_ABG (i_best - 1, j_best, &(guess[3]));
    Grid_ABG (i_best, j_best + 1, &(guess[4]));
    Grid_ABG (i_best, j_best - 1, &(guess[5]));
    Grid_ABG (i_best + 1, j_best + 1, &(guess[6]));
    Grid_ABG (i_best - 1, j_best - 1, &(guess[7]));
    Grid_ABG (i_best + 1, j_best - 1, &(guess[8]));
    Grid_ABG (i_best - 1, j_best + 1, &(guess[9]));

    qsort ((void *) guess, count, sizeof (guess_type), compare_guesses);

    if (Debug (DEBUG_BEST_GUESS))
      {
	int k;
	fprintf (stderr, "after\n");
	for (k = 0; k <= 6; k++)
	  {
	    fprintf (stderr, "%3d  ", k);
	    fprintf (stderr, "%10.5f ", guess[k].a);
	    fprintf (stderr, "%10.5f ", guess[k].b);
	    fprintf (stderr, "%10.5f ", guess[k].g);
	    fprintf (stderr, "%10.5f\n", guess[k].distance);
	  }
      }
  }



  {
    int k, kk;

    p[1][1] = a2acalc (guess[0].a);
    p[1][2] = b2bcalc (guess[0].b);

    for (k = 1; k < 7; k++)
      {
	if (guess[0].a != guess[k].a)
	  break;
      }

    p[2][1] = a2acalc (guess[k].a);
    p[2][2] = b2bcalc (guess[k].b);

    for (kk = 1; kk < 7; kk++)
      {
	if (guess[0].b != guess[kk].b && guess[k].b != guess[kk].b)
	  break;
      }
    p[3][1] = a2acalc (guess[kk].a);
    p[3][2] = b2bcalc (guess[kk].b);

    if (Debug (DEBUG_BEST_GUESS))
      {
	fprintf (stderr, "guess 1");
	fprintf (stderr, "%10.5f ", guess[0].a);
	fprintf (stderr, "%10.5f ", guess[0].b);
	fprintf (stderr, "%10.5f ", guess[0].g);
	fprintf (stderr, "%10.5f\n", guess[0].distance);
	fprintf (stderr, "guess 2");
	fprintf (stderr, "%10.5f ", guess[k].a);
	fprintf (stderr, "%10.5f ", guess[k].b);
	fprintf (stderr, "%10.5f ", guess[k].g);
	fprintf (stderr, "%10.5f\n", guess[k].distance);
	fprintf (stderr, "guess 3");
	fprintf (stderr, "%10.5f ", guess[kk].a);
	fprintf (stderr, "%10.5f ", guess[kk].b);
	fprintf (stderr, "%10.5f ", guess[kk].g);
	fprintf (stderr, "%10.5f\n", guess[kk].distance);
      }
  }




  for (i = 1; i <= 3; i++)
    {
      x[1] = p[i][1];
      x[2] = p[i][2];
      y[i] = Find_AB_fn (x);
    }


  amoeba (p, y, 2, r->tolerance, Find_AB_fn, &r->iterations);

  r->final_distance = 10;
  for (i = 1; i <= 3; i++)
    {
      if (y[i] < r->final_distance)
	{
	  r->slab.a = acalc2a (p[i][1]);
	  r->slab.b = bcalc2b (p[i][2]);
	  r->final_distance = y[i];
	}
    }




  free_dvector (x, 1, 2);
  free_dvector (y, 1, 3);
  free_dmatrix (p, 1, 3, 1, 2);



  r->a = r->slab.a;
  r->b = r->slab.b;
  r->g = r->slab.g;
  r->found = (r->tolerance <= r->final_distance);


}




void
U_Find_BG (struct measure_type m, struct invert_type *r)
{

  int i, i_best, j_best;
  double *x, *y, **p;

  x = dvector (1, 2);
  y = dvector (1, 3);
  p = dmatrix (1, 3, 1, 2);



  if (Debug (DEBUG_SEARCH))
    {
      fprintf (stderr, "In U_Find_BG");
      fprintf (stderr, " (mu=%6.4f)", r->slab.cos_angle);
      if (r->default_a != UNINITIALIZED)
	fprintf (stderr, "  default_a = %8.5f", r->default_a);
      fprintf (stderr, "\n");
    }

  r->slab.a = (r->default_a == UNINITIALIZED) ? 0 : r->default_a;
  Set_Calc_State (m, *r);


  {

    size_t count = NUMBER_OF_GUESSES;

    abg_distance (r->slab.a, r->slab.b, r->slab.g, &(guess[0]));

    if (!Valid_Grid (m, r->search))
      Fill_Grid (m, *r);


    Near_Grid_Points (m.m_r, m.m_t, r->search, &i_best, &j_best);
    Grid_ABG (i_best, j_best, &(guess[1]));
    Grid_ABG (i_best + 1, j_best, &(guess[2]));
    Grid_ABG (i_best - 1, j_best, &(guess[3]));
    Grid_ABG (i_best, j_best + 1, &(guess[4]));
    Grid_ABG (i_best, j_best - 1, &(guess[5]));
    Grid_ABG (i_best + 1, j_best + 1, &(guess[6]));
    Grid_ABG (i_best - 1, j_best - 1, &(guess[7]));
    Grid_ABG (i_best + 1, j_best - 1, &(guess[8]));
    Grid_ABG (i_best - 1, j_best + 1, &(guess[9]));

    qsort ((void *) guess, count, sizeof (guess_type), compare_guesses);

    if (Debug (DEBUG_BEST_GUESS))
      {
	int k;
	fprintf (stderr, "after\n");
	for (k = 0; k <= 6; k++)
	  {
	    fprintf (stderr, "%3d  ", k);
	    fprintf (stderr, "%10.5f ", guess[k].a);
	    fprintf (stderr, "%10.5f ", guess[k].b);
	    fprintf (stderr, "%10.5f ", guess[k].g);
	    fprintf (stderr, "%10.5f\n", guess[k].distance);
	  }
      }
  }



  {
    int k, kk;

    p[1][1] = b2bcalc (guess[0].b);
    p[1][2] = g2gcalc (guess[0].g);

    for (k = 1; k < 7; k++)
      {
	if (guess[0].b != guess[k].b)
	  break;
      }

    p[2][1] = b2bcalc (guess[k].b);
    p[2][2] = g2gcalc (guess[k].g);

    for (kk = 1; kk < 7; kk++)
      {
	if (guess[0].g != guess[kk].g && guess[k].g != guess[kk].g)
	  break;
      }
    p[3][1] = b2bcalc (guess[kk].b);
    p[3][2] = g2gcalc (guess[kk].g);

    if (Debug (DEBUG_BEST_GUESS))
      {
	fprintf (stderr, "guess 1");
	fprintf (stderr, "%10.5f ", guess[0].a);
	fprintf (stderr, "%10.5f ", guess[0].b);
	fprintf (stderr, "%10.5f ", guess[0].g);
	fprintf (stderr, "%10.5f\n", guess[0].distance);
	fprintf (stderr, "guess 2");
	fprintf (stderr, "%10.5f ", guess[k].a);
	fprintf (stderr, "%10.5f ", guess[k].b);
	fprintf (stderr, "%10.5f ", guess[k].g);
	fprintf (stderr, "%10.5f\n", guess[k].distance);
	fprintf (stderr, "guess 3");
	fprintf (stderr, "%10.5f ", guess[kk].a);
	fprintf (stderr, "%10.5f ", guess[kk].b);
	fprintf (stderr, "%10.5f ", guess[kk].g);
	fprintf (stderr, "%10.5f\n", guess[kk].distance);
      }
  }




  for (i = 1; i <= 3; i++)
    {
      x[1] = p[i][1];
      x[2] = p[i][2];
      y[i] = Find_BG_fn (x);
    }


  amoeba (p, y, 2, r->tolerance, Find_BG_fn, &r->iterations);

  r->final_distance = 10;
  for (i = 1; i <= 3; i++)
    {
      if (y[i] < r->final_distance)
	{
	  r->slab.b = bcalc2b (p[i][1]);
	  r->slab.g = gcalc2g (p[i][2]);
	  r->final_distance = y[i];
	}
    }




  free_dvector (x, 1, 2);
  free_dvector (y, 1, 3);
  free_dmatrix (p, 1, 3, 1, 2);



  r->a = r->slab.a;
  r->b = r->slab.b;
  r->g = r->slab.g;
  r->found = (r->tolerance <= r->final_distance);


}




void
U_Find_BaG (struct measure_type m, struct invert_type *r)
{

  int i, i_best, j_best;
  double *x, *y, **p;

  x = dvector (1, 2);
  y = dvector (1, 3);
  p = dmatrix (1, 3, 1, 2);


  Set_Calc_State (m, *r);

  {

    size_t count = NUMBER_OF_GUESSES;

    abg_distance (r->slab.a, r->slab.b, r->slab.g, &(guess[0]));

    if (!Valid_Grid (m, r->search))
      Fill_Grid (m, *r);


    Near_Grid_Points (m.m_r, m.m_t, r->search, &i_best, &j_best);
    Grid_ABG (i_best, j_best, &(guess[1]));
    Grid_ABG (i_best + 1, j_best, &(guess[2]));
    Grid_ABG (i_best - 1, j_best, &(guess[3]));
    Grid_ABG (i_best, j_best + 1, &(guess[4]));
    Grid_ABG (i_best, j_best - 1, &(guess[5]));
    Grid_ABG (i_best + 1, j_best + 1, &(guess[6]));
    Grid_ABG (i_best - 1, j_best - 1, &(guess[7]));
    Grid_ABG (i_best + 1, j_best - 1, &(guess[8]));
    Grid_ABG (i_best - 1, j_best + 1, &(guess[9]));

    qsort ((void *) guess, count, sizeof (guess_type), compare_guesses);

    if (Debug (DEBUG_BEST_GUESS))
      {
	int k;
	fprintf (stderr, "after\n");
	for (k = 0; k <= 6; k++)
	  {
	    fprintf (stderr, "%3d  ", k);
	    fprintf (stderr, "%10.5f ", guess[k].a);
	    fprintf (stderr, "%10.5f ", guess[k].b);
	    fprintf (stderr, "%10.5f ", guess[k].g);
	    fprintf (stderr, "%10.5f\n", guess[k].distance);
	  }
      }
  }




  if (guess[0].b > r->default_bs)
    {
      p[1][1] = b2bcalc (guess[0].b - r->default_bs);
      p[2][1] = b2bcalc (2 * (guess[0].b - r->default_bs));
      p[3][1] = p[1][1];
    }
  else
    {
      p[1][1] = b2bcalc (0.0001);
      p[2][1] = b2bcalc (0.001);
      p[3][1] = p[1][1];
    }

  p[1][2] = g2gcalc (guess[0].g);
  p[2][2] = p[1][2];
  p[3][2] = g2gcalc (0.9 * guess[0].g + 0.05);




  for (i = 1; i <= 3; i++)
    {
      x[1] = p[i][1];
      x[2] = p[i][2];
      y[i] = Find_BaG_fn (x);
    }


  amoeba (p, y, 2, r->tolerance, Find_BaG_fn, &r->iterations);

  r->final_distance = 10;
  for (i = 1; i <= 3; i++)
    {
      if (y[i] < r->final_distance)
	{
	  r->slab.b = bcalc2b (p[i][1]) + r->default_bs;
	  r->slab.a = r->default_bs / r->slab.b;
	  r->slab.g = gcalc2g (p[i][2]);
	  r->final_distance = y[i];
	}
    }




  free_dvector (x, 1, 2);
  free_dvector (y, 1, 3);
  free_dmatrix (p, 1, 3, 1, 2);



  r->a = r->slab.a;
  r->b = r->slab.b;
  r->g = r->slab.g;
  r->found = (r->tolerance <= r->final_distance);


}




void
U_Find_BsG (struct measure_type m, struct invert_type *r)
{

  int i, i_best, j_best;
  double *x, *y, **p;

  x = dvector (1, 2);
  y = dvector (1, 3);
  p = dmatrix (1, 3, 1, 2);



  if (Debug (DEBUG_SEARCH))
    {
      fprintf (stderr, "In U_Find_BsG");
      fprintf (stderr, " (mu=%6.4f)", r->slab.cos_angle);
      if (r->default_ba != UNINITIALIZED)
	fprintf (stderr, "  default_ba = %8.5f", r->default_ba);
      fprintf (stderr, "\n");
    }

  Set_Calc_State (m, *r);

  {

    size_t count = NUMBER_OF_GUESSES;

    abg_distance (r->slab.a, r->slab.b, r->slab.g, &(guess[0]));

    if (!Valid_Grid (m, r->search))
      Fill_Grid (m, *r);


    Near_Grid_Points (m.m_r, m.m_t, r->search, &i_best, &j_best);
    Grid_ABG (i_best, j_best, &(guess[1]));
    Grid_ABG (i_best + 1, j_best, &(guess[2]));
    Grid_ABG (i_best - 1, j_best, &(guess[3]));
    Grid_ABG (i_best, j_best + 1, &(guess[4]));
    Grid_ABG (i_best, j_best - 1, &(guess[5]));
    Grid_ABG (i_best + 1, j_best + 1, &(guess[6]));
    Grid_ABG (i_best - 1, j_best - 1, &(guess[7]));
    Grid_ABG (i_best + 1, j_best - 1, &(guess[8]));
    Grid_ABG (i_best - 1, j_best + 1, &(guess[9]));

    qsort ((void *) guess, count, sizeof (guess_type), compare_guesses);

    if (Debug (DEBUG_BEST_GUESS))
      {
	int k;
	fprintf (stderr, "after\n");
	for (k = 0; k <= 6; k++)
	  {
	    fprintf (stderr, "%3d  ", k);
	    fprintf (stderr, "%10.5f ", guess[k].a);
	    fprintf (stderr, "%10.5f ", guess[k].b);
	    fprintf (stderr, "%10.5f ", guess[k].g);
	    fprintf (stderr, "%10.5f\n", guess[k].distance);
	  }
      }
  }




  p[1][1] = b2bcalc (guess[0].b - r->default_ba);
  p[1][2] = g2gcalc (guess[0].g);

  p[2][1] = b2bcalc (2 * guess[0].b - 2 * r->default_ba);
  p[2][2] = p[1][2];

  p[3][1] = p[1][1];
  p[3][2] = g2gcalc (0.9 * guess[0].g + 0.05);





  for (i = 1; i <= 3; i++)
    {
      x[1] = p[i][1];
      x[2] = p[i][2];
      y[i] = Find_BsG_fn (x);
    }


  amoeba (p, y, 2, r->tolerance, Find_BsG_fn, &r->iterations);

  r->final_distance = 10;
  for (i = 1; i <= 3; i++)
    {
      if (y[i] < r->final_distance)
	{
	  r->slab.b = bcalc2b (p[i][1]) + r->default_ba;
	  r->slab.a = 1 - r->default_ba / r->slab.b;
	  r->slab.g = gcalc2g (p[i][2]);
	  r->final_distance = y[i];
	}
    }


  free_dvector (x, 1, 2);
  free_dvector (y, 1, 3);
  free_dmatrix (p, 1, 3, 1, 2);



  r->a = r->slab.a;
  r->b = r->slab.b;
  r->g = r->slab.g;
  r->found = (r->tolerance <= r->final_distance);


}
