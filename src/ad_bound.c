

#include <math.h>
#include <stdio.h>
#include "ad_globl.h"
#include "ad_bound.h"
#include "ad_frsnl.h"
#include "ad_matrx.h"
#include "nr_util.h"


static void A_Add_Slide (int n, double **R12, double **R21, double **T12,
			 double **T21, double *R10, double *T01, double **R20,
			 double **T02, double **atemp, double **btemp);

static void B_Add_Slide (int n, double **R12, double **T21,
			 double *R01, double *R10, double *T01, double *T10,
			 double **R02, double **T20,
			 double **atemp, double **btemp);



void
Init_Boundary (struct AD_slab_type slab, int n,
	       double *R01, double *R10, double *T01, double *T10,
	       char boundary)
{
  if (boundary == TOP_BOUNDARY)
    {
      Boundary_RT (1.0, slab.n_top_slide, slab.n_slab, n, slab.b_top_slide,
		   R01, T01);
      Boundary_RT (slab.n_slab, slab.n_top_slide, 1.0, n, slab.b_top_slide,
		   R10, T10);
    }
  else
    {
      Boundary_RT (1.0, slab.n_bottom_slide, slab.n_slab, n,
		   slab.b_bottom_slide, R10, T10);
      Boundary_RT (slab.n_slab, slab.n_bottom_slide, 1.0, n,
		   slab.b_bottom_slide, R01, T01);
    }
}





void
Boundary_RT (double n_i, double n_g, double n_t, int n, double b,
	     double *R, double *T)
{
  int i;
  double refl, trans;
  double mu;

  for (i = 1; i <= n; i++)
    {
      if (n_i == 1.0)
	mu = Cos_Snell (n_t, angle[i], n_i);
      else
	mu = angle[i];

      Absorbing_Glass_RT (n_i, n_g, n_t, mu, b, &refl, &trans);
      R[i] = refl * twoaw[i];
      T[i] = trans;
    }

}




void
Add_Top (int n, double *R01, double *R10, double *T01, double *T10,
	 double **R12, double **R21, double **T12, double **T21,
	 double **R02, double **R20, double **T02, double **T20,
	 double **atemp, double **btemp)
{
  A_Add_Slide (n, R12, R21, T12, T21, R10, T01, R20, T02, atemp, btemp);
  B_Add_Slide (n, R12, T21, R01, R10, T01, T10, R02, T20, atemp, btemp);
}




void
Add_Bottom (int n, double **R01, double **R10, double **T01, double **T10,
	    double *R12, double *R21, double *T12, double *T21,
	    double **R02, double **R20, double **T02, double **T20,
	    double **atemp, double **btemp)
{
  A_Add_Slide (n, R10, R01, T10, T01, R12, T21, R02, T20, atemp, btemp);
  B_Add_Slide (n, R10, T01, R21, R12, T21, T12, R20, T02, atemp, btemp);
}




static void
A_Add_Slide (int n, double **R12, double **R21, double **T12, double **T21,
	     double *R10, double *T01, double **R20, double **T02,
	     double **atemp, double **btemp)
{
  double **ctemp;

  ctemp = R20;
  Left_Diagonal_Multiply (n, R10, R12, atemp);
  One_Minus (n, atemp);
  Left_Inverse_Multiply (n, atemp, T12, ctemp);
  Right_Diagonal_Multiply (n, ctemp, T01, T02);
  Right_Diagonal_Multiply (n, ctemp, R10, btemp);
  Matrix_Multiply (n, btemp, T21, atemp);
  Matrix_Sum (n, R21, atemp, R20);
}




static void
B_Add_Slide (int n, double **R12, double **T21,
	     double *R01, double *R10, double *T01, double *T10,
	     double **R02, double **T20, double **atemp, double **btemp)
{
  double **ctemp;
  int i;
  ctemp = R02;

  Right_Diagonal_Multiply (n, R12, R10, atemp);
  One_Minus (n, atemp);
  Diagonal_To_Matrix (n, T10, btemp);
  Left_Inverse_Multiply (n, atemp, btemp, ctemp);
  Matrix_Multiply (n, ctemp, T21, T20);
  Matrix_Multiply (n, ctemp, R12, btemp);
  Right_Diagonal_Multiply (n, btemp, T01, R02);
  for (i = 1; i <= n; i++)
    R02[i][i] += R01[i] / twoaw[i] / twoaw[i];
}




void
Add_Slides (int n, double *R01, double *R10, double *T01, double *T10,
	    double **R, double **T,
	    double **R_total, double **T_total,
	    double **atemp, double **btemp)
{
  int i;
  double **R12, **R21, **T12, **T21;
  double temp;

  R12 = R;
  R21 = R;
  T21 = T;
  T12 = T;
  Left_Diagonal_Multiply (n, R10, R12, atemp);
  One_Minus (n, atemp);
  Left_Inverse_Multiply (n, atemp, T12, T_total);
  Right_Diagonal_Multiply (n, T_total, R10, btemp);
  Matrix_Multiply (n, btemp, T21, R_total);
  Matrix_Sum (n, R_total, R21, R_total);

  Right_Diagonal_Multiply (n, R_total, R10, atemp);
  One_Minus (n, atemp);
  Matrix_Inverse (n, atemp, btemp);
  Left_Diagonal_Multiply (n, T10, btemp, atemp);
  Matrix_Multiply (n, atemp, T_total, btemp);
  Right_Diagonal_Multiply (n, btemp, T01, T_total);
  Matrix_Multiply (n, atemp, R_total, btemp);
  Right_Diagonal_Multiply (n, btemp, T01, R_total);
  for (i = 1; i <= n; i++)
    {
      temp = twoaw[i];
      R_total[i][i] += R01[i] / (temp * temp);
    }
}




void
Sp_RT (int n, struct AD_slab_type slab, double *ur1, double *ut1, double *uru,
       double *utu)
{
  double mu_outside, r, t;
  int i;

  *uru = 0;
  *utu = 0;

  for (i = 1; i <= n; i++)
    {
      mu_outside = Cos_Snell (slab.n_slab, angle[i], 1.0);
      if (mu_outside != 0)
	{
	  Sp_mu_RT (slab.n_top_slide, slab.n_slab, slab.n_bottom_slide,
		    slab.b_top_slide, slab.b, slab.b_bottom_slide, mu_outside,
		    &r, &t);
	  *uru += twoaw[i] * r;
	  *utu += twoaw[i] * t;
	}
    }

  Sp_mu_RT (slab.n_top_slide, slab.n_slab, slab.n_bottom_slide,
	    slab.b_top_slide, slab.b, slab.b_bottom_slide, slab.cos_angle,
	    ur1, ut1);

  *uru *= slab.n_slab * slab.n_slab;
  *utu *= slab.n_slab * slab.n_slab;
}
