
#include <math.h>
#include <float.h>
#include <stdio.h>
#include "nr_util.h"
#include "ad_globl.h"
#include "ad_matrx.h"
#include "ad_bound.h"
#include "ad_doubl.h"
#include "ad_start.h"



void
RT_Cone (int n,
	 struct AD_slab_type *slab,
	 int use_cone, double *UR1, double *UT1, double *URU, double *UTU)
{

  struct AD_method_type method;
  double *R01, *R10, *T01, *T10;
  double *R23, *R32, *T23, *T32;
  double **R12, **T12;
  double **R02, **T02, **T20, **R20;
  double **R03, **T03, **T30, **R30;
  double **atemp, **btemp;
  double d;
  *UR1 = -1;
  *URU = -1;
  *UT1 = -1;
  *UTU = -1;



  if (slab->n_slab < 0)
    return;
  if (slab->n_top_slide < 0)
    return;
  if (slab->n_bottom_slide < 0)
    return;
  if (slab->a < 0 || slab->a > 1)
    return;
  if (slab->g < -1 || slab->g > 1)
    return;
  if (slab->b < 0)
    return;
  if (slab->cos_angle < 0 || slab->cos_angle > 1)
    return;


  n = 12 * (n / 12);
  if (n < 12)
    n = 12;
  method.quad_pts = n;



  R12 = dmatrix (1, n, 1, n);
  T12 = dmatrix (1, n, 1, n);
  R02 = dmatrix (1, n, 1, n);
  T02 = dmatrix (1, n, 1, n);
  R20 = dmatrix (1, n, 1, n);
  T20 = dmatrix (1, n, 1, n);
  R03 = dmatrix (1, n, 1, n);
  T03 = dmatrix (1, n, 1, n);
  R30 = dmatrix (1, n, 1, n);
  T30 = dmatrix (1, n, 1, n);
  atemp = dmatrix (1, n, 1, n);
  btemp = dmatrix (1, n, 1, n);




  Choose_Cone_Method (slab, &method);

  if (slab->b <= 0)
    {
      Zero_Layer (n, R12, T12);
      return;
    }

  n = method.quad_pts;
  Init_Layer (*slab, method, R12, T12);

  d = 1.0;
  if (slab->b != HUGE_VAL)
    d = method.b_thinnest * slab->b / method.b_calc;

  Double_Until (n, R12, T12, d, slab->b);



  R01 = dvector (1, n);
  R10 = dvector (1, n);
  T01 = dvector (1, n);
  T10 = dvector (1, n);
  Init_Boundary (*slab, n, R01, R10, T01, T10, TOP_BOUNDARY);

  R23 = dvector (1, n);
  R32 = dvector (1, n);
  T23 = dvector (1, n);
  T32 = dvector (1, n);
  Init_Boundary (*slab, n, R23, R32, T23, T32, BOTTOM_BOUNDARY);





  Add_Top (n, R01, R10, T01, T10, R12, R12, T12, T12, R02, R20, T02, T20,
	   atemp, btemp);
  Add_Bottom (n, R02, R20, T02, T20, R23, R32, T23, T32, R03, R30, T03, T30,
	      atemp, btemp);

  if (use_cone == CONE)
    {
      URU_and_UR1_Cone (n, slab->n_slab, slab->cos_angle, R03, URU, UR1);
      Transpose_Matrix (n, T03);
      URU_and_UR1_Cone (n, slab->n_slab, slab->cos_angle, T03, UTU, UT1);
    }
  else
    {
      if (use_cone != OBLIQUE)
	fprintf (stderr,
		 "Unknown type for use_cone.  Assuming oblique incidence.\n");
      URU_and_URx_Cone (n, slab->n_slab, slab->cos_angle, R03, URU, UR1);
      Transpose_Matrix (n, T03);
      URU_and_URx_Cone (n, slab->n_slab, slab->cos_angle, T03, UTU, UT1);
    }



  free_dvector (R01, 1, n);
  free_dvector (R10, 1, n);
  free_dvector (T01, 1, n);
  free_dvector (T10, 1, n);

  free_dmatrix (R12, 1, n, 1, n);
  free_dmatrix (T12, 1, n, 1, n);

  free_dmatrix (R03, 1, n, 1, n);
  free_dmatrix (R30, 1, n, 1, n);
  free_dmatrix (T03, 1, n, 1, n);
  free_dmatrix (T30, 1, n, 1, n);

  free_dmatrix (R02, 1, n, 1, n);
  free_dmatrix (R20, 1, n, 1, n);
  free_dmatrix (T02, 1, n, 1, n);
  free_dmatrix (T20, 1, n, 1, n);

  free_dmatrix (atemp, 1, n, 1, n);
  free_dmatrix (btemp, 1, n, 1, n);

  free_dvector (R32, 1, n);
  free_dvector (R23, 1, n);
  free_dvector (T32, 1, n);
  free_dvector (T23, 1, n);


}




void
ez_RT_Cone (int n,
	    double nslab,
	    double ntopslide,
	    double nbottomslide,
	    double a,
	    double b,
	    double g,
	    double cos_cone_angle,
	    double *UR1, double *UT1, double *URU, double *UTU)
{
  struct AD_slab_type slab;

  slab.n_slab = nslab;
  slab.n_top_slide = ntopslide;
  slab.n_bottom_slide = nbottomslide;
  slab.b_top_slide = 0;
  slab.b_bottom_slide = 0;
  slab.a = a;
  slab.b = b;
  slab.g = g;
  slab.cos_angle = cos_cone_angle;
  slab.phase_function = HENYEY_GREENSTEIN;

  RT_Cone (n, &slab, CONE, UR1, UT1, URU, UTU);
}




void
ez_RT_Oblique (int n,
	       double nslab,
	       double ntopslide,
	       double nbottomslide,
	       double a,
	       double b,
	       double g,
	       double cos_oblique_angle,
	       double *URx, double *UTx, double *URU, double *UTU)
{
  struct AD_slab_type slab;

  slab.n_slab = nslab;
  slab.n_top_slide = ntopslide;
  slab.n_bottom_slide = nbottomslide;
  slab.b_top_slide = 0;
  slab.b_bottom_slide = 0;
  slab.a = a;
  slab.b = b;
  slab.g = g;
  slab.cos_angle = cos_oblique_angle;
  slab.phase_function = HENYEY_GREENSTEIN;

  RT_Cone (n, &slab, OBLIQUE, URx, UTx, URU, UTU);
}
