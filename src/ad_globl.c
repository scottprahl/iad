
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ad_globl.h"
#include "ad_frsnl.h"


#define AD_GLOBAL_SOURCE
double angle[MAX_QUAD_PTS + 1];
double weight[MAX_QUAD_PTS + 1];
double twoaw[MAX_QUAD_PTS + 1];
int Martin_Hammer = 0;




void
Zero_Layer (int n, double **r, double **t)
{
  int i, j;

  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      {
	t[i][j] = 0.0;
	r[i][j] = 0.0;
      }

  for (i = 1; i <= n; i++)
    t[i][i] = 1 / twoaw[i];
}




void
AD_error (char error_text[])
{
  fprintf (stderr, "Adding-Doubling error\n");
  fprintf (stderr, "%s\n", error_text);
  fprintf (stderr, "...now exiting to system...\n");
  exit (1);
}




void
URU_and_UR1 (int n, double n_slab, double **R, double *URU, double *UR1)
{
  URU_and_UR1_Cone (n, n_slab, 0.0, R, URU, UR1);
}





void
URU_and_UR1_Cone (int n, double n_slab, double mu, double **R, double *URU,
		  double *UR1)
{
  int i, j, last_j;
  double mu_slab;
  double temp = 0.0;

  if (n_slab == 1)
    mu_slab = mu;
  else
    mu_slab = sqrt (n_slab * n_slab - 1 + mu * mu) / n_slab;

  last_j = 1;
  while (angle[last_j] <= mu_slab)
    last_j++;

  *URU = 0.0;
  for (i = 1; i <= n; i++)
    {
      temp = 0.0;
      for (j = last_j; j <= n; j++)
	temp += R[i][j] * twoaw[j];
      *URU += temp * twoaw[i];
    }
  *UR1 = temp;
  *URU *= n_slab * n_slab / (1 - mu * mu);
}




void
URU_and_URx_Cone (int n, double n_slab, double mu, double **R, double *URU,
		  double *URx)
{
  int i, j, cone_index;
  double mu_slab, urx, delta, closest_delta;
  double degrees = 180.0 / 3.1415926535;

  mu_slab = sqrt (n_slab * n_slab - 1 + mu * mu) / n_slab;

  closest_delta = 1;
  cone_index = n;

  for (i = n; i >= 1; i--)
    {
      delta = fabs (angle[i] - mu_slab);
      if (delta < closest_delta)
	{
	  closest_delta = delta;
	  cone_index = i;
	}
    }

  if (fabs (angle[cone_index] - mu_slab) > 1e-5)
    {
      fprintf (stderr, "Something is wrong with the quadrature\n");
      fprintf (stderr, "theta_i = %5.2f degrees or ", acos (mu) * degrees);
      fprintf (stderr, "cos(theta_i) = %8.5f\n", mu);
      fprintf (stderr, "theta_t = %5.2f degrees or ",
	       acos (mu_slab) * degrees);
      fprintf (stderr, "cos(theta_t) = %8.5f\n", mu_slab);
      fprintf (stderr, " index  degrees cosine\n");
      for (i = n; i >= 1; i--)
	{
	  fprintf (stderr, " %5d   %5.2f ", i, acos (angle[i]) * degrees);
	  fprintf (stderr, " %8.5f\n", angle[i]);
	}

      fprintf (stderr, "Closest quadrature angle is i=%5d ", cone_index);
      fprintf (stderr, "or cos(theta)=%8.5f\n", angle[cone_index]);
      fprintf (stderr, "Assuming normal incidence\n");
    }

  *URU = 0.0;
  for (i = 1; i <= n; i++)
    {
      urx = 0.0;
      for (j = 1; j <= n; j++)
	urx += R[i][j] * twoaw[j];

      *URU += urx * twoaw[i];
      if (i == cone_index)
	*URx = urx;
    }
  *URU *= n_slab * n_slab;
}




void
UFU_and_UF1 (int n, double n_slab,
	     double **Lup, double **Ldown, double *UFU, double *UF1)
{
  int i, j;
  double temp = 0.0;

  *UFU = 0.0;
  for (j = 1; j <= n; j++)
    {
      temp = 0.0;
      for (i = 1; i <= n; i++)
	temp += (Lup[i][j] + Ldown[i][j]) * 2 * weight[i];
      *UFU += twoaw[j] * temp;
    }
  *UF1 = temp * n_slab * n_slab;
  *UFU *= n_slab * n_slab / 2;
}





void
wrmatrix (int n, double **a)
{
  int i, j;
  double tflux, flux;

  printf ("%9.5f", 0.0);
  for (i = 1; i <= n; i++)
    printf ("%9.5f", angle[i]);

  printf ("     flux\n");

  tflux = 0.0;
  for (i = 1; i <= n; i++)
    {
      printf ("%9.5f", angle[i]);
      for (j = 1; j <= n; j++)
	if ((a[i][j] > 10) || (a[i][j] < -10))
	  printf ("    *****");
	else
	  printf ("%9.5f", a[i][j]);
      flux = 0.0;
      for (j = 1; j <= n; j++)
	if ((a[i][j] < 10) && (a[i][j] > -10))
	  flux += a[i][j] * twoaw[j];
      printf ("%9.5f\n", flux);
      tflux += flux * twoaw[i];
    }

  printf ("%9s", "flux   ");
  for (i = 1; i <= n; i++)
    {
      flux = 0.0;
      for (j = 1; j <= n; j++)
	if ((a[j][i] < 10) && (a[j][i] > -10))
	  flux += a[j][i] * twoaw[j];
      printf ("%9.5f", flux);
    }
  printf ("%9.5f\n", tflux);
  for (i = 1; i <= (n + 2); i++)
    printf ("*********");
  printf ("\n\n");
}





void
wrarray (int n, double *a)
{
  int i;
  double sum;

  for (i = 1; i <= n; i++)
    printf ("%9.5f", angle[i]);
  printf ("%9s\n", " angles");

  sum = 0.0;
  for (i = 1; i <= n; i++)
    {
      if (a[i] > 10 || a[i] < -10)
	printf ("    *****");
      else
	printf ("%9.5f", a[i]);
      if (a[i] < 10 && a[i] < -10)
	sum += a[i];
    }
  printf ("%9.5f", sum);
  printf ("%9s\n", " (natural)");

  sum = 0.0;
  for (i = 1; i <= n; i++)
    {
      if (a[i] > 10 || a[i] < -10)
	printf ("    *****");
      else
	printf ("%9.5f", a[i] / twoaw[i]);
      if (a[i] < 10 && a[i] < -10)
	sum += a[i];
    }
  printf ("%9.5f", sum);
  printf ("%9s\n", "*2aw");
  for (i = 1; i <= (n + 2); i++)
    printf ("*********");
  printf ("\n\n");
}
