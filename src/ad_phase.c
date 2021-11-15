
#include <stdlib.h>
#include <math.h>
#include "nr_util.h"
#include "ad_globl.h"
#include "ad_phase.h"



void
Get_Phi (int n, int phase_function, double g, double **h)
{

  int i, j, k;
  double g2M, gk, x;
  double *chi;
  double **p;



  if (g != 0 && phase_function != HENYEY_GREENSTEIN)
    AD_error
      ("Only the Henyey-Greenstein phase function has been implemented\n");

  if (fabs (g) >= 1)
    AD_error ("Get_Phi was called with a bad g_calc value");



  for (i = -n; i <= n; i++)
    for (j = -n; j <= n; j++)
      h[i][j] = 1;


  for (i = -n; i <= n; i++)
    {
      h[i][0] = 0.0;
      h[0][i] = 0.0;
    }



  if (g == 0)
    return;



  chi = dvector (1, n);
  g2M = pow (g, n);
  gk = 1.0;
  for (k = 1; k < n; k++)
    {
      gk *= g;
      chi[k] = (2 * k + 1) * (gk - g2M) / (1 - g2M);
    }




  p = dmatrix (0, n, -n, n);



  for (j = 1; j <= n; j++)
    {
      p[0][j] = 1;
      x = angle[j];
      p[1][j] = x;
      for (k = 1; k < n; k++)
	p[k + 1][j] = ((2 * k + 1) * x * p[k][j] - k * p[k - 1][j]) / (k + 1);
    }



  for (j = 1; j <= n; j++)
    for (k = 1; k < n; k++)
      {
	p[k][-j] = -p[k][j];
	k++;
	p[k][-j] = p[k][j];
      }





  for (i = 1; i <= n; i++)
    {
      for (j = i; j <= n; j++)
	{
	  for (k = 1; k < n; k++)
	    {
	      h[i][j] += chi[k] * p[k][i] * p[k][j];
	      h[-i][j] += chi[k] * p[k][-i] * p[k][j];
	    }
	}
    }



  for (i = n; i >= 2; i--)
    for (j = 1; j < i; j++)
      {
	h[-i][j] = h[-j][i];
	h[-i][-j] = h[j][i];
      }

  for (i = 1; i <= n; i++)
    h[-i][-i] = h[i][i];

  for (i = -n; i <= n; i++)
    for (j = i + 1; j <= n; j++)
      h[j][i] = h[i][j];



  free_dmatrix (p, 0, n, -n, n);
  free_dvector (chi, 1, n);
}
