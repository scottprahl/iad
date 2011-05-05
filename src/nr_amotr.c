#define NRANSI
#include "nr_amotr.h"
#include "nr_util.h"

double 
amotry (double **p, double y[], double psum[], int ndim,
	double (*funk) (double[]), int ihi, double fac)
{
  int j;
  double fac1, fac2, ytry, *ptry;

  ptry = dvector (1, ndim);
  fac1 = (1.0 - fac) / ndim;
  fac2 = fac1 - fac;
  for (j = 1; j <= ndim; j++)
    ptry[j] = psum[j] * fac1 - p[ihi][j] * fac2;
  ytry = (*funk) (ptry);
  if (ytry < y[ihi])
    {
      y[ihi] = ytry;
      for (j = 1; j <= ndim; j++)
	{
	  psum[j] += ptry[j] - p[ihi][j];
	  p[ihi][j] = ptry[j];
	}
    }
  free_dvector (ptry, 1, ndim);
  return ytry;
}

#undef NRANSI
