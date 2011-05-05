#include "nr_zbrak.h"

void 
zbrak (double (*fx) (double), double x1, double x2, int n,
       double xb1[], double xb2[], int *nb)
{
  int nbb, i;
  double x, fp, fc, dx;

  nbb = 0;
  dx = (x2 - x1) / n;
  fp = (*fx) (x = x1);
  for (i = 1; i <= n; i++)
    {
      fc = (*fx) (x += dx);
      if (fc * fp < 0.0)
	{
	  xb1[++nbb] = x - dx;
	  xb2[nbb] = x;
	  if (*nb == nbb)
	    return;
	}
      fp = fc;
    }
  *nb = nbb;
}
