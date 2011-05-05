#include <math.h>
#include "nr_gaulg.h"
#define EPS 3.0e-11
#define Pi              3.14159265358979323

void 
gauleg (double x1, double x2, double x[], double w[], int n)
{
  int m, j, i;
  double z1, z, xm, xl, pp, p3, p2, p1;

  m = (n + 1) / 2;
  xm = 0.5 * (x2 + x1);
  xl = 0.5 * (x2 - x1);
  for (i = 1; i <= m; i++)
    {
      z = cos (Pi * (i - 0.25) / (n + 0.5));
      do
	{
	  p1 = 1.0;
	  p2 = 0.0;
	  for (j = 1; j <= n; j++)
	    {
	      p3 = p2;
	      p2 = p1;
	      p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
	    }
	  pp = n * (z * p1 - p2) / (z * z - 1.0);
	  z1 = z;
	  z = z1 - p1 / pp;
	}
      while (fabs (z - z1) > EPS);
      x[i] = xm - xl * z;
      x[n + 1 - i] = xm + xl * z;
      w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
      w[n + 1 - i] = w[i];
    }
}
