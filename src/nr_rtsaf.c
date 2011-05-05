#include <math.h>
#include "nr_rtsaf.h"
/* #define MAXIT 100 */

double 
rtsafe (void (*funcd) (double, double *, double *), double x1, double x2, double xacc)
{
  void nrerror (char error_text[]);
  double df, dx, dxold, f, fh, fl;
  double temp, xh, xl, rts;
  double temp1, temp2;
  int j, MAXIT = 100;

  (*funcd) (x1, &fl, &df);
  (*funcd) (x2, &fh, &df);
  if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
    nrerror ("Root must be bracketed in rtsafe");

  if (fl == 0.0)
    return x1;
  if (fh == 0.0)
    return x2;

  if (fl < 0.0)
    {
      xl = x1;
      xh = x2;
    }
  else
    {
      xh = x1;
      xl = x2;
    }
  rts = 0.5 * (x1 + x2);
  dxold = fabs (x2 - x1);
  dx = dxold;
  (*funcd) (rts, &f, &df);

  for (j = 1; j <= MAXIT; j++)
    {

      temp1 = (rts - xh) * df - f;
      temp2 = (rts - xl) * df - f;

      if ((temp1 * temp2 >= 0.0) || (fabs (2.0 * f) > fabs (dxold * df)))
	{
	  dxold = dx;
	  dx = 0.5 * (xh - xl);
	  rts = xl + dx;
	  if (xl == rts)
	    return rts;
	}
      else
	{
	  dxold = dx;
	  dx = f / df;
	  temp = rts;
	  rts -= dx;
	  if (temp == rts)
	    return rts;
	}
      if (fabs (dx) < xacc)
	return rts;

      (*funcd) (rts, &f, &df);

      if (f < 0.0)
	xl = rts;
      else
	xh = rts;
    }
  nrerror ("Maximum number of iterations exceeded in rtsafe");
  return 0.0;			/* never get here */
}
