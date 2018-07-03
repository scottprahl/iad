#include <math.h>
#include "nr_amoeb.h"
#include "nr_amotr.h"
#include "stdio.h"
#include "ad_globl.h"
#include "iad_type.h"

#define NRANSI
#include "nr_util.h"
#define GET_PSUM \
for (j=1;j<=ndim;j++) {\
for (sum=0.0,i=1;i<=mpts;i++) sum += p[i][j];\
psum[j]=sum;}

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
/* p[1..ndim+1][1..ndim] 
   for a 2D problem this gives three rows, each with a point
   thus      (x1,y1)         funk(x1,y1)
         p = (x2,y2)     y = funk(x2,y2)
             (x3,y3)         funk(x3,y3)
    */

void 
amoeba (double **p, double y[], int ndim, double ftol,
    double (*funk) (double[]), int *nfunk)
{
    int i, ihi, ilo, inhi, j, mpts = ndim + 1;
    double sum, swap, ysave, ytry, *psum;
/*    double rtol;*/
   
    psum = dvector (1, ndim);
    *nfunk = 0;
    GET_PSUM
    for (;;) {
        ilo = 1;
        ihi = y[1] > y[2] ? (inhi = 2, 1) : (inhi = 1, 2);
        
        for (i = 1; i <= mpts; i++) {
            if (y[i] <= y[ilo])
            ilo = i;
            if (y[i] > y[ihi]) {
                inhi = ihi;
                ihi = i;
            }
            else if (y[i] > y[inhi] && i != ihi)
                inhi = i;
        }
        
/*        rtol = 2.0 * fabs (y[ihi] - y[ilo]) / (fabs (y[ihi]) + fabs (y[ilo])); */
        if (y[ilo] < ftol) {
            SWAP (y[1], y[ilo])
            for (i = 1; i <= ndim; i++)
                SWAP (p[1][i], p[ilo][i])
            break;
        }

        if (*nfunk >= IAD_MAX_ITERATIONS) {   
            /* fprintf (stderr, "Max number of iterations exceeded ... \n"); */
            free_dvector (psum, 1, ndim);
            return;
        }
        
        *nfunk += 2;
        ytry = amotry (p, y, psum, ndim, funk, ihi, -1.0);
        if (ytry <= y[ilo])
            ytry = amotry (p, y, psum, ndim, funk, ihi, 2.0);

        else if (ytry >= y[inhi]) {
            ysave = y[ihi];
            ytry = amotry (p, y, psum, ndim, funk, ihi, 0.5);
            if (ytry >= ysave) {

                for (i = 1; i <= mpts; i++) {
                    if (i != ilo) {
                        for (j = 1; j <= ndim; j++)
                            p[i][j] = psum[j] = 0.5 * (p[i][j] + p[ilo][j]);
                        y[i] = (*funk) (psum);
                    }
                }
                *nfunk += ndim;
                GET_PSUM
            }
        }
        else
            --(*nfunk);
    }
    free_dvector (psum, 1, ndim);
/*  fprintf(stderr, "<%5d>\n", *nfunk); */
}

#undef SWAP
#undef GET_PSUM
#undef NRANSI
