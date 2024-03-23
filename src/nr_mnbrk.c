#include <math.h>
#define NRANSI
#include "nr_util.h"
#include "nr_mnbrk.h"
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
static double   maxarg1, maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

void
mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
    double (*func) (double), int *iterations)
{
    double          ulim, u, r, q, fu, dum;

    *fa = (*func) (*ax);
    *fb = (*func) (*bx);
    if (*fb > *fa)
    {
        SHFT(dum, *ax, *bx, dum)
        SHFT(dum, *fb, *fa, dum)
    }
    *cx = (*bx) + GOLD * (*bx - *ax);
    *fc = (*func) (*cx);
    *iterations=3;
    while (*fb > *fc)
    {
        r = (*bx - *ax) * (*fb - *fc);
        q = (*bx - *cx) * (*fb - *fa);
        u = (*bx) - ((*bx - *cx) * q - (*bx - *ax) * r) /
            (2.0 * SIGN(FMAX(fabs(q - r), TINY), q - r));
        ulim = (*bx) + GLIMIT * (*cx - *bx);
        if ((*bx - u) * (u - *cx) > 0.0)
        {
            fu = (*func) (u);
            (*iterations)++;
            if (fu < *fc)
            {
                *ax = (*bx);
                *bx = u;
                *fa = (*fb);
                *fb = fu;
                return;
            }
            else if (fu > *fb)
            {
                *cx = u;
                *fc = fu;
                return;
            }
            u = (*cx) + GOLD * (*cx - *bx);
            fu = (*func) (u);
            (*iterations)++;
        }
        else if ((*cx - u) * (u - ulim) > 0.0)
        {
            fu = (*func) (u);
            (*iterations)++;
            if (fu < *fc)
            {
                SHFT(*bx, *cx, u, *cx + GOLD * (*cx - *bx))
                SHFT(*fb, *fc, fu, (*func) (u))
                (*iterations)++;
            }
        }
        else if ((u - ulim) * (ulim - *cx) >= 0.0)
        {
            u = ulim;
            fu = (*func) (u);
            (*iterations)++;
        }
        else
        {
            u = (*cx) + GOLD * (*cx - *bx);
            fu = (*func) (u);
            (*iterations)++;
        }
        SHFT(*ax, *bx, *cx, u)
        SHFT(*fa, *fb, *fc, fu)
    }
}

#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
#undef NRANSI
