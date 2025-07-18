#include <math.h>
#define NRANSI
#include "nr_util.h"
#include "nr_brent.h"
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double brent(double ax, double bx, double cx, double (*f)(double), double tol, double *xmin, int *iterations)
{
    int iter;
    double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
    double e = 0.0;

    *iterations = 3;
    a = (ax < cx ? ax : cx);
    b = (ax > cx ? ax : cx);
    d = 1;
    x = w = v = bx;
    fw = fv = fx = (*f) (x);
    (*iterations)++;
    for (iter = 1; iter <= ITMAX; iter++) {
        xm = 0.5 * (a + b);
        tol2 = 2.0 * (tol1 = tol * fabs(x) + ZEPS);
        if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) {
            *xmin = x;
            return fx;
        }
        if (fabs(e) > tol1) {
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0)
                p = -p;
            q = fabs(q);
            etemp = e;
            e = d;
            if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x)
                || p >= q * (b - x))
                d = CGOLD * (e = (x >= xm ? a - x : b - x));
            else {
                d = p / q;
                u = x + d;
                if (u - a < tol2 || b - u < tol2)
                    d = SIGN(tol1, xm - x);
            }
        }
        else {
            d = CGOLD * (e = (x >= xm ? a - x : b - x));
        }
        u = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
        fu = (*f) (u);
        (*iterations)++;
        if (fu <= fx) {
            if (u >= x)
                a = x;
            else
                b = x;
        SHFT(v, w, x, u) SHFT(fv, fw, fx, fu)}
        else {
            if (u < x)
                a = u;
            else
                b = u;
            if (fu <= fw || w == x) {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            }
            else if (fu <= fv || v == x || v == w) {
                v = u;
                fv = fu;
            }
        }
    }
    *xmin = x;
    return fx;
}

#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef NRANSI
