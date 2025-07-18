/* adaptive_quadtree.c  ────────────────────────────────────────────────────
   Adaptive quadtree that stores every (x,y,f) triple in the generic cache.
   Requires cache.c / cache.h  (a,b,g,ur1,ut1,uru,utu)
   ---------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cache.h"

/* same Sample struct as before, used only for the final flat array */
typedef struct { double x, y, f; } Sample;

/* ---------------- helper: memoised evaluation via Cache ---------------- */
static double eval_f_cached(Cache *C,
                            double (*f)(double,double),
                            double a, double b, double g)
{
    if (cache_get(C, a, b, g, &ur1, &ut1, NULL, NULL))   /* hit */
        return val;

    f(a, b, g, &ur1, &ut1, &uru, &utu); 
    cache_put(C, a, b, g, ur1, ut1, uru, utu);
    return val;
}

/* ---------------- public API (unchanged signature) --------------------- */
int adaptive_quadtree_sample(
        double (*f)(double, double),
        double x0, double x1, double y0, double y1,
        double tol, int max_depth,
        int min_depth, int init_grid,
        Sample **out, size_t *out_n)
{
    Cache C;
    if (cache_init(&C) != 0) return -1;

    /* optional uniform base grid */
    if (init_grid > 1) {
        for (int ix = 0; ix < init_grid; ++ix) {
            double x = x0 + (x1 - x0) * ix / (init_grid - 1);
            for (int iy = 0; iy < init_grid; ++iy) {
                double y = y0 + (y1 - y0) * iy / (init_grid - 1);
                eval_f_cached(&C, f, x, y);
            }
        }
    }

    /* cell stack */
    typedef struct { double x0,x1,y0,y1; int depth; } Cell;
    size_t cap = 128, top = 0;
    Cell *stack = malloc(cap * sizeof *stack);
    if (!stack) { cache_free(&C); return -1; }
    stack[top++] = (Cell){ x0,x1,y0,y1,0 };

    /* adaptive loop */
    while (top) {
        Cell c = stack[--top];
        double xm = 0.5*(c.x0+c.x1), ym = 0.5*(c.y0+c.y1);

        double f00 = eval_f_cached(&C,f,c.x0,c.y0);
        double f10 = eval_f_cached(&C,f,c.x1,c.y0);
        double f01 = eval_f_cached(&C,f,c.x0,c.y1);
        double f11 = eval_f_cached(&C,f,c.x1,c.y1);
        double fc  = eval_f_cached(&C,f,xm,  ym);

        double interp = 0.25*(f00+f10+f01+f11);
        double err = fabs(fc - interp);

        int split = (c.depth < min_depth) || (err > tol && c.depth < max_depth);
        if (split) {
            if (top+4 > cap) { cap*=2; stack = realloc(stack,cap*sizeof*stack); }
            int d = c.depth+1;
            stack[top++] = (Cell){ c.x0,xm, c.y0,ym, d };
            stack[top++] = (Cell){ xm, c.x1,c.y0,ym, d };
            stack[top++] = (Cell){ c.x0,xm, ym, c.y1,d };
            stack[top++] = (Cell){ xm, c.x1,ym, c.y1,d };
        }
    }
    free(stack);

    /* dump cache → flat array ----------------------------------------- */
    size_t n = cache_size(&C);
    Sample *pts = malloc(n * sizeof *pts);
    if (!pts) { cache_free(&C); return -1; }

    for (size_t i = 0; i < n; ++i) {
        double x,y,g,fv;
        cache_get_by_index(&C, i, &x, &y, &g, &fv, NULL, NULL, NULL);
        pts[i] = (Sample){ x, y, fv };
    }

    cache_free(&C);
    *out = pts;
    *out_n = n;
    return 0;
}

/* ------------- optional stand-alone demo (compile with -DDEMO) --------- */
#ifdef DEMO
static double test_fun(double x,double y)
{
    return sin(3*M_PI*x)*cos(3*M_PI*y)*exp(-((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5))/0.1);
}
#include <stdio.h>
int main(void)
{
    Sample *pts; size_t n;
    if (adaptive_quadtree_sample(test_fun,0,1,0,1,
                                 0.002,7,3,10,&pts,&n)) return 1;
    printf("# %zu samples\n",n);
    for (size_t i=0;i<n;++i) printf("%.10f %.10f %.10f\n",pts[i].x,pts[i].y,pts[i].f);
    free(pts);
}
#endif
