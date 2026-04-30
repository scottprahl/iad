/* Nelder-Mead simplex optimiser — scipy-compatible implementation.
 *
 * Replaces the Numerical Recipes amoeba() with the algorithm used by
 * scipy.optimize.minimize(..., method='Nelder-Mead').  Key differences
 * from the NR version:
 *
 *  1. Convergence criterion:
 *       OLD (NR):  y[ilo] < ftol   — stops as soon as the best function value
 *                                   is near zero, even while the simplex is
 *                                   still broad.
 *       NEW (scipy):  max|f(xi)-f(x_best)| <= ftol  AND
 *                     max|xi[j]-x_best[j]| <= ftol
 *                   — stops once the simplex has contracted sufficiently,
 *                     regardless of the absolute minimum value.
 *
 *  2. Both outside AND inside contraction steps are performed (the NR code
 *     only did outside contraction).
 *
 *  3. Shrinkage (multi-point contraction) step is included.
 *
 * Parameters (all 1-indexed, NR/IAD convention):
 *   p[1..ndim+1][1..ndim]  simplex vertices (modified in place)
 *   y[1..ndim+1]           function values at simplex vertices (modified in place)
 *   ndim                   number of parameters (2 for find_ab/find_ag/find_bg)
 *   ftol                   convergence tolerance for both f-spread and x-spread
 *   funk                   objective function:  funk(x[1..ndim]) -> double
 *   nfunk                  output: total function evaluations used
 */

#include <math.h>
#include <stdlib.h>
#include "nr_amoeb.h"
#include "nr_util.h"
#include "ad_globl.h"
#include "iad_type.h"

/* scipy default Nelder-Mead coefficients */
#define NM_RHO   1.0            /* reflection   */
#define NM_CHI   2.0            /* expansion    */
#define NM_PSI   0.5            /* contraction  */
#define NM_SIGMA 0.5            /* shrinkage    */

/* Clip point x[1..ndim] to [lo[j], hi[j]] when bounds are provided. */
static void clip_to_bounds(double x[], int ndim, double lo[], double hi[])
{
    int j;
    if (lo == NULL)
        return;
    for (j = 1; j <= ndim; j++) {
        if (x[j] < lo[j])
            x[j] = lo[j];
        if (x[j] > hi[j])
            x[j] = hi[j];
    }
}

static void sort_simplex(double **p, double y[], int mpts, int ndim)
{
    int i, j, k;
    double *xtmp = dvector(1, ndim);
    for (i = 2; i <= mpts; i++) {
        double ytmp = y[i];
        for (k = 1; k <= ndim; k++)
            xtmp[k] = p[i][k];
        for (j = i - 1; j >= 1 && y[j] > ytmp; j--) {
            y[j + 1] = y[j];
            for (k = 1; k <= ndim; k++)
                p[j + 1][k] = p[j][k];
        }
        y[j + 1] = ytmp;
        for (k = 1; k <= ndim; k++)
            p[j + 1][k] = xtmp[k];
    }
    free_dvector(xtmp, 1, ndim);
}

void amoeba(double **p, double y[], int ndim,
    double ftol, double (*funk)(double[]), int *nfunk, double lo[], double hi[])
{
    int i, j, mpts = ndim + 1;
    double *xbar, *xr, *xe, *xc;
    double fxr, fxe, fxc;
    double fatol, xatol;
    int do_shrink;

    xbar = dvector(1, ndim);
    xr = dvector(1, ndim);
    xe = dvector(1, ndim);
    xc = dvector(1, ndim);

    /* Count the mpts initial evaluations that the caller already did. */
    *nfunk = mpts;
    sort_simplex(p, y, mpts, ndim);

    for (;;) {
        /* ---- convergence check ---- */
        /* Full scipy-style check: both function-value spread and parameter
           spread must be small. */
        fatol = 0.0;
        for (i = 1; i <= mpts; i++) {
            double df = fabs(y[i] - y[1]);
            if (df > fatol)
                fatol = df;
        }
        xatol = 0.0;
        for (i = 2; i <= mpts; i++) {
            for (j = 1; j <= ndim; j++) {
                double dx = fabs(p[i][j] - p[1][j]);
                if (dx > xatol)
                    xatol = dx;
            }
        }
        if (fatol <= ftol && xatol <= ftol)
            break;
        if (*nfunk >= IAD_MAX_ITERATIONS)
            break;

        /* ---- centroid of all vertices except worst ---- */
        for (j = 1; j <= ndim; j++) {
            double sum = 0.0;
            for (i = 1; i <= mpts; i++)
                if (i != mpts)
                    sum += p[i][j];
            xbar[j] = sum / ndim;
        }

        /* ---- reflection ---- */
        for (j = 1; j <= ndim; j++)
            xr[j] = (1.0 + NM_RHO) * xbar[j] - NM_RHO * p[mpts][j];
        clip_to_bounds(xr, ndim, lo, hi);
        fxr = (*funk) (xr);
        (*nfunk)++;

        do_shrink = 0;

        if (fxr < y[1]) {
            /* ---- expansion ---- */
            for (j = 1; j <= ndim; j++)
                xe[j] = (1.0 + NM_RHO * NM_CHI) * xbar[j]
                    - NM_RHO * NM_CHI * p[mpts][j];
            clip_to_bounds(xe, ndim, lo, hi);
            fxe = (*funk) (xe);
            (*nfunk)++;
            if (fxe < fxr) {
                for (j = 1; j <= ndim; j++)
                    p[mpts][j] = xe[j];
                y[mpts] = fxe;
            }
            else {
                for (j = 1; j <= ndim; j++)
                    p[mpts][j] = xr[j];
                y[mpts] = fxr;
            }

        }
        else if (fxr < y[mpts - 1]) {
            /* ---- accept reflection ---- */
            for (j = 1; j <= ndim; j++)
                p[mpts][j] = xr[j];
            y[mpts] = fxr;

        }
        else {
            /* ---- contraction ---- */
            if (fxr < y[mpts]) {
                /* outside contraction */
                for (j = 1; j <= ndim; j++)
                    xc[j] = (1.0 + NM_PSI * NM_RHO) * xbar[j]
                        - NM_PSI * NM_RHO * p[mpts][j];
                clip_to_bounds(xc, ndim, lo, hi);
                fxc = (*funk) (xc);
                (*nfunk)++;
                if (fxc <= fxr) {
                    for (j = 1; j <= ndim; j++)
                        p[mpts][j] = xc[j];
                    y[mpts] = fxc;
                }
                else {
                    do_shrink = 1;
                }
            }
            else {
                /* inside contraction */
                for (j = 1; j <= ndim; j++)
                    xc[j] = (1.0 - NM_PSI) * xbar[j]
                        + NM_PSI * p[mpts][j];
                clip_to_bounds(xc, ndim, lo, hi);
                fxc = (*funk) (xc);
                (*nfunk)++;
                if (fxc < y[mpts]) {
                    for (j = 1; j <= ndim; j++)
                        p[mpts][j] = xc[j];
                    y[mpts] = fxc;
                }
                else {
                    do_shrink = 1;
                }
            }
        }

        if (do_shrink) {
            /* ---- shrinkage: contract all vertices towards best ---- */
            for (i = 2; i <= mpts; i++) {
                for (j = 1; j <= ndim; j++)
                    p[i][j] = p[1][j]
                        + NM_SIGMA * (p[i][j] - p[1][j]);
                clip_to_bounds(p[i], ndim, lo, hi);
                y[i] = (*funk) (p[i]);
                (*nfunk)++;
            }
        }

        sort_simplex(p, y, mpts, ndim);
    }

    free_dvector(xbar, 1, ndim);
    free_dvector(xr, 1, ndim);
    free_dvector(xe, 1, ndim);
    free_dvector(xc, 1, ndim);
}
