/* Autogenerated v3-16-3 from https://github.com/scottprahl/iad */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "ad_globl.h"
#include "ad_prime.h"
#include "ad_start.h"
#include "ad_phase.h"
#include "ad_doubl.h"
#include "ad_bound.h"
#include "nr_util.h"

int main(int argc, char **argv)
{

#define N_QUAD_PT 4

    struct AD_slab_type slab;
    struct AD_method_type method;
    double **h, **hpart, **r, **t;
    int i, j;
    double *R01, *R10, *T10, *T01;

    slab.a = 0.9;
    slab.b = 1;
    slab.g = 0.9;
    slab.phase_function = HENYEY_GREENSTEIN;
    slab.n_slab = 1.5;
    slab.n_top_slide = 1.0;
    slab.n_bottom_slide = 1.0;
    slab.b_top_slide = 0.0;
    slab.b_bottom_slide = 0.0;

    method.quad_pts = N_QUAD_PT;
    h = dmatrix(-N_QUAD_PT, N_QUAD_PT, -N_QUAD_PT, N_QUAD_PT);
    hpart = dmatrix(1, N_QUAD_PT, 1, N_QUAD_PT);
    r = dmatrix(1, N_QUAD_PT, 1, N_QUAD_PT);
    t = dmatrix(1, N_QUAD_PT, 1, N_QUAD_PT);

    Choose_Method(&slab, &method);
    printf("The quadrature angles are (cf 5.63)\n");
    wrarray(N_QUAD_PT, angle);
    printf("The quadrature weights are (cf 5.63)\n");
    wrarray(N_QUAD_PT, weight);

    Get_Phi(N_QUAD_PT, HENYEY_GREENSTEIN, method.g_calc, h);

    for (i = 1; i <= N_QUAD_PT; i++) {
        for (j = 1; j <= N_QUAD_PT; j++) {
            hpart[i][j] = h[-i][-j];
        }
    }
    printf("The h-- redistribution matrix is (cf 5.64)\n");
    wrmatrix(N_QUAD_PT, hpart);

    for (i = 1; i <= N_QUAD_PT; i++) {
        for (j = 1; j <= N_QUAD_PT; j++) {
            hpart[i][j] = h[i][-j];
        }
    }
    printf("The h+- redistribution matrix is (cf 5.64)\n");
    wrmatrix(N_QUAD_PT, hpart);

    for (i = 1; i <= N_QUAD_PT; i++) {
        for (j = 1; j <= N_QUAD_PT; j++) {
            hpart[i][j] = h[-i][j];
        }
    }
    printf("The h-+ redistribution matrix is (cf 5.64)\n");
    wrmatrix(N_QUAD_PT, hpart);

    for (i = 1; i <= N_QUAD_PT; i++) {
        for (j = 1; j <= N_QUAD_PT; j++) {
            hpart[i][j] = h[i][j];
        }
    }
    printf("The h++ redistribution matrix is (cf 5.64)\n");
    wrmatrix(N_QUAD_PT, hpart);
    free_matrix(h, -N_QUAD_PT, N_QUAD_PT, -N_QUAD_PT, N_QUAD_PT);
    free_matrix(hpart, 1, N_QUAD_PT, 1, N_QUAD_PT);

    {
        double **temp;

        temp = dmatrix(1, N_QUAD_PT, 1, N_QUAD_PT);

        Martin_Hammer = 1;
        Init_Layer(slab, method, r, t);
        Martin_Hammer = 0;

        for (i = 1; i <= N_QUAD_PT; i++) {
            for (j = 1; j <= N_QUAD_PT; j++) {
                temp[i][j] = r[i][j] * twoaw[j];
            }
        }

        printf("R for dtau=0.25 before dividing by 2mu[i]w[i] (eq 5.54)\n");
        wrmatrix(N_QUAD_PT, temp);

        printf("R for dtau=0.25 after incorporating 2mu[i]w[i] (cf 5.65)\n");
        wrmatrix(N_QUAD_PT, r);

        for (i = 1; i <= N_QUAD_PT; i++) {
            for (j = 1; j <= N_QUAD_PT; j++) {
                temp[i][j] = t[i][j] * twoaw[j];
            }
        }

        printf("T for dtau=0.25 before dividing by 2mu[i]w[i] (eq 5.54)\n");
        wrmatrix(N_QUAD_PT, temp);

        printf("T for dtau=0.25 after incorporating 2mu[i]w[i] (cf 5.66)\n");
        wrmatrix(N_QUAD_PT, t);
        free_matrix(temp, 1, N_QUAD_PT, 1, N_QUAD_PT);
    }

    Double_Once(N_QUAD_PT, r, t);
    printf("R for dtau=0.5 after incorporating 2mu[i]w[i] (cf 5.67)\n");
    wrmatrix(N_QUAD_PT, r);

    printf("T for dtau=0.5 after incorporating 2mu[i]w[i] (cf 5.68)\n");
    wrmatrix(N_QUAD_PT, t);

    Double_Once(N_QUAD_PT, r, t);
    printf("R for dtau=1.0 after incorporating 2mu[i]w[i] (cf 5.69)\n");
    wrmatrix(N_QUAD_PT, r);

    printf("T for dtau=1.0 after incorporating 2mu[i]w[i] (cf 5.70)\n");
    wrmatrix(N_QUAD_PT, t);

    {
        double *temp;

        R01 = dvector(1, N_QUAD_PT);
        R10 = dvector(1, N_QUAD_PT);
        T01 = dvector(1, N_QUAD_PT);
        T10 = dvector(1, N_QUAD_PT);
        temp = dvector(1, N_QUAD_PT);
        Init_Boundary(slab, N_QUAD_PT, R01, R10, T01, T10, TOP_BOUNDARY);

        printf("R01 for the boundary is (cf 5.71)\n");
        wrarray(N_QUAD_PT, R01);

        printf("R10 should be the same\n");
        wrarray(N_QUAD_PT, R10);

        for (i = 1; i <= N_QUAD_PT; i++)
            temp[i] = T01[i] * twoaw[i];

        printf("T01 for the boundary is (cf 5.72)\n");
        wrarray(N_QUAD_PT, temp);

        for (i = 1; i <= N_QUAD_PT; i++)
            temp[i] = T01[i] / twoaw[i];

        printf("T01 for the boundary is (cf 5.72)\n");
        wrarray(N_QUAD_PT, temp);

        printf("My program's T01 is divided by 2 nu w\n");
        wrarray(N_QUAD_PT, T01);

    }

    {
        double **r2, **t2, **a, **b;

        r2 = dmatrix(1, N_QUAD_PT, 1, N_QUAD_PT);
        t2 = dmatrix(1, N_QUAD_PT, 1, N_QUAD_PT);
        a = dmatrix(1, N_QUAD_PT, 1, N_QUAD_PT);
        b = dmatrix(1, N_QUAD_PT, 1, N_QUAD_PT);
        Add_Slides(N_QUAD_PT, R01, R10, T01, T10, r, t, r2, t2, a, b);

        printf("R for dtau=1.0 after adding boundaries (cf 5.73)\n");
        wrmatrix(N_QUAD_PT, r2);

        printf("T for dtau=1.0 after adding boundaries  (cf 5.74)\n");
        wrmatrix(N_QUAD_PT, t2);
    }
    exit(EXIT_SUCCESS);
}
