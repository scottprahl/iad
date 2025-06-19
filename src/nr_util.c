#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "nr_util.h"
#define NR_END 1
#define FREE_ARG char*

void nrerror(char error_text[])
{
    fprintf(stderr, "Numerical Recipes run-time error...\n");
    fprintf(stderr, "%s\n", error_text);
    fprintf(stderr, "...now exiting to system...\n");
    exit(1);
}

double *vector(long nl, long nh)
{
    double *v;

    v = (double *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(double)));
    if (!v)
        nrerror("allocation failure in vector()");
    return v - nl + NR_END;
}

int *ivector(long nl, long nh)
{
    int *v;

    v = (int *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(int)));
    if (!v)
        nrerror("allocation failure in ivector()");
    return v - nl + NR_END;
}

unsigned char *cvector(long nl, long nh)
{
    unsigned char *v;

    v = (unsigned char *)
        malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(unsigned char)));
    if (!v)
        nrerror("allocation failure in cvector()");
    return v - nl + NR_END;
}

unsigned long *lvector(long nl, long nh)
{
    unsigned long *v;

    v = (unsigned long *)
        malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(long)));
    if (!v)
        nrerror("allocation failure in lvector()");
    return v - nl + NR_END;
}

double *dvector(long nl, long nh)
{
    double *v;

    v = (double *) malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(double)));
    if (!v)
        nrerror("allocation failure in dvector()");
    return v - nl + NR_END;
}

double **matrix(long nrl, long nrh, long ncl, long nch)
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    double **m;

    m = (double **) malloc((size_t) ((nrow + NR_END) * sizeof(double *)));
    if (!m)
        nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;

    m[nrl] = (double *) malloc((size_t) ((nrow * ncol + NR_END) * sizeof(double)));
    if (!m[nrl])
        nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    double **m;

    m = (double **) malloc((size_t) ((nrow + NR_END) * sizeof(double *)));
    if (!m)
        nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;

    m[nrl] = (double *) malloc((size_t) ((nrow * ncol + NR_END) * sizeof(double)));
    if (!m[nrl])
        nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    int **m;

    m = (int **) malloc((size_t) ((nrow + NR_END) * sizeof(int *)));
    if (!m)
        nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;

    m[nrl] = (int *) malloc((size_t) ((nrow * ncol + NR_END) * sizeof(int)));
    if (!m[nrl])
        nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    return m;
}

double **submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch, long newrl, long newcl)
{
    long i, j, nrow = oldrh - oldrl + 1, ncol = oldcl - newcl;
    double **m;
    (void) oldch;               // explicit cast to avoid unused error

    m = (double **) malloc((size_t) ((nrow + NR_END) * sizeof(double *)));
    if (!m)
        nrerror("allocation failure in submatrix()");
    m += NR_END;
    m -= newrl;

    for (i = oldrl, j = newrl; i <= oldrh; i++, j++)
        m[j] = a[i] + ncol;

    return m;
}

double **convert_matrix(double *a, long nrl, long nrh, long ncl, long nch)
/* declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
    long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    double **m;

    m = (double **) malloc((size_t) ((nrow + NR_END) * sizeof(double *)));
    if (!m)
        nrerror("allocation failure in convert_matrix()");
    m += NR_END;
    m -= nrl;

    m[nrl] = a - ncl;
    for (i = 1, j = nrl + 1; i < nrow; i++, j++)
        m[j] = m[j - 1] + ncol;
    return m;
}

double ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
    long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
    double ***t;

    t = (double ***) malloc((size_t) ((nrow + NR_END) * sizeof(double **)));
    if (!t)
        nrerror("allocation failure 1 in f3tensor()");
    t += NR_END;
    t -= nrl;

    t[nrl] = (double **)
        malloc((size_t) ((nrow * ncol + NR_END) * sizeof(double *)));
    if (!t[nrl])
        nrerror("allocation failure 2 in f3tensor()");
    t[nrl] += NR_END;
    t[nrl] -= ncl;

    t[nrl][ncl] = (double *)
        malloc((size_t) ((nrow * ncol * ndep + NR_END) * sizeof(double)));
    if (!t[nrl][ncl])
        nrerror("allocation failure 3 in f3tensor()");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;

    for (j = ncl + 1; j <= nch; j++)
        t[nrl][j] = t[nrl][j - 1] + ndep;
    for (i = nrl + 1; i <= nrh; i++) {
        t[i] = t[i - 1] + ncol;
        t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
        for (j = ncl + 1; j <= nch; j++)
            t[i][j] = t[i][j - 1] + ndep;
    }

    return t;
}

void free_vector(double *v, long nl, long nh)
{
    (void) nh;                  // explicit cast to avoid unused error
    free((FREE_ARG) (v + nl - NR_END));
}

void free_ivector(int *v, long nl, long nh)
{
    (void) nh;                  // explicit cast to avoid unused error
    free((FREE_ARG) (v + nl - NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
{
    (void) nh;                  // explicit cast to avoid unused error
    free((FREE_ARG) (v + nl - NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
{
    (void) nh;                  // explicit cast to avoid unused error
    free((FREE_ARG) (v + nl - NR_END));
}

void free_dvector(double *v, long nl, long nh)
{
    (void) nh;                  // explicit cast to avoid unused error
    free((FREE_ARG) (v + nl - NR_END));
}

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
{
    (void) nrh;                 // explicit cast to avoid unused error
    (void) nch;                 // explicit cast to avoid unused error
    free((FREE_ARG) (m[nrl] + ncl - NR_END));
    free((FREE_ARG) (m + nrl - NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
    (void) nrh;                 // explicit cast to avoid unused error
    (void) nch;                 // explicit cast to avoid unused error
    free((FREE_ARG) (m[nrl] + ncl - NR_END));
    free((FREE_ARG) (m + nrl - NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
{
    (void) nrh;                 // explicit cast to avoid unused error
    (void) nch;                 // explicit cast to avoid unused error
    free((FREE_ARG) (m[nrl] + ncl - NR_END));
    free((FREE_ARG) (m + nrl - NR_END));
}

void free_submatrix(double **b, long nrl, long nrh, long ncl, long nch)
{
    (void) nrh;                 // explicit cast to avoid unused error
    (void) nch;                 // explicit cast to avoid unused error
    (void) ncl;                 // explicit cast to avoid unused error
    free((FREE_ARG) (b + nrl - NR_END));
}

void free_convert_matrix(double **b, long nrl, long nrh, long ncl, long nch)
{
    (void) nrh;                 // explicit cast to avoid unused error
    (void) nch;                 // explicit cast to avoid unused error
    (void) ncl;                 // explicit cast to avoid unused error
    free((FREE_ARG) (b + nrl - NR_END));
}

void free_f3tensor(double ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
    (void) nrh;                 // explicit cast to avoid unused error
    (void) nch;                 // explicit cast to avoid unused error
    (void) ndh;                 // explicit cast to avoid unused error
    free((FREE_ARG) (t[nrl][ncl] + ndl - NR_END));
    free((FREE_ARG) (t[nrl] + ncl - NR_END));
    free((FREE_ARG) (t + nrl - NR_END));
}
