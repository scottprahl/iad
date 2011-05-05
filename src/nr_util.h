#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI)	/* ANSI */

void nrerror (char error_text[]);
double *vector (long nl, long nh);
int *ivector (long nl, long nh);
unsigned char *cvector (long nl, long nh);
unsigned long *lvector (long nl, long nh);
double *dvector (long nl, long nh);
double **matrix (long nrl, long nrh, long ncl, long nch);
double **dmatrix (long nrl, long nrh, long ncl, long nch);
int **imatrix (long nrl, long nrh, long ncl, long nch);
double **submatrix (double **a, long oldrl, long oldrh, long oldcl, long oldch,
		    long newrl, long newcl);
double **convert_matrix (double *a, long nrl, long nrh, long ncl, long nch);
double ***f3tensor (long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_vector (double *v, long nl, long nh);
void free_ivector (int *v, long nl, long nh);
void free_cvector (unsigned char *v, long nl, long nh);
void free_lvector (unsigned long *v, long nl, long nh);
void free_dvector (double *v, long nl, long nh);
void free_matrix (double **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix (double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix (int **m, long nrl, long nrh, long ncl, long nch);
void free_submatrix (double **b, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix (double **b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor (double ***t, long nrl, long nrh, long ncl, long nch,
		    long ndl, long ndh);

#else /* ANSI */
/* traditional - K&R */

void nrerror ();
double *vector ();
double **matrix ();
double **submatrix ();
double **convert_matrix ();
double ***f3tensor ();
double *dvector ();
double **dmatrix ();
int *ivector ();
int **imatrix ();
unsigned char *cvector ();
unsigned long *lvector ();
void free_vector ();
void free_dvector ();
void free_ivector ();
void free_cvector ();
void free_lvector ();
void free_matrix ();
void free_submatrix ();
void free_convert_matrix ();
void free_dmatrix ();
void free_imatrix ();
void free_f3tensor ();

#endif /* ANSI */

#endif /* _NR_UTILS_H_ */
