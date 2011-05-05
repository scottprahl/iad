void amoeba (double **p, double y[], int ndim, double ftol,
	     double (*funk) (double[]), int *nfunk);

double amotry (double **p, double y[], double psum[], int ndim,
	       double (*funk) (double[]), int ihi, double fac);
