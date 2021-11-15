
#define TOP_BOUNDARY 0
#define BOTTOM_BOUNDARY 1 \




void Init_Boundary (struct AD_slab_type slab, int n,
		    double *R01, double *R10, double *T01, double *T10,
		    char boundary);

void Boundary_RT (double n_i, double n_g, double n_t, int n, double b,
		  double *R, double *T);

void Add_Top (int n, double *R01, double *R10, double *T01, double *T10,
	      double **R12, double **R21, double **T12, double **T21,
	      double **R02, double **R20, double **T02, double **T20,
	      double **atemp, double **btemp);

void Add_Bottom (int n, double **R01, double **R10, double **T01,
		 double **T10, double *R12, double *R21, double *T12,
		 double *T21, double **R02, double **R20, double **T02,
		 double **T20, double **atemp, double **btemp);

void Add_Slides (int n, double *R01, double *R10, double *T01, double *T10,
		 double **R, double **T,
		 double **R_total, double **T_total,
		 double **atemp, double **btemp);

void Sp_RT (int n, struct AD_slab_type slab, double *ur1, double *ut1,
	    double *uru, double *utu);
