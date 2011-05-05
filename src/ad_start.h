

double Get_Start_Depth (double mu, double d);

void Choose_Method (struct AD_slab_type *slab, struct AD_method_type *method);

void Choose_Cone_Method (struct AD_slab_type *slab,
			 struct AD_method_type *method);

void Init_Layer (struct AD_slab_type slab, struct AD_method_type method,
		 double **R, double **T);

void Quadrature (int n, double n_slab, double *x, double *w);
