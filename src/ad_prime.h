
#define MAX_FLUENCE_INTERVALS 200 \
 \




void RT_Matrices (int n, struct AD_slab_type *slab,
		  struct AD_method_type *method, double **R, double **T);

void RT (int n, struct AD_slab_type *slab, double *UR1, double *UT1,
	 double *URU, double *UTU);

void ez_RT (int n, double nslab,
	    double ntopslide,
	    double nbottomslide,
	    double a,
	    double b,
	    double g, double *UR1, double *UT1, double *URU, double *UTU);

void RTabs (int n, struct AD_slab_type *slab, double *UR1, double *UT1,
	    double *URU, double *UTU);

void Flux_Fluence (int n, struct AD_slab_type *slab, double zmin, double zmax,
		   int intervals, double *UF1_array, double *UFU_array,
		   double *flux_up, double *flux_down);

void ez_RT_unscattered (int n,
			double nslab,
			double ntopslide,
			double nbottomslide,
			double a,
			double b,
			double g,
			double *UR1, double *UT1, double *URU, double *UTU);
