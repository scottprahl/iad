



void RT_Cone (int n,
	      struct AD_slab_type *slab,
	      int use_cone,
	      double *UR1, double *UT1, double *URU, double *UTU);

void ez_RT_Cone (int n,
		 double nslab,
		 double ntopslide,
		 double nbottomslide,
		 double a,
		 double b,
		 double g,
		 double cos_cone_angle,
		 double *UR1, double *UT1, double *URU, double *UTU);

void ez_RT_Oblique (int n,
		    double nslab,
		    double ntopslide,
		    double nbottomslide,
		    double a,
		    double b,
		    double g,
		    double cos_oblique_angle,
		    double *URx, double *UTx, double *URU, double *UTU);
