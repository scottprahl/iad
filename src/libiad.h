

void 
ez_RT(int n, double nslab,
      double ntopslide,
      double nbottomslide,
      double a,
      double b,
      double g,
      double *UR1, double *UT1, double *URU, double *UTU);

void 
ez_RT_unscattered(int n,
		  double nslab,
		  double ntopslide,
		  double nbottomslide,
		  double a,
		  double b,
		  double g,
		  double *UR1, double *UT1, double *URU, double *UTU);



void 
ez_Inverse_RT(double n, double nslide, double UR1, double UT1, double Tc,
	      double *a, double *b, double *g, int *error);

void 
Spheres_Inverse_RT(double *setup,
		   double *analysis,
		   double *sphere_r,
		   double *sphere_t,
		   double *measurements,
		   double *results);





void 
ez_RT_Cone(int n,
	   double nslab,
	   double ntopslide,
	   double nbottomslide,
	   double a,
	   double b,
	   double g,
	   double cos_cone_angle,
	   double *UR1, double *UT1, double *URU, double *UTU);

void 
ez_RT_Oblique(int n,
	      double nslab,
	      double ntopslide,
	      double nbottomslide,
	      double a,
	      double b,
	      double g,
	      double cos_oblique_angle,
	      double *URx, double *UTx, double *URU, double *UTU);





void 
RT_Layers(int n,
	  double nslab,
	  double ntopslide,
	  double nbottomslide,
	  int nlayers,
	  double a[],
	  double b[],
	  double g[],
	  double *UR1, double *UT1, double *URU, double *UTU);

	void		RT_Layers_All(int n,
			    		double	nslab ,
			    		double	ntopslide,
			    		double	nbottomslide,
			    		int		nlayers ,
			    		double	a      [],
			    		double	b      [],
			    		double	g      [],
   		double       *dUR1, double *dUT1, double *dURU, double *dUTU,
			    		double       *uUR1, double *uUT1, double *uURU, double *uUTU);
