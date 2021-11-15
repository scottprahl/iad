



void RT_Layers (int n,
		double nslab,
		double ntopslide,
		double nbottomslide,
		int nlayers,
		double a[],
		double b[],
		double g[],
		double *UR1, double *UT1, double *URU, double *UTU);

void RT_Layers_All (int n,
		    double nslab,
		    double ntopslide,
		    double nbottomslide,
		    int nlayers,
		    double a[],
		    double b[],
		    double g[],
		    double *dUR1, double *dUT1, double *dURU, double *dUTU,
		    double *uUR1, double *uUT1, double *uURU, double *uUTU);
