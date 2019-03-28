

void Inverse_RT (struct measure_type m, struct invert_type *r);

int measure_OK (struct measure_type m, struct invert_type r);

search_type determine_search (struct measure_type m, struct invert_type r);

void Initialize_Result (struct measure_type m, struct invert_type *r);

void ez_Inverse_RT (double n, double nslide, double UR1, double UT1,
		    double Tc, double *a, double *b, double *g, int *error);

void Initialize_Measure (struct measure_type *m);

void Calculate_MR_MT (struct measure_type m,
		      struct invert_type r,
		      int include_MC, double *M_R, double *M_T);

int MinMax_MR_MT (struct measure_type m, struct invert_type r);

void Calculate_Minimum_MR (struct measure_type m,
			   struct invert_type r, double *mr, double *mt);

void Spheres_Inverse_RT2 (double *sample,
			  double *illumination,
			  double *sphere_r,
			  double *sphere_t,
			  double *analysis,
			  double *measurement,
			  double *a, double *b, double *g);
