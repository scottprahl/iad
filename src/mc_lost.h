void MC_Lost(struct measure_type m, struct invert_type r,  long n_photons,
             double *ur1,      double *ut1,      double *uru,      double *utu, 
             double *ur1_lost, double *ut1_lost, double *uru_lost, double *utu_lost);

void MC_RT(struct AD_slab_type s, long n_photons, double *UR1, double *UT1, double *URU, double *UTU);
