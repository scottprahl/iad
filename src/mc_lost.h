void MC_Lost(struct measure_type m, struct invert_type r,  long n_photons,
             double *ur1,      double *ut1,      double *uru,      double *utu, 
             double *ur1_lost, double *ut1_lost, double *uru_lost, double *utu_lost);

void MC_RT(struct AD_slab_type s, long n_photons, double t_sample, double t_slide,
           double *UR1, double *UT1, double *URU, double *UTU);

void MC_Radial(long photons, double a, double b, double g, double n_sample,
               double n_slide, int collimated, double cos_incidence,
               double t_sample, double t_slide, double mua_slide, double dr_port, double dt_port,
               double d_beam, double *r_total, double *t_total, double *r_lost, double *t_lost);
