/*2:*/
#line 57 "mc_lost.w"

/*13:*/
#line 142 "mc_lost.w"

void MC_Set_Seed(unsigned long seed)

/*:13*/
#line 58 "mc_lost.w"
;
/*46:*/
#line 718 "mc_lost.w"

void MC_Lost(struct measure_type m,struct invert_type r,long n_photons,
double*ur1,double*ut1,double*uru,double*utu,
double*ur1_lost,double*ut1_lost,double*uru_lost,double*utu_lost)

/*:46*/
#line 59 "mc_lost.w"
;
/*52:*/
#line 855 "mc_lost.w"

void MC_RT(struct AD_slab_type s,long n_photons,double t_sample,
double t_top_slide,double t_bottom_slide,
double*UR1,double*UT1,double*URU,double*UTU)

/*:52*/
#line 60 "mc_lost.w"
;
/*41:*/
#line 507 "mc_lost.w"

void MC_Radial(long photons,double a,double b,double g,double n_sample,
double n_top_slide,double n_bottom_slide,
double cos_cone_angle,double cos_incidence,
double t_sample,double t_top_slide,double t_bottom_slide,
double b_top_slide,double b_bottom_slide,
double dr_port,double dt_port,double d_beam,double*r_total,double*t_total,double*r_lost,double*t_lost)

/*:41*/
#line 61 "mc_lost.w"
;
/*54:*/
#line 880 "mc_lost.w"

void MC_Print_RT_Arrays(int status)

/*:54*/
#line 62 "mc_lost.w"
;

/*:2*/
