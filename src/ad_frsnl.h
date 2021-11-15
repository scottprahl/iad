/*2:*/
#line 25 "ad_frsnl.w"

/*3:*/
#line 57 "ad_frsnl.w"

double Cos_Critical_Angle(double ni,double nt)

/*:3*/
#line 26 "ad_frsnl.w"
;
/*5:*/
#line 110 "ad_frsnl.w"

double Cos_Snell(double n_i,double mu_i,double n_t)

/*:5*/
#line 27 "ad_frsnl.w"
;
/*12:*/
#line 312 "ad_frsnl.w"

void Absorbing_Glass_RT(double n_i,double n_g,double n_t,double mu_i,double b,
double*r,double*t)

/*:12*/
#line 28 "ad_frsnl.w"
;
/*17:*/
#line 373 "ad_frsnl.w"

void Sp_mu_RT(double n_top,double n_slab,double n_bottom,
double tau_top,double tau_slab,double tau_bottom,double mu,
double*r,double*t)

/*:17*/
#line 29 "ad_frsnl.w"
;
/*15:*/
#line 349 "ad_frsnl.w"

void Sp_mu_RT_Flip(int flip,double n_top,double n_slab,double n_bottom,
double tau_top,double tau_slab,double tau_bottom,double mu,
double*r,double*t)

/*:15*/
#line 30 "ad_frsnl.w"
;
/*24:*/
#line 514 "ad_frsnl.w"

double Diffuse_Glass_R(double nair,double nslide,double nslab)

/*:24*/
#line 31 "ad_frsnl.w"
;
/*10:*/
#line 250 "ad_frsnl.w"

double Glass(double n_i,double n_g,double n_t,double mu_i)

/*:10*/
#line 32 "ad_frsnl.w"
;

/*:2*/
