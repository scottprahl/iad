/*2:*/
#line 27 "ad_prime.w"

#define MAX_FLUENCE_INTERVALS 200 \
 \


#line 28 "ad_prime.w"

/*4:*/
#line 47 "ad_prime.w"

void RT_Matrices(int n,struct AD_slab_type*slab,struct AD_method_type*method,
double**R,double**T)

/*:4*/
#line 29 "ad_prime.w"
;
/*6:*/
#line 105 "ad_prime.w"

void RT(int n,struct AD_slab_type*slab,double*UR1,double*UT1,double*URU,double*UTU)

/*:6*/
#line 30 "ad_prime.w"
;
/*21:*/
#line 262 "ad_prime.w"

void ez_RT(int n,double nslab,
double ntopslide,
double nbottomslide,
double a,
double b,
double g,
double*UR1,double*UT1,double*URU,double*UTU)

/*:21*/
#line 31 "ad_prime.w"
;
/*25:*/
#line 336 "ad_prime.w"

void RTabs(int n,struct AD_slab_type*slab,double*UR1,double*UT1,double*URU,double*UTU)

/*:25*/
#line 32 "ad_prime.w"
;
/*35:*/
#line 443 "ad_prime.w"

void Flux_Fluence(int n,struct AD_slab_type*slab,double zmin,double zmax,int intervals,
double*UF1_array,double*UFU_array,double*flux_up,double*flux_down)

/*:35*/
#line 33 "ad_prime.w"
;
/*23:*/
#line 296 "ad_prime.w"

void ez_RT_unscattered(int n,
double nslab,
double ntopslide,
double nbottomslide,
double a,
double b,
double g,
double*UR1,double*UT1,double*URU,double*UTU)

/*:23*/
#line 34 "ad_prime.w"
;

/*:2*/
