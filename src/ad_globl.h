/*2:*/
#line 29 "ad_globl.w"

#define MAX_QUAD_PTS 128
#define DEFAULT_QUAD_PTS 4 \

#define ISOTROPIC 0
#define HENYEY_GREENSTEIN 1 \

#define DIAMOND 0
#define INFINITESIMAL_GENERATOR 1 \

#define MARTIN_HAMMER 1 \

#define CONE 1
#define OBLIQUE 0 \


#line 30 "ad_globl.w"

/*9:*/
#line 108 "ad_globl.w"


typedef struct AD_slab_type{
double a;
double b;
double g;
int phase_function;
double n_slab;
double n_top_slide;
double n_bottom_slide;
double b_top_slide;
double b_bottom_slide;
double cos_angle;
}slab_type;

/*:9*//*10:*/
#line 123 "ad_globl.w"

typedef struct AD_method_type{
int quad_pts;
double a_calc,b_calc,g_calc,b_thinnest;
}method_type;

/*:10*/
#line 31 "ad_globl.w"

/*12:*/
#line 140 "ad_globl.w"

#ifndef AD_GLOBAL_SOURCE
extern double angle[MAX_QUAD_PTS+1];
extern double weight[MAX_QUAD_PTS+1];
extern double twoaw[MAX_QUAD_PTS+1];
extern int Martin_Hammer;

#endif

/*:12*/
#line 32 "ad_globl.w"

/*15:*/
#line 163 "ad_globl.w"

void Zero_Layer(int n,double**r,double**t)

/*:15*/
#line 33 "ad_globl.w"
;
/*13:*/
#line 151 "ad_globl.w"

void AD_error(char error_text[])

/*:13*/
#line 34 "ad_globl.w"
;
/*21:*/
#line 299 "ad_globl.w"

void URU_and_UR1(int n,double n_slab,double**R,double*URU,double*UR1)

/*:21*/
#line 35 "ad_globl.w"
;
/*17:*/
#line 196 "ad_globl.w"

void URU_and_UR1_Cone(int n,double n_slab,double mu,double**R,double*URU,double*UR1)

/*:17*/
#line 36 "ad_globl.w"
;
/*19:*/
#line 242 "ad_globl.w"

void URU_and_URx_Cone(int n,double n_slab,double mu,double**R,double*URU,double*URx)

/*:19*/
#line 37 "ad_globl.w"
;
/*23:*/
#line 309 "ad_globl.w"

void UFU_and_UF1(int n,double n_slab,
double**Lup,double**Ldown,double*UFU,double*UF1)

/*:23*/
#line 38 "ad_globl.w"
;
/*25:*/
#line 331 "ad_globl.w"

void wrmatrix(int n,double**a)

/*:25*/
#line 39 "ad_globl.w"
;
/*27:*/
#line 375 "ad_globl.w"

void wrarray(int n,double*a)

/*:27*/
#line 40 "ad_globl.w"
;

/*:2*/
