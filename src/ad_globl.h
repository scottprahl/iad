
#define MAX_QUAD_PTS 128
#define DEFAULT_QUAD_PTS 4 \

#define ISOTROPIC 0
#define HENYEY_GREENSTEIN 1 \

#define DIAMOND 0
#define INFINITESIMAL_GENERATOR 1 \

#define MARTIN_HAMMER 1 \

#define CONE 1
#define OBLIQUE 0 \





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


typedef struct AD_method_type{
int quad_pts;
double a_calc,b_calc,g_calc,b_thinnest;
}method_type;



#ifndef AD_GLOBAL_SOURCE
extern double angle[MAX_QUAD_PTS+1];
extern double weight[MAX_QUAD_PTS+1];
extern double twoaw[MAX_QUAD_PTS+1];
extern int Martin_Hammer;

#endif



void Zero_Layer(int n,double**r,double**t)

;

void AD_error(char error_text[])

;

void URU_and_UR1(int n,double n_slab,double**R,double*URU,double*UR1)

;

void URU_and_UR1_Cone(int n,double n_slab,double mu,double**R,double*URU,double*UR1)

;

void URU_and_URx_Cone(int n,double n_slab,double mu,double**R,double*URU,double*URx)

;

void UFU_and_UF1(int n,double n_slab,
double**Lup,double**Ldown,double*UFU,double*UF1)

;

void wrmatrix(int n,double**a)

;

void wrarray(int n,double*a)

;

