/* Autogenerated v3-16-3 from https://github.com/scottprahl/iad */

#define FIND_A 0
#define FIND_B 1
#define FIND_AB 2
#define FIND_AG 3
#define FIND_AUTO 4
#define FIND_BG 5
#define FIND_BaG 6
#define FIND_BsG 7
#define FIND_Ba 8
#define FIND_Bs 9
#define FIND_G 10
#define FIND_B_WITH_NO_ABSORPTION 11
#define FIND_B_WITH_NO_SCATTERING 12 \

#define RELATIVE 0
#define ABSOLUTE 1 \

#define COLLIMATED 1.0
#define DIFFUSE 0.0 \

#define FALSE 0
#define TRUE 1 \

#define IAD_MAX_ITERATIONS 200
#define BIG_A_CALC_VALUE 999999.0
#define SMALL_A_CALC_VALUE 0.00001
#define MAX_ABS_G 0.999999
#define GRID_SIZE 101 \

#define IAD_NO_ERROR 0 \

#define IAD_TOO_MANY_ITERATIONS 1 \

#define IAD_AS_NOT_VALID 16
#define IAD_AE_NOT_VALID 17
#define IAD_AD_NOT_VALID 18
#define IAD_RW_NOT_VALID 19
#define IAD_RD_NOT_VALID 20
#define IAD_RSTD_NOT_VALID 21 \

#define IAD_GAMMA_NOT_VALID 22
#define IAD_F_NOT_VALID 23
#define IAD_BAD_PHASE_FUNCTION 24 \

#define IAD_QUAD_PTS_NOT_VALID 25
#define IAD_BAD_G_VALUE 26
#define IAD_TOO_MANY_LAYERS 27 \

#define IAD_MEMORY_ERROR 28
#define IAD_FILE_ERROR 29 \

#define IAD_EXCESSIVE_LIGHT_LOSS 30
#define IAD_RT_LT_MINIMUM 31 \

#define IAD_MR_TOO_SMALL 32
#define IAD_MR_TOO_BIG 33
#define IAD_MT_TOO_SMALL 34
#define IAD_MT_TOO_BIG 35
#define IAD_MU_TOO_SMALL 36
#define IAD_MU_TOO_BIG 37
#define IAD_TOO_MUCH_LIGHT 38
#define IAD_TSTD_NOT_VALID 39 \

#define UNINITIALIZED -99 \

#define DEBUG_A_LITTLE 1
#define DEBUG_GRID 2
#define DEBUG_ITERATIONS 4
#define DEBUG_LOST_LIGHT 8
#define DEBUG_BEST_GUESS 16
#define DEBUG_SEARCH 32
#define DEBUG_GRID_CALC 64
#define DEBUG_SPHERE_GAIN 128
#define DEBUG_EVERY_CALC 256
#define DEBUG_MC 512
#define DEBUG_ANY 0xFFFFFFFF \

#define UNKNOWN 0
#define COMPARISON 1
#define SUBSTITUTION 2 \

#define MC_NONE 0
#define MC_USE_EXISTING 1
#define MC_REDO 2 \
 \

typedef struct measure_type {
    double slab_index;
    double slab_thickness;

    double slab_top_slide_index;
    double slab_top_slide_b;
    double slab_top_slide_thickness;

    double slab_bottom_slide_index;
    double slab_bottom_slide_b;
    double slab_bottom_slide_thickness;
    double slab_cos_angle;

    int num_spheres;
    int num_measures;
    int method;
    int flip_sample;
    int baffle_r, baffle_t;

    double d_beam;
    double fraction_of_ru_in_mr;
    double fraction_of_tu_in_mt;

    double m_r, m_t, m_u;

    double lambda;
    double as_r, ad_r, at_r, aw_r, rd_r, rw_r, rstd_r, f_r;
    double as_t, ad_t, at_t, aw_t, rd_t, rw_t, rstd_t;
    double ur1_lost, uru_lost, ut1_lost, utu_lost;
    double d_sphere_r, d_sphere_t;
} IAD_measure_type;

typedef struct invert_type {
    double a;
    double b;
    double g;

    int found;
    int search;
    int metric;
    double tolerance;
    double MC_tolerance;
    double final_distance;
    int error;

    struct AD_slab_type slab;
    struct AD_method_type method;

    int AD_iterations;
    int MC_iterations;

    double default_a;
    double default_b;
    double default_g;
    double default_ba;
    double default_bs;
    double default_mua;
    double default_mus;

} IAD_invert_type;

typedef int search_type;

typedef int boolean_type;

typedef int illumination_type;

typedef struct guess_t {
    double distance;
    double a;
    double b;
    double g;
    double ur1_lost;
    double ut1_lost;
    double uru_lost;
    double utu_lost;
} guess_type;

extern double FRACTION;
