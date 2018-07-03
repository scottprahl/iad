@*1 IAD Types.
This file has no routines.  It is responsible for creating
the header file \.{iad\_type.h} and nothing else.
Altered 3/3/95 to change the version number below.
Change June 95 to improve cross referencing using CTwill.
Change August 97 to add root finding with known absorption

@ 
These are the various optical properties that can
be found with this program.  |FIND_AUTO| allows one
to let the computer figure out what it should be looking
for.

These determine what metric is used in the minimization
process. 

These give the two different types of illumination
allowed.  

Finally, for convenience I create a Boolean type. 

@(iad_type.h@>=
#undef FALSE
#undef TRUE

@h

@<Structs to export from IAD Types@>@;

@
@d FIND_A 0
@d FIND_B 1
@d FIND_AB 2
@d FIND_AG 3
@d FIND_AUTO 4
@d FIND_BG 5
@d FIND_BaG 6
@d FIND_BsG 7
@d FIND_Ba 8
@d FIND_Bs 9
@d FIND_G 10
@d FIND_B_WITH_NO_ABSORPTION 11
@d FIND_B_WITH_NO_SCATTERING 12

@d RELATIVE 0
@d ABSOLUTE 1

@d COLLIMATED 0
@d DIFFUSE 1

@d FALSE 0
@d TRUE 1

@d IAD_MAX_ITERATIONS          500

@ Need error codes for this silly program

@d IAD_NO_ERROR                0

@d IAD_TOO_MANY_ITERATIONS     1

@d IAD_AS_NOT_VALID            16
@d IAD_AE_NOT_VALID            17
@d IAD_AD_NOT_VALID            18
@d IAD_RW_NOT_VALID            19
@d IAD_RD_NOT_VALID            20
@d IAD_RSTD_NOT_VALID          21

@d IAD_GAMMA_NOT_VALID         22
@d IAD_F_NOT_VALID             23
@d IAD_BAD_PHASE_FUNCTION      24

@d IAD_QUAD_PTS_NOT_VALID      25
@d IAD_BAD_G_VALUE             26
@d IAD_TOO_MANY_LAYERS         27

@d IAD_MEMORY_ERROR            28
@d IAD_FILE_ERROR              29

@d IAD_EXCESSIVE_LIGHT_LOSS    30
@d IAD_RT_LT_MINIMUM           31

@d IAD_MR_TOO_SMALL            32
@d IAD_MR_TOO_BIG              33
@d IAD_MT_TOO_SMALL            34
@d IAD_MT_TOO_BIG              35
@d IAD_MU_TOO_SMALL            36
@d IAD_MU_TOO_BIG              37
@d IAD_TOO_MUCH_LIGHT          38
@d IAD_TSTD_NOT_VALID          39

@d UNINITIALIZED              -99 

@d DEBUG_A_LITTLE            1
@d DEBUG_GRID                2
@d DEBUG_ITERATIONS          4
@d DEBUG_LOST_LIGHT          8
@d DEBUG_SPHERE_EFFECTS     16
@d DEBUG_BEST_GUESS         32
@d DEBUG_EVERY_CALC         64
@d DEBUG_SEARCH            128
@d DEBUG_RD_ONLY           256
@d DEBUG_GRID_CALC         512
@d DEBUG_ANY               0xFFFFFFFF

@d UNKNOWN                   0
@d COMPARISON                1
@d SUBSTITUTION                2

@ The idea of the structure |measure_type| is collect
all the information regarding a single measurement together
in one spot.  No information regarding how the inversion
procedure is supposed to be done is contained in this 
structure, unlike in previous incarnations of this program.

@<Structs to export from IAD Types@>=
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

  double d_beam;
  double fraction_of_rc_in_mr;
  double fraction_of_tc_in_mt;

  double m_r, m_t, m_u;

  double lambda;
  double as_r, ad_r, ae_r, aw_r, rd_r, rw_r, rstd_r, f_r;
  double as_t, ad_t, ae_t, aw_t, rd_t, rw_t, rstd_t, f_t;
  double ur1_lost, uru_lost, ut1_lost, utu_lost;
  double d_sphere_r, d_sphere_t;
} IAD_measure_type;

@ This describes how the inversion process should proceed
and also contains the results of that inversion process.

@<Structs to export from IAD Types@>=
typedef struct invert_type {
  double a;		/* the calculated albedo */
  double b;		/* the calculated optical depth */
  double g;		/* the calculated anisotropy */
  
  int found;
  int search;
  int metric;
  double tolerance;
  double MC_tolerance;
  double final_distance;
  int iterations;
  int error;
  
  struct AD_slab_type slab;
  struct AD_method_type method;

  double default_a;
  double default_b;
  double default_g;
  double default_ba;
  double default_bs;
  double default_mua;
  double default_mus;
  
  } IAD_invert_type;

@ A few types that used to be enum's are now int's.

@<Structs to export from IAD Types@>=

typedef int search_type;

typedef int boolean_type;

typedef int illumination_type;

typedef struct guess_t {
	double distance;
	double a;
	double b;
	double g;
} guess_type;

extern double FRACTION;
