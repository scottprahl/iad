@** IAD Calculation.

\def\rdirect{r_s^{\hbox{\sevenrm{} direct}}}
\def\tdirect{t_s^{\hbox{\sevenrm{} direct}}}
\def\rdiffuse{r_s^{\hbox{\sevenrm{}}}}
\def\tdiffuse{t_s^{\hbox{\sevenrm{}}}}
\def\std{{\hbox{\sevenrm{}std}}}
\def\third{{\hbox{\sevenrm{}third}}}
\def\first{{\hbox{\sevenrm{}first}}}
\def\lost{{\hbox{\sevenrm{}lost}}}
\def\scat{{\hbox{\sevenrm{}scat}}}
\def\unscat{{\hbox{\sevenrm{}unscat}}}
\def\miss{{\hbox{\sevenrm{}miss}}}


@(iad_calc.c@>=

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "nr_util.h"
#include "nr_zbrent.h"
#include "ad_globl.h"
#include "ad_frsnl.h"
#include "ad_prime.h"
#include "iad_type.h"
#include "iad_util.h"
#include "iad_calc.h"

#define ABIT 1e-6
#define A_COLUMN 1
#define B_COLUMN 2
#define G_COLUMN 3
#define URU_COLUMN 4
#define UTU_COLUMN 5
#define UR1_COLUMN 6
#define UT1_COLUMN 7
#define REFLECTION_SPHERE 1
#define TRANSMISSION_SPHERE 0
#define T_TRUST_FACTOR 1

#define SWAP(a,b) {double swap=(a);(a)=(b);(b)=swap;}

static int CALCULATING_GRID = 0;
static struct measure_type MM;
static struct invert_type RR;
static struct measure_type MGRID;
static struct invert_type RGRID;
static double ** The_Grid=NULL;
static double GG_a;
static double GG_b;
static double GG_g;
static double GG_bs;
static double GG_ba;
static boolean_type The_Grid_Initialized = FALSE;
static boolean_type The_Grid_Search = -1;

    @<Definition for |Set_Calc_State|@>@;
    @<Definition for |Get_Calc_State|@>@;
    @<Definition for |Same_Calc_State|@>@;

    @<Prototype for |Fill_AB_Grid|@>;
    @<Prototype for |Fill_AG_Grid|@>;

    @<Definition for |RT_Flip|@>@;
    @<Definition for |Allocate_Grid|@>@;
    @<Definition for |Valid_Grid|@>@;
    @<Definition for |fill_grid_entry|@>@;
    @<Definition for |Fill_Grid|@>@;
    @<Definition for |Near_Grid_Points|@>@;
    @<Definition for |Fill_AB_Grid|@>@;
    @<Definition for |Fill_AG_Grid|@>@;
    @<Definition for |Fill_BG_Grid|@>@;
    @<Definition for |Fill_BaG_Grid|@>@;
    @<Definition for |Fill_BsG_Grid|@>@;
    @<Definition for |Grid_ABG|@>@;
    @<Definition for |Gain|@>@;
    @<Definition for |Gain_11|@>@;
    @<Definition for |Gain_22|@>@;
    @<Definition for |Two_Sphere_R|@>@;
    @<Definition for |Two_Sphere_T|@>@;

    @<Definition for |Calculate_Distance_With_Corrections|@>@;
    @<Definition for |Calculate_Grid_Distance|@>@;
    @<Definition for |Calculate_Distance|@>@;
    @<Definition for |abg_distance|@>@;

    @<Definition for |Find_AG_fn|@>@;
    @<Definition for |Find_AB_fn|@>@;
    @<Definition for |Find_Ba_fn|@>@;
    @<Definition for |Find_Bs_fn|@>@;
    @<Definition for |Find_A_fn|@>@;
    @<Definition for |Find_B_fn|@>@;
    @<Definition for |Find_G_fn|@>@;
    @<Definition for |Find_BG_fn|@>@;
    @<Definition for |Find_BaG_fn|@>@;
    @<Definition for |Find_BsG_fn|@>@;
    @<Definition for |maxloss|@>@;
    @<Definition for |Max_Light_Loss|@>@;

@

@(iad_calc.h@>=
    @<Prototype for |Gain|@>;
    @<Prototype for |Gain_11|@>;
    @<Prototype for |Gain_22|@>;
    @<Prototype for |Two_Sphere_R|@>;
    @<Prototype for |Two_Sphere_T|@>;
    @<Prototype for |Set_Calc_State|@>;
    @<Prototype for |Get_Calc_State|@>;
    @<Prototype for |Same_Calc_State|@>;
    @<Prototype for |Valid_Grid|@>;
    @<Prototype for |Allocate_Grid|@>;
    @<Prototype for |Fill_Grid|@>;
    @<Prototype for |Near_Grid_Points|@>;
    @<Prototype for |Grid_ABG|@>;
    @<Prototype for |Find_AG_fn|@>;
    @<Prototype for |Find_AB_fn|@>;
    @<Prototype for |Find_Ba_fn|@>;
    @<Prototype for |Find_Bs_fn|@>;
    @<Prototype for |Find_A_fn|@>;
    @<Prototype for |Find_B_fn|@>;
    @<Prototype for |Find_G_fn|@>;
    @<Prototype for |Find_BG_fn|@>;
    @<Prototype for |Find_BsG_fn|@>;
    @<Prototype for |Find_BaG_fn|@>;
    @<Prototype for |Fill_BG_Grid|@>;
    @<Prototype for |Fill_BsG_Grid|@>;
    @<Prototype for |Fill_BaG_Grid|@>;
    @<Prototype for |Calculate_Distance_With_Corrections|@>;
    @<Prototype for |Calculate_Distance|@>;
    @<Prototype for |Calculate_Grid_Distance|@>;
    @<Prototype for |abg_distance|@>;
    @<Prototype for |maxloss|@>;
    @<Prototype for |Max_Light_Loss|@>;
    @<Prototype for |RT_Flip|@>;


@*1 Initialization.

The functions in this file assume that the local variables |MM| and |RR|
have been initialized appropriately.  The variable |MM| contains all the
information about how a particular experiment was done.  The structure
|RR| contains the data structure that is passed to the adding-doubling routines
as well as the number of quadrature points.

@*1 Gain.

Assume that a sphere is illuminated with diffuse light having a power
$P$. This light will undergo multiple reflections in the sphere walls that
will increase the power falling on the detector.
The gain on the detector due to integrating sphere effects varies with
the presence of a baffle between the sample and the detector.  If a baffle is
present then
$$
G_{\rm no\ baffle}(\rdiffuse, r_t) = {1 \over 1 - a_w r_w - a_d r_d - a_s \rdiffuse - a_s r_t}
$$
or with a baffle as
$$
G_{\rm baffle}(\rdiffuse, r_t) = {1 \over 1- a_w r_w' - r_w' (1-a_t) (a_d r_d + a_s \rdiffuse)}
$$
where
$$
r_w' = r_w + (a_t/a_w) r_t
$$
For a black sphere the gain does not depend on the diffuse reflectivity of the sample
and is unity.  $G(\rdiffuse) = 1$, which is easily verified by setting $r_w=0$.

The value |uru_sample| is the total reflectance from the sample for diffuse incident light
and |uru_third| is the total reflectance from the third port for diffuse incident light.
For a reflection sphere, the third port is the entrance port and |uru_third=0|.

@ @<Prototype for |Gain|@>=
double Gain(int sphere, struct measure_type m, double uru_sample, double uru_third)

@ @<Definition for |Gain|@>=
    @<Prototype for |Gain|@>
{
double inv_gain;

if (sphere == REFLECTION_SPHERE) {
    if (m.baffle_r) {
        inv_gain = m.rw_r + (m.at_r / m.aw_r) * uru_third;
        inv_gain *= m.aw_r + (1 - m.at_r) * (m.ad_r * m.rd_r + m.as_r * uru_sample);
        inv_gain = 1.0 - inv_gain;
    } else {
        inv_gain = 1.0 - m.aw_r * m.rw_r - m.ad_r * m.rd_r - m.as_r * uru_sample - m.at_r * uru_third;
    }
} else
    if (m.baffle_t) {
        inv_gain = m.rw_t + (m.at_t / m.aw_t) * uru_third;
        inv_gain *= m.aw_t + (1 - m.at_t) * (m.ad_t * m.rd_t + m.as_t * uru_sample);
        inv_gain = 1.0 - inv_gain;
    } else {
        inv_gain = 1.0 - m.aw_t * m.rw_t - m.ad_t * m.rd_t - m.as_t * uru_sample - m.at_t * uru_third;
    }
return 1.0 / inv_gain;
}

@ The gain for light on the detector in the first sphere for diffuse light starting
in that same sphere is defined as
$$
G_{1\rightarrow1}(r_s,t_s) \equiv
{P_{1\rightarrow1}(r_s,t_s)/A_d\over P/A}
$$
then the full expression for the gain is
$$
G_{1\rightarrow1}(r_s,t_s) =
{G(r_s) \over 1-a_s a_s' r_w r_w' (1-a_t)(1-a_t') G(r_s) G'(r_s)t_s^2  }
$$


@<Prototype for |Gain_11|@>=
double Gain_11(struct measure_type m, double URU, double tdiffuse)

@ @<Definition for |Gain_11|@>=
    @<Prototype for |Gain_11|@>
{
    double G, GP, G11;

    G        = Gain(REFLECTION_SPHERE, m, URU, 0);
    GP       = Gain(TRANSMISSION_SPHERE, m, URU, 0);

    G11 = G / (1-m.as_r * m.as_t * m.aw_r * m.aw_t * (1-m.at_r) * (1-m.at_t)
               * G * GP * tdiffuse * tdiffuse);

    return G11;
}

@ Similarly, when the light starts in the second sphere, the gain for light
on the detector in the second sphere $G_{2\rightarrow2}$ is found by switching
all primed variables to unprimed.  Thus $G_{2\rightarrow1}(r_s,t_s)$ is
$$
G_{2\rightarrow2}(r_s,t_s) = {G'(r_s) \over 1-a_s a_s' r_w r_w'
                              (1-a_t)(1-a_t') G(r_s) G'(r_s)t_s^2  }
$$

@<Prototype for |Gain_22|@>=
double Gain_22(struct measure_type m, double URU, double tdiffuse)

@ @<Definition for |Gain_22|@>=
    @<Prototype for |Gain_22|@>
{
    double G, GP, G22;

    G        = Gain(REFLECTION_SPHERE, m, URU, 0);
    GP       = Gain(TRANSMISSION_SPHERE, m, URU, 0);

    G22 = GP / (1-m.as_r * m.as_t * m.aw_r * m.aw_t * (1-m.at_r) * (1-m.at_t)
               * G * GP * tdiffuse * tdiffuse);

    return G22;
}

@ The reflected power for two spheres makes use of the formulas for
|Gain_11| above.

The light
on the detector in the reflection (first) sphere arises from three
sources: the fraction of light directly reflected off the sphere wall $f
r_w^2 (1-a_t) P$, the fraction of light reflected by the sample $(1-f)
\rdirect r_w^2 (1-a_t) P$, and the light transmitted through the sample
$(1-f) \tdirect r_w' (1-a_t') P$,
$$
\eqalign{
R(\rdirect,\rdiffuse,\tdirect,\tdiffuse)
&= G_{1\rightarrow1}(\rdiffuse,\tdiffuse) \cdot a_d (1-a_t) r_w^2 f  P \cr
&+ G_{1\rightarrow1}(\rdiffuse,\tdiffuse) \cdot a_d (1-a_t) r_w (1-f) \rdirect  P \cr
&+ G_{2\rightarrow1}(\rdiffuse,\tdiffuse) \cdot a_d (1-a_t') r_w' (1-f) \tdirect  P \cr
}
$$
which simplifies slightly to
$$
\eqalign{
R(\rdirect,\rdiffuse,\tdirect,\tdiffuse)
&= a_d (1-a_t) r_w P \cdot G_{1\rightarrow1}(\rdiffuse,\tdiffuse) \cr
&\times \bigg[(1-f)\rdirect + f r_w +
(1-f)a_s'(1-a_t')r_w'\tdirect\tdiffuse G'(\rdiffuse)\bigg] \cr
}
$$

@<Prototype for |Two_Sphere_R|@>=
double Two_Sphere_R(struct measure_type m,
                    double UR1, double URU, double UT1, double UTU)

@ @<Definition for |Two_Sphere_R|@>=
    @<Prototype for |Two_Sphere_R|@>
{
    double x, GP;
    GP = Gain(TRANSMISSION_SPHERE, m, URU, 0);

    x = m.ad_r*(1-m.at_r)*m.rw_r*Gain_11(m,URU,UTU);
    x *= (1-m.f_r)*UR1+m.rw_r*m.f_r+(1-m.f_r)*m.as_t*(1-m.at_t)*m.rw_t*UT1*UTU*GP;
    return x;
}

@ For the power on the detector in the transmission (second) sphere we
have the same three sources.  The only difference is that the subscripts
on the gain terms now indicate that the light ends up in the second
sphere
$$
\eqalign{
T(\rdirect,\rdiffuse,\tdirect,\tdiffuse)
&= G_{1\rightarrow2}(\rdiffuse,\tdiffuse)  \cdot a_d' (1-a_t) r_w^2 f P \cr
&+ G_{1\rightarrow2}(\rdiffuse,\tdiffuse) \cdot a_d' (1-a_t) r_w (1-f) \rdirect P  \cr
&+ G_{2\rightarrow2}(\rdiffuse,\tdiffuse) \cdot a_d' (1-a_t') r_w' (1-f) \tdirect  P \cr
}
$$
or
$$
\eqalign{
T(\rdirect,\rdiffuse,\tdirect,\tdiffuse)
&= a_d' (1-a_t') r_w' P\cdot G_{2\rightarrow2}(\rdiffuse,\tdiffuse) \cr
&\times \bigg[(1-f)\tdirect
 + (1-a_t)r_w a_s \tdiffuse (f r_w+(1-f) \rdirect)G(\rdiffuse) \bigg] \cr
}
$$

@<Prototype for |Two_Sphere_T|@>=
double Two_Sphere_T(struct measure_type m,
                    double UR1, double URU, double UT1, double UTU)

@ @<Definition for |Two_Sphere_T|@>=
    @<Prototype for |Two_Sphere_T|@>
{
    double x, G;
    G = Gain(REFLECTION_SPHERE, m, URU, 0);
    x = m.ad_t*(1-m.at_t)*m.rw_t*Gain_22(m,URU,UTU);
    x *= (1-m.f_r)*UT1+(1-m.at_r)*m.rw_r*m.as_r*UTU*(m.f_r*m.rw_r+(1-m.f_r)*UR1)*G;
    return x;
}

@*1 Grid Routines.
There is a long story associated with these routines.  I spent a lot of time
trying to find an empirical function to allow a guess at a starting value for
the inversion routine.  Basically nothing worked very well.  There were
too many special cases and what not.  So I decided to calculate a whole bunch
of reflection and transmission values and keep their associated optical
properties linked nearby.

I did the very simplest thing.  I just allocate a matrix that is five columns wide.
Then I fill every row with a calculated set of optical properties and observables.
The distribution of values that I use could certainly use some work, but
they currently work.

SO... how does this thing work anyway?   There are two possible grids one for
calculations requiring the program to find the albedo and the optical depth ($a$
and $b$) and one to find the albedo and anisotropy ($a$ and $g$).  These grids
must be allocated and initialized before use.

@  This is a pretty important routine that should have some explanation.  The
reason that it exists, is that we need some `out-of-band' information during
the minimization process.  Since the light transport calculation depends on
all sorts of stuff (e.g., the sphere parameters) and the minimization routines
just vary one or two parameters this information needs to be put somewhere.

I chose the global variables |MM| and |RR| to save things in.

The bottom line is that you cannot do a light transport calculation without
calling this routine first.

@<Prototype for |Set_Calc_State|@>=
void Set_Calc_State(struct measure_type m, struct invert_type r)

@ @<Definition for |Set_Calc_State|@>=
    @<Prototype for |Set_Calc_State|@>
{
    memcpy(&MM, &m, sizeof(struct measure_type));
    memcpy(&RR, &r, sizeof(struct invert_type));
}

@ The inverse of the previous routine.  Note that you must have space for
the parameters |m| and |r| already allocated.

@<Prototype for |Get_Calc_State|@>=
void Get_Calc_State(struct measure_type *m, struct invert_type *r)

@ @<Definition for |Get_Calc_State|@>=
    @<Prototype for |Get_Calc_State|@>
{
    memcpy(m, &MM, sizeof(struct measure_type));
    memcpy(r, &RR, sizeof(struct invert_type));
}

@ The inverse of the previous routine.  Note that you must have space for
the parameters |m| and |r| already allocated.

@<Prototype for |Same_Calc_State|@>=
boolean_type Same_Calc_State(struct measure_type m, struct invert_type r)

@ @<Definition for |Same_Calc_State|@>=
    @<Prototype for |Same_Calc_State|@>
{
    if (The_Grid==NULL)        return FALSE;
    if (!The_Grid_Initialized) return FALSE ;

    if (r.search              != RR.search)              return FALSE;
    if (r.method.quad_pts     != RR.method.quad_pts)     return FALSE;
    if (r.slab.a              != RR.slab.a)              return FALSE;
    if (r.slab.b              != RR.slab.b)              return FALSE;
    if (r.slab.g              != RR.slab.g)              return FALSE;
    if (r.slab.phase_function != RR.slab.phase_function) return FALSE;
    if (r.slab.n_slab         != RR.slab.n_slab)         return FALSE;
    if (r.slab.n_top_slide    != RR.slab.n_top_slide)    return FALSE;
    if (r.slab.n_bottom_slide != RR.slab.n_bottom_slide) return FALSE;
    if (r.slab.b_top_slide    != RR.slab.b_top_slide)    return FALSE;
    if (r.slab.b_bottom_slide != RR.slab.b_bottom_slide) return FALSE;
    if (r.slab.cos_angle      != RR.slab.cos_angle)      return FALSE;
    if ((m.num_measures==3) && (m.m_u!=MGRID.m_u)) return (FALSE);
    return TRUE;
}

@ @<Prototype for |Allocate_Grid|@>=
void Allocate_Grid(search_type s)

@ @<Definition for |Allocate_Grid|@>=
    @<Prototype for |Allocate_Grid|@>
{
    (void) s;
    The_Grid = dmatrix(0,GRID_SIZE*GRID_SIZE,1,7);
    if (The_Grid==NULL) AD_error("unable to allocate the grid matrix");
    The_Grid_Initialized = FALSE;
}

@ This routine will return the |a|, |b|, and |g| values for a particular
row in the grid.

@<Prototype for |Grid_ABG|@>=
void Grid_ABG(int i, int j, guess_type *guess)

@ @<Definition for |Grid_ABG|@>=
    @<Prototype for |Grid_ABG|@>
{
    if (0<=i && i<GRID_SIZE && 0<=j && j<GRID_SIZE) {
        guess->a = The_Grid[GRID_SIZE*i+j][A_COLUMN];
        guess->b = The_Grid[GRID_SIZE*i+j][B_COLUMN];
        guess->g = The_Grid[GRID_SIZE*i+j][G_COLUMN];
        guess->distance = Calculate_Grid_Distance(i,j);
    } else {
        guess->a = 0.5;
        guess->b = 0.5;
        guess->g = 0.5;
        guess->distance = 999;
    }
}

@ This routine is used to figure out if the current grid is valid.
This can fail for several reasons.  First the grid may not have
been allocated.  Or it may not have been initialized.
The boundary conditions may have changed.
The number or values of the sphere parameters may have changed.
It is tedious, but straightforward to check these cases out.

If this routine returns true, then it is a pretty good bet that the values
in the current grid can be used to guess the next starting set of values.

@<Prototype for |Valid_Grid|@>=
boolean_type Valid_Grid(struct measure_type m, struct invert_type r)

@ @<Definition for |Valid_Grid|@>=
    @<Prototype for |Valid_Grid|@>
{
    int s = r.search;
    @<Tests for invalid grid@>@;

    return(TRUE);
}

@ First check are to test if the grid has ever been filled

@<Tests for invalid grid@>=
    if (The_Grid==NULL) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr,"GRID: Fill because NULL\n");
        return(FALSE);
    }
    if (!The_Grid_Initialized) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr,"GRID: Fill because not initialized\n");
        return(FALSE);
    }

@ If the type of search has changed then report the grid as invalid

@<Tests for invalid grid@>=
    if (The_Grid_Search != s) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr,"GRID: Fill because search type changed\n");
        return(FALSE);
    }

@ Compare the |m.m_u| value only if there are three measurements

@<Tests for invalid grid@>=

    if ((m.num_measures==3) && (m.m_u!=MGRID.m_u)) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr,"GRID: Fill because unscattered light changed\n");
        return (FALSE);
    }

@ Make sure that the boundary conditions have not changed.

@<Tests for invalid grid@>=
    if (m.slab_index              != MGRID.slab_index) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr,"GRID: Fill because slab refractive index changed\n");
        return(FALSE);
    }
    if (m.slab_cos_angle          != MGRID.slab_cos_angle) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr,"GRID: Fill because light angle changed\n");
        return(FALSE);
    }

    if (m.slab_top_slide_index    != MGRID.slab_top_slide_index) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr,"GRID: Fill because top slide index changed\n");
        return(FALSE);
    }

    if (m.slab_bottom_slide_index != MGRID.slab_bottom_slide_index) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr,"GRID: Fill because bottom slide index changed\n");
        return(FALSE);
    }

    if (s == FIND_AB && r.slab.g != RGRID.slab.g) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr,"GRID: Fill because anisotropy changed\n");
        return(FALSE);
    }

    if (s == FIND_AG && r.slab.b != RGRID.slab.b) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr,"GRID: Fill because optical depth changed\n");
        return(FALSE);
    }

    if (s == FIND_BsG && r.default_ba != RGRID.default_ba) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr,"GRID: Fill because mu_a changed\n");
        return(FALSE);
    }

    if (s == FIND_BaG && r.default_bs != RGRID.default_bs) {
        if (Debug(DEBUG_GRID))
            fprintf(stderr,"GRID: Fill because mu_s changed\n");
        return(FALSE);
    }


@ Routine to just figure out the distance to a particular a, b, g point

@<Prototype for |abg_distance|@>=
void abg_distance(double a, double b, double g, guess_type *guess)

@ @<Definition for |abg_distance|@>=
    @<Prototype for |abg_distance|@>
{
    double m_r, m_t, distance;
    struct measure_type old_mm;
    struct invert_type old_rr;

    Get_Calc_State(&old_mm, &old_rr);

    RR.slab.a = a;
    RR.slab.b = b;
    RR.slab.g = g;

    Calculate_Distance(&m_r,&m_t,&distance);

    Set_Calc_State(old_mm, old_rr);

    guess->a = a;
    guess->b = b;
    guess->g = g;
    guess->distance = distance;
}

@ This just searches through the grid to find the minimum entry and returns the
optical properties of that entry.  The
smallest, the next smallest, and the third smallest values are returned.

This has been rewritten to use |Calculate_Distance_With_Corrections| so that changes in
sphere parameters won't necessitate recalculating the grid.

@<Prototype for |Near_Grid_Points|@>=
void Near_Grid_Points(double r, double t, search_type s, int *i_min, int *j_min)

@ @<Definition for |Near_Grid_Points|@>=
    @<Prototype for |Near_Grid_Points|@>
{
    int i,j;
    double  fval;
    double smallest=10.0;
    struct measure_type old_mm;
    struct invert_type old_rr;
    (void) r;
    (void) t;
    (void) s;

    if (Debug(DEBUG_GRID))
        fprintf(stderr,"GRID: Finding best grid points\n");
    Get_Calc_State(&old_mm, &old_rr);

    *i_min = 0;
    *j_min = 0;
    for(i=0; i<GRID_SIZE; i++){
        for(j=0; j<GRID_SIZE; j++){

            CALCULATING_GRID = 1;
            fval = Calculate_Grid_Distance(i,j);
            CALCULATING_GRID = 0;

            if (fval<smallest){
                *i_min = i;
                *j_min = j;
                smallest=fval;
            }
        }
    }

    Set_Calc_State(old_mm, old_rr);
}

@ Routine to incorporate flipping of sample if needed.  This is pretty
simple.  The assumption is that flipping is handled relative to the
reflection side of the sphere.  Thus even when flipping is needed,
the usual call to |RT()| will result in the correct values for the
reflectances.  The transmission values can then be calculated by
swapping the top and bottom slides.

Technically, the value of slab should be |const| but it is not so that
we don't pay a copying overhead whenever |flip| is false (the usual case).

@<Prototype for |RT_Flip|@>=
void RT_Flip(int flip, int n, struct AD_slab_type * slab, double *UR1, double *UT1,
             double *URU, double *UTU)

@ @<Definition for |RT_Flip|@>=
    @<Prototype for |RT_Flip|@>
{
    double correct_UR1, correct_URU;
    RT(n, slab, UR1, UT1, URU, UTU);
    if (flip) {
        correct_UR1 = *UR1;
        correct_URU = *URU;
        SWAP(slab->n_top_slide, slab->n_bottom_slide)
        SWAP(slab->b_top_slide, slab->b_bottom_slide)
        RT(n, slab, UR1, UT1, URU, UTU);
        SWAP(slab->n_top_slide, slab->n_bottom_slide)
        SWAP(slab->b_top_slide, slab->b_bottom_slide)
        *UR1 = correct_UR1;
        *URU = correct_URU;
    }
}

@ Simple routine to put values into the grid

Presumes that |RR.slab| is properly set up.

@<Definition for |fill_grid_entry|@>=
static void fill_grid_entry(int i, int j)
{
    double ur1,ut1,uru,utu;

    if (RR.slab.b <= 1e-6 ) RR.slab.b = 1e-6;

    if (Debug(DEBUG_GRID_CALC) && i==0 && j==0) {
        fprintf(stderr, "+   i   j ");
        fprintf(stderr, "      a         b          g     |");
        fprintf(stderr, "     M_R        grid  |");
        fprintf(stderr, "     M_T        grid\n");
    }

    if (Debug(DEBUG_EVERY_CALC) ) {
        if (!CALCULATING_GRID)
            fprintf(stderr, "a=%8.5f b=%10.5f g=%8.5f ", RR.slab.a, RR.slab.b, RR.slab.g);
        else {
            if (j==0) fprintf(stderr, ".");
            if (i+1 == GRID_SIZE && j==0) fprintf(stderr, "\n");
        }
    }

    RT_Flip(MM.flip_sample, RR.method.quad_pts, &RR.slab, &ur1, &ut1, &uru, &utu);

    if (Debug(DEBUG_EVERY_CALC) && !CALCULATING_GRID)
        fprintf(stderr, "ur1=%8.5f ut1=%8.5f\n", ur1, ut1);

    The_Grid[GRID_SIZE*i+j][A_COLUMN]=RR.slab.a;
    The_Grid[GRID_SIZE*i+j][B_COLUMN]=RR.slab.b;
    The_Grid[GRID_SIZE*i+j][G_COLUMN]=RR.slab.g;
    The_Grid[GRID_SIZE*i+j][UR1_COLUMN]=ur1;
    The_Grid[GRID_SIZE*i+j][UT1_COLUMN]=ut1;
    The_Grid[GRID_SIZE*i+j][URU_COLUMN]=uru;
    The_Grid[GRID_SIZE*i+j][UTU_COLUMN]=utu;

    if (Debug(DEBUG_GRID_CALC)) {
        fprintf(stderr, "+ %3d %3d ",i,j);
        fprintf(stderr, "%10.5f %10.5f %10.5f |", RR.slab.a, RR.slab.b, RR.slab.g);
        fprintf(stderr, "%10.5f %10.5f |", MM.m_r, uru);
        fprintf(stderr, "%10.5f %10.5f \n", MM.m_t, utu);
    }
}

@ This routine fills the grid with a proper set of values.  With a little work, this
routine could be made much faster by (1) only generating the phase function
matrix once, (2) Making only one pass through the array for each albedo value,
i.e., using the matrix left over from $b=1$ to generate the solution for $b=2$.
Unfortunately this would require a complete revision of the |Calculate_Distance|
routine.  Fortunately, this routine should only need to be calculated once at the
beginning of each run.

@<Prototype for |Fill_AB_Grid|@>=
void Fill_AB_Grid(struct measure_type m, struct invert_type r)

@ @<Definition for |Fill_AB_Grid|@>=
    @<Prototype for |Fill_AB_Grid|@>
{
    int i,j;
    double min_log_b = -8;  /* exp(-10) is smallest thickness */
    double max_log_b = +8;   /* exp(+8) is greatest thickness */

    if (Debug(DEBUG_GRID))
        fprintf(stderr, "GRID: Filling AB grid (g=%.5f)\n", RR.slab.g);

    if (The_Grid==NULL) Allocate_Grid(r.search);
    @<Zero \\{GG}@>@;

    Set_Calc_State(m,r);

    GG_g = RR.slab.g;
    for(i=0; i<GRID_SIZE; i++){
        double x = (double) i/(GRID_SIZE-1.0);
        RR.slab.b = exp(min_log_b + (max_log_b-min_log_b) *x);
        for(j=0; j<GRID_SIZE; j++){
            @<Generate next albedo using j@>@;
            fill_grid_entry(i,j);
        }
    }

The_Grid_Initialized=TRUE;
The_Grid_Search = FIND_AB;
}

@ Now it seems that we must be a bit more subtle in choosing the values
of albedos to use in the grid.  Spacing the points sometimes gives too
coarse of spacing as $a\rightarrow0$ or $a\rightarrow1$.  A function that
seems to work is
$$
a = 1 - (1-x)^2\cdot (1+2x)\qquad{\rm where}\qquad x={j-1\over n-1}
$$
where $j = 0,1,\ldots,n$.  Thus instead of getting
$$
0.0, 0.1,0.2,\ldots,0.9,1.0
$$
we get
$$
0.000, 0.028, 0.104, 0.216, 0.352, 0.500, 0.648, 0.784, 0.896, 0.972, 1.000
$$
which is has symmetric spacing at both ends of the range of $a$.

@ Here is heuristic that seems to work well

@<Generate next albedo using j@>=
    {
    double x = (double) j/(GRID_SIZE-1.0);
    RR.slab.a = 1.0 - (1.0 - x) * (1.0 - x) * (1.0 + 2.0 * x);
    }

@ Let's do the same thing for $g$ grid points with the exception that points
are contrained between |-MAX_ABS_G| and |+MAX_ABS_G|

@<Generate next anisotropy using i@>=
    {
    double x = (double) i/(GRID_SIZE-1.0);
    double xx = (1.0 - x) * (1.0 - x) * (1.0 + 2.0 * x);
    RR.slab.g = (1.0 - 2.0 * xx) * MAX_ABS_G;
    }

@ @<Generate next anisotropy using j@>=
    {
    double x = (double) j/(GRID_SIZE-1.0);
    double xx = (1.0 - x) * (1.0 - x) * (1.0 + 2.0 * x);
    RR.slab.g = (1.0 - 2.0 * xx) * MAX_ABS_G;
    }

@ This is quite similar to |Fill_AB_Grid|, with the exception of the
little shuffle I do at the beginning to figure out the optical thickness to
use.  The problem is that the optical thickness may not be known.  If it is
known then the only way that we could have gotten here is if the user
dictated |FIND_AG| and specified |b| and only provided two measurements.
Otherwise, the user must have made three measurements and the optical
depth can be figured out from |m.m_u|.

This routine could also be improved by not recalculating the anisotropy
matrix for every point.  But this would only end up being a minor performance
enhancement if it were fixed.

@<Prototype for |Fill_AG_Grid|@>=
void Fill_AG_Grid(struct measure_type m, struct invert_type r)

@ @<Definition for |Fill_AG_Grid|@>=
    @<Prototype for |Fill_AG_Grid|@>
{
    int i,j;
    double max_a = -10;
    double min_a = 10;

    if (Debug(DEBUG_GRID))
        fprintf(stderr, "GRID: Filling AG grid\n");

    if (The_Grid==NULL) Allocate_Grid(r.search);
    @<Zero \\{GG}@>@;

    Set_Calc_State(m,r);
    GG_b=r.slab.b;
    for(i=0; i<GRID_SIZE; i++){
        @<Generate next anisotropy using i@>@;
        for(j=0; j<GRID_SIZE; j++){
            @<Generate next albedo using j@>@;
            fill_grid_entry(i,j);
            if (RR.slab.a < min_a) min_a = RR.slab.a;
            if (RR.slab.a > max_a) max_a = RR.slab.a;
        }
    }

    if (Debug(DEBUG_GRID)) {
        fprintf(stderr, "GRID: a        = %9.7f to %9.7f \n", min_a, max_a);
        fprintf(stderr, "GRID: b        = %9.5f \n", r.slab.b);
        fprintf(stderr, "GRID: g  range = %9.6f to %9.6f \n", -MAX_ABS_G, MAX_ABS_G);
    }

    The_Grid_Initialized=TRUE;
    The_Grid_Search = FIND_AG;
}

@   @<Zero \\{GG}@>=
GG_a = 0.0;
GG_b = 0.0;
GG_g = 0.0;
GG_bs = 0.0;
GG_ba = 0.0;


@ This is quite similar to |Fill_AB_Grid|, with the exception of the
that the albedo is held fixed while $b$ and $g$ are varied.

This routine could also be improved by not recalculating the anisotropy
matrix for every point.  But this would only end up being a minor performance
enhancement if it were fixed.

@<Prototype for |Fill_BG_Grid|@>=
void Fill_BG_Grid(struct measure_type m, struct invert_type r)

@ @<Definition for |Fill_BG_Grid|@>=
    @<Prototype for |Fill_BG_Grid|@>
{
    int i,j;
    double min_log_b = -8;  /* exp(-10) is smallest thickness */
    double max_log_b = +10;  /* exp(+8)  is greatest thickness */

    if (The_Grid==NULL) Allocate_Grid(r.search);
    @<Zero \\{GG}@>@;

    if (Debug(DEBUG_GRID))
        fprintf(stderr, "GRID: Filling BG grid\n");

    Set_Calc_State(m,r);
    RR.slab.a = RR.default_a;
    GG_a = RR.slab.a;

    for(i=0; i<GRID_SIZE; i++){
        double x = (double) i/(GRID_SIZE-1.0);
        RR.slab.b = exp(min_log_b + (max_log_b-min_log_b) *x);
        for(j=0; j<GRID_SIZE; j++){
            @<Generate next anisotropy using j@>@;
            fill_grid_entry(i,j);
        }
    }

    if (Debug(DEBUG_GRID)) {
        fprintf(stderr, "GRID: a        = %9.7f \n", RR.default_a);
        fprintf(stderr, "GRID: b  range = %9.5f to %9.3f \n", exp(min_log_b), exp(max_log_b));
        fprintf(stderr, "GRID: g  range = %9.6f to %9.6f \n", -MAX_ABS_G, MAX_ABS_G);
    }

    The_Grid_Initialized=TRUE;
    The_Grid_Search = FIND_BG;
}

@ This is quite similar to |Fill_BG_Grid|, with the exception of the
that the $b_s=\mu_s d$ is held fixed.  Here $b$ and $g$ are varied
on the usual grid, but the albedo is forced to take whatever value
is needed to ensure that the scattering remains fixed.

@<Prototype for |Fill_BaG_Grid|@>=
void Fill_BaG_Grid(struct measure_type m, struct invert_type r)

@ @<Definition for |Fill_BaG_Grid|@>=
    @<Prototype for |Fill_BaG_Grid|@>
{
    int i,j;
    double max_a = -10;
    double min_a = 10;
    double bs = r.default_bs;
    double min_log_ba = -8;  /* exp(-10) is smallest thickness */
    double max_log_ba = +10;  /* exp(+8)  is greatest thickness */

    if (The_Grid==NULL) Allocate_Grid(r.search);
    @<Zero \\{GG}@>@;

    if (Debug(DEBUG_GRID)) {
        fprintf(stderr, "GRID: Filling BaG grid\n");
        fprintf(stderr, "GRID:       bs = %9.5f\n", bs);
        fprintf(stderr, "GRID: ba range = %9.6f to %9.3f \n", exp(min_log_ba), exp(max_log_ba));
    }

    Set_Calc_State(m,r);
    GG_bs = bs;
    for(i=0; i<GRID_SIZE; i++){
        double x = (double) i/(GRID_SIZE-1.0);
        double ba = exp(min_log_ba + (max_log_ba-min_log_ba) *x);
        RR.slab.b = ba + bs;
        if (RR.slab.b > 0)
            RR.slab.a = bs / RR.slab.b;
        else
            RR.slab.a = 0;
        if (RR.slab.a < 0) RR.slab.a = 0;
        if (RR.slab.a < min_a) min_a = RR.slab.a;
        if (RR.slab.a > max_a) max_a = RR.slab.a;

        for(j=0; j<GRID_SIZE; j++){
            @<Generate next anisotropy using j@>@;
            fill_grid_entry(i,j);
        }
    }

    if (Debug(DEBUG_GRID)) {
        fprintf(stderr, "GRID: a        = %9.7f to %9.7f \n", min_a, max_a);
        fprintf(stderr, "GRID: b  range = %9.5f to %9.3f \n", exp(min_log_ba)+bs, exp(max_log_ba)+bs);
        fprintf(stderr, "GRID: g  range = %9.6f to %9.6f \n", -MAX_ABS_G, MAX_ABS_G);
    }

    The_Grid_Initialized=TRUE;
    The_Grid_Search = FIND_BaG;
}

@ Very similiar to the above routine, but holding $b_a=\mu_a d$ fixed.
Here $b$ and $g$ are varied on the usual grid, but the albedo is forced
to take whatever value is needed to ensure that the absorption remains fixed.

@<Prototype for |Fill_BsG_Grid|@>=
void Fill_BsG_Grid(struct measure_type m, struct invert_type r)

@ @<Definition for |Fill_BsG_Grid|@>=
    @<Prototype for |Fill_BsG_Grid|@>
{
    int i,j;
    double max_a = -10;
    double min_a = 10;
    double ba = r.default_ba;
    double min_log_bs = -8;  /* exp(-10) is smallest thickness */
    double max_log_bs = +10;  /* exp(+8)  is greatest thickness */

    if (The_Grid==NULL) Allocate_Grid(r.search);
    @<Zero \\{GG}@>@;

    if (Debug(DEBUG_GRID)) {
        fprintf(stderr, "GRID: Filling BsG grid\n");
        fprintf(stderr, "GRID:       ba = %9.5f\n", ba);
        fprintf(stderr, "GRID: bs range = %9.6f to %9.3f \n", exp(min_log_bs), exp(max_log_bs));
    }

    Set_Calc_State(m,r);
    GG_ba = RR.default_ba;
    for(i=0; i<GRID_SIZE; i++){
        double x = (double) i/(GRID_SIZE-1.0);
        double bs = exp(min_log_bs + (max_log_bs-min_log_bs) *x);
        RR.slab.b = ba + bs;
        if (RR.slab.b > 0)
            RR.slab.a = 1-RR.default_ba/RR.slab.b;
        else
            RR.slab.a = 0;
        if (RR.slab.a < 0) RR.slab.a = 0;
        if (RR.slab.a < min_a) min_a = RR.slab.a;
        if (RR.slab.a > max_a) max_a = RR.slab.a;

        for(j=0; j<GRID_SIZE; j++){
            @<Generate next anisotropy using j@>@;
            fill_grid_entry(i,j);
        }
    }

    if (Debug(DEBUG_GRID)) {
        fprintf(stderr, "GRID: a  range = %9.7f to %9.7f \n", min_a, max_a);
        fprintf(stderr, "GRID: b  range = %9.5f to %9.3f \n", exp(min_log_bs)+ba, exp(max_log_bs)+ba);
        fprintf(stderr, "GRID: g  range = %9.6f to %9.6f \n", -MAX_ABS_G, MAX_ABS_G);
    }

The_Grid_Initialized=TRUE;
The_Grid_Search = FIND_BsG;
}

@ @<Prototype for |Fill_Grid|@>=
void Fill_Grid(struct measure_type m, struct invert_type r, int force_new)

@ @<Definition for |Fill_Grid|@>=
    @<Prototype for |Fill_Grid|@>
{
    if (force_new || !Same_Calc_State(m,r)) {
        switch (r.search) {
            case FIND_AB:
                Fill_AB_Grid(m,r);
                break;
            case FIND_AG:
                Fill_AG_Grid(m,r);
                break;
            case FIND_BG:
                Fill_BG_Grid(m,r);
                break;
            case FIND_BaG:
                Fill_BaG_Grid(m,r);
                break;
            case FIND_BsG:
                Fill_BsG_Grid(m,r);
                break;
            default:
                AD_error("Attempt to fill grid for unknown search case.");
        }
    }

    Get_Calc_State(&MGRID, &RGRID);
}

@*1 Calculating R and T.


|Calculate_Distance| returns the distance between the measured
values in |MM| and the calculated values for the current
guess at the optical properties.  It
assumes that the everything in the local variables |MM| and |RR| have
been set appropriately.

@<Prototype for |Calculate_Distance|@>=
void Calculate_Distance(double *M_R, double *M_T, double *deviation)

@ @<Definition for |Calculate_Distance|@>=
    @<Prototype for |Calculate_Distance|@>
{
    double Ru, Tu, ur1, ut1, uru, utu;

    if (RR.slab.b <= 1e-6 )
        RR.slab.b = 1e-6;

    RT_Flip(MM.flip_sample, RR.method.quad_pts, &RR.slab, &ur1, &ut1, &uru, &utu);

    Sp_mu_RT_Flip(MM.flip_sample,
             RR.slab.n_top_slide, RR.slab.n_slab, RR.slab.n_bottom_slide,
             RR.slab.b_top_slide, RR.slab.b,      RR.slab.b_bottom_slide,
             RR.slab.cos_angle, &Ru, &Tu);

    if ((!CALCULATING_GRID && Debug(DEBUG_ITERATIONS)) ||
        ( CALCULATING_GRID && Debug(DEBUG_GRID_CALC)))
            fprintf(stderr, "        ");

    Calculate_Distance_With_Corrections(ur1,ut1,Ru,Tu,uru,utu,M_R,M_T,deviation);
}

@ @<Prototype for |Calculate_Grid_Distance|@>=
double Calculate_Grid_Distance(int i, int j)

@ @<Definition for |Calculate_Grid_Distance|@>=
    @<Prototype for |Calculate_Grid_Distance|@>
{
    double ur1,ut1,uru,utu,Ru,Tu,b,dev,LR,LT;

    if (Debug(DEBUG_GRID_CALC) && i==0 && j==0) {
        fprintf(stderr, "+   i   j ");
        fprintf(stderr, "      a         b          g     |");
        fprintf(stderr, "     M_R        grid   |");
        fprintf(stderr, "     M_T        grid   |  distance\n");
    }

    if (Debug(DEBUG_GRID_CALC))
        fprintf(stderr, "g %3d %3d ",i,j);

    b   = The_Grid[GRID_SIZE*i+j][B_COLUMN];
    ur1 = The_Grid[GRID_SIZE*i+j][UR1_COLUMN];
    ut1 = The_Grid[GRID_SIZE*i+j][UT1_COLUMN];
    uru = The_Grid[GRID_SIZE*i+j][URU_COLUMN];
    utu = The_Grid[GRID_SIZE*i+j][UTU_COLUMN];
    RR.slab.a = The_Grid[GRID_SIZE*i+j][A_COLUMN];
    RR.slab.b = The_Grid[GRID_SIZE*i+j][B_COLUMN];
    RR.slab.g = The_Grid[GRID_SIZE*i+j][G_COLUMN];

    Sp_mu_RT_Flip(MM.flip_sample,
             RR.slab.n_top_slide, RR.slab.n_slab, RR.slab.n_bottom_slide,
             RR.slab.b_top_slide, b,      RR.slab.b_bottom_slide,
             RR.slab.cos_angle, &Ru, &Tu);

    CALCULATING_GRID = 1;
    Calculate_Distance_With_Corrections(ur1,ut1,Ru,Tu,uru,utu,&LR,&LT,&dev);
    CALCULATING_GRID = 0;

    return dev;
}

@ This is the routine that actually finds the distance.  I have factored
this part out so that it can be used in the |Near_Grid_Points| routine.

|Ru| and |Tu| refer to the unscattered (collimated) reflection and transmission.

The only tricky part is to remember that the we are trying to match the
measured values.  The measured values are affected by sphere parameters
and light loss.  Since the values |UR1| and |UT1| are for an infinite slab
sample with no light loss, the light loss out the edges must be subtracted.
It is these values that are used with the sphere formulas to convert the
modified |UR1| and |UT1| to values for |*M_R| and |*M_T|.

@<Prototype for |Calculate_Distance_With_Corrections|@>=
void Calculate_Distance_With_Corrections(
                double UR1, double UT1,
                double Ru,  double Tu,
                double URU, double UTU,
                double *M_R, double *M_T, double *dev)

@ @<Definition for |Calculate_Distance_With_Corrections|@>=
    @<Prototype for |Calculate_Distance_With_Corrections|@>
{
    @<Determine calculated light to be used@>@;

    switch (MM.num_spheres) {
        case 0:
            @<Calc |M_R| and |M_T| for no spheres@>@;
            break;

        case 1:
            if (MM.method == COMPARISON) {
                @<Calc |M_R| and |M_T| for dual beam sphere@>@;
            } else {
                @<Calc |M_R| and |M_T| for single beam sphere@>@;
            }break;

        case 2:
            @<Calc |M_R| and |M_T| for two spheres@>@;
            break;

        default:
            fprintf(stderr, "Bad number of spheres = %d\n", MM.num_spheres);
            exit(EXIT_FAILURE);
    }

    @<Calculate the deviation@>@;
    @<Print diagnostics@>@;
}

@ The calculated values for |M_R| and |M_T| must be adapted to match the measurements.
The diffuse light |URU| and |UTU| are used to determine the gain from the sphere.
They're only modified by the lost light calculation.
All values can become slightly negative because the Monte Carlo is noisy.
Negative values are set to zero.

@<Determine calculated light to be used@>=
    double UR1_calc, UT1_calc, URU_calc, UTU_calc;

    URU_calc = URU - MM.uru_lost;
    if (URU_calc < 0) URU_calc = 0;

    UTU_calc = UTU - MM.utu_lost;
    if (UTU_calc < 0) UTU_calc = 0;

@ The measurements for |UR1| and |UT1| need to be modified to accommodate light that
misses the detector either because it is intentionally not collected
(unscattered light) or it leaks out (lost light).  Since none of the light
that leaks out could be unscattered light, these two are independent
of one another.

@<Determine calculated light to be used@>=
    UR1_calc = UR1 - (1.0 - MM.fraction_of_ru_in_mr) * Ru - MM.ur1_lost;
    if (UR1_calc < 0) UR1_calc = 0;

    UT1_calc = UT1 - (1.0 - MM.fraction_of_tu_in_mt) * Tu - MM.ut1_lost;
    if (UT1_calc < 0) UT1_calc = 0;

@ When no spheres are used, then no corrections can or need to
be made.  The lost light estimates in |MM.ur1_lost| and |MM.ut1_lost|
should be zero and so the values for |UR1_calc| and |UT1_calc|
properly account for the presence or absence of unscattered light.

@<Calc |M_R| and |M_T| for no spheres@>=
{
    *M_R = UR1_calc;
    *M_T = UT1_calc;
}

@*2 Reflectance measurement for one sphere.

The total reflection from a slab may be broken down into light that has and has
not been scattered.  The total reflectance for normal incidence on an infinite slab is
denoted by van de Hulst by {\tt UR1}.  To track integrating sphere corrections, this
is separated into scattered and  unscattered parts.
$$
{\tt UR1} = R_\unscat + R_\scat
$$
{\tt UR1} is calculated from the current guess of optical properties, but we can also
calculate $R_\unscat$ by setting $\mu_s=0$ which means that the albedo is zero. Assuming
the incident power is $P_i$ then the scattered light in the sphere $P_{ss}$ and the
unscattered light in the sphere $P_{su}$ start as
$$
P_{ss} = ({\tt UR1} - R_\unscat) P_i
    \qquad\hbox{\it and}\qquad
P_{su} =  R_\unscat P_i
$$

@ Light that is lost.

In an experiment, the scattered light in the sphere will be reduced by any light that
leaks out.  This is determined by doing a Monte Carlo simulation to determine ${\tt UR1}_\lost$
or how much light is scattered within the sample and escapes outside the sphere. Since
the only way that light may be lost is through scattering, this does not affect the
unscattered fraction.

The fraction of unscattered light collected by the sphere $f_\unscat$ is
determined by the experimentalist. The unscattered reflected light can aligned so that
it bounces off the sample and exits completely through the entrance port so that
$f_\unscat=0$.  Alternatively, the beam may be incident on the sample at an
angle so that the unscattered light will bounce and remain completely within
the sphere so that $f_\unscat=1$.
$$
P_{ss} = ({\tt UR1} - {\tt UR1}_\unscat-R_\lost) P_i
    \qquad\hbox{\it and}\qquad
P_{su} =  f_\unscat R_\unscat P_i
$$

@ Incident light that misses the sample.

In an experiment, a fraction $f_\miss$ of the incident beam might miss the sample
and hit the sphere wall instead.  In this case, both $R_{ss}$ and $R_{su}$ are affected.
Typically the beam is much smaller than the sample and $f_\miss=0$, however sometimes
the beam diverges too much and some of the beam hits the wall. After hitting the sphere
wall then $f_\miss r_w$ will be added to the scattered light in the sphere. However,
now only $(1-f_\miss)$ hits the sample directly, both the scattered and
unscattered light in the sphere must be adjusted accordingly.  So in the
unusual case of non-zero $f_\miss$, we have
$$
R_{ss} = (1-f_\miss) ({\tt UR1} - R_\unscat - {\tt UR1}_\lost)P_i + f_\miss r_w P_i
    \qquad\hbox{\it and}\qquad
R_{su} = (1-f_\miss) f_\unscat R_\unscat P_i
$$

@ Reflectance with baffle.

When a baffle blocks light from passing directly between the sample to the detector then
the light reflected by the sample must bounce once.  Some of the light will be
reflected by the sphere walls, but some may be reflected by the third port. The
weighted reflectance of the first bounce is
$$
r_\first=(1-a_t) r_w + a_t r_t
$$
We can safely assume that the fraction of light generated by $f_\miss$ will
originate close enough to the sample that a baffle, if present, will prevent light
from directly reaching the detector.  The scattered light in the sphere after
the first bounce will be
$$
P_{ss} = r_\first \left[(1-f_\miss) ({\tt UR1} - R_\unscat - R_\lost) + f_\miss r_w\right] P_i
$$
All unscattered reflectance that is collected must hit the sphere wall (otherwise it
would exit through the entrance port and not be collected!).  This means
that the correction in $r_\first$ for the entrance port is not needed and the
unreflected light can just be multiplied by $r_w$
$$
P_{su} =  r_w (1-f_\miss) f_\unscat R_\unscat P_i
$$
The last step is to account for the sphere gain. The sample is held in the sample port
and the entrance port is empty.  The total reflection for diffuse illumination of the
sample is {\tt URU}.  This quantity must also be corrected for light that is not
collected by the sphere ${\tt UR1}_\lost$.  The gain for such a sphere is
$G({\tt URU}-{\tt URU}_\lost,0)$.  Finally,
$$
P_s = [a_d (1-r_d)] \cdot G({\tt URU}-{\tt URU}_\lost,0)\left[P_{ss}+P_{su}\right]
$$

@ The reflection standard.

We let ${\tt UR1}={\tt URU}=r_\std$, $R_\unscat=0$, $R_\lost=0$ to get
$$
P_\std = [a_d (1-r_d)] \cdot r_\first \left[(1-f_\miss) r_\std + f_\miss r_w\right] G(r_\std,0) P_i
$$

@ The open port.

We let ${\tt UR1}={\tt URU}=0$, $R_\unscat=0$, $R_\lost=0$ to get
$$
P_0 = [a_d (1-r_d)] \cdot r_\first f_\miss r_w G(0,0) P_i
$$

@ The unbaffled reflectance sphere.

In this case, light can reach the detector from the sample.  The first bounce is
not needed which can be accommodated by letting $r_\first=1$.  We, of course,
assume that the gain is calculated assuming that no baffle is present.
The incident power $P_i$ and the quantities in square brackets are identical for
$P_s$, $P_\std$, and $P_0$ and cancel in the normalized reflection fraction
$$
M_R = r_\std\cdot{P_s-P_0 \over P_\std - P_0}
$$
In addition, the entrance port is empty and therefore $r_t=0$ and can be omitted
from the $r_\first$ calculation.  This leads to the following code for |M_R|

@<Calc |M_R| and |M_T| for single beam sphere@>=
    double P_std, P, P_0, G, G_0, G_std, r_first, P_ss, P_su;
    r_first = 1;

    if (MM.baffle_r) r_first = MM.rw_r * (1 - MM.at_r);

    UR1_calc = UR1 - Ru - MM.ur1_lost;
    if (UR1_calc < 0) UR1_calc = 0;

    G_0      = Gain(REFLECTION_SPHERE, MM, 0.0, 0.0);
    G        = Gain(REFLECTION_SPHERE, MM, URU_calc, 0.0);
    G_std    = Gain(REFLECTION_SPHERE, MM, MM.rstd_r, 0.0);

    P_std    = G_std * (MM.rstd_r * (1-MM.f_r) + MM.f_r * MM.rw_r);
    P_0      = G_0   * (                         MM.f_r * MM.rw_r);
    P_ss = r_first * (UR1_calc * (1-MM.f_r) + MM.f_r * MM.rw_r);
    P_su = MM.rw_r * (1-MM.f_r) * MM.fraction_of_ru_in_mr * Ru;
    P = G * (P_ss + P_su);

    *M_R     = MM.rstd_r * (P - P_0)/(P_std - P_0);

    if (Debug(DEBUG_SPHERE_GAIN) && !CALCULATING_GRID) {
        fprintf(stderr, "SPHERE: REFLECTION\n");
        fprintf(stderr, "SPHERE:      G0 = %7.3f      G  = %7.3f G_std = %7.3f\n", G_0, G, G_std);
        fprintf(stderr, "SPHERE:      P0 = %7.3f      P  = %7.3f P_std = %7.3f\n", P_0, P, P_std);
        fprintf(stderr, "SPHERE:     UR1 = %7.3f UR1calc = %7.3f   M_R = %7.3f\n", UR1, UR1_calc, *M_R);
    }

@*2 Transmittance measurement for one sphere.

Like in the reflection case, the total transmission can be split into unscattered transmission and
scattered transmission,
$$
UT1=T_u + T_s
$$
We define $P_{ss}$ as the scattered light in the sphere and $P_{su}$ as the unscattered
light in the sphere for an incident power $P_i$
$$
P_{ss} = ({\tt UT1} - T_u)P_i
    \qquad\hbox{\it and}\qquad
P_{su} = f_\unscat T_u P_i
$$

@ Transmitted light not collected.

The scattered light will be reduced by {\tt UT1}$_\lost$. The unscattered light
will be affected by the fraction of unscattered light collected by the sphere.
When transmission measurements are made, the third port (opposite the sample port
in the sphere) is typically filled with a known standard. However, the third port
might also be left open so that all the scattered light might exit and only scattered
light is collected.  In the former case $f_\unscat=1$ and in the latter case $f_\unscat=0$.
So including these effects gives
$$
P_{ss} = ({\tt UT1} - T_u - T_\lost)P_i
    \qquad\hbox{\it and}\qquad
P_{su} = f_\unscat T_u
$$

@ The baffling case of transmission.

With a baffle, then scattered light from the sample cannot reach
the detector without a bounce.  This weighted reflection is given by
$$
r_\first= (1-a_t) r_w + a_t r_t
$$
The unscattered light will be reflected by the standard $r_t$ in the third port
and so
$$
P_{ss} = r_\first ({\tt UT1} - T_u - T_\lost) P_i
    \qquad\hbox{\it and}\qquad
P_{su} = r_t f_\unscat T_u
$$
The last step is to account for the sphere gain. The sample is held in the sample port
and the third port reflects $r_t$.  The total reflection for diffuse illumination of the
sample is {\tt URU}.  This quantity must also be corrected for light that is not
collected by the sphere ${\tt UR1}_\lost$.  The gain for such a sphere is
$G({\tt URU}-{\tt URU}_\lost,r_t)$.  Finally,
$$
P_s = G({\tt URU}- {\tt URU}_\lost, r_t) (P_{ss}+P_{su})
$$

@ No baffle.

Here $r_\first=1$ and the gain should be calculated assuming no baffle.

@ The standard measurement.

When transmission measurements are made, typically the third port (the one that
let the light into the sphere for the reflection measurement) is filled with a known
standard. In this case, the natural way to make the 100\% transmission measurement
is to shine the beam through the empty sample port onto the known standard.

We let ${\tt UT1}=T_u=1$, $T_\lost=0$, ${\tt URU} =0$, $r_t=r_\std$, and $f_\unscat=1$ so
$$
P_\std = G(0, r_\std) r_\std P_i
$$

The estimate for the measured transmittance is
$$
M_T = r_\std {P_s \over P_\std}
= {P_{ss}+P_{su}\over P_i} \cdot {G({\tt URU}- {\tt URU}_\lost, r_\third)\over G(0, r_\std)}
$$

@<Calc |M_R| and |M_T| for single beam sphere@>=
{
    double r_first = 1;
    double r_third = MM.rstd_t;

    if (MM.fraction_of_tu_in_mt == 0) r_third = 0;

    if (MM.baffle_t) r_first = MM.rw_t * (1-MM.at_t) + MM.rstd_t * MM.at_t;

    UT1_calc = UT1 - Tu - MM.ut1_lost;
    if (UT1_calc < 0) UT1_calc = 0;

    G = Gain(TRANSMISSION_SPHERE, MM, URU_calc, r_third);
    G_std = Gain(TRANSMISSION_SPHERE, MM, 0, MM.rstd_t);

    *M_T  = (r_third * Tu * MM.fraction_of_tu_in_mt + r_first * UT1_calc) * G/G_std;

    if (Debug(DEBUG_SPHERE_GAIN) && !CALCULATING_GRID) {
        fprintf(stderr, "SPHERE: TRANSMISSION\n");
        fprintf(stderr, "SPHERE:      G  = %7.3f   G_std = %7.3f\n", G, G_std);
        fprintf(stderr, "SPHERE:     UT1 = %7.3f UT1calc = %7.3f T_c = %7.3f\n", UT1, UT1_calc, Tu);
        fprintf(stderr, "SPHERE:     M_T = %7.3f\n", *M_T);
        fprintf(stderr, "\n");
    }
}

@*2 Dual beam case for one sphere.

@ The dual beam case is different because the sphere efficiency
is equivalent for measurement of light hitting the sample first or
hitting the reference standard first.  The dual beam measurement
should report the ratio of these two reflectance measurements, thereby
eliminating the need to calculate the sphere gain.

The only correction that needs to be made have already been made, namely
subtracting the |UR1| or |UT1| lost light and also accounting for whether
or not unscattered light is collected.

Originally, I had a bunch of calculations trying to account for light
that hits the sphere wall first.  Since the exact details of how a
dual beam spectrometer reports its measurements is unknown, it makes
no sense to try and account for it.

@<Calc |M_R| and |M_T| for dual beam sphere@>=
{
    *M_R     = UR1_calc;
    *M_T     = UT1_calc;
}

@*2 Double integrating spheres.

@ When two integrating spheres are present then the
double integrating sphere formulas are slightly more complicated.

The normalized sphere measurements for two spheres are

$$
M_R = {R(\rdirect,\rdiffuse,\tdirect,\tdiffuse) - R(0,0,0,0)
\over R(r_\std,r_\std,0,0) - R(0,0,0,0)}
$$
and
$$
M_T = {T(\rdirect,\rdiffuse,\tdirect,\tdiffuse) - T(0,0,0,0)
                   \over T(0,0,1,1) - T(0,0,0,0)}
$$

Note that |R_0| and |T_0| will be zero unless one has explicitly
set the fraction |m.f_r| or |m.f_t| to be non-zero.

@<Calc |M_R| and |M_T| for two spheres@>=
{
double R_0, T_0;
R_0 = Two_Sphere_R(MM, 0, 0, 0, 0);
T_0 = Two_Sphere_T(MM, 0, 0, 0, 0);

*M_R = MM.rstd_r * (Two_Sphere_R(MM, UR1_calc, URU_calc, UT1_calc, UTU_calc) - R_0)/
                   (Two_Sphere_R(MM, MM.rstd_r, MM.rstd_r, 0, 0) - R_0);
*M_T =  (Two_Sphere_T(MM, UR1_calc, URU_calc, UT1_calc, UTU_calc) - T_0)/
                   (Two_Sphere_T(MM, 0, 0, 1, 1) - T_0);
}

@ There are at least three things that need to be considered here.
First, the number of measurements.  Second, is the metric is relative or absolute.
And third, is the albedo fixed at zero which means that the transmission
measurement should be used instead of the reflection measurement.

@<Calculate the deviation@>=

if (RR.search==FIND_A  || RR.search==FIND_G || RR.search==FIND_B ||
    RR.search==FIND_Bs || RR.search == FIND_Ba) {
        @<One parameter deviation@>@;
} else {
        @<Two parameter deviation@>@;
}

@ This part was slightly tricky.  The crux of the problem was to
decide if the transmission or the reflection was trustworthy.  After
looking a bunches of measurements, I decided that the transmission
measurement was almost always more reliable.  So when there is just
a single measurement known, then use the total transmission if it
exists.

@<One parameter deviation@>=

if ( MM.m_t > 0 ){
    if (RR.metric == RELATIVE)
        *dev = fabs(MM.m_t - *M_T) / (MM.m_t + ABIT);
    else
        *dev = fabs(MM.m_t - *M_T) ;
} else {
    if (RR.metric == RELATIVE)
        *dev = fabs(MM.m_r - *M_R) / (MM.m_r + ABIT);
    else
        *dev = fabs(MM.m_r - *M_R) ;
}

@ This stuff happens when we are doing two parameter searches.
In these cases there should be information in both R and T.
The distance should be calculated using the deviation from
both.  The albedo stuff might be able to be take out.  We'll see.

@<Two parameter deviation@>=

    if (RR.metric == RELATIVE) {
        if (MM.m_t > ABIT)
            *dev = T_TRUST_FACTOR* fabs(MM.m_t - *M_T) / (UTU_calc + ABIT);
        if ( RR.default_a != 0 ) {
            *dev += fabs(MM.m_r - *M_R) / (URU_calc + ABIT);
        }
    } else {
        *dev = T_TRUST_FACTOR * fabs(MM.m_t - *M_T);
        if ( RR.default_a != 0 )
            *dev += fabs(MM.m_r - *M_R);
    }


@ This is here so that I can figure out why the program is not converging.
  This is a little convoluted so that the global constants at the top of
  this file interact properly.

@<Print diagnostics@>=
if ((Debug(DEBUG_ITERATIONS) && !CALCULATING_GRID) || @|
    (Debug(DEBUG_GRID_CALC)  &&  CALCULATING_GRID)) {
    fprintf(stderr, "%10.5f %10.4f %10.5f |", RR.slab.a, RR.slab.b, RR.slab.g);
    fprintf(stderr, " %10.5f %10.5f |", MM.m_r, *M_R);
    fprintf(stderr, " %10.5f %10.5f |", MM.m_t, *M_T);
    fprintf(stderr, "%10.3f\n", *dev);
}

@ @<Prototype for |Find_AG_fn|@>=
double Find_AG_fn(double x[])

@ @<Definition for |Find_AG_fn|@>=
    @<Prototype for |Find_AG_fn|@>
{
    double m_r,m_t,deviation;
    RR.slab.a = acalc2a(x[1]);
    RR.slab.g = gcalc2g(x[2]);
    Calculate_Distance(&m_r,&m_t,&deviation);
    return deviation;
}

@ @<Prototype for |Find_AB_fn|@>=
double Find_AB_fn(double x[])

@ @<Definition for |Find_AB_fn|@>=
    @<Prototype for |Find_AB_fn|@>
{
    double m_r,m_t,deviation;
    RR.slab.a = acalc2a(x[1]);
    RR.slab.b = bcalc2b(x[2]);
    Calculate_Distance(&m_r,&m_t,&deviation);
    return deviation;
}

@ @<Prototype for |Find_Ba_fn|@>=
double Find_Ba_fn(double x)

@ This is tricky only because the value in |RR.slab.b| is used to hold the
value of |bs| or $d \cdot \mu_s$.  It must be switched to the correct value
for the optical thickness and then switched back at the end of the routine.

@<Definition for |Find_Ba_fn|@>=
    @<Prototype for |Find_Ba_fn|@>
{
    double m_r,m_t,deviation,ba,bs;

    bs        = RR.slab.b;
    ba        = bcalc2b(x);
    RR.slab.b = ba+bs;        /* unswindle */
    RR.slab.a = bs/(ba+bs);

    Calculate_Distance(&m_r,&m_t,&deviation);

    RR.slab.b = bs;          /* swindle */
    return deviation;
}

@ See the comments for the |Find_Ba_fn| routine above.  Play the same trick
but use |ba|.

@<Prototype for |Find_Bs_fn|@>=
double Find_Bs_fn(double x)

@ @<Definition for |Find_Bs_fn|@>=
    @<Prototype for |Find_Bs_fn|@>
{
    double m_r,m_t,deviation,ba,bs;

    ba        = RR.slab.b;  /* unswindle */
    bs        = bcalc2b(x);
    RR.slab.b = ba+bs;
    RR.slab.a = bs/(ba+bs);

    Calculate_Distance(&m_r,&m_t,&deviation);

    RR.slab.b = ba;         /* swindle */
    return deviation;
}

@ @<Prototype for |Find_A_fn|@>=
double Find_A_fn(double x)

@ @<Definition for |Find_A_fn|@>=
    @<Prototype for |Find_A_fn|@>
{
    double m_r,m_t,deviation;
    RR.slab.a = acalc2a(x);
    Calculate_Distance(&m_r,&m_t,&deviation);
    return deviation;
}

@ @<Prototype for |Find_B_fn|@>=
double Find_B_fn(double x)

@ @<Definition for |Find_B_fn|@>=
    @<Prototype for |Find_B_fn|@>
{
    double m_r,m_t,deviation;
    RR.slab.b = bcalc2b(x);
    Calculate_Distance(&m_r,&m_t,&deviation);
    return deviation;
}

@ @<Prototype for |Find_G_fn|@>=
double Find_G_fn(double x)

@ @<Definition for |Find_G_fn|@>=
    @<Prototype for |Find_G_fn|@>
{
    double m_r,m_t,deviation;
    RR.slab.g = gcalc2g(x);
    Calculate_Distance(&m_r,&m_t,&deviation);
    return deviation;
}

@ @<Prototype for |Find_BG_fn|@>=
double Find_BG_fn(double x[])

@ @<Definition for |Find_BG_fn|@>=
    @<Prototype for |Find_BG_fn|@>
{
    double m_r,m_t,deviation;
    RR.slab.b = bcalc2b(x[1]);
    RR.slab.g = gcalc2g(x[2]);
    RR.slab.a = RR.default_a;
    Calculate_Distance(&m_r,&m_t,&deviation);
    return deviation;
}

@ For this function the first term |x[1]| will contain the value of
$\mu_s d$, the second term will contain the anisotropy.  Of course the
first term is in the bizarre calculation space and needs to be
translated back into normal terms before use.  We just at the scattering
back on and voil{\'a} we have a useable value for the optical depth.

@<Prototype for |Find_BaG_fn|@>=
double Find_BaG_fn(double x[])

@ @<Definition for |Find_BaG_fn|@>=
    @<Prototype for |Find_BaG_fn|@>
{
    double m_r,m_t,deviation;

    RR.slab.b = bcalc2b(x[1])+RR.default_bs;
    if (RR.slab.b<=0 )
        RR.slab.a = 0;
    else
        RR.slab.a = RR.default_bs/RR.slab.b;

    RR.slab.g = gcalc2g(x[2]);

    Calculate_Distance(&m_r,&m_t,&deviation);
    return deviation;
}

@ @<Prototype for |Find_BsG_fn|@>=
double Find_BsG_fn(double x[])

@ @<Definition for |Find_BsG_fn|@>=
    @<Prototype for |Find_BsG_fn|@>
{
    double m_r,m_t,deviation;

    RR.slab.b = bcalc2b(x[1])+RR.default_ba;
    if (RR.slab.b<=0 )
        RR.slab.a = 0;
    else
        RR.slab.a = 1.0 - RR.default_ba/RR.slab.b;

    RR.slab.g = gcalc2g(x[2]);
    Calculate_Distance(&m_r,&m_t,&deviation);
    return deviation;
}

@ Routine to figure out if the light loss exceeds what is
physically possible.  Returns the descrepancy betweent the
current values and the maximum possible values for the the
measurements |m_r| and |m_t|.

@<Prototype for |maxloss|@>=
double maxloss(double f)

@ @<Definition for |maxloss|@>=
    @<Prototype for |maxloss|@>
{
    struct measure_type m_old;
    struct invert_type r_old;
    double m_r, m_t, deviation;

    Get_Calc_State(&m_old, &r_old);

    RR.slab.a = 1.0;
    MM.ur1_lost *= f;
    MM.ut1_lost *= f;

    Calculate_Distance(&m_r,&m_t,&deviation);

    Set_Calc_State(m_old, r_old);
    deviation = ( (MM.m_r + MM.m_t) - ( m_r + m_t ) );

    return deviation;
}

@ This checks the two light loss values |ur1_loss| and |ut1_loss|
to see if they exceed what is physically possible.  If they do, then
these values are replaced by a couple that are the maximum possible
for the current values in |m| and |r|.

@<Prototype for |Max_Light_Loss|@>=
void Max_Light_Loss(struct measure_type m, struct invert_type r,
                    double *ur1_loss, double *ut1_loss)

@ @<Definition for |Max_Light_Loss|@>=
    @<Prototype for |Max_Light_Loss|@>
{
    struct measure_type m_old;
    struct invert_type r_old;

    *ur1_loss = m.ur1_lost;
    *ut1_loss = m.ut1_lost;

    if (Debug(DEBUG_LOST_LIGHT))
        fprintf(stderr, "\nlost before ur1=%7.5f, ut1=%7.5f\n", *ur1_loss, *ut1_loss);

    Get_Calc_State(&m_old, &r_old);

    Set_Calc_State(m, r);

    if (maxloss(1.0) * maxloss(0.0) < 0) {
        double frac;
        frac = zbrent(maxloss,0.00,1.0,0.001);

        *ur1_loss = m.ur1_lost * frac;
        *ut1_loss = m.ut1_lost * frac;
    }

    Set_Calc_State(m_old, r_old);
    if (Debug(DEBUG_LOST_LIGHT))
        fprintf(stderr, "lost after  ur1=%7.5f, ut1=%7.5f\n", *ur1_loss, *ut1_loss);
}

@ this is currently unused

@<Unused diffusion fragment@>=
typedef struct {
  double f;
  double aprime;
  double bprime;
  double gprime;
  double boundary_method;
  double n_top;
  double n_bottom;
  double slide_top;
  double slide_bottom;
  double F0;
  double depth;
  double Exact_coll_flag;
} slabtype;
@#
static void DE_RT(int nfluxes, AD_slab_type slab,
double *UR1, double *UT1, double *URU, double *UTU)
{
    slabtype s;
    double rp, tp, rs, ts;

    s.f = slab.g * slab.g;
    s.gprime = slab.g / (1 + slab.g);
    s.aprime = (1 - s.f) * slab.a / (1 - slab.a * s.f);
    s.bprime = (1 - slab.a * s.f) * slab.b;
    s.boundary_method = Egan;
    s.n_top = slab.n_slab;
    s.n_bottom = slab.n_slab;
    s.slide_top = slab.n_top_slide;
    s.slide_bottom = slab.n_bottom_slide;
    s.F0 = 1 / M_PI;
    s.depth = 0.0;
    s.Exact_coll_flag = FALSE;
    if (MM.illumination == collimated) {
        compute_R_and_T(&s, 1.0, &rp, &rs, &tp, &ts);
        *UR1 = rp + rs;
        *UT1 = tp + ts;
        *URU = 0.0;
        *UTU = 0.0;
        return;
    }
    quad_Dif_Calc_R_and_T(&s, &rp, &rs, &tp, &ts);
    *URU = rp + rs;
    *UTU = tp + ts;
    *UR1 = 0.0;
    *UT1 = 0.0;
}
