@*1 IAD Calculation.

\def\rdirect{r_s^{\hbox{\sevenrm{} direct}}}
\def\tdirect{t_s^{\hbox{\sevenrm{} direct}}}
\def\rdiffuse{r_s^{\hbox{\sevenrm{}}}}
\def\tdiffuse{t_s^{\hbox{\sevenrm{}}}}
\def\std{{\hbox{\sevenrm{}std}}}


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
#define GRID_SIZE 101
#define T_TRUST_FACTOR 2

static int CALCULATING_GRID = 1;
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


@*2 Initialization.

The functions in this file assume that the local variables |MM| and |RR|
have been initialized appropriately.  The variable |MM| contains all the 
information about how a particular experiment was done.  The structure
|RR| contains the data structure that is passed to the adding-doubling routines
as well as the number of quadrature points.  

history
6/8/94 changed error output to |stderr|.

@*2 Gain.

Assume that a sphere is illuminated with diffuse light having a power
$P$. 
This light can reach all parts of sphere --- specifically, light from
this source is not blocked by a baffle.  Multiple reflections in the
sphere will increase the power falling on non-white areas in the sphere
(e.g., the sample, detector, and entrance)  To find the total we need to
sum all the total of all incident light at a point.  The first incidence
is
$$
P_w^{(1)} = a_w P, \qquad P_s^{(1)} = a_s P, \qquad P_d^{(1)} = a_d P
$$
The light from the detector and sample is multiplied by $(1-a_e)$ and
not by $a_w$ because the light from the detector (and sample) is not
allowed to hit either the detector or sample.  The second incidence on
the wall is
$$
P_w^{(2)} = a_w r_w P_w^{(1)} + (1-a_e) r_d P_d^{(1)} + (1-a_e)\rdiffuse
 P_s^{(1)}
$$
The light that hits the walls after $k$ bounces has the same form as
above 
$$
P_w^{(k)} = a_w r_w P_w^{(k-1)} + (1-a_e) r_d  P_d^{(k-1)} + (1-a_e)
\rdiffuse P_s^{(k-1)}
$$
Since the light falling on the sample and detector must come from the
wall
$$
P_s^{(k)} = a_s r_w P_w^{(k-1)} \qquad\hbox{and}\qquad  P_d^{(k)} = a_d
r_w P_w^{(k-1)},
$$
Therefore,
$$
P_w^{(k)} = a_w r_w P_w^{(k-1)} + (1-a_e) r_w (a_d r_d + a_s \rdiffuse)
P_w^{(k-2)}
$$
The total power falling on the walls is just
$$
P_w = \sum_{k=1}^\infty P_w^{(k)} = {a_w + (1-a_e) ( a_d r_d +a_s
\rdiffuse) 
\over 1- a_w r_w - (1-a_e)r_w(a_d r_d + a_s \rdiffuse)} P
$$
The total power falling the detector is
$$
P_d = a_d P + \sum_{k=2}^\infty a_d r_w P_w^{(k-1)} = a_d P + a_d r_w
P_w
$$

The gain $G(\rdiffuse)$ on the irradiance on the detector (relative to a
black sphere),
$$
G(\rdiffuse) \equiv {P_d /A_d\over P/A} 
$$
in terms of the sphere parameters
$$
G(\rdiffuse)  = 1+ {1\over a_w} \cdot {a_w r_w + (1-a_e) r_w (a_d r_d
+a_s \rdiffuse)  \over 1- a_w r_w - (1-a_e)r_w(a_d r_d + a_s \rdiffuse)}
$$
The gain for a detector in a transmission sphere is similar, but with
primed parameters to designate a second potential sphere that is used.
For a black sphere the gain $G(0) = 1$, which is easily verified by
setting $r_w=0$, $\rdiffuse=0$, and $r_d=0$.  Conversely, when the
sphere walls and sample are perfectly white, the irradiance at the
entrance port, the sample port, and the detector port must increase so
that the total power leaving via these ports is equal to the incident
diffuse power $P$.  
Thus the gain should be the ratio of the sphere wall area over the area
of the ports through which light leaves or $G(1)=A/(A_e+A_d)$ which
follows immediately from the gain formula with $r_w=1$, $\rdiffuse=1$,
and $r_d=0$.


@ The gain $G(r_s)$ on the irradiance on the detector (relative to a black sphere),
$$
G(r_s) \equiv {P_d /A_d\over P/A} 
$$
in terms of the sphere parameters
$$
G(r_s)  = 1+ {a_w r_w + (1-a_e) r_w (a_d r_d +a_s r_s) 
\over 1- a_w r_w - (1-a_e)r_w(a_d r_d + a_s r_s)}
$$

@<Prototype for |Gain|@>=
double Gain(int sphere, struct measure_type m, double URU)

@ @<Definition for |Gain|@>=
    @<Prototype for |Gain|@>
{
double G, tmp;

if (sphere==REFLECTION_SPHERE) {
    tmp = m.rw_r*(m.aw_r + (1 - m.ae_r) * (m.ad_r * m.rd_r + m.as_r * URU));
    if (tmp == 1.0) 
        G = 1;
    else
        G = 1.0 + tmp / (1.0 - tmp);

} else {
    tmp = m.rw_t*(m.aw_t + (1 - m.ae_t) * (m.ad_t * m.rd_t + m.as_t * URU));
    if (tmp == 1.0) 
        G = 1;
    else
        G = 1.0 + tmp / (1.0 - tmp);
}

return G;
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
{G(r_s) \over 1-a_s a_s' r_w r_w' (1-a_e)(1-a_e') G(r_s) G'(r_s)t_s^2  }
$$


@<Prototype for |Gain_11|@>=
double Gain_11(struct measure_type m, double URU, double tdiffuse)

@ @<Definition for |Gain_11|@>=
    @<Prototype for |Gain_11|@>
{
    double G, GP, G11;

    G        = Gain(REFLECTION_SPHERE, m, URU );
    GP       = Gain(TRANSMISSION_SPHERE, m, URU );

    G11 = G / (1-m.as_r * m.as_t * m.aw_r * m.aw_t * (1-m.ae_r) * (1-m.ae_t)
               * G * GP * tdiffuse * tdiffuse);
    
    return G11;
}

@ Similarly, when the light starts in the second sphere, the gain for light 
on the detector in the second sphere $G_{2\rightarrow2}$ is found by switching 
all primed variables to unprimed.  Thus $G_{2\rightarrow1}(r_s,t_s)$ is
$$
G_{2\rightarrow2}(r_s,t_s) = {G'(r_s) \over 1-a_s a_s' r_w r_w' 
                              (1-a_e)(1-a_e') G(r_s) G'(r_s)t_s^2  }
$$

@<Prototype for |Gain_22|@>=
double Gain_22(struct measure_type m, double URU, double tdiffuse)

@ @<Definition for |Gain_22|@>=
    @<Prototype for |Gain_22|@>
{
    double G, GP, G22;

    G        = Gain(REFLECTION_SPHERE, m, URU );
    GP       = Gain(TRANSMISSION_SPHERE, m, URU );

    G22 = GP / (1-m.as_r * m.as_t * m.aw_r * m.aw_t * (1-m.ae_r) * (1-m.ae_t)
               * G * GP * tdiffuse * tdiffuse);
    
    return G22;
}

@ The reflected power for two spheres makes use of the formulas for
|Gain_11| above.

The light
on the detector in the reflection (first) sphere arises from three
sources: the fraction of light directly reflected off the sphere wall $f
r_w^2 (1-a_e) P$, the fraction of light reflected by the sample $(1-f)
\rdirect r_w^2 (1-a_e) P$, and the light transmitted through the sample
$(1-f) \tdirect r_w' (1-a_e') P$,
$$
\eqalign{
R(\rdirect,\rdiffuse,\tdirect,\tdiffuse) 
&= G_{1\rightarrow1}(\rdiffuse,\tdiffuse) \cdot a_d (1-a_e) r_w^2 f  P \cr
&+ G_{1\rightarrow1}(\rdiffuse,\tdiffuse) \cdot a_d (1-a_e) r_w (1-f) \rdirect  P \cr
&+ G_{2\rightarrow1}(\rdiffuse,\tdiffuse) \cdot a_d (1-a_e') r_w' (1-f) \tdirect  P \cr
}
$$
which simplifies slightly to
$$
\eqalign{
R(\rdirect,\rdiffuse,\tdirect,\tdiffuse) 
&= a_d (1-a_e) r_w P \cdot G_{1\rightarrow1}(\rdiffuse,\tdiffuse) \cr
&\times \bigg[(1-f)\rdirect + f r_w + 
(1-f)a_s'(1-a_e')r_w'\tdirect\tdiffuse G'(\rdiffuse)\bigg] \cr
}
$$

@<Prototype for |Two_Sphere_R|@>=
double Two_Sphere_R(struct measure_type m, 
                    double UR1, double URU, double UT1, double UTU)

@ @<Definition for |Two_Sphere_R|@>=
    @<Prototype for |Two_Sphere_R|@>
{
    double x, GP;
    GP = Gain(TRANSMISSION_SPHERE, m, URU );

    x = m.ad_r*(1-m.ae_r)*m.rw_r*Gain_11(m,URU,UTU);
    x *= (1-m.f_r)*UR1+m.rw_r*m.f_r+(1-m.f_r)*m.as_t*(1-m.ae_t)*m.rw_t*UT1*UTU*GP;
    return x;
}

@ For the power on the detector in the transmission (second) sphere we
have the same three sources.  The only difference is that the subscripts
on the gain terms now indicate that the light ends up in the second
sphere
$$
\eqalign{
T(\rdirect,\rdiffuse,\tdirect,\tdiffuse) 
&= G_{1\rightarrow2}(\rdiffuse,\tdiffuse)  \cdot a_d' (1-a_e) r_w^2 f P \cr
&+ G_{1\rightarrow2}(\rdiffuse,\tdiffuse) \cdot a_d' (1-a_e) r_w (1-f) \rdirect P  \cr
&+ G_{2\rightarrow2}(\rdiffuse,\tdiffuse) \cdot a_d' (1-a_e') r_w' (1-f) \tdirect  P \cr
}
$$
or 
$$
\eqalign{
T(\rdirect,\rdiffuse,\tdirect,\tdiffuse) 
&= a_d' (1-a_e') r_w' P\cdot G_{2\rightarrow2}(\rdiffuse,\tdiffuse) \cr
&\times \bigg[(1-f)\tdirect 
 + (1-a_e)r_w a_s \tdiffuse (f r_w+(1-f) \rdirect)G(\rdiffuse) \bigg] \cr
}
$$

@<Prototype for |Two_Sphere_T|@>=
double Two_Sphere_T(struct measure_type m, 
                    double UR1, double URU, double UT1, double UTU)

@ @<Definition for |Two_Sphere_T|@>=
    @<Prototype for |Two_Sphere_T|@>
{
    double x, G;
    G = Gain(REFLECTION_SPHERE, m, URU );
    x = m.ad_t*(1-m.ae_t)*m.rw_t*Gain_22(m,URU,UTU);
    x *= (1-m.f_r)*UT1+(1-m.ae_r)*m.rw_r*m.as_r*UTU*(m.f_r*m.rw_r+(1-m.f_r)*UR1)*G;
    return x;
}

@*2 Grid Routines.
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
    if ( Debug(DEBUG_ITERATIONS) && !CALCULATING_GRID) {
        fprintf(stderr,"UR1 loss=%g, UT1 loss=%g  ", m.ur1_lost, m.ut1_lost);
        fprintf(stderr,"URU loss=%g, UTU loss=%g\n", m.uru_lost, m.utu_lost);
    }
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
boolean_type Valid_Grid(struct measure_type m, search_type s)

@ @<Definition for |Valid_Grid|@>=
    @<Prototype for |Valid_Grid|@>
{   
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
            fprintf(stderr,"GRID: Fill slab refractive index changed\n");
        return(FALSE);
    }
    if (m.slab_cos_angle          != MGRID.slab_cos_angle) {
        if (Debug(DEBUG_GRID)) 
            fprintf(stderr,"GRID: Fill incident light changed\n");
        return(FALSE);
    }

    if (m.slab_top_slide_index    != MGRID.slab_top_slide_index) {
        if (Debug(DEBUG_GRID)) 
            fprintf(stderr,"GRID: Fill top slide refractive index changed\n");
        return(FALSE);
    }

    if (m.slab_bottom_slide_index != MGRID.slab_bottom_slide_index) {
        if (Debug(DEBUG_GRID)) 
            fprintf(stderr,"GRID: Fill bottom slide refractive index changed\n");
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
    double swap, correct_UR1, correct_URU;
    RT(n, slab, UR1, UT1, URU, UTU);
    if (flip) {
        correct_UR1 = *UR1;
        correct_URU = *URU;
        swap=slab->n_top_slide;
        slab->n_top_slide = slab->n_bottom_slide;
        slab->n_bottom_slide = swap;
        swap=slab->b_top_slide;
        slab->b_top_slide = slab->b_bottom_slide;
        slab->b_bottom_slide = swap;
        RT(n, slab, UR1, UT1, URU, UTU);
        swap=slab->n_top_slide;
        slab->n_top_slide = slab->n_bottom_slide;
        slab->n_bottom_slide = swap;
        swap=slab->b_top_slide;
        slab->b_top_slide = slab->b_bottom_slide;
        slab->b_bottom_slide = swap;
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
        fprintf(stderr, "+ %2d %2d ",i,j);
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
    double a;
    double min_b = -8;  /* exp(-10) is smallest thickness */
    double max_b = +8;   /* exp(+8) is greatest thickness */

    if (Debug(Debug(DEBUG_GRID))) 
        fprintf(stderr, "Filling AB grid\n");
        
    if (The_Grid==NULL) Allocate_Grid(r.search);
    @<Zero |GG|@>@;

    Set_Calc_State(m,r);
        
    GG_g = RR.slab.g;
    for(i=0; i<GRID_SIZE; i++){
        double x = (double) i/(GRID_SIZE-1.0);
        
        RR.slab.b = exp(min_b + (max_b-min_b) *x);
        for(j=0; j<GRID_SIZE; j++){
            @<Generate next albedo using j@>@;
            
            fill_grid_entry(i,j);
        }
    }

The_Grid_Initialized=TRUE;
The_Grid_Search = FIND_AB;
}

@ Now it seems that I must be a bit more subtle in choosing the range
of albedos to use in the grid.  Originally I just spaced them according
to
$$
a = 1 - \left[ {j-1\over n-1} \right]^3
$$
where $1\le j\le n$.  Long ago it seems that I based things only on the
square of the bracketed term, but I seem to remember that I was forced to
change it from a square to a cube to get more global convergence.

So why am I rewriting this?  Well, because it works very poorly for samples
with small albedos.  For example, when $n=11$ then the values chosen for
|a| are (1, .999, .992, .973, .936, .875, .784, .657, .488, .271, 0).
Clearly very skewed towards high albedos.  

I am considering a two part division.  I'm not too sure how it should go.
Let the first half be uniformly divided and the last half follow the
cubic scheme given above.  The list of values should then be
(1, .996, .968, .892, 0.744, .5, .4, .3, .2, .1, 0).

Maybe it would be best if I just went back to a quadratic term.  
Who knows?

In the |if| statement below, note that it could read |j>=k| and still
generate the same results.

@<Nonworking code@>=
            k = floor((GRID_SIZE-1) / 2);
            if (j > k) {
                a = 0.5 *(1-(j-k-1)/(GRID_SIZE-k-1));
                RR.slab.a = a;
            } else {
                a = (j-1.0)/(GRID_SIZE-k-1);
                RR.slab.a = 1.0-a*a*a/2;
            }

@ Well, the above code did not work well.  So I futzed around and 
sort of empirically ended up using the very simple method 
below.  The only real difference from the previous method what
that the method is now quadratic and not cubic.

@<Generate next albedo using j@>=
            a = (double) j/(GRID_SIZE-1.0);
            if (a < 0.25)
                RR.slab.a = 1.0-a*a;
            else if (a > 0.75)              
                RR.slab.a = (1.0 - a) * (1.0 -a);
            else 
                RR.slab.a = 1 - a;

@ Well, the above code has gaps.  Here is an attempt to eliminate the gaps

@<Generate next albedo using j@>=
            a = (double) j/(GRID_SIZE-1.0);
            RR.slab.a = (1.0-a*a)*(1.0-a) + (1.0-a)*(1.0-a)*a;

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
double a;

    if (Debug(Debug(DEBUG_GRID)))
        fprintf(stderr, "Filling AG grid\n");

    if (The_Grid==NULL) Allocate_Grid(r.search);
    @<Zero |GG|@>@;
        
    Set_Calc_State(m,r);
    GG_b=r.slab.b;
    for(i=0; i<GRID_SIZE; i++){
        RR.slab.g =0.9999*(2.0*i/(GRID_SIZE-1.0)-1.0);
        for(j=0; j<GRID_SIZE; j++){
        
            @<Generate next albedo using j@>@;
            fill_grid_entry(i,j);
        }
    }

The_Grid_Initialized=TRUE;
The_Grid_Search = FIND_AG;
}

@   @<Zero |GG|@>=
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

    if (The_Grid==NULL) Allocate_Grid(r.search);
    @<Zero |GG|@>@;
    
    if (Debug(Debug(DEBUG_GRID))) 
        fprintf(stderr, "Filling BG grid\n");

    Set_Calc_State(m,r);
    RR.slab.b = 1.0/32.0;
    RR.slab.a = RR.default_a;
    GG_a = RR.slab.a;
    
    for(i=0; i<GRID_SIZE; i++){
        RR.slab.b *=2;
        for(j=0; j<GRID_SIZE; j++){
            RR.slab.g =0.9999*(2.0*j/(GRID_SIZE-1.0)-1.0);
            fill_grid_entry(i,j);
        }
    }

The_Grid_Initialized=TRUE;
The_Grid_Search = FIND_BG;
}

@ This is quite similar to |Fill_BG_Grid|, with the exception of the 
that the $b_s=\mu_s d$ is held fixed.  Here $b$ and $g$ are varied
on the usual grid, but the albedo is forced to take whatever value
is needed to ensure that the scattering constant remains fixed.

@<Prototype for |Fill_BaG_Grid|@>=
void Fill_BaG_Grid(struct measure_type m, struct invert_type r)

@ @<Definition for |Fill_BaG_Grid|@>=
    @<Prototype for |Fill_BaG_Grid|@>
{ 
int i,j;
double bs, ba;

    if (The_Grid==NULL) Allocate_Grid(r.search);
    @<Zero |GG|@>@;

    if (Debug(Debug(DEBUG_GRID))) 
        fprintf(stderr, "Filling BaG grid\n");

    Set_Calc_State(m,r);
    ba = 1.0/32.0;
    bs = RR.default_bs;
    GG_bs = bs;
    for(i=0; i<GRID_SIZE; i++){
        ba *=2;
        ba = exp((double) i/(GRID_SIZE-1.0)*log(1024.0))/16.0;
        RR.slab.b = ba+bs;
        if (RR.slab.b>0) 
            RR.slab.a = bs/RR.slab.b;
        else
            RR.slab.a = 0;
        for(j=0; j<GRID_SIZE; j++){
            RR.slab.g = 0.9999*(2.0*j/(GRID_SIZE-1.0)-1.0);
            fill_grid_entry(i,j);
        }
    }

The_Grid_Initialized=TRUE;
The_Grid_Search = FIND_BaG;
}

@ Very similiar to the above routine.  The value of $b_a=\mu_a d$ is held constant.

@<Prototype for |Fill_BsG_Grid|@>=
void Fill_BsG_Grid(struct measure_type m, struct invert_type r)

@ @<Definition for |Fill_BsG_Grid|@>=
    @<Prototype for |Fill_BsG_Grid|@>
{ 
int i,j;
double bs, ba;

    if (The_Grid==NULL) Allocate_Grid(r.search);
    @<Zero |GG|@>@;
    
    Set_Calc_State(m,r);
    bs = 1.0/32.0;
    ba = RR.default_ba;
    GG_ba = ba;
    for(i=0; i<GRID_SIZE; i++){
        bs *=2;
        RR.slab.b = ba+bs;
        if (RR.slab.b>0) 
            RR.slab.a = bs/RR.slab.b;
        else
            RR.slab.a = 0;
        for(j=0; j<GRID_SIZE; j++){
            RR.slab.g = 0.9999*(2.0*j/(GRID_SIZE-1.0)-1.0);
            fill_grid_entry(i,j);
        }
    }

The_Grid_Initialized=TRUE;
The_Grid_Search = FIND_BsG;
}

@ @<Prototype for |Fill_Grid|@>=
void Fill_Grid(struct measure_type m, struct invert_type r)

@ @<Definition for |Fill_Grid|@>=
    @<Prototype for |Fill_Grid|@>
{ 
    if (!Same_Calc_State(m,r)) {
        switch (r.search) {
            case FIND_AB:
                if (Debug(DEBUG_SEARCH)) fprintf(stderr,"filling AB Grid\n");
                Fill_AB_Grid(m,r);
                break;
            case FIND_AG: 
                if (Debug(DEBUG_SEARCH)) fprintf(stderr,"filling AG Grid\n");
                Fill_AG_Grid(m,r);
                break;
            case FIND_BG: 
                if (Debug(DEBUG_SEARCH)) fprintf(stderr,"filling BG Grid\n");
                Fill_BG_Grid(m,r);
                break;
            case FIND_BaG: 
                if (Debug(DEBUG_SEARCH)) fprintf(stderr,"filling BaG Grid\n");
                Fill_BaG_Grid(m,r);
                break;
            case FIND_BsG: 
                if (Debug(DEBUG_SEARCH)) fprintf(stderr,"filling BsG Grid\n");
                Fill_BsG_Grid(m,r);
                break;
            default:
                AD_error("Attempt to fill grid for unusual search case.");
        }
    }

    Get_Calc_State(&MGRID, &RGRID);
}

@*2 Calculating R and T.


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
    double Rc, Tc, ur1, ut1, uru, utu;

    if (RR.slab.b <= 1e-6 ) RR.slab.b = 1e-6;
    if (Debug(DEBUG_EVERY_CALC)) 
        fprintf(stderr, "a=%8.5f b=%10.5f g=%8.5f ", RR.slab.a, RR.slab.b, RR.slab.g);
                
    RT_Flip(MM.flip_sample, RR.method.quad_pts, &RR.slab, &ur1, &ut1, &uru, &utu);
    
    if (Debug(DEBUG_EVERY_CALC)) 
        fprintf(stderr, "ur1=%8.5f ut1=%8.5f (not M_R and M_T!)\n", ur1, ut1);

    Sp_mu_RT_Flip(MM.flip_sample, 
             RR.slab.n_top_slide, RR.slab.n_slab, RR.slab.n_bottom_slide, 
             RR.slab.b_top_slide, RR.slab.b,      RR.slab.b_bottom_slide, 
             RR.slab.cos_angle, &Rc, &Tc);

    if ((!CALCULATING_GRID && Debug(DEBUG_ITERATIONS)) ||
        ( CALCULATING_GRID && Debug(DEBUG_GRID_CALC)))
            fprintf(stderr, "        ");
            
    Calculate_Distance_With_Corrections(ur1,ut1,Rc,Tc,uru,utu,M_R,M_T,deviation);
}

@ @<Prototype for |Calculate_Grid_Distance|@>=
double Calculate_Grid_Distance(int i, int j)

@ @<Definition for |Calculate_Grid_Distance|@>=
    @<Prototype for |Calculate_Grid_Distance|@>
{
    double ur1,ut1,uru,utu,Rc,Tc,b,dev,LR,LT;
    
    if (Debug(DEBUG_GRID_CALC)) 
        fprintf(stderr, "g %2d %2d ",i,j);
        
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
             RR.slab.cos_angle, &Rc, &Tc);

    CALCULATING_GRID = 1;
    Calculate_Distance_With_Corrections(ur1,ut1,Rc,Tc,uru,utu,&LR,&LT,&dev);
    CALCULATING_GRID = 0;
    
    return dev;
}

@ This is the routine that actually finds the distance.  I have factored
this part out so that it can be used in the |Near_Grid_Point| routine.

|Rc| and |Tc| refer to the unscattered (collimated) reflection and transmission.  

The only tricky part is to remember that the we are trying to match the
measured values.  The measured values are affected by sphere parameters
and light loss.  Since the values |UR1| and |UT1| are for an infinite slab
sample with no light loss, the light loss out the edges must be subtracted.
It is these values that are used with the sphere formulas to convert the
modified |UR1| and |UT1| to values for |*M_R| and |*M_T|.

@<Prototype for |Calculate_Distance_With_Corrections|@>=
void Calculate_Distance_With_Corrections(
                double UR1, double UT1, 
                double Rc,  double Tc, 
                double URU, double UTU, 
                double *M_R, double *M_T, double *dev)

@ @<Definition for |Calculate_Distance_With_Corrections|@>=
    @<Prototype for |Calculate_Distance_With_Corrections|@>
{
    double R_direct, T_direct, R_diffuse, T_diffuse;

    R_diffuse = URU - MM.uru_lost;
    T_diffuse = UTU - MM.utu_lost;
    
    R_direct = UR1 - MM.ur1_lost - (1.0 - MM.fraction_of_rc_in_mr) * Rc;
    T_direct = UT1 - MM.ut1_lost - (1.0 - MM.fraction_of_tc_in_mt) * Tc;

    switch (MM.num_spheres) {
        case 0:
            @<Calc |M_R| and |M_T| for no spheres@>@;
            break;

        case 1:
        case -2:
            if (MM.method == COMPARISON)
                @<Calc |M_R| and |M_T| for dual beam sphere@>@;
            else
                @<Calc |M_R| and |M_T| for single beam sphere@>@;
            break;

        case 2:
            @<Calc |M_R| and |M_T| for two spheres@>@;
            break;
    }

    @<Calculate the deviation@>@;
    @<Print diagnostics@>@;
}

@ If no spheres were used in the measurement, then presumably the
measured values are the reflection and transmission.  Consequently,
we just acertain what the irradiance was and whether the
specular reflection ports were blocked and proceed accordingly.
Note that blocking the ports does not have much meaning unless
the light is collimated, and therefore the reflection and
transmission is only modified for collimated irradiance.

@<Calc |M_R| and |M_T| for no spheres@>=
    *M_R = R_direct;
    *M_T = T_direct;
        
@ The direct incident power is $(1-f) P$. The reflected power will 
be $(1-f)\rdirect P$.  Since baffles ensure that the light cannot
reach the detector, we must bounce the light off the sphere walls to
use to above gain formulas.  The contribution will then be 
$(1-f) \rdirect (1-a_e) r_w P$.  The measured power will be 
$$
P_d = a_d (1-a_e) r_w [(1-f) \rdirect + f r_w] P \cdot G(r_s)
$$
Similarly the power falling on the detector measuring transmitted light is
$$
P_d'= a_d' \tdirect r_w'  (1-a_e')  P \cdot G'(r_s)
$$
when the `entrance' port in the transmission sphere is closed, $a_e'=0$.

The normalized sphere measurements are
$$
M_R = r_\std\cdot{R(\rdirect,r_s)-R(0,0) \over R(r_\std,r_\std)-R(0,0)}
$$
and
$$
M_T = t_\std\cdot{T(\tdirect,r_s)-T(0,0) \over T(t_\std,r_\std)-T(0,0)}
$$

@<Calc |M_R| and |M_T| for single beam sphere@>=
{
    double P_std, P_d, P_0;
    double G, G_0, G_std, GP_std, GP;
    
    G_0      = Gain(REFLECTION_SPHERE, MM, 0.0);
    G        = Gain(REFLECTION_SPHERE, MM, R_diffuse);
    G_std    = Gain(REFLECTION_SPHERE, MM, MM.rstd_r);
    
    P_d      = G     * (R_direct  * (1-MM.f_r) + MM.f_r*MM.rw_r);
    P_std    = G_std * (MM.rstd_r * (1-MM.f_r) + MM.f_r*MM.rw_r);
    P_0      = G_0   * (                         MM.f_r*MM.rw_r);   
    *M_R     = MM.rstd_r * (P_d - P_0)/(P_std - P_0);

    GP       = Gain(TRANSMISSION_SPHERE, MM, R_diffuse);
    GP_std   = Gain(TRANSMISSION_SPHERE, MM, 0.0);
    *M_T     = T_direct * GP / GP_std;
}

@ The dual beam case is different because the sphere efficiency
is equivalent for measurement of light hitting the sample first or
hitting the reference standard first.  The dual beam measurement
should report the ratio of these two reflectance measurements, thereby
eliminating the need to calculate the gain completely.   The same
holds when no sample is present.

The normalized reflectance measurement (the difference between
dual beam measurement for a port with the sample and with nothing)
is
$$
M_R = r_\std\cdot{(1-f)\rdirect  + f r_w\over (1-f')r_\std -f' r_w}
    - r_\std\cdot{(1-f) (0)  + f r_w\over (1-f')r_\std -f' r_w}
$$
or
$$
M_R = {(1-f)\rdirect \over (1-f') -f' r_w/r_\std}
$$
When $f=f'=1$, then $M_R=1$ no matter what the reflectance is.
(Leave it in this form to avoid division by zero when $f=1$.)

The normalized transmittance is simply $\tdirect$.

When $f=0$ then this result is essentially the same as the
no spheres result (because no sphere corrections are needed).
However if the number of spheres is zero, then no lost light
calculations are made and therefore that is a potential error.

@<Calc |M_R| and |M_T| for dual beam sphere@>=
{
    *M_R     = (1-MM.f_r) * R_direct/((1-MM.f_r) + MM.f_r*MM.rw_r/MM.rstd_r);
    *M_T     = T_direct;
}

@ When two integrating spheres are present then the
double integrating sphere formulas are slightly more complicated.

I am not sure what it means when |rstd_t| is not unity.

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
set the fraction |m.f_r| ore |m.f_t| to be non-zero.

@<Calc |M_R| and |M_T| for two spheres@>=
{
double R_0, T_0;
R_0 = Two_Sphere_R(MM, 0, 0, 0, 0);
T_0 = Two_Sphere_T(MM, 0, 0, 0, 0);

*M_R = MM.rstd_r * (Two_Sphere_R(MM, R_direct, R_diffuse, T_direct, T_diffuse) - R_0)/
                   (Two_Sphere_R(MM, MM.rstd_r, MM.rstd_r, 0, 0) - R_0);
*M_T =  (Two_Sphere_T(MM, R_direct, R_diffuse, T_direct, T_diffuse) - T_0)/
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
        *dev = 0;
        if (MM.m_t > ABIT)
            *dev = T_TRUST_FACTOR* fabs(MM.m_t - *M_T) / (MM.m_t + ABIT);
        if ( RR.default_a != 0 )
            *dev += fabs(MM.m_r - *M_R) / (MM.m_r + ABIT);          
    } else {
        *dev = T_TRUST_FACTOR * fabs(MM.m_t - *M_T);
        if ( RR.default_a != 0 )
            *dev += fabs(MM.m_r - *M_R);            
    } 


@ This is here so that I can figure out why the program is not converging.
  This is a little convoluted so that the global constants at the top of
  this file interact properly.  
  
@<Print diagnostics@>=
if ((Debug(DEBUG_ITERATIONS) && !CALCULATING_GRID) ||
    (Debug(DEBUG_GRID_CALC)  &&  CALCULATING_GRID)) {
    static int once = 0;
        
    if (once==0) {
            fprintf(stderr, "%10s %10s %10s |%10s %10s |%10s %10s |%10s\n",
                    "a", "b", "g", "m_r", "fit", "m_t", "fit", "delta");
            once = 1;
        }
    
    fprintf(stderr, "%10.5f %10.5f %10.5f |", RR.slab.a, RR.slab.b, RR.slab.g);
    fprintf(stderr, "%10.5f %10.5f |", MM.m_r, *M_R);
    fprintf(stderr, "%10.5f %10.5f |", MM.m_t, *M_T);
    fprintf(stderr, "%10.5f \n", *dev);
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
    s.F0 = 1 / pi;
    s.depth = 0.0;
    s.Exact_coll_flag = false;
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
