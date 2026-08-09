@** AD Cone.
This file provides routines to obtain reflection and transmission values 
returning within a cone assuming normal illumination.

@(ad_cone.c@>=
#include <math.h>
#include <float.h>
#include <stdio.h>
#include "nr_util.h"
#include "ad_globl.h"
#include "ad_matrx.h"
#include "ad_bound.h"
#include "ad_doubl.h"
#include "ad_start.h"
#include "ad_frsnl.h"

@<Definition for |RT_Cone|@>@;
@<Definition for |ez_RT_Cone|@>@;
@<Definition for |ez_RT_Oblique|@>@;

@ @(ad_cone.h@>=
    @h
    @<Prototype for |RT_Cone|@>;
    @<Prototype for |ez_RT_Cone|@>;
    @<Prototype for |ez_RT_Oblique|@>;

@ @(ad_cone_ez.h@>=
    @h
    @<Prototype for |ez_RT_Cone|@>;
    @<Prototype for |ez_RT_Oblique|@>;

@*1 RT Cone.
Sometimes you just need to know the total reflection and transmission from a
target within a specified cone of angles.  For example, you might want to test
a Monte Carlo implementation of fiber illumination.  The way that this works is
to divide the integration over angles into two or three pieces.  A separate 
quadrature is done over each integration range.  For example if $\nu_{\hbox{cone}}$
is the cosine of the cone angle and there are no index of refraction changes that
need to accounted for, then 
$$
\int_0^1A(\nu,\nu') B(\nu',\nu'')\,d\nu'  =
          \int_0^{\nu_{\rm cone}}A(\nu,\nu') B(\nu',\nu'')\,d\nu' +
             \int_{\nu_{\rm cone}}^1 A(\nu,\nu') B(\nu',\nu'')\,d\nu' .
$$
otherwise one needs to include the critical angle as a special point in the
integration and the integration becomes
$$
\eqalign{\int_0^1A(\nu,\nu') B(\nu',\nu'')\,d\nu'  &=
          \int_0^{\nu_{\rm crit}}A(\nu,\nu') B(\nu',\nu'')\,d\nu' \cr
          &+
          \int_{\nu_{\rm crit}}^{\nu_{\rm cone}}A(\nu,\nu') B(\nu',\nu'')\,d\nu' +
             \int_{\nu_{\rm cone}}^1 A(\nu,\nu') B(\nu',\nu'')\,d\nu' .\cr}
$$
Radau quadrature is chosen for the integration range from $\nu_{\hbox{cone}}$ to
1.  The other two use Gaussian quadrature.  

@<Prototype for |RT_Cone|@>=
void RT_Cone(int n,     
                struct AD_slab_type * slab,
                int use_cone,
                double *UR1, double *UT1, double *URU, double *UTU)

@ @<Definition for |RT_Cone|@>=
    @<Prototype for |RT_Cone|@>
{
    @<|RT_Cone| Declare variables@>@;
    @<|RT_Cone| Check inputs@>@;
    Choose_Cone_Method(slab, &method);
    @<|RT_Cone| Handle a slab of zero thickness@>@;
    @<|RT_Cone| Allocate slab memory@>@;
    @<|RT_Cone| Initialize homogeneous layer@>@;
    @<|RT_Cone| Allocate and generate top and bottom boundaries@>@;
    @<|RT_Cone| Add top and bottom boundaries@>@;
    @<|RT_Cone| Free memory@>@;
}

@ @<|RT_Cone| Declare variables@>=
    struct AD_method_type method;
    double *R01, *R10, *T01, *T10;
    double *R23, *R32, *T23, *T32;
    double **R12, **T12;
    double **R02, **T02, **T20, **R20;
    double **R03, **T03, **T30, **R30;
    double **atemp, **btemp;
    double d;
    *UR1=-1;
    *URU=-1;
    *UT1=-1;
    *UTU=-1;

@       @<|RT_Cone| Check inputs@>=
    if (slab->n_slab<0) return;
    if (slab->n_top_slide<0) return;
    if (slab->n_bottom_slide<0) return;
    if (slab->a<0  || slab->a>1) return;
    if (slab->g<-1 || slab->g>1) return;
    if (slab->b<0) return;
    if (slab->cos_angle<0 || slab->cos_angle>1) return;

@ The number of quadrature points must be fixed before
starting to allocate memory.  We want the number of points
to be at least twelve so that each of the three integrals 
will have four quadrature points.  

@<|RT_Cone| Check inputs@>=
    n = 12 * (n / 12);
    if (n < 12) n = 12;
    method.quad_pts = n;

@ @<|RT_Cone| Allocate slab memory@>=
    R12 = dmatrix(1, n, 1, n);
    T12 = dmatrix(1, n, 1, n);  
    R02 = dmatrix(1, n, 1, n);
    T02 = dmatrix(1, n, 1, n);  
    R20 = dmatrix(1, n, 1, n);
    T20 = dmatrix(1, n, 1, n);  
    R03 = dmatrix(1, n, 1, n);
    T03 = dmatrix(1, n, 1, n);  
    R30 = dmatrix(1, n, 1, n);
    T30 = dmatrix(1, n, 1, n);  
    atemp = dmatrix(1, n, 1, n);
    btemp = dmatrix(1, n, 1, n);

@ The homogeneous layer initially has 0\% reflection and 100\% transmission.
We cannot fob the details on how this layer is created to |RT_Matrices| 
because we need to (1) set the quadrature angles to a multiple of three,
and (2) explicitly make a call to |Choose_Cone_Method| so that the quadrature angles
will get chosen appropriately.

This code is directly lifted from the |RT_Matrices| routine.

@<|RT_Cone| Initialize homogeneous layer@>=

    n = method.quad_pts;
    Init_Layer(*slab, method, R12, T12);

    d= 1.0;
    if (slab->b != HUGE_VAL)
        d = method.b_thinnest * slab->b / method.b_calc;

    Double_Until(n, R12, T12, d, slab->b);

@ A slab of zero optical thickness is perfectly legitimate --- a bare slide,
or a window whose index matches what surrounds it.  The adding machinery
cannot express one, though.  A layer with no thickness reflects and transmits
specularly, and a specular layer is a delta function in angle, which no
quadrature matrix represents.  Pushing it through anyway produced a singular
matrix and killed the program inside |Left_Inverse_Multiply| the moment the
slab and its boundaries had different indices.

This is why |RT| never sends such a slab through the matrices at all; it hands
the case to |Sp_RT|.  The code here once tried to build a |Zero_Layer| and
|return|, which was a fragment lifted from |RT_Matrices|.  That |return| is
right in |RT_Matrices|, where filling |R| and |T| is all the routine ever
does, and wrong here, where the boundaries have still to be applied: it threw
the layer away the instant it was built and left the four fluxes holding the
$-1$ that marks refused input.  So a slab that merely had no thickness was
reported as one that could not be understood --- and every matrix leaked on
the way out.

Oblique incidence needs nothing new.  There |cos_angle| is the direction the
light arrives from rather than the width of a cone, which is precisely what
|Sp_RT| already computes.

For a cone the delta function is what makes the answer easy.  Normally
incident light leaves a specular slab along the normal, so the whole beam is
caught by any cone of nonzero width and none of it by a cone of no width:
|UR1| is all or nothing.  The diffuse quantity follows from the definition
$$
\hbox{URU} \equiv {n^2\over1-\mu^2} \int_\mu^1\!\!\int_0^1
R(\nu',\nu'')\,2\nu'd\nu'\,2\nu''d\nu''
$$
by collapsing one integral against the delta, which leaves a single sum over
the quadrature angles that lie inside the cone.  Angles past the critical
angle are trapped and contribute nothing, which is why |Cos_Snell| returning
zero is skipped exactly as |Sp_RT| skips it.  At |mu| of zero the cone opens
to the whole hemisphere, every untrapped angle is included, $1-\mu^2$ is one,
and the sum reduces term by term to |Sp_RT| --- so the full cone agrees with
|RT| by construction rather than by luck.

@<|RT_Cone| Handle a slab of zero thickness@>=

    if (slab->b <= 0) {
        double mu = slab->cos_angle;
        double mu_slab, mu_outside, r, t;
        int i;

        if (use_cone != CONE) {
            Sp_RT(n, *slab, UR1, UT1, URU, UTU);
            return;
        }

        if (mu >= 1) {
            *UR1 = 0;
            *UT1 = 0;
            *URU = 0;
            *UTU = 0;
            return;
        }

        Sp_mu_RT(slab->n_top_slide, slab->n_slab, slab->n_bottom_slide,
                 slab->b_top_slide, 0.0, slab->b_bottom_slide, 1.0, UR1, UT1);

        if (slab->n_slab == 1)
            mu_slab = mu;
        else
            mu_slab = sqrt(slab->n_slab*slab->n_slab - 1 + mu*mu)/slab->n_slab;

        *URU = 0.0;
        *UTU = 0.0;
        for (i = 1; i <= n; i++) {
            if (angle[i] <= mu_slab)
                continue;
            mu_outside = Cos_Snell(slab->n_slab, angle[i], 1.0);
            if (mu_outside == 0)
                continue;
            Sp_mu_RT(slab->n_top_slide, slab->n_slab, slab->n_bottom_slide,
                     slab->b_top_slide, 0.0, slab->b_bottom_slide, mu_outside, &r, &t);
            *URU += twoaw[i] * r;
            *UTU += twoaw[i] * t;
        }

        *URU *= slab->n_slab*slab->n_slab/(1 - mu*mu);
        *UTU *= slab->n_slab*slab->n_slab/(1 - mu*mu);
        return;
    }

@ Create the matrices needed for the top and bottom
@<|RT_Cone| Allocate and generate top and bottom boundaries@>=
    R01 = dvector(1, n);
    R10 = dvector(1, n);
    T01 = dvector(1, n);
    T10 = dvector(1, n);
    Init_Boundary(*slab, n, R01, R10, T01, T10, TOP_BOUNDARY);

    R23 = dvector(1, n);
    R32 = dvector(1, n);
    T23 = dvector(1, n);
    T32 = dvector(1, n);
    Init_Boundary(*slab, n, R23, R32, T23, T32, BOTTOM_BOUNDARY);


@ Here the layer numbering is pretty consistent.  The top slide is 01, the
scattering layer is 12, and the bottom slide is 23.  Light going from the top
of the slide to the bottom of the scattering layer is 02 and similarly light
going all the way through is 03.

The only tricky part is that the definitions of |UR1| and |URU| have changed
from their usual definitions.  When |use_cone==OBLIQUE| then |UR1| refers to 
the light reflected back
into the specified cone for normal irradiance and |URU| is for light reflected
back into the cone for light incident uniformly at all angles within that cone.
Otherwise, assume that the incidence is oblique.  |UR1| then refers to the 
total amount of light reflected back for light incident only at the cone angle.

@<|RT_Cone| Add top and bottom boundaries@>=

    Add_Top   (n, R01, R10, T01, T10, R12, R12, T12, T12, R02, R20, T02, T20, atemp, btemp);
    Add_Bottom(n, R02, R20, T02, T20, R23, R32, T23, T32, R03, R30, T03, T30, atemp, btemp);
    
    if (use_cone==CONE) {
        URU_and_UR1_Cone(n, slab->n_slab, slab->cos_angle, R03, URU, UR1);
        Transpose_Matrix(n,T03);
        URU_and_UR1_Cone(n, slab->n_slab, slab->cos_angle, T03, UTU, UT1);
    } else {
        {
        double unused;
        if (use_cone!=OBLIQUE) 
    		fprintf(stderr,"Unknown type for use_cone.  Assuming oblique incidence.\n");
        URU_and_URx_Cone(n, slab->n_slab, slab->cos_angle, R03, URU, UR1);
        URU_and_UR1(n, slab->n_slab, R03, URU, &unused);
        Transpose_Matrix(n,T03);
        URU_and_URx_Cone(n, slab->n_slab, slab->cos_angle, T03, UTU, UT1);
        URU_and_UR1(n, slab->n_slab, T03, UTU, &unused);
        }
        
    } 

@ @<|RT_Cone| Free memory@>=
    free_dvector(R01, 1, n);
    free_dvector(R10, 1, n);
    free_dvector(T01, 1, n);
    free_dvector(T10, 1, n);

    free_dmatrix(R12, 1, n, 1, n);
    free_dmatrix(T12, 1, n, 1, n);

    free_dmatrix(R03, 1, n, 1, n);
    free_dmatrix(R30, 1, n, 1, n);
    free_dmatrix(T03, 1, n, 1, n);
    free_dmatrix(T30, 1, n, 1, n);

    free_dmatrix(R02, 1, n, 1, n);
    free_dmatrix(R20, 1, n, 1, n);
    free_dmatrix(T02, 1, n, 1, n);
    free_dmatrix(T20, 1, n, 1, n);

    free_dmatrix(atemp, 1, n, 1, n);
    free_dmatrix(btemp, 1, n, 1, n);

    free_dvector(R32, 1, n);
    free_dvector(R23, 1, n);
    free_dvector(T32, 1, n);
    free_dvector(T23, 1, n);

@ Simple wrapper that avoids data structures.  As with the other |ez_RT|
routines the slides refract but do not absorb: no slide optical depth can be
passed and both are set to zero below.  Call |RT_Cone| with a filled-in
|AD_slab_type| when the slides absorb.

@<Prototype for |ez_RT_Cone|@>=
void ez_RT_Cone(int n,  
                double nslab, 
                double ntopslide, 
                double nbottomslide, 
                double a,
                double b,
                double g,
                double cos_cone_angle,
                double *UR1, double *UT1, double *URU, double *UTU)


@ @<Definition for |ez_RT_Cone|@>=
    @<Prototype for |ez_RT_Cone|@>
{
struct AD_slab_type slab;

    slab.n_slab = nslab;
    slab.n_top_slide = ntopslide;
    slab.n_bottom_slide = nbottomslide;
    slab.b_top_slide = 0;
    slab.b_bottom_slide = 0;
    slab.a = a;
    slab.b = b;
    slab.g = g;
    slab.cos_angle = cos_cone_angle;
    slab.phase_function = HENYEY_GREENSTEIN;

    RT_Cone(n, &slab, CONE, UR1, UT1, URU, UTU);
}

@ This routine calculates reflection and transmission
for oblique incidence.  |URx| and |UTx| are the total light
reflected and transmitted for light incident at |cos_oblique_angle|.
|URU| and |UTU| are the same thing for diffuse incident light.

The slides refract but do not absorb here too; both slide optical depths are
set to zero below and there is no argument for them.

@<Prototype for |ez_RT_Oblique|@>=
void ez_RT_Oblique(int n,   
                double nslab, 
                double ntopslide, 
                double nbottomslide, 
                double a,
                double b,
                double g,
                double cos_oblique_angle,
                double *URx, double *UTx, double *URU, double *UTU)

@ @<Definition for |ez_RT_Oblique|@>=
    @<Prototype for |ez_RT_Oblique|@>
{
struct AD_slab_type slab;

    slab.n_slab = nslab;
    slab.n_top_slide = ntopslide;
    slab.n_bottom_slide = nbottomslide;
    slab.b_top_slide = 0;
    slab.b_bottom_slide = 0;
    slab.a = a;
    slab.b = b;
    slab.g = g;
    slab.cos_angle = cos_oblique_angle;
    slab.phase_function = HENYEY_GREENSTEIN;

    RT_Cone(n, &slab, OBLIQUE, URx, UTx, URU, UTU);
}
