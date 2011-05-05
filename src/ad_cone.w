@*1 AD Cone.
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

@*2 RT Cone.
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

    Choose_Cone_Method(slab, &method);

    if (slab->b <= 0) {
        Zero_Layer(n, R12, T12);
        return;
    }

    n = method.quad_pts;
    Init_Layer(*slab, method, R12, T12);

    d= 1.0;
    if (slab->b != HUGE_VAL) 
        d = method.b_thinnest * slab->b / method.b_calc;

    Double_Until(n, R12, T12, d, slab->b);

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
        if (use_cone!=OBLIQUE) 
    		fprintf(stderr,"Unknown type for use_cone.  Assuming oblique incidence.\n");
        URU_and_URx_Cone(n, slab->n_slab, slab->cos_angle, R03, URU, UR1);
        Transpose_Matrix(n,T03);
        URU_and_URx_Cone(n, slab->n_slab, slab->cos_angle, T03, UTU, UT1);
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

@ Simple wrapper that avoids data structures

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
reflected and transmitted for light incident at at |cos_oblique_angle|. 
|URU| and |UTU| are the same thing for diffuse incident light.

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
