@*1 AD Prime.
This has the rather stupid name prime because I was at a loss for another.
Currently this is poorly commented.  The fluence routine has not even been
checked.  There may or may not be errors associated with the $n^2$-law 
in there.  It just needs to be checked.

@(ad_prime.c@>=
#include <math.h>
#include <float.h>
#include <stdio.h>
#include "nr_util.h"
#include "ad_globl.h"
#include "ad_bound.h"
#include "ad_start.h"
#include "ad_doubl.h"
#include "ad_prime.h"
#include "ad_matrx.h"
#include "ad_cone.h"

@<Definition for |RT_Matrices|@>@;
@<Definition for |RT|@>@;
@<Definition for |ez_RT|@>@;
@<Definition for |RTabs|@>@;
@<Definition for |Flux_Fluence|@>@;
@<Definition for |ez_RT_unscattered|@>@;

@ @(ad_prime.h@>=
    @h
    @<Prototype for |RT_Matrices|@>;
    @<Prototype for |RT|@>;
    @<Prototype for |ez_RT|@>;
    @<Prototype for |RTabs|@>;
    @<Prototype for |Flux_Fluence|@>;
    @<Prototype for |ez_RT_unscattered|@>;

@ @(lib_ad.h@>=
    @<Prototype for |ez_RT|@>;
    @<Prototype for |ez_RT_unscattered|@>;

@*2 R and T Matrix routines.
This section contains the routine to calculate the reflection and transmission
matrix for a scattering and absorbing slab.  Basically you just need to set the
number of quadrature points |method->quad_pts| and the optical properties
(the albedo, anisotropy, optical thickness, and choice of phase function) in |slab|.
Call this routine and get back matrices filled with cool numbers.

@<Prototype for |RT_Matrices|@>=
void RT_Matrices(int n, struct AD_slab_type * slab, struct AD_method_type * method, 
double **R, double **T)

@ @<Definition for |RT_Matrices|@>=
    @<Prototype for |RT_Matrices|@>
{
    double d;

    if (n < 3)
        method->quad_pts = DEFAULT_QUAD_PTS;
    else
        if (n > MAX_QUAD_PTS)
            method->quad_pts = MAX_QUAD_PTS;
    else
        if ((n & 1) == 1)
            method->quad_pts = n / 2 * 2;
    else
        method->quad_pts = n;

    Choose_Method(slab, method);

    if (slab->b <= 0) {
        Zero_Layer(n, R, T);
        return;
    }

    n = method->quad_pts;
    Init_Layer(*slab, *method, R, T);

    if (slab->b == HUGE_VAL) 
        d = 1.0; /* Ignored ... just set it something. */
    else
        d = method->b_thinnest * slab->b / method->b_calc;
    
    Double_Until(n, R, T, d, slab->b);
}


@*2 Total reflection and transmission.

|RT| is the top level routine for accessing the adding-doubling algorithm.
By passing the optical paramters characteristic of the slab, this routine will
do what it must to return the total reflection and transmission for collimated
and diffuse irradiance.

This routine has three different components based on if zero, one, or two boundary
layers must be included.  If the index of refraction of the slab and the top and bottom
slides are all one, then no boundaries need to be included.  If the top and bottom slides
are identical, then some simplifications can be made and some time saved as a consequence.
If the top and bottom slides are different, then the full red carpet treatment is required.

Since the calculation time increases for each of these cases we test 
for matched boundaries first.  If the boundaries are matched then don't bother
with boundaries for the top and bottom.  Just calculate the integrated reflection
and transmission.   Similarly, if the top and bottom slides are similar, then quickly
calculate these.

@<Prototype for |RT|@>=
void RT(int n, struct AD_slab_type * slab, double *UR1, double *UT1, double *URU, double *UTU)

@ @<Definition for |RT|@>=
    @<Prototype for |RT|@>
{
    @<Declare variables for |RT|@>@;

	if (slab->cos_angle != 1.0) {
		RT_Cone(n,slab,OBLIQUE,UR1,UT1,URU,UTU);
		return;
	}
	
    @<Validate input parameters@>@;
    
    @<Allocate and calculate R and T for homogeneous slab@>@;

    if (slab->b == 0) {
        Sp_RT(n, *slab, UR1, UT1, URU, UTU);

    } else if (slab->n_slab == 1 && slab->n_top_slide == 1 && slab->n_bottom_slide == 1
                          && slab->b_top_slide == 0 && slab->b_bottom_slide == 0) {
        @<Do slab with no boundaries@>@;
    } else if (slab->n_top_slide == slab->n_bottom_slide &&
               slab->b_top_slide == 0 && slab->b_bottom_slide == 0) {
        @<Allocate and generate top boundary@>@;
        @<Do slab with matched top and bottom boundaries@>@;
        @<Free top boundary@>@;
    } else {
        @<Allocate and generate top boundary@>@;
        @<Allocate and generate bottom boundary@>@;
        @<Allocate misc matrices@>@;
        @<Do slab with mismatched boundaries@>@;
        @<Free misc matrices@>@;
        @<Free bottom boundary@>@;
        @<Free top boundary@>@;
    }

    @<Free R and T@>@;
}

@ @<Declare variables for |RT|@>=
    double **R, **T, **R2, **T2;
    double *R01, *R10, *T01, *T10;
    double *R23, *R32, *T23, *T32;
    double **R02, **R20, **T02, **T20;
    double **R03, **R30, **T03, **T30;
    double **atemp, **btemp;
    struct AD_method_type method;
    *UR1=-1;
    *URU=-1;
    *UT1=-1;
    *UTU=-1;

@   @<Validate input parameters@>=
    if (slab->n_slab<0) return;
    if (slab->n_top_slide<0) return;
    if (slab->n_bottom_slide<0) return;
    if (slab->a<0  || slab->a>1) return;
    if (slab->g<-1 || slab->g>1) return;
    if (slab->b<0) return;
    
@ Find the R and T for a homogeneous slab without boundaries 
@<Allocate and calculate R and T for homogeneous slab@>=
    R = dmatrix(1, n, 1, n);
    T = dmatrix(1, n, 1, n);
    RT_Matrices(n, slab, &method, R, T);

@ @<Do slab with no boundaries@>=
    URU_and_UR1(n, slab->n_slab, R, URU, UR1);
    URU_and_UR1(n, slab->n_slab, T, UTU, UT1);
 
@ @<Allocate and generate top boundary@>=
    R01 = dvector(1, n);
    R10 = dvector(1, n);
    T01 = dvector(1, n);
    T10 = dvector(1, n);
    Init_Boundary(*slab, method.quad_pts, R01, R10, T01, T10, TOP_BOUNDARY);

@ @<Do slab with matched top and bottom boundaries@>=
    atemp = dmatrix(1, n, 1, n);
    btemp = dmatrix(1, n, 1, n);
    R2 = dmatrix(1, n, 1, n);
    T2 = dmatrix(1, n, 1, n);
    Add_Slides(n, R01, R10, T01, T10, R, T, R2, T2, atemp, btemp);
    URU_and_UR1(n, slab->n_slab, R2, URU, UR1);
    URU_and_UR1(n, slab->n_slab, T2, UTU, UT1);
    free_dmatrix(atemp, 1, n, 1, n);
    free_dmatrix(btemp, 1, n, 1, n);
    free_dmatrix(R2, 1, n, 1, n);
    free_dmatrix(T2, 1, n, 1, n);
    
@ @<Free top boundary@>=
    free_dvector(R01, 1, n);
    free_dvector(R10, 1, n);
    free_dvector(T01, 1, n);
    free_dvector(T10, 1, n);

@ @<Allocate and generate bottom boundary@>=
    R23 = dvector(1, n);
    R32 = dvector(1, n);
    T23 = dvector(1, n);
    T32 = dvector(1, n);
    Init_Boundary(*slab, method.quad_pts, R23, R32, T23, T32, BOTTOM_BOUNDARY);

@ @<Allocate misc matrices@>=

    R02 = dmatrix(1, n, 1, n);
    R20 = dmatrix(1, n, 1, n);
    T02 = dmatrix(1, n, 1, n);
    T20 = dmatrix(1, n, 1, n);
    R03 = dmatrix(1, n, 1, n);
    R30 = dmatrix(1, n, 1, n);
    T03 = dmatrix(1, n, 1, n);
    T30 = dmatrix(1, n, 1, n);
    atemp = dmatrix(1, n, 1, n);
    btemp = dmatrix(1, n, 1, n);

@ @<Do slab with mismatched boundaries@>=
    Add_Top(n, R01, R10, T01, T10, R, R, T, T, R02, R20, T02, T20, atemp, btemp);
    Add_Bottom(n, R02, R20, T02, T20, R23, R32, T23, T32, R03, R30, T03, T30, atemp, btemp);
    URU_and_UR1(n, slab->n_slab, R03, URU, UR1);
    Transpose_Matrix(n,T03);
    URU_and_UR1(n, slab->n_slab, T03, UTU, UT1);

@ @<Free misc matrices@>=
    free_dmatrix(R02, 1, n, 1, n);
    free_dmatrix(R20, 1, n, 1, n);
    free_dmatrix(T02, 1, n, 1, n);
    free_dmatrix(T20, 1, n, 1, n);

    free_dmatrix(R03, 1, n, 1, n);
    free_dmatrix(R30, 1, n, 1, n);
    free_dmatrix(T03, 1, n, 1, n);
    free_dmatrix(T30, 1, n, 1, n);

    free_dmatrix(atemp, 1, n, 1, n);
    free_dmatrix(btemp, 1, n, 1, n);

@ @<Free bottom boundary@>=
    free_dvector(R23, 1, n);
    free_dvector(R32, 1, n);
    free_dvector(T23, 1, n);
    free_dvector(T32, 1, n);


@ @<Free R and T@>=
    free_dmatrix(R, 1, n, 1, n);
    free_dmatrix(T, 1, n, 1, n);

@*2 Simple interfaces for Perl, Python, or Mathematica.

|ez_RT| is a top level routine for accessing the adding-doubling algorithm.
This routine was originally created so that I could make a Perl .xs module.  
Since I did not know how to mess around with passing structures, I changed 
the interface to avoid using structures.

@<Prototype for |ez_RT|@>=
void ez_RT(int n,   double nslab, 
                    double ntopslide, 
                    double nbottomslide, 
                    double a,
                    double b,
                    double g,
                    double *UR1, double *UT1, double *URU, double *UTU)

@ @<Definition for |ez_RT|@>=
    @<Prototype for |ez_RT|@>
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
    slab.phase_function = HENYEY_GREENSTEIN;
    slab.cos_angle = 1.0;
    RT(n, &slab, UR1, UT1, URU, UTU);
}

@*2 Unscattered reflection and transmission. 
                
|ez_RT_unscattered| is a top level routine for accessing the adding-doubling algorithm.
This routine was created so that I could make a Perl module.  Since I did
not know how to mess around with passing structures, I changed the interface
to avoid using structures.  
        
@<Prototype for |ez_RT_unscattered|@>=
void ez_RT_unscattered(int n,       
                       double nslab,
                       double ntopslide,
                       double nbottomslide,
                       double a,
                       double b,
                       double g,
                       double *UR1, double *UT1, double *URU, double *UTU)
        
@ @<Definition for |ez_RT_unscattered|@>=
        @<Prototype for |ez_RT_unscattered|@> 
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
        slab.phase_function = HENYEY_GREENSTEIN;
        slab.cos_angle = 1.0;
        Sp_RT(n, slab, UR1, UT1, URU, UTU);
}       

@*2 Including absorbing slides.

The idea is to create a function that includes absorption in the top
and bottom slides.  This is done by creating two extra layers, finding
the full reflection and transmission matrices for these layers and adding
them to the slab.  Of course this only works when all the indices of refraction
are the same.  Yikes!

This routine returns |UR1| and |UT1| for light incident from the top of
the slab.  The values for light incident from the bottom will be different
when the slides on the top and bottom are different.  {\it Caveat emptor!}

@<Prototype for |RTabs|@>=
void RTabs(int n, struct AD_slab_type * slab, double *UR1, double *UT1, double *URU, double *UTU)

@ @<Definition for |RTabs|@>=
    @<Prototype for |RTabs|@>
{
    @<Declare variables for |RTabs|@>@;
    double **Rtop, **Ttop, **Rbottom, **Tbottom;
    struct AD_slab_type slab1;
    double btop, bbottom;
    
    @<Allocate and calculate R and T for homogeneous slab@>@;
    @<Allocate and calculate top absorbing slide@>@;
    @<Allocate and calculate bottom absorbing slide@>@;
    @<Allocate misc matrices@>@;

    @<Allocate and calculate top non-absorbing boundary@>@;
    @<Allocate and calculate bottom non-absorbing boundary@>@;
    
    @<Add all the stuff together@>@;

    @<Free misc matrices@>@;
    @<Free bottom boundary@>@;
    @<Free top boundary@>@;
    @<Free R and T@>@;
    @<Free matrices for the top and bottom absorbing slides@>@;
}

@ @<Declare variables for |RTabs|@>=
    double **R, **T;
    double *R01, *R10, *T01, *T10;
    double *R23, *R32, *T23, *T32;
    double **R02, **R20, **T02, **T20;
    double **R03, **R30, **T03, **T30;
    double **atemp, **btemp;
    struct AD_method_type method;

@ @<Allocate and calculate top absorbing slide@>=
  
  slab1.b = slab->b_top_slide;
  slab1.cos_angle = slab->cos_angle;
  slab1.a = 0;
  slab1.g = 0;
  slab1.phase_function = HENYEY_GREENSTEIN;
  slab1.n_slab = slab->n_slab;
  slab1.n_top_slide = 1.0;
  slab1.n_bottom_slide = 1.0;
  slab1.b_top_slide = 0.0;
  slab1.b_bottom_slide = 0.0;

  Rtop = dmatrix(1, n, 1, n);
  Ttop = dmatrix(1, n, 1, n);
  
  RT_Matrices(n, &slab1, &method, Rtop, Ttop);

@ @<Allocate and calculate bottom absorbing slide@>=
  slab1.b = slab->b_bottom_slide;
  slab1.cos_angle = slab->cos_angle;

  Rbottom = dmatrix(1, n, 1, n);
  Tbottom = dmatrix(1, n, 1, n);
  RT_Matrices(n, &slab1, &method, Rbottom, Tbottom);
 
 @      @<Allocate and calculate top non-absorbing boundary@>=
    btop=slab->b_top_slide;
    slab->b_top_slide = 0;
    @<Allocate and generate top boundary@>@;
    slab->b_top_slide = btop;

@       @<Allocate and calculate bottom non-absorbing boundary@>=
    bbottom=slab->b_bottom_slide;
    slab->b_bottom_slide = 0;
    @<Allocate and generate bottom boundary@>@;
    slab->b_bottom_slide = bbottom;

@   @<Add all the stuff together@>=
    Add(n, Rtop, Rtop, Ttop, Ttop, R, R, T, T, R02, R20, T02, T20);
    Add(n, R02, R20, T02, T20, Rbottom, Rbottom, Tbottom, Tbottom, R03, R30, T03, T30);
    Add_Top(n, R01, R10, T01, T10, R03, R30, T03, T30, R02, R20, T02, T20, atemp, btemp);
    Add_Bottom(n, R02, R20, T02, T20, R23, R32, T23, T32, R03, R30, T03, T30, atemp, btemp);
    URU_and_UR1(n, slab->n_slab, R03, URU, UR1);
    Transpose_Matrix(n,T03);
    URU_and_UR1(n, slab->n_slab, T03, UTU, UT1);

@   @<Free matrices for the top and bottom absorbing slides@>=
    free_dmatrix(Rtop, 1, n, 1, n);
    free_dmatrix(Ttop, 1, n, 1, n);
    free_dmatrix(Rbottom, 1, n, 1, n);
    free_dmatrix(Tbottom, 1, n, 1, n);

@*2 Flux and Fluence.

Calculates the flux and fluence at various depths between the
optical depths |zmin| and |zmax| for a slab.  The number of values
is |intervals+1| times...i.e. it calculates at |zmin|, |zmin +(zmax-zmin)/intervals|, ... , |zmax|

The fluence and fluxes
at |0| and |slab.b| are calculated just inside the boundary, i.e. beneath any
existing glass slide or just below a mismatched boundary.

This routine could be improved dramatically.  I just have not had the need so far.

This has not been adequately tested.

@d MAX_FLUENCE_INTERVALS  200


@ @<Prototype for |Flux_Fluence|@>=
    void Flux_Fluence(int n, struct AD_slab_type * slab, double zmin, double zmax, int intervals,
                        double *UF1_array, double *UFU_array, double *flux_up, double *flux_down)

@ @<Definition for |Flux_Fluence|@>=
    @<Prototype for |Flux_Fluence|@>
{
    @<Declare variables for |Flux_Fluence|@>@;

    if (intervals > MAX_FLUENCE_INTERVALS)
        AD_error("too many intervals requested.  increase the const max_fluence_intervals\n");

    @<Find the 02 matrix for the slab above all layers@>@;
    @<Find the 46 matrix for the slab below all layers@>@;
    @<Allocate intermediate matrices@>@;

    for (i = 0; i <= intervals; i++) {

        @<Find radiance at each depth@>@;
        
        @<Calculate Fluence and Flux@>@;
    }

    @<Free all those intermediate matrices@>@;
}

@ @<Declare variables for |Flux_Fluence|@>=

    double *R01, *R10, *T01, *T10;
    double *R56, *R65, *T56, *T65;
    double **R12, **T12;
    double **R23, **T23;
    double **R34, **T34;
    double **R45, **T45;
    double **R02, **R20, **T02, **T20;
    double **R46, **R64, **T46, **T64;
    double **R03, **R30, **T03, **T30;
    double **R36, **R63, **T36, **T63;
    double **Lup, **Ldown;
    double **a, **b;
    double flx_down, flx_up, UFU, UF1;
    double slab_thickness;

    struct AD_method_type method;
    int i, j;

@   @<Find the 02 matrix for the slab above all layers@>=
    slab_thickness = slab->b; /* save it for later */
    slab->b = zmin;
    R12 = dmatrix(1, n, 1, n);
    T12 = dmatrix(1, n, 1, n);
    RT_Matrices(n, slab, &method, R12, T12);

    R01 = dvector(1, n);
    R10 = dvector(1, n);
    T01 = dvector(1, n);
    T10 = dvector(1, n);
    Init_Boundary(*slab, method.quad_pts, R01, R10, T01, T10, TOP_BOUNDARY);

    R20 = dmatrix(1, n, 1, n);
    T20 = dmatrix(1, n, 1, n);
    R02 = dmatrix(1, n, 1, n);
    T02 = dmatrix(1, n, 1, n);
    a = dmatrix(1, n, 1, n);
    b = dmatrix(1, n, 1, n);
    Add_Top(n, R01, R10, T01, T10, R12, R12, T12, T12, R02, R20, T02, T20, a, b);

    free_dmatrix(R12, 1, n, 1, n);
    free_dmatrix(T12, 1, n, 1, n);
    free_dvector(R01, 1, n);
    free_dvector(R10, 1, n);
    free_dvector(T01, 1, n);
    free_dvector(T10, 1, n);

@ @<Find the 46 matrix for the slab below all layers@>=
    slab->b = slab_thickness - zmax;
    R45 = dmatrix(1, n, 1, n);
    T45 = dmatrix(1, n, 1, n);
    RT_Matrices(n, slab, &method, R45, T45);
    R56 = dvector(1, n);
    R65 = dvector(1, n);
    T56 = dvector(1, n);
    T65 = dvector(1, n);
    Init_Boundary(*slab, method.quad_pts, R56, R65, T56, T65, BOTTOM_BOUNDARY);
    R46 = dmatrix(1, n, 1, n);
    T46 = dmatrix(1, n, 1, n);
    R64 = dmatrix(1, n, 1, n);
    T64 = dmatrix(1, n, 1, n);
    Add_Bottom(n, R45, R45, T45, T45, R56, R65, T56, T65, R46, R64, T46, T64, a, b);
    free_dmatrix(R45, 1, n, 1, n);
    free_dmatrix(T45, 1, n, 1, n);
    free_dvector(R56, 1, n);
    free_dvector(R65, 1, n);
    free_dvector(T56, 1, n);
    free_dvector(T65, 1, n);
    free_dmatrix(a, 1, n, 1, n);
    free_dmatrix(b, 1, n, 1, n);

@ @<Allocate intermediate matrices@>=
    R23 = dmatrix(1, n, 1, n);
    T23 = dmatrix(1, n, 1, n);
    R03 = dmatrix(1, n, 1, n);
    T03 = dmatrix(1, n, 1, n);
    R30 = dmatrix(1, n, 1, n);
    T30 = dmatrix(1, n, 1, n);

    R34 = dmatrix(1, n, 1, n);
    T34 = dmatrix(1, n, 1, n);
    R63 = dmatrix(1, n, 1, n);
    T63 = dmatrix(1, n, 1, n);
    R36 = dmatrix(1, n, 1, n);
    T36 = dmatrix(1, n, 1, n);

    Lup = dmatrix(1, n, 1, n);
    Ldown = dmatrix(1, n, 1, n);

@ @<Find radiance at each depth@>=
    slab->b = (zmax - zmin) / intervals * i;
    RT_Matrices(n, slab, &method, R23, T23);
    Add(n, R02, R20, T02, T20, R23, R23, T23, T23, R03, R30, T03, T30);

    slab->b = (zmax - zmin) - slab->b;
    RT_Matrices(n, slab, &method, R34, T34);
    Add(n, R34, R34, T34, T34, R46, R64, T46, T64, R36, R63, T36, T63);

    Between(n, R03, R30, T03, T30, R36, R63, T36, T63, Lup, Ldown);

@ @<Calculate Fluence and Flux@>=
    UFU_and_UF1(n, slab->n_slab, Lup, Ldown, &UFU, &UF1);
    UF1_array[i] = UF1;
    UFU_array[i] = UFU;

    flx_down = 0.0;
    flx_up = 0.0;
    for (j = 1; j <= n; j++) {
        flx_down += twoaw[j] * Ldown[j][n];
        flx_up += twoaw[j] * Lup[j][n];
    }
    flux_down[i] = flx_down * slab->n_slab * slab->n_slab;
    flux_up[i] = flx_up * slab->n_slab * slab->n_slab;
        
@ @<Free all those intermediate matrices@>=
    free_dmatrix(R02, 1, n, 1, n);
    free_dmatrix(T02, 1, n, 1, n);
    free_dmatrix(R20, 1, n, 1, n);
    free_dmatrix(T20, 1, n, 1, n);

    free_dmatrix(R23, 1, n, 1, n);
    free_dmatrix(T23, 1, n, 1, n);

    free_dmatrix(R03, 1, n, 1, n);
    free_dmatrix(T03, 1, n, 1, n);
    free_dmatrix(R30, 1, n, 1, n);
    free_dmatrix(T30, 1, n, 1, n);

    free_dmatrix(R34, 1, n, 1, n);
    free_dmatrix(T34, 1, n, 1, n);

    free_dmatrix(R63, 1, n, 1, n);
    free_dmatrix(T63, 1, n, 1, n);
    free_dmatrix(R36, 1, n, 1, n);
    free_dmatrix(T36, 1, n, 1, n);

    free_dmatrix(R64, 1, n, 1, n);
    free_dmatrix(T64, 1, n, 1, n);
    free_dmatrix(R46, 1, n, 1, n);
    free_dmatrix(T46, 1, n, 1, n);

    free_dmatrix(Lup, 1, n, 1, n);
    free_dmatrix(Ldown, 1, n, 1, n);

