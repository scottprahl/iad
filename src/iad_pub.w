@*1 IAD Public.

This contains the routine |Inverse_RT| that should generally be
the basic entry point into this whole mess.  Call this routine with
the proper values and true happiness is bound to be yours.

Altered accuracy of the
standard method of root finding from 0.001 to 0.00001.
Note, it really doesn't help to change the method from
|ABSOLUTE| to |RELATIVE|, but I did anyway.  (3/3/95)

@(iad_pub.c@>=
#include <stdio.h>
#include <math.h>
#include "nr_util.h"
#include "ad_globl.h"
#include "ad_frsnl.h"
#include "iad_type.h"
#include "iad_util.h"
#include "iad_calc.h"
#include "iad_find.h"
#include "iad_pub.h"
#include "iad_io.h"
#include "mc_lost.h"

@<Definition for |Inverse_RT|@>@;
    @<Definition for |measure_OK|@>@;
    @<Definition for |determine_search|@>@;
    @<Definition for |Initialize_Result|@>@;
    @<Definition for |Initialize_Measure|@>@;
    @<Definition for |ez_Inverse_RT|@>@;
    @<Definition for |Spheres_Inverse_RT|@>@;
    @<Definition for |Spheres_Inverse_RT2|@>@;
    @<Definition for |Calculate_MR_MT|@>@;
    @<Definition for |MinMax_MR_MT|@>@;
    @<Definition for |Calculate_Minimum_MR|@>@;

    
@ All the information that needs to be written to
the header file \.{iad\_pub.h}.  This eliminates the need to maintain a set of 
header files as well.

@(iad_pub.h@>=
    @<Prototype for |Inverse_RT|@>;
    @<Prototype for |measure_OK|@>;
    @<Prototype for |determine_search|@>;
    @<Prototype for |Initialize_Result|@>;
    @<Prototype for |ez_Inverse_RT|@>;
    @<Prototype for |Initialize_Measure|@>;
    @<Prototype for |Calculate_MR_MT|@>;
    @<Prototype for |MinMax_MR_MT|@>;
    @<Prototype for |Calculate_Minimum_MR|@>;
    @<Prototype for |Spheres_Inverse_RT2|@>;


@ Here is the header file needed to access one interesting routine in the 
\.{libiad.so} library.

@(lib_iad.h@>=
    @<Prototype for |ez_Inverse_RT|@>;
    @<Prototype for |Spheres_Inverse_RT|@>;
    @<Prototype for |Spheres_Inverse_RT2|@>;

@*2 Inverse RT.
|Inverse_RT| is the main function in this whole package.
You pass the variable |m| containing
your experimentally measured values to the function
|Inverse_RT|.  It hopefully returns the optical properties in |r| 
that are appropriate for your experiment.  

@<Prototype for |Inverse_RT|@>=
    void Inverse_RT(struct measure_type m, struct invert_type *r)

@ @<Definition for |Inverse_RT|@>=
    @<Prototype for |Inverse_RT|@>
{
    if (0 && Debug(DEBUG_LOST_LIGHT)) {
        fprintf(stderr, "** Inverse_RT (%d spheres) **\n", m.num_spheres);
        fprintf(stderr, "    M_R      = %8.5f, MT       = %8.5f\n", m.m_r, m.m_t);
        fprintf(stderr, "    UR1 lost = %8.5f, UT1 lost = %8.5f\n", m.ur1_lost, m.ut1_lost);
    }

    r->found = FALSE;

    @<Exit with bad input data@>@;
    
    r->search = determine_search(m,*r);
    
    if (r->search == FIND_B_WITH_NO_ABSORPTION) {
        r->default_a = 1;
        r->search = FIND_B;
    }
    
    if (r->search == FIND_B_WITH_NO_SCATTERING) {
        r->default_a = 0;
        r->search = FIND_B;
    }

    @<Find the optical properties@>@;
    if ( r->final_distance <= r->tolerance) r->found=TRUE;
}

@ There is no sense going to all the trouble to try a multivariable
minimization if the input data is bogus.  So I wrote a 
single routine |measure_OK| to do just this.

@<Exit with bad input data@>=
    r->error=measure_OK(m,*r);
        
    if (r->method.quad_pts<4) 
        r->error = IAD_QUAD_PTS_NOT_VALID;
        
    if (0 && (r->error != IAD_NO_ERROR)) 
        return;
  
@ Now I fob the real work off to the unconstrained minimization
routines.  Ultimately, I would like to replace all these by constrained
minimization routines.  Actually the first five already are constrained.
The real work will be improving the last five because these are 2-D
minimization routines.

@<Find the optical properties@>=
switch (r->search){
    case FIND_A:  U_Find_A(m,r);
                  break;
    case FIND_B:  U_Find_B(m,r);
                  break;
    case FIND_G:  U_Find_G(m,r);
                  break;
    case FIND_Ba: U_Find_Ba(m,r);
                  break;
    case FIND_Bs: U_Find_Bs(m,r);
                  break;

    case FIND_AB: U_Find_AB(m,r);
                  break;
    case FIND_AG: U_Find_AG(m,r);
                  break;
    case FIND_BG: U_Find_BG(m,r);
                  break;
    case FIND_BsG: U_Find_BsG(m,r);
                  break;
    case FIND_BaG: U_Find_BaG(m,r);
                  break;
}
if (r->iterations==IAD_MAX_ITERATIONS)
    r->error=IAD_TOO_MANY_ITERATIONS;

@*2 Validation.

@ Now the question is --- just what is bad data?  Here's the prototype.

@<Prototype for |measure_OK|@>=
int measure_OK(struct measure_type m, struct invert_type r)

@ It would just be nice to stop computing with bad data.  This does not
work in practice becasue it turns out that there is often bogus data in a full wavelength
scan.  Often the reflectance is too low for short wavelengths and at
long wavelengths the detector (photomultiplier tube) does not work worth
a damn.

The two sphere checks are more complicated.  For example, we can no longer
categorically state that the transmittance is less than one or that the sum
of the reflectance and transmittance is less than one.  Instead we use
the transmittance to bound the values for the reflectance --- see the 
routine |MinMax_MR_MT| below.

@<Definition for |measure_OK|@>=
@<Prototype for |measure_OK|@>
{
    double ru, tu;

    if (m.num_spheres != 2) {  
        @<Check MR for zero or one spheres@>@;
        @<Check MT for zero or one spheres@>@;
    } else {
        int error = MinMax_MR_MT(m,r);
        if (error != IAD_NO_ERROR) return error;
    }

    @<Check MU@>@;

    if (m.num_spheres != 0) {
        @<Check sphere parameters@>@;
    }

    return IAD_NO_ERROR;
}

@ The reflectance is constrained by the index of refraction of the material
and the transmission.  The upper bound for the reflectance is just one minus the
transmittance.  The specular (unscattered) reflectance from the boundaries
imposes minimum for the reflectance. Obviously, the reflected light cannot be
less than that from the first boundary.  This might be calculated by assuming an
infinite layer thickness.  But we can do better.

There is a definite bound on the minimum reflectance from a sample. If you
have a sample with a given transmittance |m_t|, the minimum reflectance possible
is found by assuming that the sample does not scatter any light.

Knowledge of the indicies of refraction makes it a relatively simple matter
to determine the optical thickness |b=mu_a*d| of the slab. The minimum
reflection is obtained by including all the specular reflectances from all the
surfaces.

If the default albedo has been specified as zero, then there is really
no need to check |MR| because it is ignored.

@<Check MR for zero or one spheres@>=

    if (r.default_a == UNINITIALIZED || r.default_a > 0) {
        double mr,mt;
        Calculate_Minimum_MR(m,r,&mr,&mt);
        if (m.m_r < mr)
            return IAD_MR_TOO_SMALL;
    }

@ The transmittance is also constrained by the index of refraction of the
material.  The minimum transmittance is zero, but the maximum transmittance
cannot exceed the total light passing through the sample when there is no
scattering or absorption.  This is calculated by assuming an infinitely thin (to
eliminate any scattering or absorption effects).

There is a problem when spheres are present.  The estimated values
for the transmittance using |Sp_mu_RT| are not actually limiting cases.
This will require a bit of fixing, but for now that test is omitted
if the number of spheres is more than zero.

@<Check MT for zero or one spheres@>=

    if (m.m_t < 0)
        return IAD_MT_TOO_SMALL;      

    Sp_mu_RT_Flip(m.flip_sample, r.slab.n_top_slide, r.slab.n_slab, r.slab.n_bottom_slide, 
             r.slab.b_top_slide, 0, r.slab.b_bottom_slide, r.slab.cos_angle, &ru, &tu);
    
    if (m.num_spheres == 0 && m.m_t > tu) {
fprintf(stderr,"ntop=%7.5f, nslab=%7.5f, nbottom=%7.5f\n", 
r.slab.n_top_slide,r.slab.n_slab,r.slab.n_bottom_slide);
        fprintf(stderr,"tu_max=%7.5f, m_t=%7.5f, t_std=%7.5f\n", tu, m.m_t, m.rstd_t);
        return IAD_MT_TOO_BIG;
    }

@  The unscattered transmission is now always included in the total 
transmittance.  Therefore the unscattered transmittance must fall betwee
zero and |M_T|

@<Check MU @>=

    if (m.m_u < 0)
        return IAD_MU_TOO_SMALL;
    
    if (m.m_u > m.m_t)
        return IAD_MU_TOO_BIG;
    

@ Make sure that reflection sphere parameters are reasonable

@<Check sphere parameters@>=
  
    if (m.as_r < 0 || m.as_r >= 0.2) 
        return IAD_AS_NOT_VALID;
        
    if (m.ad_r < 0 || m.ad_r >= 0.2) 
        return IAD_AD_NOT_VALID;
    
    if (m.ae_r < 0 || m.ae_r >= 0.2) 
        return IAD_AE_NOT_VALID;

    if (m.rw_r < 0 || m.rw_r > 1.0) 
        return IAD_RW_NOT_VALID;

    if (m.rd_r < 0 || m.rd_r > 1.0) 
        return IAD_RD_NOT_VALID;

    if (m.rstd_r < 0 || m.rstd_r > 1.0) 
        return IAD_RSTD_NOT_VALID;

    if (m.rstd_t < 0 || m.rstd_t > 1.0) 
        return IAD_TSTD_NOT_VALID;

    if (m.f_r < 0 || m.f_r > 1) 
        return IAD_F_NOT_VALID;
  
@ Make sure that transmission sphere parameters are reasonable

@<Check sphere parameters@>=
  
    if (m.as_t < 0 || m.as_t >= 0.2) 
        return IAD_AS_NOT_VALID;
        
    if (m.ad_t < 0 || m.ad_t >= 0.2) 
        return IAD_AD_NOT_VALID;
    
    if (m.ae_t < 0 || m.ae_t >= 0.2) 
        return IAD_AE_NOT_VALID;

    if (m.rw_t < 0 || m.rw_r > 1.0) 
        return IAD_RW_NOT_VALID;

    if (m.rd_t < 0 || m.rd_t > 1.0) 
        return IAD_RD_NOT_VALID;

    if (m.rstd_t < 0 || m.rstd_t > 1.0) 
        return IAD_TSTD_NOT_VALID;

    if (m.f_t < 0 || m.f_t > 1) 
        return IAD_F_NOT_VALID;

@*2 Searching Method.  

The original idea was that this routine would automatically determine
what optical parameters could be figured out from the input data.  This
worked fine for a long while, but I discovered that often it was convenient
to constrain the optical properties in various ways.  Consequently, this
routine got more and more complicated.  

What should be done is to figure out whether the search will be 1D or 2D
and split this routine into two parts.

It would be nice to enable the user to constrain two parameters, but
the infrastructure is missing at this point.

@<Prototype for |determine_search|@>=
search_type determine_search(struct measure_type m, struct invert_type r)

@ This routine is responsible for selecting the appropriate 
optical properties to determine.  

@<Definition for |determine_search|@>=
@<Prototype for |determine_search|@>
{
    double rt, tt, rd, td, tc, rc;
    int search=0;
    int independent = m.num_measures;

    if (Debug(DEBUG_SEARCH)) {
        fprintf(stderr,"\n*** Determine_Search()\n");
        fprintf(stderr,"    starting with %d measurement(s)\n",m.num_measures);
        fprintf(stderr,"    m_r=%.5f\n",m.m_r);
        fprintf(stderr,"    m_t=%.5f\n",m.m_t);
    }

    Estimate_RT(m, r, &rt, &tt, &rd, &rc, &td, &tc);

    if (m.m_u==0 && independent == 3) {
        if (Debug(DEBUG_SEARCH)) fprintf(stderr,"    no information in tc\n");
        independent--;
    }

    if (rd==0 && independent == 2) {
        if (Debug(DEBUG_SEARCH)) fprintf(stderr,"    no information in rd\n");
        independent--;
    }
        
    if (td==0 && independent == 2) {
        if (Debug(DEBUG_SEARCH)) fprintf(stderr,"    no information in td\n");
        independent--;
    }
        
    if (independent == 1) {
        @<One parameter search@>@;
    }
    
    else if (independent == 2) {
        @<Two parameter search@>@;
    }

    /* three real parameters with information! */   
    else {
        search = FIND_AG;
    }
    
    if (Debug(DEBUG_SEARCH)) {
        fprintf(stderr,"    independent measurements = %3d\n",independent);
        fprintf(stderr,"    m_r=%8.5f m_t=%8.5f (rd = %8.5f td=%8.5f)\n",m.m_r, m.m_t, rd,td);
        if (search==FIND_A)    fprintf(stderr,"    search = FIND_A\n");
        if (search==FIND_B)    fprintf(stderr,"    search = FIND_B\n");
        if (search==FIND_AB)   fprintf(stderr,"    search = FIND_AB\n");
        if (search==FIND_AG)   fprintf(stderr,"    search = FIND_AG\n");
        if (search==FIND_AUTO) fprintf(stderr,"    search = FIND_AUTO\n");
        if (search==FIND_BG)   fprintf(stderr,"    search = FIND_BG\n");
        if (search==FIND_BaG)  fprintf(stderr,"    search = FIND_BaG\n");
        if (search==FIND_BsG)  fprintf(stderr,"    search = FIND_BsG\n");
        if (search==FIND_Ba)   fprintf(stderr,"    search = FIND_Ba\n");
        if (search==FIND_Bs)   fprintf(stderr,"    search = FIND_Bs\n");
        if (search==FIND_G)    fprintf(stderr,"    search = FIND_G\n");
        if (search==FIND_B_WITH_NO_ABSORPTION)
                                fprintf(stderr,"    search = FIND_B_WITH_NO_ABSORPTION\n");
        if (search==FIND_B_WITH_NO_SCATTERING)
                                fprintf(stderr,"    search = FIND_B_WITH_NO_SCATTERING\n");
    }
    
    return search;
}


@ The fastest inverse problems are those in which just one
measurement is known.  This corresponds to a simple one-dimensional
minimization problem.  The only complexity is deciding exactly what
should be allowed to vary.  The basic assumption is that the 
anisotropy has been specified or will be assumed to be zero.

If the anistropy is assumed known, then one other assumption will
allow us to figure out the last parameter to solve for.

Ultimately, if no default values are given, then we look at the
value of the total transmittance.  If this is zero, then we
assume that the optical thickness is infinite and solve for
the albedo.  Otherwise we will just make a stab at solving
for the optical thickness assuming the albedo is one.  

@<One parameter search@>=
    if (r.default_a  != UNINITIALIZED) {
        if (r.default_a == 0)
            search = FIND_B_WITH_NO_SCATTERING;
        else if (r.default_a == 1)
            search = FIND_B_WITH_NO_ABSORPTION;
        else if (tt == 0) 
            search = FIND_G;
        else
            search = FIND_B;
    }
    else if (r.default_b  != UNINITIALIZED)
        search = FIND_A;

    else if (r.default_bs != UNINITIALIZED)
        search = FIND_Ba;

    else if (r.default_ba != UNINITIALIZED)
        search = FIND_Bs;

    else if (td == 0)
        search = FIND_A;

    else if (rd == 0)
        search = FIND_B_WITH_NO_SCATTERING;

    else 
        search = FIND_B_WITH_NO_ABSORPTION;

@ If the absorption depth $\mu_a d$ is constrained return |FIND_BsG|.
Recall that I use the bizarre mnemonic $bs=\mu_s d$ here and so this
means that the program will search over various values of $\mu_s d$
and $g$.

If there are just two measurements then I assume that the 
anisotropy is not of interest and the only thing to calculate
is the reduced albedo and optical thickness based on an assumed
anisotropy.  

@<Two parameter search@>=
    if (r.default_a != UNINITIALIZED) {
        
        if (r.default_a  == 0)
            search =  FIND_B;
        else if (r.default_g  != UNINITIALIZED)
            search =  FIND_B;
        else
            search =  FIND_BG;
    
    } else if (r.default_b  != UNINITIALIZED) {
    
        if (r.default_g  != UNINITIALIZED)
            search =  FIND_A;
        else
            search =  FIND_AG;

    } else if (r.default_ba != UNINITIALIZED) {
    
        if (r.default_g  != UNINITIALIZED)
            search =  FIND_Bs;
        else
            search =  FIND_BsG; 

    } else if (r.default_bs != UNINITIALIZED) {
    
        if (r.default_g  != UNINITIALIZED)
            search =  FIND_Ba;
        else
            search =  FIND_BaG; 
    
    } else if (rt + tt > 1 && 0 && m.num_spheres != 2) 
        search =  FIND_B_WITH_NO_ABSORPTION;
        
    else
        search = FIND_AB;


@ This little routine just stuffs reasonable values into the 
structure we use to return the solution.  This does not replace 
the values for |r.default_g| nor for |r.method.quad_pts|.  Presumably
these have been set correctly elsewhere.

@<Prototype for |Initialize_Result|@>=
void Initialize_Result(struct measure_type m, struct invert_type *r)

@ @<Definition for |Initialize_Result|@>=
    @<Prototype for |Initialize_Result|@>
{
@<Fill |r| with reasonable values@>@;
}

@ Start with the optical properties.
@<Fill |r| with reasonable values@>=
    r->a = 0.0;
    r->b = 0.0;
    r->g = 0.0;

@ Continue with other useful stuff.
@<Fill |r| with reasonable values@>=
    r->found = FALSE;
    r->tolerance = 0.0001;
    r->MC_tolerance = 0.01;  /* percent */
    r->search = FIND_AUTO;
    r->metric = RELATIVE;
    r->final_distance = 10;
    r->iterations =0;
    r->error = IAD_NO_ERROR;

@ The defaults might be handy

@<Fill |r| with reasonable values@>=
    r->default_a=UNINITIALIZED;
    r->default_b=UNINITIALIZED;
    r->default_g=UNINITIALIZED;
    r->default_ba=UNINITIALIZED;
    r->default_bs=UNINITIALIZED;
    r->default_mua=UNINITIALIZED;
    r->default_mus=UNINITIALIZED;


@ It is necessary to set up the slab correctly so, I stuff reasonable
values into this record as well.
@<Fill |r| with reasonable values@>=
    
    r->slab.a = 0.5;
    r->slab.b = 1.0;
    r->slab.g = 0;
    r->slab.phase_function = HENYEY_GREENSTEIN;
    r->slab.n_slab = m.slab_index;
    r->slab.n_top_slide = m.slab_top_slide_index;
    r->slab.n_bottom_slide = m.slab_bottom_slide_index;
    r->slab.b_top_slide = m.slab_top_slide_b;
    r->slab.b_bottom_slide = m.slab_bottom_slide_b;
    r->slab.cos_angle = m.slab_cos_angle;
    
    r->method.a_calc=0.5;
    r->method.b_calc=1;
    r->method.g_calc=0.5;
    r->method.quad_pts=8;
    r->method.b_thinnest = 1.0/32.0;

@*2 EZ Inverse RT.
|ez_Inverse_RT| is a simple interface to the main function |Inverse_RT|
in this package.  It eliminates the need for complicated data structures
so that the command line interface (as well as those to Perl and Mathematica)
will be simpler.  This function assumes
that the reflection and transmission include specular reflection and that
the transmission also include unscattered transmission.

Other assumptions are that the top and bottom slides have the same
index of refraction, that the illumination is collimated.  Of course
no sphere parameters are included.

@<Prototype for |ez_Inverse_RT|@>=
    void ez_Inverse_RT(double n, double nslide, double UR1, double UT1, double Tc, 
                       double *a, double *b, double *g, int *error)

@ @<Definition for |ez_Inverse_RT|@>=
    @<Prototype for |ez_Inverse_RT|@>
{
  struct measure_type m;
  struct invert_type r;
  *a = 0;
  *b = 0;
  *g = 0;

  Initialize_Measure(&m);
  
  m.slab_index = n;
  m.slab_top_slide_index=nslide;
  m.slab_bottom_slide_index=nslide;
  m.slab_cos_angle=1.0;
  
  m.num_measures=3;
  if (UT1 == 0) m.num_measures--;
  if (Tc  == 0) m.num_measures--;

  m.m_r = UR1;
  m.m_t = UT1;
  m.m_u = Tc;

  Initialize_Result(m,&r);
  r.method.quad_pts=8;
  
  Inverse_RT (m, &r);

  *error = r.error;
  if (r.error == IAD_NO_ERROR) {
    *a = r.a;
    *b = r.b;
    *g = r.g;
  } 
}

@ @<Prototype for |Initialize_Measure|@>=
void Initialize_Measure(struct measure_type *m)

@ @<Definition for |Initialize_Measure|@>=
    @<Prototype for |Initialize_Measure|@>
{
    double default_sphere_d   = 8.0 * 25.4;
    double default_sample_d   = 0.0 * 25.4;
    double default_detector_d = 0.1 * 25.4;
    double default_entrance_d = 0.5 * 25.4;
    double sphere = default_sphere_d * default_sphere_d;
    
    m->slab_index=1.0;
    m->slab_top_slide_index=1.0;
    m->slab_top_slide_b=0.0;
    m->slab_top_slide_thickness=0.0;
    m->slab_bottom_slide_index=1.0;
    m->slab_bottom_slide_b=0.0;
    m->slab_bottom_slide_thickness=0.0;
    m->slab_thickness=1.0;
    m->slab_cos_angle=1.0;
    
    m->num_spheres=0;
    m->num_measures=1;
    m->method = UNKNOWN;
    
    m->fraction_of_rc_in_mr=1.0;
    m->fraction_of_tc_in_mt=1.0;

    m->flip_sample = 0;
    
    m->m_r=0.0;
    m->m_t=0.0;
    m->m_u=0.0;
    
    m->d_sphere_r = default_sphere_d;
    m->as_r = default_sample_d * default_sample_d / sphere;
    m->ad_r = default_detector_d * default_detector_d / sphere;
    m->ae_r = default_entrance_d * default_entrance_d / sphere;
    m->aw_r = 1.0 - m->as_r - m->ad_r -m->ae_r;
    m->rd_r = 0.0;
    m->rw_r = 1.0;
    m->rstd_r = 1.0;
    m->f_r = 0.0;
    
    m->d_sphere_t = default_sphere_d;
    m->as_t = m->as_r;
    m->ad_t = m->ad_r;
    m->ae_t = m->ae_r;
    m->aw_t = m->aw_r;
    m->rd_t = 0.0;
    m->rw_t = 1.0;
    m->rstd_t = 1.0;
    m->f_t = 0.0;

    m->lambda = 0.0;
    m->d_beam = 0.0;
    m->ur1_lost = 0;
    m->uru_lost = 0;
    m->ut1_lost = 0;
    m->utu_lost = 0;
}

@ To avoid interfacing with C-structures it is necessary to pass the information
as arrays.  Here I have divided the experiment into (1) setup, (2) reflection 
sphere coefficients, (3) transmission sphere coefficients, (4) measurements, and (5) 
results.

@<Prototype for |Spheres_Inverse_RT|@>=
    void Spheres_Inverse_RT(double *setup, 
                            double *analysis, 
                            double *sphere_r, 
                            double *sphere_t,
                            double *measurements,
                            double *results)
    
@ @<Definition for |Spheres_Inverse_RT|@>=
    @<Prototype for |Spheres_Inverse_RT|@>
{
    struct measure_type m;
    struct invert_type r;
    long num_photons;
    double ur1,ut1,uru,utu;
    int i, mc_runs = 1;
    
    Initialize_Measure(&m);
    
    @<handle setup @>@;
    @<handle reflection sphere @>@;
    @<handle transmission sphere @>@;
    @<handle measurement @>@;

    Initialize_Result(m,&r);
    results[0]=0;
    results[1]=0;
    results[2]=0;

    @<handle analysis @>@;
    
    Inverse_RT (m, &r);
    for (i=0; i<mc_runs; i++) {
        MC_Lost(m, r, num_photons, &ur1, &ut1, &uru, &utu, 
                     &m.ur1_lost, &m.ut1_lost, &m.uru_lost, &m.utu_lost);   
        Inverse_RT (m, &r);
    }
    
    if (r.error == IAD_NO_ERROR) {
        results[0]=(1-r.a)*r.b/m.slab_thickness;
        results[1]=(r.a  )*r.b/m.slab_thickness;
        results[2]=r.g;
    } 
    
    results[3]=r.error;
}

@ These are in exactly the same order as the parameters in the .rxt header
@<handle setup @>=
{  
    double d_sample_r, d_entrance_r, d_detector_r;
    double d_sample_t, d_entrance_t, d_detector_t;
    
    m.slab_index               = setup[0];
    m.slab_top_slide_index     = setup[1];
    m.slab_thickness           = setup[2];
    m.slab_top_slide_thickness = setup[3];
    m.d_beam                   = setup[4];
    m.rstd_r                   = setup[5];
    m.num_spheres              = (int) setup[6];  
    
    m.d_sphere_r               = setup[7];
    d_sample_r                 = setup[8];
    d_entrance_r               = setup[9];
    d_detector_r               = setup[10];
    m.rw_r                     = setup[11];
    
    m.d_sphere_t               = setup[12];
    d_sample_t                 = setup[13];
    d_entrance_t               = setup[14];
    d_detector_t               = setup[15];
    m.rw_t                     = setup[16];
    
    r.default_g                = setup[17];
    num_photons                = (long) setup[18];
    
    m.as_r = (d_sample_r   / m.d_sphere_r) * (d_sample_r   / m.d_sphere_r);
    m.ae_r = (d_entrance_r / m.d_sphere_r) * (d_entrance_r / m.d_sphere_r);
    m.ad_r = (d_detector_r / m.d_sphere_r) * (d_detector_r / m.d_sphere_r);
    m.aw_r = 1.0 - m.as_r - m.ae_r - m.ad_r;
    m.as_t = (d_sample_t   / m.d_sphere_t) * (d_sample_t   / m.d_sphere_t);
    m.ae_t = (d_entrance_t / m.d_sphere_t) * (d_entrance_t / m.d_sphere_t);
    m.ad_t = (d_detector_t / m.d_sphere_t) * (d_detector_t / m.d_sphere_t);
    m.aw_t = 1.0 - m.as_t - m.ae_t - m.ad_t;

    m.slab_bottom_slide_index = m.slab_top_slide_index;
    m.slab_bottom_slide_thickness = m.slab_top_slide_thickness;
    
    fprintf(stderr,"**** executing FIXME ****/n");
    m.slab_cos_angle = 1.0; /* FIXME */
    
}

@ @<handle analysis@>=

  r.method.quad_pts        = (int) analysis[0];
  mc_runs                  = (int) analysis[1];
  
@  @<handle measurement @>=
  m.m_r = measurements[0];
  m.m_t = measurements[1];
  m.m_u = measurements[2];
 
  m.num_measures=3;
  if (m.m_t == 0) m.num_measures--;
  if (m.m_u == 0) m.num_measures--;

@  @<handle reflection sphere @>=
  m.as_r    = sphere_r[0];
  m.ae_r    = sphere_r[1];
  m.ad_r    = sphere_r[2];
  m.rw_r    = sphere_r[3];
  m.rd_r    = sphere_r[4];
  m.rstd_r  = sphere_r[5];
  m.f_r     = sphere_r[7];
  
@  @<handle transmission sphere @>=
  m.as_t    = sphere_t[0];
  m.ae_t    = sphere_t[1];
  m.ad_t    = sphere_t[2];
  m.rw_t    = sphere_t[3];
  m.rd_t    = sphere_t[4];
  m.rstd_t  = sphere_t[5];
  m.f_t     = sphere_t[7];

@ I needed a routine that would calculate the values of |M_R|
and |M_T| without doing the whole inversion process.  It
seems odd that this does not exist yet.

The values for the lost light |m.uru_lost| etc., should be
calculated before calling this routine.

@<Prototype for |Calculate_MR_MT|@>=
void Calculate_MR_MT(struct measure_type m, 
                     struct invert_type r, 
                     int include_MC,
                     double *M_R, 
                     double *M_T)

@ @<Definition for |Calculate_MR_MT|@>=
    @<Prototype for |Calculate_MR_MT|@>
{
    double distance, ur1, ut1, uru, utu;
    struct measure_type old_mm;
    struct invert_type old_rr;
    
    if (include_MC && m.num_spheres > 0) 
        MC_Lost(m, r, -2000, &ur1, &ut1, &uru, &utu, 
            &(m.ur1_lost), &(m.ut1_lost), &(m.uru_lost), &(m.utu_lost));   
        
    Get_Calc_State(&old_mm, &old_rr);
    Set_Calc_State(m, r);

    Calculate_Distance(M_R, M_T, &distance);

    Set_Calc_State(old_mm, old_rr);
}

@ So, it turns out that the minimum measured |M_R| can
be less than four percent for black glass!  This is because 
the sphere efficiency is much worse for the glass than for
the white standard.

@<Prototype for |Calculate_Minimum_MR|@>=
void Calculate_Minimum_MR(struct measure_type m, 
                     struct invert_type r, double *mr, double *mt)

@ @<Definition for |Calculate_Minimum_MR|@>=
    @<Prototype for |Calculate_Minimum_MR|@>
{
    if (r.default_b == UNINITIALIZED)
        r.slab.b = 9999;
    else
        r.slab.b = r.default_b;

    if (r.default_a == UNINITIALIZED)
        r.slab.a = 0;
    else
        r.slab.a = r.default_a;

    if (r.default_g == UNINITIALIZED)
        r.slab.g = 0.99;
    else
        r.slab.g = r.default_g;

    r.a = r.slab.a;
    r.b = r.slab.b;
    r.g = r.slab.g;
    
    Calculate_MR_MT(m,r,0,mr,mt);
}

@ The minimum possible value of |MR| for a given |MT| will be when
the albedo is zero and the maximum value will be when the albedo
is one.  In the first case there will be no light loss and in the
second we will assume that any light loss is neglible (to maximize |MR|).

The second case is perhaps over-simplified.  Obviously for a fixed thickness
as the albedo increases, the reflectance will increase.  So how does |U_Find_B()|
work when the albedo is set to 1?  

The problem is that to calculate these values one must know the 
optical thickness.  Fortunately with the recent addition of 
constrained minimization, we can do exactly this.

The only thing that remains is to sort out the light lost effect.

@<Prototype for |MinMax_MR_MT|@>=
int MinMax_MR_MT(struct measure_type m, 
                  struct invert_type r)

@ @<Definition for |MinMax_MR_MT|@>=
    @<Prototype for |MinMax_MR_MT|@>
{
    double distance, measured_m_r, min_possible_m_r, max_possible_m_r, temp_m_t;
        
    if (m.m_r < 0)
        return IAD_MR_TOO_SMALL;
    if (m.m_r*m.rstd_r > 1)
        return IAD_MR_TOO_BIG;
    if (m.m_t < 0)
        return IAD_MT_TOO_SMALL;
    if (m.m_t == 0) 
        return IAD_NO_ERROR;
        
    measured_m_r = m.m_r;
    
    m.m_r = 0;
    r.search = FIND_B;
    
    r.default_a = 0;
    U_Find_B(m, &r);
    Calculate_Distance(&min_possible_m_r, &temp_m_t, &distance);
    if (measured_m_r < min_possible_m_r) 
        return IAD_MR_TOO_SMALL;

    r.default_a = 1.0;
    U_Find_B(m, &r);
    Calculate_Distance(&max_possible_m_r, &temp_m_t, &distance);
    if (measured_m_r > max_possible_m_r) 
        return IAD_MR_TOO_BIG;
    
    return IAD_NO_ERROR;
}

@ @<Prototype for |Spheres_Inverse_RT2|@>=
    void Spheres_Inverse_RT2(double *sample, 
                            double *illumination, 
                            double *sphere_r, 
                            double *sphere_t,
                            double *analysis,
                            double *measurement,
                            double *a,
                            double *b,
                            double *g)
    
@ @<Definition for |Spheres_Inverse_RT2|@>=
    @<Prototype for |Spheres_Inverse_RT2|@>
{
    struct measure_type m;
    struct invert_type r;
    long num_photons;
    double ur1,ut1,uru,utu;
    int i, mc_runs = 1;
    
    Initialize_Measure(&m);
    
    @<handle2 sample @>@;
    @<handle2 illumination @>@;
    @<handle2 reflection sphere @>@;
    @<handle2 transmission sphere @>@;
    @<handle2 analysis @>@;
    @<handle2 measurement @>@;

    Initialize_Result(m,&r);
    
    Inverse_RT (m, &r);
    for (i=0; i<mc_runs; i++) {
        MC_Lost(m, r, num_photons, &ur1, &ut1, &uru, &utu, 
                     &m.ur1_lost, &m.ut1_lost, &m.uru_lost, &m.utu_lost);   
        Inverse_RT (m, &r);
    }
    
    if (r.error == IAD_NO_ERROR) {
        *a = r.a;
        *b = r.b;
        *g = r.g;
    } 
}

@ Just move the values from the sample array into the right places
@<handle2 sample @>=
    m.slab_index                  = sample[0];
    m.slab_top_slide_index        = sample[1];
    m.slab_bottom_slide_index     = sample[2];
    m.slab_thickness              = sample[3];
    m.slab_top_slide_thickness    = sample[4];
    m.slab_bottom_slide_thickness = sample[5];
    m.slab_top_slide_thickness    = 0;
    m.slab_bottom_slide_thickness = 0;

@ Just move the values from the illumination array into the right places.  Need
to spend time to figure out how to integrate items 2, 3, and 4
@<handle2 illumination @>=
    m.d_beam                       = illumination[0];
/*  m.lambda                       = illumination[1]; */
/*  m.specular-reflection-excluded = illumination[2]; */
/*  m.direct-transmission-excluded = illumination[3]; */
/*  m.diffuse-illumination         = illumination[4]; */
    m.num_spheres                  = illumination[5];

@  @<handle2 reflection sphere @>=
{
    double d_sample_r, d_entrance_r, d_detector_r;
    
    m.d_sphere_r = sphere_r[0];
    d_sample_r   = sphere_r[1];
    d_entrance_r = sphere_r[2];
    d_detector_r = sphere_r[3];
    m.rw_r       = sphere_r[4];
    m.rd_r       = sphere_r[5];

    m.as_r = (d_sample_r   / m.d_sphere_r) * (d_sample_r   / m.d_sphere_r);
    m.ae_r = (d_entrance_r / m.d_sphere_r) * (d_entrance_r / m.d_sphere_r);
    m.ad_r = (d_detector_r / m.d_sphere_r) * (d_detector_r / m.d_sphere_r);
    m.aw_r = 1.0 - m.as_r - m.ae_r - m.ad_r;
}

@  @<handle2 transmission sphere @>=
{
    double d_sample_t, d_entrance_t, d_detector_t;
    
    m.d_sphere_t = sphere_t[0];
    d_sample_t   = sphere_t[1];
    d_entrance_t = sphere_t[2];
    d_detector_t = sphere_t[3];
    m.rw_t       = sphere_t[4];
    m.rd_t       = sphere_t[5];

    m.as_t = (d_sample_t   / m.d_sphere_t) * (d_sample_t   / m.d_sphere_t);
    m.ae_t = (d_entrance_t / m.d_sphere_t) * (d_entrance_t / m.d_sphere_t);
    m.ad_t = (d_detector_t / m.d_sphere_t) * (d_detector_t / m.d_sphere_t);
    m.aw_t = 1.0 - m.as_t - m.ae_t - m.ad_t;
}

@  @<handle2 analysis @>=
    r.method.quad_pts        = (int) analysis[0];
    mc_runs                  = (int) analysis[1];
    num_photons              = (long) analysis[2];

  
@  @<handle2 measurement @>=
  m.rstd_r = measurement[0];
  m.m_r    = measurement[1];
  m.m_t    = measurement[2];
  m.m_u    = measurement[3];
 
  m.num_measures=3;
  if (m.m_t == 0) m.num_measures--;
  if (m.m_u == 0) m.num_measures--;
