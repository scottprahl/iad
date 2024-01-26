@** IAD Input Output.

The special define below is to get Visual C to suppress silly warnings.

@(iad_io.c@>=
#define _CRT_SECURE_NO_WARNINGS
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "ad_globl.h"
#include "iad_type.h"
#include "iad_io.h"
#include "iad_pub.h"
#include "version.h"

@<Definition for |skip_white|@>@;
@<Definition for |read_number|@>@;
@<Definition for |check_magic|@>@;
@<Definition for |Read_Header|@>@;
@<Definition for |Write_Header|@>@;
@<Definition for |Read_Data_Line|@>@;

@ @(iad_io.h@>=
@<Prototype for |Read_Header|@>;
@<Prototype for |Write_Header|@>;
@<Prototype for |Read_Data_Line|@>;

@*1 Reading the file header.

@ @<Prototype for |Read_Header|@>=
int Read_Header(FILE *fp, struct measure_type *m, int *params)
        
@ Pretty straightforward stuff.  The only thing that needs to be commented
on is that only one slide thickness/index is specified in the file.  This
must be applied to both the top and bottom slides.  Finally, to specify no
slide, then either setting the slide index to 1.0 or the thickness to 0.0
should do the trick.

@<Definition for |Read_Header|@>=
@<Prototype for |Read_Header|@>
{
    double x;
    Initialize_Measure(m);
    if (check_magic(fp))                               return 1;
    if (read_number(fp,&m->slab_index))                return 1;
    if (read_number(fp,&m->slab_top_slide_index))      return 1;
    if (read_number(fp,&m->slab_thickness))            return 1;
    if (read_number(fp,&m->slab_top_slide_thickness))  return 1;
    if (read_number(fp,&m->d_beam))                    return 1;

    if (m->slab_top_slide_thickness == 0.0) m->slab_top_slide_index = 1.0;
    if (m->slab_top_slide_index == 1.0) m->slab_top_slide_thickness = 0.0;
    if (m->slab_top_slide_index == 0.0) {
        m->slab_top_slide_thickness = 0.0;
        m->slab_top_slide_index = 1.0;
    }
    
    m->slab_bottom_slide_index = m->slab_top_slide_index;
    m->slab_bottom_slide_thickness = m->slab_top_slide_thickness;
    
    if (read_number(fp,&m->rstd_r))                    return 1;
    
    if (read_number(fp,&x))                            return 1;
    m->num_spheres = (int) x;
    m->method = SUBSTITUTION;
    
    @<Read coefficients for reflection sphere@>@;

    @<Read coefficients for transmission sphere@>@;
    
    fprintf(stderr, "here\n");
    if (read_number(fp,&x))                            return 1;
    fprintf(stderr, "%d\n", (int) x);
    *params = (int) x;
    m->num_measures = (*params >= 3) ? 3 : *params;

    return 0;
}

@ @<Read coefficients for reflection sphere@>=
{ 
    double d_sample_r, d_entrance_r, d_detector_r;
    if (read_number(fp,&m->d_sphere_r))   return 1;
    if (read_number(fp,&d_sample_r))      return 1;
    if (read_number(fp,&d_entrance_r))    return 1;
    if (read_number(fp,&d_detector_r))    return 1;
    if (read_number(fp,&m->rw_r))         return 1;
    
    m->as_r = (d_sample_r   / m->d_sphere_r) * (d_sample_r   / m->d_sphere_r)/4.0;
    m->ae_r = (d_entrance_r / m->d_sphere_r) * (d_entrance_r / m->d_sphere_r)/4.0;
    m->ad_r = (d_detector_r / m->d_sphere_r) * (d_detector_r / m->d_sphere_r)/4.0;
    m->aw_r = 1.0 - m->as_r - m->ae_r - m->ad_r;
}

@ @<Read coefficients for transmission sphere@>=
{    
    double d_sample_t, d_entrance_t, d_detector_t;
    if (read_number(fp,&m->d_sphere_t))   return 1;
    if (read_number(fp,&d_sample_t))     return 1;
    if (read_number(fp,&d_entrance_t)) return 1;
    if (read_number(fp,&d_detector_t)) return 1;
    if (read_number(fp,&m->rw_t))         return 1;
    
    m->as_t = (d_sample_t   / m->d_sphere_t) * (d_sample_t   / m->d_sphere_t)/4.0;
    m->ae_t = (d_entrance_t / m->d_sphere_t) * (d_entrance_t / m->d_sphere_t)/4.0;
    m->ad_t = (d_detector_t / m->d_sphere_t) * (d_detector_t / m->d_sphere_t)/4.0;
    m->aw_t = 1.0 - m->as_t - m->ae_t - m->ad_t;
}

@*1 Reading just one line of a data file.

This reads a line of data based on the value of |params|.  

If the first number is greater than one then it
is assumed to be the wavelength and is ignored.  
test on the first value of the line. 

A non-zero value is returned upon a failure.

@<Prototype for |Read_Data_Line|@>=
        int Read_Data_Line(FILE *fp, struct measure_type *m, int params)
        
@ @<Definition for |Read_Data_Line|@>=
        @<Prototype for |Read_Data_Line|@>
{
    if (read_number(fp,&m->m_r)) return 1;
    if (m->m_r > 1) {
        m->lambda = m->m_r;
        if (read_number(fp,&m->m_r)) return 1;
    }

    if (params == -1) {
        m->m_t = m->m_r;
        m->m_r = 0;
        return 0;
    }
    
    if (params == 1)       return 0;
    
    if (read_number(fp,&m->m_t))     return 1;
    if (params == 2)       return 0;
    
    if (read_number(fp,&m->m_u))     return 1;
    if (params == 3)       return 0;

    if (read_number(fp,&m->rw_r))   return 1;
    m->rw_t = m->rw_r;
    if (params == 4)       return 0;

    if (read_number(fp,&m->rw_t)) return 1;
    if (params == 5)       return 0;

    if (read_number(fp,&m->rstd_r)) return 1;
    if (params == 6)       return 0;

    if (read_number(fp,&m->rstd_t)) return 1;
    return 0;
}


@ Skip over white space and comments.  It is assumed that \# starts all
comments and continues to the end of a line.  This routine should work on
files with nearly any line ending CR, LF, CRLF.

Failure is indicated by a non-zero return value.

@<Prototype for |skip_white|@>=
int skip_white(FILE *fp)

@ @<Definition for |skip_white|@>=
@<Prototype for |skip_white|@>@;
{
    int c=fgetc(fp);
    
    while (!feof(fp)) {
        if (isspace(c)) 
            c=fgetc(fp);
        else if (c == '#') 
            do c=fgetc(fp); while (!feof(fp) && c!='\n' && c!='\r');
        else
            break;
    }
    
    if (feof(fp)) return 1;

    ungetc(c,fp);
    return 0;
}

@ Read a single number.  Return 0 if there are no problems, otherwise return 1.
@<Prototype for |read_number|@>=
int read_number(FILE *fp, double *x)

@ @<Definition for |read_number|@>=
@<Prototype for |read_number|@>@;
{    
    if (skip_white(fp)) 
        return 1;  
    
    if (fscanf(fp, "%lf", x))
        return 0;
    else
        return 1;
}

@ Ensure that the data file is actually in the right form.  Return 0 if 
the file has the right starting characters.  Return 1 if on a failure.

@<Prototype for |check_magic|@>=
int check_magic(FILE *fp)

@ @<Definition for |check_magic|@>=
@<Prototype for |check_magic|@>@;
{
    char magic[]="IAD1";
    int i,c;
    
    for (i=0; i<4; i++) {
        c = fgetc(fp);
        if (feof(fp) || c != magic[i]) {
            fprintf(stderr, "Sorry, but iad input files must begin with IAD1\n");
            fprintf(stderr, "       as the first four characters of the file.\n");
            fprintf(stderr, "       Perhaps you are using an old iad format?\n");
            return 1;
        }
    }
    
    return 0;
}

@*1 Formatting the header information.

@<Prototype for |Write_Header|@>=
        void Write_Header(struct measure_type m, struct invert_type r, int params)

@ @<Definition for |Write_Header|@>=
        @<Prototype for |Write_Header|@>
{
@<Write slab info@>@;
@<Write irradiation info@>@;
@<Write general sphere info@>@;
@<Write first sphere info@>@;
@<Write second sphere info@>@;
@<Write measure and inversion info@>@;
}

@ @<Write slab info@>=
        double xx;
        
        printf("# Inverse Adding-Doubling %s \n",Version);
        printf("# \n");
        printf("#                        Beam diameter = %7.1f mm\n", m.d_beam);
        printf("#                     Sample thickness = %7.3f mm\n", 
                m.slab_thickness );
        printf("#                  Top slide thickness = %7.3f mm\n", 
                        m.slab_top_slide_thickness );
        printf("#               Bottom slide thickness = %7.3f mm\n", 
                        m.slab_bottom_slide_thickness );
        printf("#           Sample index of refraction = %7.4f\n", 
                m.slab_index );
        printf("#        Top slide index of refraction = %7.4f\n", 
                m.slab_top_slide_index );
        printf("#     Bottom slide index of refraction = %7.4f\n", 
                m.slab_bottom_slide_index );

@ @<Write irradiation info@>=
        printf("# \n");

@ @<Write general sphere info@>=

        printf("#    Fraction unscattered refl. in M_R = %7.1f %%\n", 
        m.fraction_of_rc_in_mr*100);
        printf("#   Fraction unscattered trans. in M_T = %7.1f %%\n", 
        m.fraction_of_tc_in_mt*100);
        printf("# \n");

@ @<Write first sphere info@>=
        printf("# Reflection sphere\n");
        printf("#                      sphere diameter = %7.1f mm\n", 
        m.d_sphere_r );
        printf("#                 sample port diameter = %7.1f mm\n", 
        2*m.d_sphere_r*sqrt(m.as_r) );
        printf("#               entrance port diameter = %7.1f mm\n", 
        2*m.d_sphere_r*sqrt(m.ae_r) );
        printf("#               detector port diameter = %7.1f mm\n", 
        2*m.d_sphere_r*sqrt(m.ad_r) );
        printf("#                     wall reflectance = %7.1f %%\n", m.rw_r*100 );
        printf("#                 standard reflectance = %7.1f %%\n", m.rstd_r*100 );
        printf("#                 detector reflectance = %7.1f %%\n", m.rd_r*100 );
        printf("#\n");

@ @<Write second sphere info@>=
        printf("# Transmission sphere\n");
        printf("#                      sphere diameter = %7.1f mm\n", 
        m.d_sphere_t );
        printf("#                 sample port diameter = %7.1f mm\n", 
        2*m.d_sphere_r*sqrt(m.as_t) );
        printf("#               entrance port diameter = %7.1f mm\n", 
        2*m.d_sphere_r*sqrt(m.ae_t) );
        printf("#               detector port diameter = %7.1f mm\n", 
        2*m.d_sphere_r*sqrt(m.ad_t) );
        printf("#                     wall reflectance = %7.1f %%\n", m.rw_t*100 );
        printf("#               standard transmittance = %7.1f %%\n", m.rstd_t*100 );
        printf("#                 detector reflectance = %7.1f %%\n", m.rd_t*100 );

@ @<Write measure and inversion info@>=
        printf("#\n");
        switch (params) {
            case -1: 
                printf("# No M_R or M_T -- forward calculation.\n");
                break;
            case 1: 
                printf("# Just M_R was measured");
                break;
            case 2:         
                printf("# M_R and M_T were measured");
                break;
            case 3: 
                printf("# M_R, M_T, and M_U were measured");
                break;
            case 4: 
                printf("# M_R, M_T, M_U, and r_w were measured");
                break;
            case 5: 
                printf("# M_R, M_T, M_U, r_w, and t_w were measured");
                break;
            case 6: 
                printf("# M_R, M_T, M_U, r_w, t_w, and r_std were measured");
                break;
            case 7: 
                printf("# M_R, M_T, M_U, r_w, t_w, r_std and t_std were measured");
                break;
            default:
                printf("# Something went wrong ... measures should be 1 to 5!\n");
                break;
        }
        
        if (1<=params && params<=7) {
            if (m.flip_sample)
                printf(" (sample flipped) ");
                
            switch (m.method) {
                case UNKNOWN:
                    printf(" using an unknown method.\n");
                    break;
                case SUBSTITUTION:
                    printf(" using the substitution (single-beam) method.\n");
                    break;
                case COMPARISON:
                    printf(" using the comparison (dual-beam) method.\n");
            }
        }

        switch (m.num_spheres) {
        case 0:
            printf("# No sphere corrections were used");
            break;

        case 1:
            printf("# Single sphere corrections were used");
            break;

        case 2:
            printf("# Double sphere corrections were used");
            break;
        }
        
        printf(" with light incident at %d degrees from the normal", 
                   (int) (acos(m.slab_cos_angle)*57.2958));
        printf(".\n");
        
        switch (r.search) {
            case FIND_AB:   
                 printf("# The inverse routine varied the albedo and optical depth.\n");
                 printf("# \n");
                 xx = (r.default_g != UNINITIALIZED) ? r.default_g : 0;
                 printf("# Default single scattering anisotropy = %7.3f \n", xx);
                 break;
            case FIND_AG:   
                 printf("# The inverse routine varied the albedo and anisotropy.\n");
                 printf("# \n");
                 if (r.default_b != UNINITIALIZED) 
                    printf("#                     Default (mu_t*d) = %7.3g\n", r.default_b);
                 else
                    printf("# \n");
                 break;
            case FIND_AUTO: 
                 printf("# The inverse routine adapted to the input data.\n");
                 printf("# \n");
                 printf("# \n");
                 break;
            case FIND_A:    
                 printf("# The inverse routine varied only the albedo.\n");
                 printf("# \n");
                 xx = (r.default_g != UNINITIALIZED) ? r.default_g : 0;
                 printf("# Default single scattering anisotropy is %7.3f ", xx);
                 xx = (r.default_b != UNINITIALIZED) ? r.default_b : HUGE_VAL;
                 printf(" and (mu_t*d) = %7.3g\n", xx);
                 break;
            case FIND_B:    
                 printf("# The inverse routine varied only the optical depth.\n");
                 printf("# \n");
                 xx = (r.default_g != UNINITIALIZED) ? r.default_g : 0;
                 printf("# Default single scattering anisotropy is %7.3f ", xx);
                 if (r.default_a != UNINITIALIZED)
                    printf("and default albedo = %7.3g\n", r.default_a);
                 else
                    printf("\n");
                 break;
            case FIND_Ba:  
                 printf("# The inverse routine varied only the absorption.\n");
                 printf("# \n");
                 xx = (r.default_bs != UNINITIALIZED) ? r.default_bs : 0;
                 printf("#                     Default (mu_s*d) = %7.3g\n", xx);
                 break;
            case FIND_Bs:  
                 printf("# The inverse routine varied only the scattering.\n");
                 printf("# \n");
                 xx = (r.default_ba != UNINITIALIZED) ? r.default_ba : 0;
                 printf("#                     Default (mu_a*d) = %7.3g\n", xx);
                 break;
            default:  
                 printf("# \n");
                 printf("# \n");
                 printf("# \n");
                 break;
        }

        printf("#                 AD quadrature points = %3d\n", 
                r.method.quad_pts);
        printf("#             AD tolerance for success = %9.5f\n", r.tolerance );
        printf("#      MC tolerance for mu_a and mu_s' = %7.3f %%\n", r.MC_tolerance );
