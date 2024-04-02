@** IAD Input Output.

The special define below is to get Visual C to suppress silly warnings.

@(iad_io.c@>=
#define _CRT_SECURE_NO_WARNINGS
#define MAX_COLUMNS 256
char COLUMN_LABELS[MAX_COLUMNS] = "";

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include "ad_globl.h"
#include "iad_type.h"
#include "iad_io.h"
#include "iad_pub.h"
#include "version.h"

@<Definition for |get_current_line_number|@>@;
@<Definition for |skip_white|@>@;
@<Definition for |read_number|@>@;
@<Definition for |check_magic|@>@;

@<Definition for |remove_whitespace|@>@;
@<Definition for |remove_comment|@>@;
@<Definition for |remove_first_char|@>@;
@<Definition for |print_maybe|@>@;
@<Definition for |Read_Data_Legend|@>@;
@<Definition for |Read_Data_Line_Per_Labels|@>@;

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
    @<Read info about measurements@>@;

    return 0;
}

@ @<Read coefficients for reflection sphere@>=
{
    double d_sample_r, d_third_r, d_detector_r;
    if (read_number(fp,&m->d_sphere_r))   return 1;
    if (read_number(fp,&d_sample_r))      return 1;
    if (read_number(fp,&d_third_r))    return 1;
    if (read_number(fp,&d_detector_r))    return 1;
    if (read_number(fp,&m->rw_r))         return 1;

    m->as_r = (d_sample_r   / m->d_sphere_r / 2.0) * (d_sample_r   / m->d_sphere_r / 2.0);
    m->at_r = (d_third_r    / m->d_sphere_r / 2.0) * (d_third_r    / m->d_sphere_r / 2.0);
    m->ad_r = (d_detector_r / m->d_sphere_r / 2.0) * (d_detector_r / m->d_sphere_r / 2.0);
    m->aw_r = 1.0 - m->as_r - m->at_r - m->ad_r;
}

@ @<Read coefficients for transmission sphere@>=
{
    double d_sample_t, d_third_t, d_detector_t;
    if (read_number(fp,&m->d_sphere_t))   return 1;
    if (read_number(fp,&d_sample_t))     return 1;
    if (read_number(fp,&d_third_t)) return 1;
    if (read_number(fp,&d_detector_t)) return 1;
    if (read_number(fp,&m->rw_t))         return 1;

    m->as_t = (d_sample_t   / m->d_sphere_t / 2.0) * (d_sample_t   / m->d_sphere_t / 2.0);
    m->at_t = (d_third_t    / m->d_sphere_t / 2.0) * (d_third_t    / m->d_sphere_t / 2.0);
    m->ad_t = (d_detector_t / m->d_sphere_t / 2.0) * (d_detector_t / m->d_sphere_t / 2.0);
    m->aw_t = 1.0 - m->as_t - m->at_t - m->ad_t;
}

@ @<Read info about measurements@>=

    *params = Read_Data_Legend(fp);
    
    if (COLUMN_LABELS[0] != '\0') {
        m->num_measures = 0;
        if ( strchr(COLUMN_LABELS, 'r') ) m->num_measures++;
        if ( strchr(COLUMN_LABELS, 't') ) m->num_measures++;
        if ( strchr(COLUMN_LABELS, 'u') ) m->num_measures++;
        if (m->num_measures == 0) {
            fprintf(stderr, "Column labels must have at least one 'r', 't', or 'u'\n");
            fprintf(stderr, "Column labels = '%s'\n", COLUMN_LABELS);
            exit(EXIT_FAILURE);
        }
    } else 
        m->num_measures = (*params >= 3) ? 3 : *params;

@*1 Reading just one line of a data file.

This reads a line of data based on the value of |params|.

If the first number is greater than one then it
is assumed to be the wavelength and is ignored.
test on the first value of the line.

A non-zero value is returned upon a failure.

@<Prototype for |Read_Data_Line|@>=
int Read_Data_Line(FILE *fp, struct measure_type *m, struct invert_type *r, int params)

@ @<Definition for |Read_Data_Line|@>=
        @<Prototype for |Read_Data_Line|@>
{
    if (strlen(COLUMN_LABELS)>0)
        return Read_Data_Line_Per_Labels(fp, m, r, params);

    if (read_number(fp,&m->m_r)) return 1;
    if (m->m_r > 1) {
        m->lambda = m->m_r;
        if (read_number(fp,&m->m_r)) return 1;
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

@ @<Definition for |Read_Data_Line_Per_Labels|@>=
int Read_Data_Line_Per_Labels(FILE *fp, struct measure_type *m, struct invert_type *r, int params)
{
    int count=0;
    double x;
    while (count<params) {
        if (read_number(fp,&x)) return 1;
        char c = COLUMN_LABELS[count];
        if (FALSE)
            fprintf(stderr, "count = %2d, option = %c, value = %10.5f\n", count, c, x);
        switch (c) {
            case 'a':
                r->default_a = x;
                break;
            case 'A':
                r->default_mua = x;
                r->default_ba = x * m->slab_thickness;
                break;
            case 'b':
                r->default_b = x;
                break;
            case 'B':
                m->d_beam = x;
                break;
            case 'c':
                m->fraction_of_ru_in_mr = x;
                break;
            case 'C':
                m->fraction_of_tu_in_mt = x;
                break;
            case 'd':
                m->slab_thickness = x;
                break;
            case 'D':
                m->slab_top_slide_thickness = x;
                m->slab_bottom_slide_thickness = x;
                break;
            case 'e':
                r->tolerance = x;
                r->MC_tolerance = x;
                break;
            case 'E':
                m->slab_bottom_slide_b = x;
                m->slab_top_slide_b = x;
                break;
            case 'F':
                r->default_mus = x;
                r->default_bs = x * m->slab_thickness;
                break;
            case 'g':
                r->default_g = x;
                break;
            case 'L':
                m->lambda = x;
                break;
            case 'M':
                m->num_spheres = (int) x;
                break;
            case 'n':
                m->slab_index = x;
                break;
            case 'N':
                m->slab_top_slide_index = x;
                m->slab_bottom_slide_index = x;
                break;
            case 'q':
                r->method.quad_pts = (int) x;
                break;
            case 'r':
                m->m_r = x;
                break;
            case 'R':
                m->rstd_r = x;
                break;
            case 't':
                m->m_t = x;
                break;
            case 'S':
                m->num_spheres = (int) x;
                break;
            case 'T':
                m->rstd_t = x;
                break;
            case 'u':
                m->m_u = x;
                break;
            case 'w':
                m->rw_r = x;
                break;
            case 'W':
                m->rw_t = x;
                break;
            default:
                fprintf(stderr, "legend variable '%c' unimplemented", c);
                return 1;
        }
        count++;
    }
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

    if (feof(fp)) 
        return 1;

    ungetc(c,fp);
    return 0;
}

@ Read a single number.  Return 0 if there are no problems, otherwise return 1.
@<Prototype for |read_number|@>=
int read_number(FILE *fp, double *x)

@ @<Definition for |read_number|@>=
@<Prototype for |read_number|@>@;
{
    long line_no = 0;
    
    if (skip_white(fp)) 
        return 1;

    if (fscanf(fp, "%lf", x))
        return 0;

    line_no = get_current_line_number(fp);
    fprintf(stderr, "\nIAD ERROR: Bad number on line %ld in input file.", line_no);
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
void Write_Header(struct measure_type m, struct invert_type r, int params, char *cmd)

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
        printf("# %s\n", cmd);
        printf("#                        Beam diameter = ");
        print_maybe('B', "%7.1f mm\n", m.d_beam);
        printf("#                     Sample thickness = ");
        print_maybe('d', "%7.3f mm\n", m.slab_thickness);
        printf("#                  Top slide thickness = ");
        print_maybe('D', "%7.3f mm\n", m.slab_top_slide_thickness);
        printf("#               Bottom slide thickness = ");
        print_maybe('D', "%7.3f mm\n", m.slab_bottom_slide_thickness);
        printf("#           Sample index of refraction = ");
        print_maybe('n', "%7.4f mm\n", m.slab_index);
        printf("#        Top slide index of refraction = ");
        print_maybe('N', "%7.4f mm\n", m.slab_top_slide_index);
        printf("#     Bottom slide index of refraction = ");
        print_maybe('N', "%7.4f mm\n", m.slab_bottom_slide_index);

@ @<Write irradiation info@>=
        printf("# \n");

@ @<Write general sphere info@>=

        printf("#  Percentage unscattered refl. in M_R = ");
        print_maybe('c', "%7.1f %%\n", m.fraction_of_ru_in_mr*100);
        printf("# Percentage unscattered trans. in M_T = ");
        print_maybe('C', "%7.1f %%\n", m.fraction_of_tu_in_mt*100);
        printf("# \n");

@ @<Write first sphere info@>=
        printf("# Reflection sphere");
        if (m.baffle_r)
            printf(" has a baffle between sample and detector");
        else
            printf(" has no baffle between sample and detector");
        if (m.num_spheres > 0)
            printf("\n");
        else
            printf(" (ignored since no spheres used)\n");
        printf("#                      sphere diameter = %7.1f mm\n", m.d_sphere_r );
        printf("#                 sample port diameter = %7.1f mm\n",
        2*m.d_sphere_r*sqrt(m.as_r) );
        printf("#               entrance port diameter = %7.1f mm\n",
        2*m.d_sphere_r*sqrt(m.at_r) );
        printf("#               detector port diameter = %7.1f mm\n",
        2*m.d_sphere_r*sqrt(m.ad_r) );
        printf("#                 detector reflectance = %7.1f %%\n", m.rd_r*100 );
        printf("#                     wall reflectance = ");
        print_maybe('w', "%7.1f %%\n", m.rw_r*100);
        printf("#                 calibration standard = ");
        print_maybe('R', "%7.1f %%\n", m.rstd_r*100);
        printf("#\n");

@ @<Write second sphere info@>=
        printf("# Transmission sphere");
        if (m.baffle_t)
            printf(" has a baffle between sample and detector");
        else
            printf(" has no baffle between sample and detector");
        if (m.num_spheres > 0)
            printf("\n");
        else
            printf(" (ignored since no spheres used)\n");
        printf("#                      sphere diameter = %7.1f mm\n",
        m.d_sphere_t );
        printf("#                 sample port diameter = %7.1f mm\n",
        2*m.d_sphere_r*sqrt(m.as_t) );
        printf("#                  third port diameter = %7.1f mm\n",
        2*m.d_sphere_r*sqrt(m.at_t) );
        printf("#               detector port diameter = %7.1f mm\n",
        2*m.d_sphere_r*sqrt(m.ad_t) );
        printf("#                 detector reflectance = %7.1f %%\n", m.rd_t*100 );
        if (m.at_t == 0)
            printf("#    wall reflectance and cal standard = ");
        else
            printf("#                     wall reflectance = ");
        print_maybe('w', "%7.1f %%\n", m.rw_t*100);
        printf("#                 calibration standard = %7.1f %%", m.rstd_t*100 );
        if (m.at_t == 0)
            printf(" (ignored)");
        printf("\n");

@ @<Write measure and inversion info@>=
        printf("#\n");
        if (COLUMN_LABELS[0]=='\0'){
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
                    printf("# Something went wrong ... measures should be 1 to 7!\n");
                    break;
            }
        } else {
            int i;
            printf("# %d input columns with LABELS:", params);
            for (i=0; i<params; i++) {
                printf(" %c ", COLUMN_LABELS[i]);
            }
        }

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


        switch (m.num_spheres) {
        case 0:
            printf("# No sphere corrections were used");
            break;

        case 1:
            if (m.method == COMPARISON)
                printf("# No sphere corrections were needed");
            else
                printf("# Single sphere corrections were used");
            break;

        case 2:
            printf("# Double sphere corrections were used");
            break;
        }

        printf(" and light was incident at %d degrees from the normal",
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


@ Get the current line number by counting.

@<Definition for |get_current_line_number|@>=

long get_current_line_number(FILE *file) {
    long line_number = 0;
    long current_position = ftell(file); 
    fseek(file, 0, SEEK_SET);
    int c;
    while ((c = fgetc(file)) != EOF && ftell(file) < current_position) {
        if (c == '\n') {
            line_number++; 
        }
    }
    fseek(file, current_position, SEEK_SET); 
    return line_number + 1; 
}

@ Discard white space and dashes in the legend string

@<Definition for |remove_whitespace|@>=

void remove_whitespace(char *str) {
    int i, j = 0;
    for (i = 0; str[i] != '\0'; i++) {
        if (!isspace(str[i]) && str[i] != '-') {
            str[j++] = str[i];
        }
    }
    str[j] = '\0';
}

@ @<Definition for |remove_comment|@>=

void remove_comment(char *str) {
    int i;
    for (i = 0; str[i] != '\0'; i++) {
        if (str[i] == '#') {
            str[i] = '\0';
            break;
        }
    }
}

@ @<Definition for |remove_first_char|@>=
void remove_first_char(char *str) {
    int len = strlen(str);
    if (len > 0) {
        for (int i = 0; i < len; i++) {
            str[i] = str[i + 1];
        }
    }
}

@ @<Definition for |print_maybe|@>=
void print_maybe(char c, char *format, double x) {
    char *result = strchr(COLUMN_LABELS, c);
    if (result == NULL) 
        printf(format, x);
    else
        printf(" (varies with input row)\n");
}

@ @<Prototype for |Read_Data_Legend|@>=
    int Read_Data_Legend(FILE *fp)

@ @<Definition for |Read_Data_Legend|@>=
    @<Prototype for |Read_Data_Legend|@>
{
    int n=0;
    char c;

    skip_white(fp);
    if (fgets(COLUMN_LABELS, MAX_COLUMNS, fp) == NULL) {
        fprintf(stderr, "could not read Data Legend String in file\n");
        exit(EXIT_FAILURE);
    }

    remove_whitespace(COLUMN_LABELS);
    remove_comment(COLUMN_LABELS);
    c = COLUMN_LABELS[0];
    
    if (c=='1' || c=='2' || c=='3' || c=='4' || c=='5' || c=='6' || c=='7') {
        n = COLUMN_LABELS[0] - '0';
        COLUMN_LABELS[0] = '\0';
    } else {
        n = strlen(COLUMN_LABELS);
    }

    return n;
}

