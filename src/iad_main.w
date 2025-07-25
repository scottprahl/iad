@** iad command-line program.

Here is a relatively robust command-line utility that shows how
the iad and ad subroutines might be called.  It suffers because
it is written in \.{CWEB} and I used the macro expansion feature
instead of creating separate functions.  Oh well.

@ All the actual output for this web file goes into \.{iad\_main.c}

@(iad_main.c@>=

@<Include files for |main|@>@;

@<print version function@>@;
@<print usage function@>@;
@<stringdup together function@>@;
@<mystrtod function@>@;
@<seconds elapsed function@>@;
@<print error legend function@>@;
@<what\_char function@>@;
@<print long error function@>@;
@<print dot function@>@;
@<calculate coefficients function@>@;
@<parse string into array function@>@;
@<print results header function@>@;
@<Print results function@>@;

int main (int argc, char **argv)
{
    @<Declare variables for |main|@>@;

    @<Save the command line for use later@>@;
    @<Handle options@>@;

    Initialize_Measure(&m);
    @<Command-line changes to |m|@>@;

    Initialize_Result(m, &r, TRUE);

    if (cl_forward_calc != UNINITIALIZED) {
        @<Command-line changes to |r|@>@;
        @<Calculate and Print the Forward Calculation@>@;
        exit(EXIT_SUCCESS);
    }

    @<prepare file for reading@>@;

    if (process_command_line) {
        @<Count command-line measurements@>@;
        @<Calculate and write optical properties@>@;
        exit(EXIT_SUCCESS);
    }

    if (Read_Header (stdin, &m, &params) != 0)
        exit(EXIT_FAILURE);

    start_time = clock();
    while (Read_Data_Line (stdin, &m, &r, params) == 0) {
        @<Command-line changes to |m|@>@;
        @<Calculate and write optical properties@>@;
    }

    @<Generate and write grid@>@;

    if (cl_verbosity>0) fprintf(stderr,"\n\n");
    if (any_error && cl_verbosity>1) print_error_legend();
    exit(EXIT_SUCCESS);
}

@ The first two defines are to stop Visual C++ from silly complaints
@<Include files for |main|@>=
#define _CRT_SECURE_NO_WARNINGS
#define _CRT_NONSTDC_NO_WARNINGS

#define NO_SLIDES                 0
#define ONE_SLIDE_ON_TOP          1
#define TWO_IDENTICAL_SLIDES      2
#define ONE_SLIDE_ON_BOTTOM       3
#define ONE_SLIDE_NEAR_SPHERE     4
#define ONE_SLIDE_NOT_NEAR_SPHERE 5

#define MR_IS_ONLY_RD        1
#define MT_IS_ONLY_TD        2
#define NO_UNSCATTERED_LIGHT 3

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>

#include "ad_globl.h"
#include "ad_prime.h"
#include "iad_type.h"
#include "iad_pub.h"
#include "iad_io.h"
#include "iad_calc.h"
#include "iad_util.h"
#include "version.h"
#include "mc_lost.h"
#include "ad_frsnl.h"

@ @<Declare variables for |main|@>=
    struct measure_type m;
    struct invert_type r;
    char *g_out_name = NULL;
    char *g_grid_name = NULL;
    int c;

    long n_photons = 100000;
    int MAX_MC_iterations = 19;
    int any_error = 0;
    int process_command_line = 0;
    int params = 0;

    int cl_quadrature_points = UNINITIALIZED;
    int cl_verbosity = 2;

    double cl_forward_calc= UNINITIALIZED;
    double cl_grid_calc   = UNINITIALIZED;
    double cl_default_a   = UNINITIALIZED;
    double cl_default_g   = UNINITIALIZED;
    double cl_default_b   = UNINITIALIZED;
    double cl_default_mua = UNINITIALIZED;
    double cl_default_mus = UNINITIALIZED;
    double cl_default_musp= UNINITIALIZED;
    double cl_tolerance   = UNINITIALIZED;
    double cl_slide_OD    = UNINITIALIZED;

    double cl_cos_angle   = UNINITIALIZED;
    double cl_beam_d      = UNINITIALIZED;
    double cl_sample_d    = UNINITIALIZED;
    double cl_sample_n    = UNINITIALIZED;
    double cl_slide_d     = UNINITIALIZED;
    double cl_slide_n     = UNINITIALIZED;
    double cl_slides      = UNINITIALIZED;
    double cl_default_fr  = UNINITIALIZED;
    double cl_rstd_t      = UNINITIALIZED;
    double cl_rstd_r      = UNINITIALIZED;
    double cl_baffle_r    = UNINITIALIZED;
    double cl_baffle_t    = UNINITIALIZED;
    double cl_ru_fraction = UNINITIALIZED;
    double cl_tu_fraction = UNINITIALIZED;
    double cl_lambda      = UNINITIALIZED;
    double cl_rwall_r     = UNINITIALIZED;
    double cl_rwall_t     = UNINITIALIZED;

    double cl_search      = UNINITIALIZED;
    double cl_mus0        = UNINITIALIZED;
    double cl_mus0_pwr    = UNINITIALIZED;
    double cl_mus0_lambda = UNINITIALIZED;

    double cl_UR1         = UNINITIALIZED;
    double cl_UT1         = UNINITIALIZED;
    double cl_Tc          = UNINITIALIZED;

    double cl_method      = UNINITIALIZED;
    int    cl_num_spheres = UNINITIALIZED;
    double cl_sphere_one[5] = {UNINITIALIZED, UNINITIALIZED, UNINITIALIZED,
                             UNINITIALIZED, UNINITIALIZED };
    double cl_sphere_two[5] = {UNINITIALIZED, UNINITIALIZED, UNINITIALIZED,
                             UNINITIALIZED, UNINITIALIZED };
    double cl_wave_limit[2] = {UNINITIALIZED, UNINITIALIZED};

    clock_t start_time=clock();
    char command_line_options[] =
             "1:2:a:A:b:B:c:C:d:D:e:E:f:F:g:G:hH:i:j:Jl:L:M:n:N:o:p:q:r:R:s:S:t:T:u:vV:w:W:x:Xz";
    char *command_line = NULL;

@ I want to include the command line to the output file.  To do this, we need to save
the entire thing before the options get processed.  The extra |+1| in the total length
calculation is for the space character between options. Finally, we need to reset
|optind| to 1 to start |getopt()| processing from the beginning.  It should be noted
that this strips any quotes from the command line.

@<Save the command line for use later@>=
{
    size_t command_line_length = 0;
    for (int i = 0; i < argc; ++i) {
        command_line_length += strlen(argv[i]) + 3;
    }

    command_line = (char *)malloc(command_line_length);
    if (command_line == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    strcpy(command_line, "");
    for (int i = 0; i < argc; ++i) {
        if (strchr(argv[i], ' ') != NULL) {
            strcat(command_line, "'");
            strcat(command_line, argv[i]);
            strcat(command_line, "' ");
        } else {
            strcat(command_line, argv[i]);
            strcat(command_line, " ");
        }
    }

    optind = 1;
}

@*1 Handling the command-line options.

@<Handle options@>=
    while ((c = getopt(argc, argv, command_line_options)) != EOF) {
        int n;
        char cc;
        char *tmp_str = NULL;

        switch (c) {

            case '1':
                tmp_str = strdup(optarg);
                parse_string_into_array(optarg, cl_sphere_one, 5);
                if (cl_sphere_one[4] == UNINITIALIZED) {
                    fprintf(stderr, "Error in the command-line argument for -1\n");
                    fprintf(stderr, "    the current argument is '%s' but it must have 5 terms: ", tmp_str);
                    fprintf(stderr, "'d_sphere d_sample d_entrance d_detector r_wall'\n");
                    exit(EXIT_FAILURE);
                }
                break;

            case '2':
                tmp_str = strdup(optarg);
                parse_string_into_array(optarg, cl_sphere_two, 5);
                if (cl_sphere_two[4] == UNINITIALIZED) {
                    fprintf(stderr, "Error in the command-line argument for -2\n");
                    fprintf(stderr, "    the current argument is '%s' but it must have 5 terms: ", tmp_str);
                    fprintf(stderr, "'d_sphere d_sample d_third d_detector r_wall'\n");
                    exit(EXIT_FAILURE);
                }
                break;

            case 'a':
                cl_default_a = my_strtod(optarg);
                if (cl_default_a<0 || cl_default_a>1) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    albedo '-a %s'\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;

            case 'A':
                cl_default_mua = my_strtod(optarg);
                if (cl_default_mua < 0) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    absorption '-A %s'\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;

            case 'b':
                cl_default_b = my_strtod(optarg);
                if (cl_default_b < 0) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    optical thickness '-b %s'\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;

            case 'B':
                cl_beam_d = my_strtod(optarg);
                if (cl_beam_d < 0) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    beam diameter '-B %s'\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;

            case 'c':
                cl_ru_fraction = my_strtod(optarg);
                if (cl_ru_fraction<0.0 || cl_ru_fraction>1.0) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    unscattered refl fraction '-c %s'\n", optarg);
                    fprintf(stderr, "    must be between 0 and 1\n");
                    exit(EXIT_SUCCESS);
                }
                break;

            case 'C':
                cl_tu_fraction = my_strtod(optarg);
                if (cl_tu_fraction < 0.0 || cl_tu_fraction > 1.0) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    unscattered trans fraction '-C %s'\n", optarg);
                    fprintf(stderr, "    must be between 0 and 1\n");
                    exit(EXIT_SUCCESS);
                }
                break;

            case 'd':
                cl_sample_d = my_strtod(optarg);
                if (cl_sample_d < 0) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    sample thickness '-d %s'\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;

            case 'D':
                cl_slide_d = my_strtod(optarg);
                if (cl_slide_d < 0) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    slide thickness '-D %s'\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;

            case 'e':
                cl_tolerance = my_strtod(optarg);
                if (cl_tolerance < 0) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    error tolerance '-e %s'\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;

            case 'E':
                cl_slide_OD = my_strtod(optarg);
                if (cl_slide_OD < 0) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    slide optical depth '-E %s'\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;

            case 'f':
                cl_default_fr = my_strtod(optarg);
                if (cl_default_fr < 0.0 || cl_default_fr > 1.0) {
                    fprintf(stderr, "Error in the command-line argument: ");
                    fprintf(stderr, "'-f %s' The argument must be between 0 and 1.\n", optarg);
                    exit(EXIT_SUCCESS);
                }
                break;

            case 'F':
                /* initial digit means that mus is constant */
                if (isdigit(optarg[0])) {
                    cl_default_mus = my_strtod(optarg);
                    if (cl_default_mus < 0) {
                        fprintf(stderr, "Error in the command line\n");
                        fprintf(stderr, "    mus '-F %s'\n", optarg);
                        exit(EXIT_FAILURE);
                    }
                    break;
                }

                /* should be a string like 'P 1000 1.2 -1.8' */
                n=sscanf(optarg, "%c %lf %lf %lf",&cc, &cl_mus0_lambda, &cl_mus0, &cl_mus0_pwr);

                if (n != 4 || (cc != 'P' && cc != 'p')) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    bad -F option. '-F %s'\n", optarg);
                    fprintf(stderr, "    -F 1.0              for mus =1.0\n");
                    fprintf(stderr, "    -F 'P 500 1.0 -1.3' for mus =1.0*(lambda/500)^(-1.3)\n");
                    exit(EXIT_FAILURE);
                }

                break;

            case 'g':
                cl_default_g = my_strtod(optarg);
                if (cl_default_g < -1 || cl_default_g > 1) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    anisotropy '-g %s'\n", optarg);
                    exit(EXIT_FAILURE);
                }
                if (cl_default_g == -1) cl_default_g=-MAX_ABS_G;
                if (cl_default_g == 1) cl_default_g= MAX_ABS_G;
                break;

            case 'G':
                if (optarg[0]=='0')
                    cl_slides = NO_SLIDES;
                else if (optarg[0]=='2')
                    cl_slides = TWO_IDENTICAL_SLIDES;
                else if (optarg[0]=='t' || optarg[0]=='T')
                    cl_slides = ONE_SLIDE_ON_TOP;
                else if (optarg[0]=='b' || optarg[0]=='B')
                    cl_slides = ONE_SLIDE_ON_BOTTOM;
                else if (optarg[0]=='n' || optarg[0]=='N')
                    cl_slides = ONE_SLIDE_NEAR_SPHERE;
                else if (optarg[0]=='f' || optarg[0]=='F')
                    cl_slides = ONE_SLIDE_NOT_NEAR_SPHERE;
                else {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    Argument for '-G %s' must be \n", optarg);
                    fprintf(stderr, "    't' --- light always hits top slide first\n");
                    fprintf(stderr, "    'b' --- light always hits bottom slide first\n");
                    fprintf(stderr, "    'n' --- slide always closest to sphere\n");
                    fprintf(stderr, "    'f' --- slide always farthest from sphere\n");
                    exit(EXIT_FAILURE);
                }
                break;

            case 'H':
                if (optarg[0]=='0') {
                    cl_baffle_r = 0;
                    cl_baffle_t = 0;
                } else if (optarg[0]=='1') {
                    cl_baffle_r = 1;
                    cl_baffle_t = 0;
                } else if (optarg[0]=='2') {
                    cl_baffle_r = 0;
                    cl_baffle_t = 1;
                } else if (optarg[0]=='3') {
                    cl_baffle_r = 1;
                    cl_baffle_t = 1;
                } else {
                    fprintf(stderr, "Error in the command-line -H argument\n");
                    fprintf(stderr, "    argument is '%s', but ", optarg);
                    fprintf(stderr, "must be 0, 1, 2, or 3\n");
                    exit(EXIT_FAILURE);
                }

            case 'i':
                cl_cos_angle = my_strtod(optarg);
                if (cl_cos_angle<0 || cl_cos_angle>90) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    incident angle '-i %s'\n", optarg);
                    fprintf(stderr, "    must be between 0 and 90 degrees\n");
                    exit(EXIT_FAILURE);
                }
                cl_cos_angle = cos(cl_cos_angle*M_PI/180.0);
                break;

            case 'j':
                cl_default_musp = my_strtod(optarg);
                if (cl_default_musp < 0) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    musp '-j %s'\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;

            case 'J':
                cl_grid_calc = 1;
                break;

            case 'l':
                tmp_str = strdup(optarg);
                parse_string_into_array(optarg, cl_wave_limit, 2);
                break;

            case 'L':
                cl_lambda = my_strtod(optarg);
                break;

            case 'M':
                MAX_MC_iterations = (int) my_strtod(optarg);
                if (MAX_MC_iterations < 0 || MAX_MC_iterations > 50) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    MC iterations '-M %s'\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;

            case 'n':
                cl_sample_n = my_strtod(optarg);
                if (cl_sample_n < 0.1 || cl_sample_n > 10) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    slab index '-n %s'\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;

            case 'N':
                cl_slide_n = my_strtod(optarg);
                if (cl_slide_n < 0.1 || cl_slide_n > 10) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    slide index '-N %s'\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;

            case 'o':
                g_out_name = strdup(optarg);
                break;

            case 'p':
                n_photons = (long) my_strtod(optarg);
                break;

            case 'q':
                cl_quadrature_points = (int) my_strtod(optarg);
                if (cl_quadrature_points % 4 != 0) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    '-q %s'\n", optarg);
                    fprintf(stderr, "    Quadrature points must be a multiple of 4\n");
                    exit(EXIT_FAILURE);
                }
                if ((cl_cos_angle != UNINITIALIZED) && (cl_quadrature_points % 12 != 0)) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    '-q %s'\n", optarg);
                    fprintf(stderr, "    Quadrature points must be multiple of 12 for oblique incidence\n");
                    exit(EXIT_FAILURE);
                }
                break;

            case 'r':
                cl_UR1 = my_strtod(optarg);
                process_command_line = 1;
                if (cl_UR1 < 0 || cl_UR1 > 1) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    UR1 value '-r %s'\n", optarg);
                    fprintf(stderr, "    must be between 0 and 1\n");
                    exit(EXIT_FAILURE);
                }
                break;

            case 'R':
                cl_rstd_r = my_strtod(optarg);
                if (cl_rstd_r < 0 || cl_rstd_r > 1) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    Rstd value '-R %s'\n", optarg);
                    fprintf(stderr, "    must be between 0 and 1\n");
                    exit(EXIT_FAILURE);
                }
                break;

            case 's':
                cl_search = (int) my_strtod(optarg);
                break;

            case 'S':
                cl_num_spheres = (int) my_strtod(optarg);
                if (cl_num_spheres != 0 && cl_num_spheres != 1 && cl_num_spheres !=2) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    sphere number '-S %s'\n", optarg);
                    fprintf(stderr, "    must be 0, 1, or 2\n");
                    exit(EXIT_FAILURE);
                }
                break;

            case 't':
                cl_UT1 = my_strtod(optarg);
                if (cl_UT1 < 0 || cl_UT1 > 1) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    UT1 value '-t %s'\n", optarg);
                    fprintf(stderr, "    must be between 0 and 1\n");
                    exit(EXIT_FAILURE);
                }
                process_command_line=1;
                break;

            case 'T':
                cl_rstd_t = my_strtod(optarg);
                if (cl_rstd_t < 0 || cl_rstd_t > 1) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    transmission standard '-T %s'\n", optarg);
                    fprintf(stderr, "    must be between 0 and 1\n");
                    exit(EXIT_FAILURE);
                }
                break;

            case 'u':
                cl_Tc = my_strtod(optarg);
                if (cl_Tc < 0 || cl_Tc > 1) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    unscattered transmission '-u %s'\n", optarg);
                    fprintf(stderr, "    must be between 0 and 1\n");
                    exit(EXIT_FAILURE);
                }
                process_command_line=1;
                break;

            case 'v':
                print_version(cl_verbosity);
                exit(EXIT_SUCCESS);
                break;

            case 'V':
                cl_verbosity = my_strtod(optarg);
                break;

            case 'w':
                cl_rwall_r = my_strtod(optarg);
                if (cl_rwall_r < 0 || cl_rwall_r > 1) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    refl sphere wall '-w %s'\n", optarg);
                    fprintf(stderr, "    must be between 0 and 1\n");
                    exit(EXIT_FAILURE);
                }
                break;

            case 'W':
                cl_rwall_t = my_strtod(optarg);
                if (cl_rwall_t < 0 || cl_rwall_r > 1) {
                    fprintf(stderr, "Error in the command line\n");
                    fprintf(stderr, "    trans sphere wall '-w %s'\n", optarg);
                    fprintf(stderr, "    must be between 0 and 1\n");
                    exit(EXIT_FAILURE);
                }
                break;

            case 'x':
                Set_Debugging( (int) my_strtod(optarg));
                break;

            case 'X':
                cl_method=COMPARISON;
                break;

            case 'z':
                cl_forward_calc = 1;
                process_command_line=1;
                break;

            default:
                fprintf(stderr, "unknown option '%c'\n", c);
                /* fall through */

            case 'h':
                print_usage();
                exit(EXIT_SUCCESS);
        }
    }

    argc -= optind;
    argv += optind;

@*1 The forward calculation.

We are doing a forward calculation.  This takes care of setting up the
sample scattering and absorption coefficients when they are specified on
the command line.  It also deals with the case when the reduced scattering
coefficient is specified.  We carefully set up the default values (scattering
is 1/mm, absorption is 0/mm, and anisotropy is 0.0).

@<Calculate and Print the Forward Calculation@>=
    double temp_mus = 1;
    double temp_mua = 0;
    r.g = 0;

    if (cl_default_mus != UNINITIALIZED)
        temp_mus = cl_default_mus;

    if (cl_default_mua != UNINITIALIZED)
        temp_mua = cl_default_mua;

    if (cl_default_g!= UNINITIALIZED)
        r.g = cl_default_g;

    if (cl_default_musp != UNINITIALIZED)
        temp_mus = cl_default_musp / (1 - r.g);

@ We still need to set the albedo and optical depth appropriately.
Obviously when the |-a| switch is used then the albedo must be set
to equal to |cl_default_a|.  It is also important to adjust the scattering
and absorption coefficients in this case so that they match the albedo.

This is easy when both the optical thickness and the physical thickness of
the sample are both given.  If this not the case, then we just assume the
scattering coefficient is correct and adjust the absorption coefficient to
be whatever it needs to be.  However in the case of |a=0| then we set the
scattering coefficient to zero and the absorption coefficient to 1/mm.

@<Calculate and Print the Forward Calculation@>=
    if (cl_default_a != UNINITIALIZED) {
        r.a = cl_default_a;

        if (cl_default_b != UNINITIALIZED && cl_sample_d != UNINITIALIZED) {
            temp_mus = cl_default_a * cl_default_b / cl_sample_d;
            temp_mua = cl_default_b / cl_sample_d - temp_mus;
        } else {
            if (cl_default_a == 0) {
                temp_mus = 0;
                temp_mua = 1;
            } else
                temp_mua = temp_mus / cl_default_a - temp_mus;
        }
    } else
        r.a = temp_mus/(temp_mus+temp_mua);

@ By this point the scattering coefficient and scattering coefficient are
set correctly.  These values will be used unless |-b| has been used to set
the optical thickness --- |cl_default_b|.  If the sample thickness is not
specified, then the thickness is assumed infinite.

@<Calculate and Print the Forward Calculation@>=
    if (cl_default_b != UNINITIALIZED) {
        r.b = cl_default_b;
    } else {
        if (cl_sample_d == UNINITIALIZED)
            r.b = HUGE_VAL;
        else
            r.b = (temp_mus+temp_mua)*cl_sample_d;
    }

@ @<Calculate and Print the Forward Calculation@>=
    r.slab.a = r.a;
    r.slab.b = r.b;
    r.slab.g = r.g;

    {
        double mu_s, mu_sp, mu_a, m_r, m_t;
        if (MAX_MC_iterations==0 || m.num_spheres==0) {
            Calculate_MR_MT(m, r, MC_NONE, TRUE, &m_r, &m_t);
        } else {
            Calculate_MR_MT(m, r, MC_REDO, TRUE, &m_r, &m_t);
        }

        Calculate_Mua_Musp(m, r, &mu_s, &mu_sp, &mu_a);
        if (cl_verbosity>0) {
            Write_Header (m, r, -1, command_line);
            print_results_header(stdout);
        }
        print_optical_property_result(stdout,m,r,m_r,m_t,mu_a,mu_sp,0);
    }

@*1 Calculating a grid for graphing.

We will start simple.  Just vary $a'$ and $b'$.

@<Generate and write grid@>=
if (cl_grid_calc != UNINITIALIZED) {
    double m_r, m_t, aprime, bprime, g;
    double aa[] = {0, 0.8, 0.9, 0.95, 0.98, 0.99, 1.0};
    double bb[] = {0, 0.2, 0.5, 1.0,  3.0, 10.0, 100};
    int i, j;
    int count=0;

    FILE *grid;

    grid = fopen(g_grid_name, "w");
    if (grid == NULL) {
        fprintf(stderr, "Could not open grid file '%s' for output\n", g_out_name);
        exit(EXIT_FAILURE);
    }

    m.ur1_lost = 0;
    m.uru_lost = 0;
    m.ut1_lost = 0;
    m.utu_lost = 0;

    if (r.default_g != UNINITIALIZED) {
        g = r.default_g;
    } else if (r.found) {
        g = r.slab.g;
    } else {
        g = 0;
    }
    fprintf(grid, "# %s (g=%6.4f)\n", command_line, g);
    fprintf(grid, "#    a'          b'          g           M_R         M_T\n");
    fprintf(stderr, "\ndoing grid calculation\n");
    for (i=0; i<7; i++) {
        aprime = aa[i];
        for (j=0; j<7; j++) {
            fprintf(stderr, "*");
            bprime = bb[j];
            r.a = aprime / (1 - g + aprime * g);
            r.b = bprime / (1 - r.slab.a * g);
            r.g = g;
            r.slab.a = r.a;
            r.slab.b = r.b;
            r.slab.g = r.g;
            if (MAX_MC_iterations==0) {
                Calculate_MR_MT(m, r, MC_NONE, TRUE, &m_r, &m_t);
            } else {
                Calculate_MR_MT(m, r, MC_REDO, TRUE, &m_r, &m_t);
            }
            count++;
            if (count % 10 == 0)
                fprintf(stderr, " ");
            if (count % 50 == 0)
                fprintf(stderr, "\n");

            fprintf(grid, "%10.5f, %10.5f, %10.5f, %10.5f, %10.5f\n", \
                    aprime, bprime, g, m_r, m_t);
        }
    }
    fclose(grid);
    fprintf(stderr, "\n");
}

@ Make sure that the file is not named '-' and warn about too many files

@<prepare file for reading@>=
    if (argc > 1) {
        fprintf(stderr, "Only a single file can be processed at a time\n");
        fprintf(stderr, "try 'apply iad file1 file2 ... fileN'\n");
        exit(EXIT_FAILURE);
    }

    if (argc == 1 && strcmp(argv[0],"-")!=0) {  /* filename exists and != "-" */
        int n;
        char *base_name, *rt_name;
        base_name = strdup(argv[0]);
        n = (int) (strlen(base_name) - strlen(".rxt"));

        if (n>0 && strstr(base_name+n,".rxt") != NULL)
            base_name[n]='\0';

        rt_name = strdup_together(base_name,".rxt");

        if (freopen(argv[0],"r",stdin)==NULL &&
            freopen(rt_name,"r",stdin)==NULL) {
            fprintf(stderr, "Could not open either '%s' or '%s'\n",
                            argv[0], rt_name);
            exit(EXIT_FAILURE);
        }

        if (g_out_name==NULL)
            g_out_name=strdup_together(base_name,".txt");

        if (g_grid_name==NULL)
            g_grid_name=strdup_together(base_name,".grid");

        free(rt_name);
        free(base_name);
        process_command_line = 0;
    }

    if (g_out_name!=NULL) {
        if (freopen(g_out_name,"w",stdout)==NULL) {
            fprintf(stderr, "Could not open file '%s' for output\n", g_out_name);
            exit(EXIT_FAILURE);
        }
    }


@ Need to explicitly reset |r.search| each time through the loop,
because it will get altered by the calculation process.  This also allows the
command line to overwrite the reflection or transmission value from the command line
with something like {\tt -r 0} or {\tt -t 0}.

We also want to be able to let different lines have different constraints.  In particular
consider the file \.{newton.tst}.  In that file the first two rows contain
three real measurements and the last two have the collimated transmission
explicitly set to zero --- in other words there are really only two
measurements.

@<Calculate and write optical properties@>=
{
    @<Local Variables for Calculation@>@;

    if (cl_wave_limit[0] != UNINITIALIZED) {
        if (m.lambda != 0) {
            if (m.lambda < cl_wave_limit[0])
                skip = TRUE;
            if (m.lambda > cl_wave_limit[1])
                skip = TRUE;
        }
    }

    if (Debug(DEBUG_ANY) && !skip) {
        fprintf(stderr, "\n-------------------NEXT DATA POINT---------------------\n");
        if (m.lambda != 0)
            fprintf(stderr, "lambda=%6.1f ", m.lambda);
        fprintf(stderr, "MR=%8.5f MT=%8.5f\n\n", m.m_r, m.m_t);
        if (skip)
            fprintf(stderr, "skipping, wavelength out of range.\n");
    }

    if (!skip) {
        rt_total++;
        Initialize_Result(m, &r, FALSE);

        @<Command-line changes to |r|@>@;
        @<Warn and quit for bad options@>@;
        @<Write Header @>@;

        m.ur1_lost = 0;
        m.uru_lost = 0;
        m.ut1_lost = 0;
        m.utu_lost = 0;

        Inverse_RT (m, &r);

        calculate_coefficients(m,r,&LR,&LT,&mu_sp,&mu_a);
        @<Improve result using Monte Carlo@>@;

        calculate_coefficients(m,r,&LR,&LT,&mu_sp,&mu_a);
        print_optical_property_result(stdout,m,r,LR,LT,mu_a,mu_sp,rt_total);

        if (r.error != IAD_NO_ERROR)
            any_error = 1;

        if (Debug(DEBUG_ANY))
            print_long_error(r.error);
        else
            print_dot(start_time, r.error, mc_total, TRUE, cl_verbosity);
    }
}

@   @<Local Variables for Calculation@>=
    static int rt_total = 0;
    static int mc_total = 0;

    double ur1=0;
    double ut1=0;
    double uru=0;
    double utu=0;
    double mu_a=0;
    double mu_sp=0;
    double LR=0;
    double LT=0;
    int skip = FALSE;

@ @<Command-line changes to |r|@>=

    if (cl_quadrature_points != UNINITIALIZED)
        r.method.quad_pts = cl_quadrature_points;
    else
        r.method.quad_pts = 8;

    if (cl_default_a != UNINITIALIZED)
        r.default_a = cl_default_a;

    if (cl_default_mua != UNINITIALIZED) {
        r.default_mua = cl_default_mua;
        if (cl_sample_d != UNINITIALIZED)
            r.default_ba = cl_default_mua * cl_sample_d;
        else
            r.default_ba = cl_default_mua * m.slab_thickness;
    }

    if (cl_default_b != UNINITIALIZED)
        r.default_b = cl_default_b;

    if (cl_default_g != UNINITIALIZED)
        r.default_g = cl_default_g;

    if (cl_tolerance != UNINITIALIZED) {
        r.tolerance = cl_tolerance;
        r.MC_tolerance = cl_tolerance;
    }

    if (cl_mus0 != UNINITIALIZED) {
        if (m.lambda != 0) {
            cl_default_mus = cl_mus0 * pow(m.lambda/cl_mus0_lambda,cl_mus0_pwr);
        } else {
            fprintf(stderr, "Seems like you want to constrain scattering to a power law.\n");
            fprintf(stderr, "Unfortunately, there is no wavelength so this cannot be done.\n");
        }
    }

    if (cl_default_mus != UNINITIALIZED) {
        r.default_mus = cl_default_mus;
        if (cl_sample_d != UNINITIALIZED)
            r.default_bs = cl_default_mus * cl_sample_d;
        else
            r.default_bs = cl_default_mus * m.slab_thickness;
    }

    if (cl_default_musp != UNINITIALIZED) {
        if (cl_default_g != UNINITIALIZED)
            r.default_mus = cl_default_musp / (1.0 - cl_default_g);
        else
            r.default_mus = cl_default_musp;

        if (cl_sample_d != UNINITIALIZED)
            r.default_bs = r.default_mus * cl_sample_d;
        else
            r.default_bs = r.default_mus * m.slab_thickness;
    }

    if (cl_search != UNINITIALIZED)
        r.search = cl_search;

@ @<Write Header @>=

if (rt_total==1 && cl_verbosity>0) {
    Write_Header (m, r, params, command_line);
    if (MAX_MC_iterations > 0) {
        if (n_photons>=0)
           fprintf(stdout,"#  Photons used to estimate lost light =   %ld\n",n_photons);
        else
           fprintf(stdout,"#     Time used to estimate lost light =   %ld ms\n",-n_photons);
    } else
        fprintf(stdout,"#  Photons used to estimate lost light =   0\n");

    fprintf(stdout,"#\n");

    print_results_header(stdout);
}

@*1 Monte Carlo light loss.
Use Monte Carlo to figure out how much light leaks out.  We use
the sphere corrected values as the starting values and only do try
Monte Carlo when spheres are used, the albedo unknown or non-zero, and
there has been no error.  The sphere parameters must be known because
otherwise the beam size and the port size are unknown.

@<Improve result using Monte Carlo@>=

if (m.num_spheres > 0) {

    if (Debug(DEBUG_LOST_LIGHT)) {
        print_results_header(stderr);
        print_optical_property_result(stderr,m,r,LR,LT,mu_a,mu_sp,rt_total);
    }

    while (r.MC_iterations < MAX_MC_iterations) {
        double last_mu_sp, last_mu_a, last_final_distance;
        double current_ur1_lost, current_ut1_lost, current_uru_lost, current_utu_lost;
        double diff_ur1_lost, diff_ut1_lost, diff_uru_lost, diff_utu_lost;
        double factor = 0.8;
        int too_much_lost;
        double tol = r.MC_tolerance;

        calculate_coefficients(m,r,&LR,&LT,&mu_sp,&mu_a);
        last_mu_sp = mu_sp;
        last_mu_a  = mu_a;
        last_final_distance = r.final_distance;

        if (Debug(DEBUG_ITERATIONS) || Debug(DEBUG_A_LITTLE)) {
            fprintf(stderr, "\n------------- Monte Carlo Iteration %d -----------------\n", r.MC_iterations+1);
        }
        MC_Lost(m, r, n_photons, &ur1, &ut1, &uru, &utu,
                &current_ur1_lost, &current_ut1_lost, &current_uru_lost, &current_utu_lost);
        diff_ur1_lost= current_ur1_lost - m.ur1_lost;
        diff_uru_lost= current_uru_lost - m.uru_lost;
        diff_ut1_lost= current_ut1_lost - m.ut1_lost;
        diff_utu_lost= current_utu_lost - m.utu_lost;

        if (diff_ur1_lost > 0.001 || diff_ut1_lost > 0.001)
            too_much_lost = 1;
        else
            too_much_lost = 0;

        m.ur1_lost +=  factor * diff_ur1_lost;
        m.uru_lost +=  factor * diff_uru_lost;
        m.ut1_lost +=  factor * diff_ut1_lost;
        m.utu_lost +=  factor * diff_utu_lost;

        mc_total++;
        r.MC_iterations++;

        Inverse_RT (m, &r);
        calculate_coefficients(m,r,&LR,&LT,&mu_sp,&mu_a);

        if (0) {
            fprintf(stderr, "%2d %2d %2d | %7.4f %7.4f %7.4f | %7.4f %7.4f %7.4f\n",
                r.MC_iterations, too_much_lost, r.found,
                m.m_r, current_ur1_lost, m.ur1_lost,
                m.m_t, current_ut1_lost, m.ut1_lost);
        }

        if (Debug(DEBUG_LOST_LIGHT))
            print_optical_property_result(stderr,m,r,LR,LT,mu_a,mu_sp,rt_total);
        else
            print_dot(start_time, r.error, mc_total, FALSE, cl_verbosity);

        if (r.found) {
            if (fabs(last_mu_a-mu_a)>tol) {
                if (Debug(DEBUG_ITERATIONS))
                    fprintf(stderr, "Repeat MC because mua is still changing\n");
                continue;
            }

            if (fabs(last_mu_sp-mu_sp)>tol) {
                if (Debug(DEBUG_ITERATIONS))
                    fprintf(stderr, "Repeat MC because musp is still changing\n");
                continue;
            }

            if (too_much_lost){
                if (Debug(DEBUG_ITERATIONS))
                    fprintf(stderr, "Repeat MC because mua and musp are still changing\n");
                continue;
            }

            if (Debug(DEBUG_ITERATIONS))
                fprintf(stderr, "found!\n");
            break;

        } else {
            if (last_final_distance - r.final_distance < r.MC_tolerance) {
                if (Debug(DEBUG_ITERATIONS))
                    fprintf(stderr, "MC does not make things better\n");
                break;
            } else {
                if (Debug(DEBUG_ITERATIONS))
                    fprintf(stderr, "Repeat MC because distance is reduced\n");
            }
        }
    }
}

@ Stuff the command-line arguments that should be constant over the entire
inversion process into the measurement record  and
set up the result record to handle the arguments properly so that the optical
properties can be determined.

@<Command-line changes to |m|@>=

    if (cl_cos_angle != UNINITIALIZED) {
        m.slab_cos_angle = cl_cos_angle;
        if (cl_quadrature_points == UNINITIALIZED)
            cl_quadrature_points = 12;

        if (cl_quadrature_points != 12 * (cl_quadrature_points / 12)) {
            fprintf(stderr, "If you use the -i option to specify an oblique incidence angle, then\n");
            fprintf(stderr, "the number of quadrature points must be a multiple of 12\n");
            exit(EXIT_SUCCESS);
        }
    }

    if (cl_sample_n != UNINITIALIZED)
        m.slab_index = cl_sample_n;

    if (cl_slide_n != UNINITIALIZED) {
        m.slab_bottom_slide_index = cl_slide_n;
        m.slab_top_slide_index    = cl_slide_n;
    }

    if (cl_slide_OD != UNINITIALIZED) {
        m.slab_bottom_slide_b = cl_slide_OD;
        m.slab_top_slide_b    = cl_slide_OD;
    }

    if (cl_sample_d != UNINITIALIZED)
        m.slab_thickness = cl_sample_d;

    if (cl_beam_d != UNINITIALIZED)
        m.d_beam = cl_beam_d;

    if (cl_slide_d != UNINITIALIZED) {
        m.slab_bottom_slide_thickness = cl_slide_d;
        m.slab_top_slide_thickness    = cl_slide_d;
    }

    if (cl_slides == NO_SLIDES  ) {
        m.slab_bottom_slide_index     = 1.0;
        m.slab_bottom_slide_thickness = 0.0;
        m.slab_top_slide_index        = 1.0;
        m.slab_top_slide_thickness    = 0.0;
    }

    if (cl_slides == ONE_SLIDE_ON_TOP ||
        cl_slides == ONE_SLIDE_NEAR_SPHERE) {
        m.slab_bottom_slide_index     = 1.0;
        m.slab_bottom_slide_thickness = 0.0;
    }

    if (cl_slides == ONE_SLIDE_ON_BOTTOM ||
        cl_slides == ONE_SLIDE_NOT_NEAR_SPHERE) {
        m.slab_top_slide_index        = 1.0;
        m.slab_top_slide_thickness    = 0.0;
    }

    if (cl_slides == ONE_SLIDE_NEAR_SPHERE ||
        cl_slides == ONE_SLIDE_NOT_NEAR_SPHERE)
        m.flip_sample = 1;
    else
        m.flip_sample = 0;

    if (cl_method != UNINITIALIZED)
        m.method = (int) cl_method;

    if (cl_rstd_r != UNINITIALIZED) {
        m.rstd_r = cl_rstd_r;
        m.rstd_t = cl_rstd_r;
    }

    if (cl_rstd_t != UNINITIALIZED) {
        m.rstd_t = cl_rstd_t;
        if (cl_rstd_r == UNINITIALIZED)
            m.rstd_r = cl_rstd_t;
    }

    if (cl_rwall_r != UNINITIALIZED) {
        if (cl_sphere_one[0] != UNINITIALIZED) {
            fprintf(stderr, "-w is overridden by -1 option. omit.\n");
            exit(EXIT_FAILURE);
        }
        m.rw_r = cl_rwall_r;
    }

    if (cl_rwall_t != UNINITIALIZED) {
        if (cl_sphere_one[0] != UNINITIALIZED || cl_sphere_one[1] != UNINITIALIZED) {
            fprintf(stderr, "-W is overridden by -1 and -2 options. omit.");
            exit(EXIT_FAILURE);
        }
        m.rw_t = cl_rwall_t;
    }

    if (cl_sphere_one[0] != UNINITIALIZED) {
        double d_sample_r, d_third_r, d_detector_r;

        m.d_sphere_r     = cl_sphere_one[0];
        d_sample_r       = cl_sphere_one[1];
        d_third_r        = cl_sphere_one[2];
        d_detector_r     = cl_sphere_one[3];
        m.rw_r           = cl_sphere_one[4];

        m.as_r = sqr(d_sample_r   / m.d_sphere_r / 2);
        m.at_r = sqr(d_third_r    / m.d_sphere_r / 2);
        m.ad_r = sqr(d_detector_r / m.d_sphere_r / 2);

        m.aw_r = 1.0 - m.as_r - m.at_r - m.ad_r;

        m.d_sphere_t = m.d_sphere_r;
        m.as_t       = m.as_r;
        m.at_t       = m.at_r;
        m.ad_t       = m.ad_r;
        m.aw_t       = m.aw_r;
        m.rw_t       = m.rw_r;

        if (cl_num_spheres == UNINITIALIZED)
            m.num_spheres = 1;
    }

    if (cl_sphere_two[0] != UNINITIALIZED) {
        double d_sample_t, d_third_t, d_detector_t;

        m.d_sphere_t     = cl_sphere_two[0];
        d_sample_t       = cl_sphere_two[1];
        d_third_t        = cl_sphere_two[2];
        d_detector_t     = cl_sphere_two[3];
        m.rw_t           = cl_sphere_two[4];

        m.as_t = sqr(d_sample_t   / m.d_sphere_t / 2);
        m.at_t = sqr(d_third_t    / m.d_sphere_t / 2);
        m.ad_t = sqr(d_detector_t / m.d_sphere_t / 2);
        m.aw_t = 1.0 - m.as_t - m.at_t - m.ad_t;

        if (cl_num_spheres == UNINITIALIZED)
            m.num_spheres = 2;
    }

    if (cl_num_spheres != UNINITIALIZED) {
        m.num_spheres = (int) cl_num_spheres;
        if (m.num_spheres > 0 && m.method == UNKNOWN)
            m.method = SUBSTITUTION;
    }

    if (cl_ru_fraction != UNINITIALIZED)
        m.fraction_of_ru_in_mr = cl_ru_fraction;

    if (cl_tu_fraction != UNINITIALIZED)
        m.fraction_of_tu_in_mt = cl_tu_fraction;

    if (cl_UR1 != UNINITIALIZED)
        m.m_r = cl_UR1;

    if (cl_UT1 != UNINITIALIZED)
        m.m_t = cl_UT1;

    if (cl_Tc != UNINITIALIZED)
        m.m_u = cl_Tc;

    if (cl_default_fr != UNINITIALIZED)
        m.f_r = cl_default_fr;

    if (cl_baffle_r != UNINITIALIZED)
        m.baffle_r = cl_baffle_r;

    if (cl_baffle_t != UNINITIALIZED)
        m.baffle_t = cl_baffle_t;

    if (cl_lambda != UNINITIALIZED)
        m.lambda = cl_lambda;

@ @<Warn and quit for bad options@>=
    if (cl_method == COMPARISON && m.d_sphere_r != 0 && m.as_r == 0) {
        fprintf(stderr, "A dual-beam measurement is specified, but no port sizes.\n");
        fprintf(stderr, "You might forsake the -X option and use zero spheres (which gives\n");
        fprintf(stderr, "the same result except lost light is not taken into account).\n");
        fprintf(stderr, "Alternatively, bite the bullet and enter your sphere parameters,\n");
        fprintf(stderr, "with the knowledge that only the beam diameter and sample port\n");
        fprintf(stderr, "diameter will be used to estimate lost light from the edges.\n");
        exit(EXIT_SUCCESS);
    }

    if (cl_method == COMPARISON && m.num_spheres == 2) {
        fprintf(stderr, "A dual-beam measurement is specified, but a two sphere experiment\n");
        fprintf(stderr, "is specified. Since this seems impossible, I will make it\n");
        fprintf(stderr, "impossible for you unless you specify 0 or 1 sphere.\n");
        exit(EXIT_SUCCESS);
    }

    if (cl_method == COMPARISON && m.f_r != 0) {
        fprintf(stderr, "A dual-beam measurement is specified, but a fraction of light\n");
        fprintf(stderr, "is specified to hit the sphere wall first.  This situation\n");
        fprintf(stderr, "is not supported by iad.  Sorry.\n");
        exit(EXIT_SUCCESS);
    }

@ Put the command-line values for reflection and transmission into the measurement record.

@<Count command-line measurements@>=

    m.num_measures=3;
    if (m.m_r == 0) m.num_measures--;
    if (m.m_t == 0) m.num_measures--;
    if (m.m_u == 0) m.num_measures--;
    params = m.num_measures;

@ @<print version function@>=

static void print_version(int verbosity)
{
    if (verbosity == 0) {
        fprintf(stdout, "%s", VersionShort);
    } else {
        fprintf(stdout, "iad %s\n",Version);
        fprintf(stdout, "Copyright 1993-2024 Scott Prahl, scott.prahl@@oit.edu\n");
        fprintf(stdout, "          (see Applied Optics, 32:559-568, 1993)\n\n");
        fprintf(stdout, "This is free software; see the source for copying conditions.\n");
        fprintf(stdout, "There is no warranty; not even for MERCHANTABILITY or FITNESS.\n");
        fprintf(stdout, "FOR A PARTICULAR PURPOSE.\n");
    }
}

@ @<print usage function@>=
static void print_usage(void)
{
fprintf(stdout, "iad %s\n\n",Version);
fprintf(stdout, "iad finds optical properties from measurements\n\n");
fprintf(stdout, "Usage:  iad [options] input\n\n");
fprintf(stdout, "Options:\n");
fprintf(stdout, "  -1 '# # # # #'   reflection sphere parameters \n");
fprintf(stdout, "                   'd_sphere d d_sample_port d_entrance_port d_detector_port r_wall'\n");
fprintf(stdout, "  -2 '# # # # #'   transmission sphere parameters \n");
fprintf(stdout, "                   'd_sphere d d_sample_port d_third_port d_detector_port r_wall'\n");
fprintf(stdout, "  -a #             use this albedo \n");
fprintf(stdout, "  -A #             use this absorption coefficient \n");
fprintf(stdout, "  -b #             use this optical thickness \n");
fprintf(stdout, "  -B #             beam diameter \n");
fprintf(stdout, "  -c #             fraction of unscattered refl in MR\n");
fprintf(stdout, "  -C #             fraction of unscattered trans in MT\n");
fprintf(stdout, "  -d #             thickness of sample \n");
fprintf(stdout, "  -D #             thickness of slide \n");
fprintf(stdout, "  -e #             error tolerance (default 0.0001) \n");
fprintf(stdout, "  -E #             optical depth (=mua*D) for slides\n");
fprintf(stdout, "  -f #             allow a fraction 0.0-1.0 of light to hit sphere wall first\n");
fprintf(stdout, "  -F #             constrain scattering coefficient \n");
fprintf(stdout, "                   # = constant: use constant scattering coefficient \n");
fprintf(stdout, "                   # = 'P lambda0 mus0 gamma' then mus=mus0*(lambda/lambda0)^gamma\n");
fprintf(stdout, "  -g #             scattering anisotropy (default 0) \n");
fprintf(stdout, "  -G #             type of boundary '0', '2', 't', 'b', 'n', 'f' \n");
fprintf(stdout, "                   '0' or '2'                --- number of slides\n");
fprintf(stdout, "                   't' (top) or 'b' (bottom) \
--- one slide that is hit by light first\n");
fprintf(stdout, "                   'n' (near) or 'f' (far)   \
--- one slide position relative to sphere\n");
fprintf(stdout, "  -h               display help\n");
fprintf(stdout, "  -H #             # = 0, no baffles for R or T spheres\n");
fprintf(stdout, "                   # = 1, baffle for R but not for T sphere\n");
fprintf(stdout, "                   # = 2, baffle for T but not for R sphere\n");
fprintf(stdout, "                   # = 3, baffle for both R and T spheres (default)\n");
fprintf(stdout, "  -i #             incident angle in degrees\n");
fprintf(stdout, "  -j #             constrain reduced scattering coefficient \n");
fprintf(stdout, "  -J               generate grid after inverse calculation\n");
fprintf(stdout, "  -l #             wavelength limits\n");
fprintf(stdout, "  -L #             specify the wavelength lambda\n");
fprintf(stdout, "  -M #             limit number of Monte Carlo iterations\n");
fprintf(stdout, "  -n #             specify index of refraction of slab\n");
fprintf(stdout, "  -N #             specify index of refraction of slides\n");
fprintf(stdout, "  -o filename      explicitly specify filename for output\n");
fprintf(stdout, "  -p #             # of Monte Carlo photons (default 100000)\n");
fprintf(stdout, "                   a negative number is max MC time in milliseconds\n");
fprintf(stdout, "  -q #             number of quadrature points (default=8)\n");
fprintf(stdout, "  -r #             total reflection measurement\n");
fprintf(stdout, "  -R #             actual reflectance for 100%% measurement \n");
fprintf(stdout, "  -s #             specify type of search to do\n");
fprintf(stdout, "  -S #             number of spheres used\n");
fprintf(stdout, "  -t #             total transmission measurement\n");
fprintf(stdout, "  -T #             actual transmission for 100%% measurement \n");
fprintf(stdout, "  -u #             unscattered transmission measurement\n");
fprintf(stdout, "  -v               version information\n");
fprintf(stdout, "  -V 0             verbosity low --- no output to stdout\n");
fprintf(stdout, "  -V 1             verbosity moderate \n");
fprintf(stdout, "  -V 2             verbosity high\n");
fprintf(stdout, "  -w #             wall reflectivity for reflection sphere\n");
fprintf(stdout, "  -W #             wall reflectivity for transmission sphere\n");
fprintf(stdout, "  -x #             set debugging level\n");
fprintf(stdout, "  -X               dual beam configuration\n");
fprintf(stdout, "  -z               do forward calculation\n");
fprintf(stdout, "Examples:\n");
fprintf(stdout, "  iad file.rxt              Results will be put in file.txt\n");
fprintf(stdout, "  iad file                  Same as above\n");
fprintf(stdout, "  iad -c 0.9 file.rxt       \
Assume M_R includes 90%% of unscattered reflectance\n");
fprintf(stdout, "  iad -C 0.8 file.rxt       \
Assume M_T includes 80%% of unscattered transmittance\n");
fprintf(stdout, "  iad -e 0.0001 file.rxt    Better convergence to R & T values\n");
fprintf(stdout, "  iad -f 1.0 file.rxt       All light hits reflectance sphere wall first\n");
fprintf(stdout, "  iad -l '500 600' file.rxt Only do wavelengths between 500 and 600\n");
fprintf(stdout, "  iad -o out file.rxt       Calculated values in out\n");
fprintf(stdout, "  iad -r 0.3                R_total=0.3, b=inf, find albedo\n");
fprintf(stdout, "  iad -r 0.3 -t 0.4         R_total=0.3, T_total=0.4, find a,b,g\n");
fprintf(stdout, "  iad -r 0.3 -t 0.4 -n 1.5  R_total=0.3, T_total=0.4, n=1.5, find a,b\n");
fprintf(stdout, "  iad -r 0.3 -t 0.4         R_total=0.3, T_total=0.4, find a,b\n");
fprintf(stdout, "  iad -p 1000 file.rxt      Only 1000 photons\n");
fprintf(stdout, "  iad -p -100 file.rxt      Allow only 100ms per iteration\n");
fprintf(stdout, "  iad -q 4 file.rxt         Four quadrature points\n");
fprintf(stdout, "  iad -M 0 file.rxt         No MC    (iad)\n");
fprintf(stdout, "  iad -M 1 file.rxt         MC once  (iad -> MC -> iad)\n");
fprintf(stdout, "  iad -M 2 file.rxt         MC twice (iad -> MC -> iad -> MC -> iad)\n");
fprintf(stdout, "  iad -M 0 -q 4 file.rxt    Fast and crude conversion\n");
fprintf(stdout, "  iad -G t file.rxt         One top slide with properties from file.rxt\n");
fprintf(stdout, "  iad -G b -N 1.5 -D 1 file Use 1 bottom slide with n=1.5 and thickness=1\n");
fprintf(stdout, "  iad -x   1 file.rxt       Show sphere and MC effects\n");
fprintf(stdout, "  iad -x   2 file.rxt       Show grid decisions\n");
fprintf(stdout, "  iad -x   4 file.rxt       Show iterations\n");
fprintf(stdout, "  iad -x   8 file.rxt       Show lost light effects\n");
fprintf(stdout, "  iad -x  16 file.rxt       Show best grid points\n");
fprintf(stdout, "  iad -x  32 file.rxt       Show decisions for type of search\n");
fprintf(stdout, "  iad -x  64 file.rxt       Show all grid calculations\n");
fprintf(stdout, "  iad -x 128 file.rxt       Show sphere calculations\n");
fprintf(stdout, "  iad -x 256 file.rxt       DEBUG_EVERY_CALC\n");
fprintf(stdout, "  iad -x 511 file.rxt       Show all debugging output\n");
fprintf(stdout, "  iad -X -i 8 file.rxt      Dual beam spectrometer with 8 degree incidence\n\n");
fprintf(stdout, "  iad -z -a 0.9 -b 1 -i 45  Forward calc assuming 45 degree incidence\n\n");
fprintf(stdout, "  apply iad x.rxt y.rxt     Process multiple files\n\n");
fprintf(stdout, "Report bugs to <scott.prahl@@oit.edu>\n\n");
}

@ This can only be called immediately after |Inverse_RT|
You have been warned!  Notice that |Calculate_Distance|
does not pass any slab properties.

@<calculate coefficients function@>=
static void calculate_coefficients(struct measure_type m,
                            struct invert_type r,
                            double *LR, double *LT, double *musp, double *mua)
{
    double delta, mus;
    *LR = 0;
    *LT = 0;
    if (r.found || (!r.found && r.error == IAD_TOO_MANY_ITERATIONS)) {
        Calculate_Distance(LR, LT, &delta);
        Calculate_Mua_Musp(m, r, &mus, musp, mua);
    } else {
        *musp=0;
        *mua=0;
    }
}

@ @<print results header function@>=
static void print_results_header(FILE *fp)
{
    if (Debug(DEBUG_LOST_LIGHT)) {
        fprintf(fp,"#      | Meas      M_R  | Meas      M_T  |  calc   calc   calc  |");
        fprintf(fp,"  Lost   Lost   Lost   Lost  | MC   IAD  Error\n");

        fprintf(fp,"# wave |  M_R      fit  |  M_T      fit  |  mu_a   mu_s'   g    |  ");
        fprintf(fp," UR1    URU    UT1    UTU  |  #    #   Type\n");

        fprintf(fp,"#  nm  |  ---      ---  |  ---      ---  |  1/mm   1/mm    ---  |");
        fprintf(fp,"   ---    ---    ---    ---  | ---  ---  ---\n");

        fprintf(fp,"#---------------------------------------------------------");
        fprintf(fp,"--------------------------------------------------------\n");
    } else {
        fprintf(fp,"#     \tMeasured \t   M_R   \tMeasured \t   M_T   \tEstimated\tEstimated\tEstimated");
        fprintf(fp,"\n");

        fprintf(fp,"##wave\t   M_R   \t   fit   \t   M_T   \t   fit   \t  mu_a   \t  mu_s'  \t    g    ");
        fprintf(fp,"\n");

        fprintf(fp,"# [nm]\t  [---]  \t  [---]  \t  [---]  \t  [---]  \t  1/mm   \t  1/mm   \t  [---]  ");
        fprintf(fp,"\n");
    }
}

@ When debugging lost light, it is handy to see how each iteration changes
the calculated values for the optical properties.  We do that here if we are
debugging, otherwise we just print a number or something to keep the user from
wondering what is going on.

@s line x @q unreserve preprocessor identifier @>

@<Print results function@>=
void print_optical_property_result(FILE *fp,
                           struct measure_type m,
                           struct invert_type r,
                           double LR,
                           double LT,
                           double mu_a,
                           double mu_sp,
                           int line)
{
if (Debug(DEBUG_LOST_LIGHT)) {
    if (m.lambda != 0)
        fprintf(fp, "%6.1f   ", m.lambda);
    else
        fprintf(fp, "%6d   ", line);

    if (mu_a >= 200) mu_a = 199.9999;
    if (mu_sp >= 1000) mu_sp = 999.9999;

    fprintf(fp, "%6.4f % 6.4f | ", m.m_r, LR);
    fprintf(fp, "%6.4f % 6.4f | ", m.m_t, LT);
    fprintf(fp, "%6.3f ", mu_a);
    fprintf(fp, "%6.3f ", mu_sp);
    fprintf(fp, "%6.3f |", r.g);

    fprintf(fp, " %6.4f %6.4f ", m.ur1_lost, m.uru_lost);
    fprintf(fp, "%6.4f %6.4f | ", m.ut1_lost, m.utu_lost);
    fprintf(fp, "%2d  ", r.MC_iterations);
    fprintf(fp, "%3d", r.AD_iterations);

    fprintf(fp, "    %c \n",what_char(r.error));
} else {
    if (m.lambda != 0)
        fprintf(fp, "%6.1f\t", m.lambda);
    else
        fprintf(fp, "%6d\t", line);

    if (mu_a >= 200) mu_a = 199.9999;
    if (mu_sp >= 1000) mu_sp = 999.9999;

    fprintf(fp, "% 9.4f\t% 9.4f\t", m.m_r, LR);
    fprintf(fp, "% 9.4f\t% 9.4f\t", m.m_t, LT);
    fprintf(fp, "% 9.4f\t", mu_a);
    fprintf(fp, "% 9.4f\t", mu_sp);
    fprintf(fp, "% 9.4f\t", r.g);
    fprintf(fp, " %c \n",what_char(r.error));
}
    fflush(fp);
}

@ @<print error legend function@>=
static void print_error_legend(void)
{
    if (Debug(DEBUG_ANY)) return;
    fprintf(stderr, "----------------- Sorry, but ... errors encountered ---------------\n");
    fprintf(stderr, "   *  ==> Success          ");
    fprintf(stderr, "  0-9 ==> Monte Carlo Iteration\n");
    fprintf(stderr, "   R  ==> M_R is too big   ");
    fprintf(stderr, "   r  ==> M_R is too small\n");
    fprintf(stderr, "   T  ==> M_T is too big   ");
    fprintf(stderr, "   t  ==> M_T is too small\n");
    fprintf(stderr, "   U  ==> M_U is too big   ");
    fprintf(stderr, "   u  ==> M_U is too small\n");
    fprintf(stderr, "   !  ==> M_R + M_T > 1    ");
    fprintf(stderr, "   +  ==> Did not converge\n\n");
}


@ returns a new string consisting of s+t
@<stringdup together function@>=
static char *  strdup_together(char *s, char *t)
{
    char * both;

    if (s==NULL) {
        if (t==NULL) return NULL;
        return strdup(t);
    }

    if (t==NULL)
        return strdup(s);

    both = malloc(strlen(s) + strlen(t) + 1);
    if (both == NULL)
        fprintf(stderr, "Could not allocate memory for both strings.\n");

    strcpy(both, s);
    strcat(both, t);
    return both;
}

@ catch parsing errors in strtod
@<mystrtod function@>=

static double my_strtod(const char *str)
{
    char* endptr;
    errno = 0;

    double val = strtod(str, &endptr);

    if (endptr == str) {
        // No digits were found
        fprintf(stderr, "Error in the command line\n");
        fprintf(stderr, "    No conversion could be performed for `%s`.\n", str);
        exit(EXIT_FAILURE);
    }

    if (*endptr != '\0') {
        // String contains extra characters after the number
        fprintf(stderr, "Error in the command line\n");
        fprintf(stderr, "    Partial conversion of string = '%s'\n", str);
        exit(EXIT_FAILURE);
    }

    if (errno == ERANGE) {
        // The converted value is out of range of representable values by a double
        fprintf(stderr, "Error in the command line\n");
        printf("    The value '%s' is out of range of double.\n", str);
        exit(EXIT_FAILURE);
    }

    return val;
}

@ assume that start time has already been set
@<seconds elapsed function@>=

static double seconds_elapsed(clock_t start_time)
{
    clock_t finish_time = clock();
    return (double)(finish_time-start_time)/CLOCKS_PER_SEC;
}

@ given a string and an array, this fills the array with numbers
from the string.  The numbers should be separated by spaces.

Returns 0 upon successfully filling |n| entries, returns
1 for any error.

@<parse string into array function@>=

static int parse_string_into_array(char *s, double *a, int n)
{
    char *t, *last, *r;
    int i=0;
    t = s;
    last = s + strlen(s);

    while(t<last) {

        /* a space should mark the end of number */
        r = t;
        while (*r != ' ' && *r != '\0') r++;
        *r = '\0';

        /* parse the number and save it */
        if (sscanf(t, "%lf", &(a[i]) )==0) return 1;
        i++;

        /* are we done? */
        if (i==n) {
            if (i==5) { /* sphere parameters case */
                if (a[i-1] <= 0 || a[i-1] > 1) {
                    fprintf(stderr, "Sphere wall reflectivity (r_w=%g) must be a fraction less than one.\n", a[i-1]);
                    exit(EXIT_FAILURE);
                }
            }
            return 0;
        }

        /* move pointer just after last number */
        t=r+1;
    }

    return 1;
}

@ @<what\_char function@>=
static char what_char(int err)
{
    if (err == IAD_NO_ERROR)            return '*';
    if (err == IAD_TOO_MANY_ITERATIONS) return '+';
    if (err == IAD_MR_TOO_BIG)          return 'R';
    if (err == IAD_MR_TOO_SMALL)        return 'r';
    if (err == IAD_MT_TOO_BIG)          return 'T';
    if (err == IAD_MT_TOO_SMALL)        return 't';
    if (err == IAD_MU_TOO_BIG)          return 'U';
    if (err == IAD_MU_TOO_SMALL)        return 'u';
    if (err == IAD_TOO_MUCH_LIGHT)      return '!';
    return '?';
}

@ @<print long error function@>=
static void print_long_error(int err)
{
    if (err == IAD_TOO_MANY_ITERATIONS) fprintf(stderr, "Failed Search, too many iterations\n");
    if (err == IAD_MR_TOO_BIG)          fprintf(stderr, "Failed Search, M_R is too big\n");
    if (err == IAD_MR_TOO_SMALL)        fprintf(stderr, "Failed Search, M_R is too small\n");
    if (err == IAD_MT_TOO_BIG)          fprintf(stderr, "Failed Search, M_T is too big\n");
    if (err == IAD_MT_TOO_SMALL)        fprintf(stderr, "Failed Search, M_T is too small\n");
    if (err == IAD_MU_TOO_BIG)          fprintf(stderr, "Failed Search, M_U is too big\n");
    if (err == IAD_MU_TOO_SMALL)        fprintf(stderr, "Failed Search, M_U is too snall\n");
    if (err == IAD_TOO_MUCH_LIGHT)      fprintf(stderr, "Failed Search, Total light bigger than 1\n");
    if (err == IAD_NO_ERROR)            fprintf(stderr, "Successful Search\n");
    fprintf(stderr,"\n");
}

@ The idea here is to show some intermediate output while a file is
being processed.

@<print dot function@>=

static void print_dot(clock_t start_time, int err, int points, int final, int verbosity)
{
    static int counter = 0;
    counter++;

    if (verbosity == 0 || Debug(DEBUG_ANY)) return;

    if (final)
        fprintf(stderr, "%c", what_char(err));
    else {
        counter--;
        fprintf(stderr, "%1d\b", points % 10);
    }

    if (final) {
        if (counter % 50 == 0) {
            double rate = (seconds_elapsed(start_time) / counter);
            fprintf(stderr, "  %3d done (%5.2f s/pt)\n", counter, rate);
        } else if (counter % 10 == 0)
            fprintf(stderr," ");
    }

    fflush(stderr);
}
