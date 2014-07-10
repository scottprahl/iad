@*1 Main Program.

Here is a quick program that I put together on the 18th of July 1996 to calculate the change
in reflection and transmission when a small change in the absorption
coefficient is made.  Specifically, the absorption coefficient will
change from $\mu_a$ to $\mu_a+\mu_a\Delta$.  

The program reads and input file that contains the optical properties
of the slab.  The output file will have the same name, but appended
by ``.out'' and contain the change in the reflection and transmission 
calculated for normal irradiance using 8 quadrature points.

Note
that the streams get redirected so that I can use the standard
streams for reading, writing, and error messages.  This makes
interactive stuff problematic, but this whole thing is a batch
sort of problem.

All the output for this web file goes into \.{ad\_main.c}
but to simplify the Makefile, I create an empty \.{ad\_main.h}.

@(ad_main.h@>=

@ The program begins here  

@(ad_main.c@>=

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "ad_globl.h"
#include "ad_prime.h"
#include "ad_cone.h"
#include "mygetopt.h"
#include "version.h"

extern char *optarg;
extern int   optind;

@<print version function@>@;
@<print usage function@>@;
@<stringdup together function@>@;
@<validate slab function@>@;

int main (int argc, char **argv){

    @<Declare variables for |main|@>@;
    
    if (argc == 1) {
        print_usage();
        exit(0);
    }

    @<Handle options@>@;
    
    if (argc >= 1) {
        @<Prepare file for reading@>@;
        @<Prepare file for writing@>@;
        while (feof(stdin)==0) {
            slab.phase_function = HENYEY_GREENSTEIN;
            @<Read line from input file@>@;
            @<Calculate and Print the Results@>@;
        }
    } else {
        @<Put optical properties into |slab|@>@;
        @<Calculate and Print the Results@>@;
    }   
    return 0;
}

@ @<Declare variables for |main|@>=
    struct AD_slab_type slab;
    int nstreams = 24;
    double anisotropy = 0;
    double albedo = 0.5;
    double index_of_refraction = 1.0;
    double index_of_slide = 1.0;
    double optical_thickness = 100;
    char *g_out_name = NULL;
    double g_incident_cosine = 1;
    int machine_readable_output = 0;
    double R1, T1, URU, UTU;
    int failed;

@ I assume that the optical properties are in the following order ---
albedo, optical thickness, anisotropy, the index of refraction of the 
slab, the index of refraction of the top slide, the index of refraction
of the bottom slide.  The slides are assumed to have no absorption.

@<Put optical properties into |slab|@>=

    slab.phase_function = HENYEY_GREENSTEIN;
    slab.a = albedo;
    slab.b = optical_thickness;
    slab.g = anisotropy;
    slab.n_slab = index_of_refraction;
    slab.n_top_slide = index_of_slide;
    slab.n_bottom_slide = index_of_slide;
    slab.b_top_slide = 0.0;
    slab.b_bottom_slide = 0.0;
    slab.cos_angle = g_incident_cosine;

@
@<Read line from input file@>=
{
    int fileflag;
    fileflag = scanf("%lf", &slab.a);
    slab.cos_angle = g_incident_cosine;
    
    if (fileflag!=EOF) 
        fileflag=scanf("%lf", &slab.b);
    if (fileflag!=EOF) 
        fileflag=scanf("%lf", &slab.g);
    if (fileflag!=EOF) 
        fileflag=scanf("%lf", &slab.n_slab);
    if (fileflag!=EOF) 
        fileflag=scanf("%lf", &slab.n_top_slide);
    if (fileflag!=EOF) 
        fileflag=scanf("%lf", &slab.n_bottom_slide);
    if (fileflag!=EOF)
            fileflag=scanf("%lf", &slab.b_top_slide);
    if (fileflag!=EOF)
            fileflag=scanf("%lf", &slab.b_bottom_slide);
    if (fileflag!=EOF)
            fileflag=scanf("%d", &nstreams);
}

@ @<Calculate and Print the Results@>=
    
    failed = validate_slab(slab,nstreams,machine_readable_output);
    R1 = failed;
    T1 = failed;
    URU = failed;
    UTU = failed;
    
    if (!failed) 
        RT(nstreams,&slab,&R1,&T1,&URU,&UTU);

    if (machine_readable_output)
        printf("%9.5f \t%9.5f \t%9.5f \t%9.5f\n", R1,T1,URU,UTU);
    else if (!failed) {
        printf("UR1 = Total Reflection   for Normal  Illumination\n");
        printf("UT1 = Total Transmission for Normal  Illumination\n");
        printf("URU = Total Reflection   for Diffuse Illumination\n");
        printf("UTU = Total Transmission for Diffuse Illumination\n\n");
        printf("   UR1    \t   UT1    \t   URU    \t   UTU\n");
        printf("%9.5f \t%9.5f \t%9.5f \t%9.5f\n", R1,T1,URU,UTU);
    }

@ use the |mygetop| to process options.  We only handle help at the moment

@<Handle options@>=
{
    char c;
    double x;
    
    while ((c = my_getopt(argc, argv, "h?vma:b:g:i:n:o:q:s:")) != EOF) {
        switch (c) {

            case 'i':
                x = strtod(optarg, NULL);
                if (x<0 || x>90) 
                    fprintf(stderr, "Incident angle must be between 0 and 90 degrees\n");
                else
                    g_incident_cosine = cos(x*3.1415926535/180.0);
                break;

            case 'o':
                g_out_name = strdup(optarg);
                break;
    
            case 'n':
                index_of_refraction = strtod(optarg, NULL);
                break;

            case 's':
                index_of_slide = strtod(optarg, NULL);
                break;

            case 'm':
                machine_readable_output=1;
                break;
    
            case 'q':
                nstreams = (int) strtod(optarg, NULL);
                break;

            case 'a':
                albedo = strtod(optarg, NULL);
                break;

            case 'b':
                optical_thickness = strtod(optarg, NULL);
                break;

            case 'g':
                anisotropy = strtod(optarg, NULL);
                break;

            case 'v':
                print_version();
                break;
                
            default:
            case 'h':
            case '?':
                print_usage();
                break;
        }
    }
    
    argc -= optind;
    argv += optind;
}



@ Make sure that the file is not named '-' and warn about too many files

@<Prepare file for reading@>=
    if (argc > 1) {
        fprintf(stderr, "Only a single file can be processed at a time\n");
        fprintf(stderr, "try 'apply ad file1 file2 ... fileN'\n");
        exit(1);
    }

    if (argc == 1 && strcmp(argv[0],"-")!=0) {  /* filename exists and != "-" */

        if (freopen(argv[0],"r",stdin)==NULL) {
            fprintf(stderr, "Could not open file '%s'\n", argv[0]);
            exit(1);
        }

        if (g_out_name==NULL)
            g_out_name=strdup_together(argv[0],".rt");
    }
    
@ Take care of all the output files

@<Prepare file for writing@>=

    if (g_out_name!=NULL) {
        if (freopen(g_out_name,"w",stdout)==NULL) {
            fprintf(stderr, "Could not open file <%s> for output", g_out_name);
            exit(1);
        }
    }

@ @<print version function@>=

static void print_version(void)
{
    fprintf(stderr, "ad %s\n\n",Version);
    fprintf(stderr, "Copyright (C) 2014 Free Software Foundation, Inc.\n");
    fprintf(stderr, "This is free software; see the source for copying conditions.\n");
    fprintf(stderr, "There is no warranty; not even for MERCHANTABILITY or FITNESS.\n");
    fprintf(stderr, "FOR A PARTICULAR PURPOSE.\n");
    fprintf(stderr, "\nWritten by Scott Prahl\n");
    exit(0);
}

@ @<print usage function@>=
static void print_usage(void)
{
    fprintf(stderr, "ad %s\n\n",Version);
    fprintf(stderr, "ad finds the reflection and transmission from optical properties\n\n");
    fprintf(stderr, "Usage:  ad [options] input\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h               display help\n");
    fprintf(stderr, "  -m               machine readable output\n");
    fprintf(stderr, "  -o filename      explicitly specify filename for output\n");
    fprintf(stderr, "  -a #             albedo (0-1)\n");
    fprintf(stderr, "  -b #             optical thickness (>0)\n");
    fprintf(stderr, "  -g #             scattering anisotropy (-1 to 1)\n");
    fprintf(stderr, "  -i theta         oblique incidence at angle theta\n");
    fprintf(stderr, "  -n #             specify index of refraction of slab\n");
    fprintf(stderr, "  -s #             specify index of refraction of slide\n");
    fprintf(stderr, "  -q #             quadrature points 4, 8, 16, 32\n");
    fprintf(stderr, "  -v               version information\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "  ad data                        UR1, UT1, URU, UTU in data.rt\n");
    fprintf(stderr, "  ad -m data                     data.rt in machine readable format\n");
    fprintf(stderr, "  ad data -o out.txt             out.txt is the \n");
    fprintf(stderr, "  ad -a 0.3                      a=0.3, b=inf, g=0.0, n=1.0\n");
    fprintf(stderr, "  ad -a 0.3 -b 0.4               a=0.3, b=0.4, g=0.0, n=1.0\n");
    fprintf(stderr, "  ad -a 0.3 -b 0.4 -g 0.5        a=0.3, b=0.4, g=0.5, n=1.0\n");
    fprintf(stderr, "  ad -a 0.3 -b 0.4 -n 1.5        a=0.3, b=0.4, g=0.0, n=1.5\n\n");
    fprintf(stderr, "inputfile has lines of the form:\n");
    fprintf(stderr, "    a b g nslab ntopslide nbottomlslide btopslide bbottomslide q\n");
    fprintf(stderr, "where:\n");
    fprintf(stderr, "    1) a = albedo\n");
    fprintf(stderr, "    2) b = optical thickness\n");
    fprintf(stderr, "    3) g = anisotropy\n");
    fprintf(stderr, "    4) nslab = index of refraction of slab\n");
    fprintf(stderr, "    5) ntopslide = index of refraction of glass slide on top\n");
    fprintf(stderr, "    6) nbottomslide = index of refraction of glass slide on bottom\n");
    fprintf(stderr, "    7) btopslide = optical depth of top slide (for IR)\n");
    fprintf(stderr, "    8) bbottomslide = optical depth of bottom slide (for IR)\n");
    fprintf(stderr, "    9) q = number of quadrature points\n\n");
    fprintf(stderr, "Report bugs to <scott.prahl@@oit.edu>\n\n");
    exit(0);
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

@ Make sure that the input values are correct
@<validate slab function@>=
static int validate_slab(struct AD_slab_type slab, int nstreams, int machine)
{
    if (slab.a<0 || slab.a>1) {
        if (!machine)
            printf("Bad Albedo a=%f\n",slab.a);
        return(-1);
    }

    if (slab.b<0 ) {
        if (!machine)
            printf("Bad Optical Thickness b=%f\n",slab.b);
        return(-2);
    }

    if (slab.g <= -1 || slab.g >= 1) {
        if (!machine)
            printf("Bad Anisotropy g=%f\n",slab.g);
        return(-3);
    }

    if (slab.n_slab<0 || slab.n_slab>10) {
        if (!machine)
            printf("Bad Slab Index n=%f\n",slab.n_slab);
        return(-4);
    }

    if (slab.n_top_slide<1 || slab.n_top_slide>10) {
        if (!machine)
            printf("Bad Top Slide Index n=%f\n",slab.n_top_slide);
        return(-5);
    }

    if (slab.n_bottom_slide<1 || slab.n_bottom_slide>10) {
        if (!machine)
            printf("Bad Top Slide Index n=%f\n",slab.n_bottom_slide);
        return(-6);
    }
    
    if (slab.b_top_slide<0 || slab.b_top_slide>10) {
        if (!machine)
            printf("Bad Top Slide Optical Thickness b=%f\n",slab.b_top_slide);
        return(-7);
    }

    if (slab.b_bottom_slide<0 || slab.b_bottom_slide>10) {
        if (!machine)
            printf("Bad Bottom Slide Optical Thickness b=%f\n",slab.b_bottom_slide);
        return(-8);
    }
    
    if (nstreams<4 || nstreams % 4 != 0){
        if (!machine) {
            printf("Bad Number of Quadrature Points npts=%d\n",nstreams);
            printf("Should be a multiple of four!\n");
        }
        return(-9);
    }
    
    return 0;
}
