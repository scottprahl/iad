#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "ad_globl.h"
#include "ad_prime.h"
#include "iad_type.h"
#include "mc_lost.h"

static void print_usage(void)
{
    fprintf(stderr, "lost1 \n\n");
    fprintf(stderr, "lost finds the reflection and transmission from optical properties\n\n");
    fprintf(stderr, "Usage:  lost [options] input\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h     display help\n");
    fprintf(stderr, "  -a #   albedo (0-1)\n");
    fprintf(stderr, "  -b #   optical thickness (>0)\n");
    fprintf(stderr, "  -g #   scattering anisotropy (-1 to 1)\n");
    fprintf(stderr, "  -n #   specify index of refraction of slab\n");
    fprintf(stderr, "  -s #   specify index of refraction of slide\n");
    fprintf(stderr, "  -N #   number of photons\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "  lost -m data                     data.rt in machine readable format\n");
    fprintf(stderr, "  ad data -o out.txt             out.txt is the \n");
    fprintf(stderr, "  lost -a 0.3                    a=0.3, b=2.0, g=0.0, n=1.0\n");
    exit(0);
}

int main(int argc, char** argv)
{
    char c;
    double UR1, UT1, URU, UTU;
    double ur1_lost, ut1_lost, uru_lost, utu_lost;
    double ur1, ut1, uru, utu;
    int collimated = 1;
    int diffuse = 0;
    double t_sample = 1.0;
    double d_port   = 10.0;
    double d_beam   = 5.0;
    double t_slide  = 1.0;
    double a = 0.5;
    double b = 2;
    double g = 0.0;
    double n_slab = 1.0;
    double n_slide = 1.0;
    double mu0 = 1;
    double mua_slide = 0;
    int n_photons = 1024*1024;

    while ((c = getopt(argc, argv, "h?va:b:g:n:N:s:q:")) != -1) {
        switch (c) {

            case 'n':
                n_slab = strtod(optarg, NULL);
                break;

            case 'N':
                n_photons = (int) strtod(optarg, NULL);
                break;

            case 'a':
                a = strtod(optarg, NULL);
                break;

            case 'b':
                b = strtod(optarg, NULL);
                break;

            case 'g':
                g = strtod(optarg, NULL);
                break;

            case 's':
                n_slide = strtod(optarg, NULL);
                break;

            default:
            case 'h':
                print_usage();
                break;
        }
    }

    MC_Radial(n_photons, a, b, g, n_slab, n_slide, collimated, mu0, t_sample, t_slide, mua_slide, d_port, d_beam, &ur1, &ut1, &ur1_lost, &ut1_lost);
    MC_Radial(n_photons, a, b, g, n_slab, n_slide, diffuse,    mu0, t_sample, t_slide, mua_slide, d_port, d_beam, &uru, &utu, &uru_lost, &utu_lost);

    printf("   UR1    \t   UT1    \t   URU    \t   UTU\n");
    printf("%9.5f \t%9.5f \t%9.5f \t%9.5f\n", ur1, ut1, uru, utu);
    printf("%9.5f \t%9.5f \t%9.5f \t%9.5f\n", ur1_lost, ut1_lost, uru_lost, utu_lost);

    ez_RT(16, n_slab, n_slide, n_slide, a, b, g, &UR1, &UT1, &URU, &UTU);
    printf("%9.5f \t%9.5f \t%9.5f \t%9.5f\n", UR1, UT1, URU, UTU);

    return 0;
}
