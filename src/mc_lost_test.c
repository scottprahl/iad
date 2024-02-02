#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ad_globl.h"
#include "ad_prime.h"
#include "iad_type.h"
#include "mc_lost.h"

static void print_usage(void)
{
    fprintf(stderr, "lost1 \n\n");
    fprintf(stderr,
        "lost finds the reflection and transmission from optical properties\n\n");
    fprintf(stderr, "Usage:  lost [options] input\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -a #   albedo (0-1)\n");
    fprintf(stderr, "  -b #   optical thickness (>0)\n");
    fprintf(stderr, "  -g #   scattering anisotropy (-1 to 1)\n");
    fprintf(stderr, "  -h     display help\n");
    fprintf(stderr, "  -i #   light is incident at this angle in degrees\n");
    fprintf(stderr, "  -n #   specify index of refraction of slab\n");
    fprintf(stderr, "  -N #   specify index of refraction of slide\n");
    fprintf(stderr, "  -p #   number of photons\n");
    fprintf(stderr, "  -v     display help\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr,
        "  mc_lost -a 0.3 -b 2.0 -g 0.8 -n 1.4 -N 1.5 -p 1000000\n");
    exit(0);
}

int main(int argc, char **argv)
{
    char c;
    double UR1, UT1, URU, UTU, x;
    double ur1_lost, ut1_lost, uru_lost, utu_lost;
    double ur1, ut1, uru, utu;
    int collimated = 1;
    int diffuse = 0;
    double t_sample = 1.0;
    double d_port = 10.0;
    double d_beam = 5.0;
    double t_slide = 1.0;
    double a = 0.9;
    double b = 1;
    double g = 0.8;
    double n_slab = 1.4;
    double n_slide = 1.5;
    double mu0 = 1;
    double mua_slide = 0;
    long n_photons = 1024 * 1024;

    while ((c = getopt(argc, argv, "hva:b:g:i:n:N:p:")) != -1) {
        switch (c) {

        case 'a':
            a = strtod(optarg, NULL);
            break;

        case 'b':
            b = strtod(optarg, NULL);
            break;

        case 'g':
            g = strtod(optarg, NULL);
            break;

        case 'i':
            x = strtod(optarg, NULL);
            if (x < 0 || x > 90) {
                fprintf(stderr,
                    "Incident angle must be between 0 and 90 degrees\n");
                return 1;
            }
            else
                mu0 = cos(x * M_PI / 180.0);
            break;

        case 'n':
            n_slab = strtod(optarg, NULL);
            break;

        case 'N':
            n_slide = strtod(optarg, NULL);
            break;

        case 'p':
            n_photons = (long) strtod(optarg, NULL);
            break;

        default:
        case 'h':
        case 'v':
            print_usage();
            break;
        }
    }

    MC_Radial(n_photons, a, b, g, n_slab, n_slide, collimated, mu0, t_sample,
        t_slide, mua_slide, d_port, d_beam, &ur1, &ut1, &ur1_lost, &ut1_lost);
    MC_Radial(n_photons, a, b, g, n_slab, n_slide, diffuse, mu0, t_sample,
        t_slide, mua_slide, d_port, d_beam, &uru, &utu, &uru_lost, &utu_lost);

    printf("   UR1    \t   UT1    \t   URU    \t   UTU\n");
    printf("%9.5f \t%9.5f \t%9.5f \t%9.5f \tMC Calc\n", ur1, ut1, uru, utu);

    ez_RT(16, n_slab, n_slide, n_slide, a, b, g, &UR1, &UT1, &URU, &UTU);
    printf("%9.5f \t%9.5f \t%9.5f \t%9.5f \tAD Calc\n", UR1, UT1, URU, UTU);
    printf
        ("-----------------------------------------------------------------------\n");
    printf("%9.5f \t%9.5f \t%9.5f \t%9.5f \tMC Loss\n", ur1_lost, ut1_lost,
        uru_lost, utu_lost);

    return 0;
}
