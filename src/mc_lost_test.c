#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ad_globl.h"
#include "ad_prime.h"
#include "ad_cone.h"
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
    fprintf(stderr, "  -B #   beam diameter\n");
    fprintf(stderr, "  -g #   scattering anisotropy (-1 to 1)\n");
    fprintf(stderr, "  -h     display help\n");
    fprintf(stderr, "  -i #   light is incident at this angle in degrees\n");
    fprintf(stderr, "  -n #   specify index of refraction of slab\n");
    fprintf(stderr, "  -N #   specify index of refraction of slide\n");
    fprintf(stderr, "  -p #   number of photons\n");
    fprintf(stderr, "  -P #   sample port diameter\n");
    fprintf(stderr, "  -v     display help\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr,
        "  mc_lost -a 0.3 -b 2.0 -g 0.8 -n 1.4 -N 1.5 -p 1000000\n");
    exit(EXIT_SUCCESS);
}

int main(int argc, char **argv)
{
    char c;
    double URx, UTx, URU, UTU, x;
    double mc_ur1_lost, mc_ut1_lost, mc_uru_lost, mc_utu_lost;
    double mc_ur1, mc_ut1, mc_uru, mc_utu;
    int collimated = 1;
    int diffuse = 0;
    double t_sample = 1.0;
    double dr_port = 10.0;
    double dt_port = 10.0;
    double d_beam = 5.0;
    double t_slide = 1.0;
    double a = 0.99;
    double b = 1;
    double g = 0.0;
    double n_slab = 1.0;
    double n_slide = 1.0;
    double mu0 = 1;
    double mua_slide = 0;
    long n_photons = 1024 * 1024;
    int machine_output = 0;

    while ((c = getopt(argc, argv, "hva:b:B:g:i:n:mN:p:P:")) != -1) {
        switch (c) {

        case 'a':
            a = strtod(optarg, NULL);
            break;

        case 'b':
            b = strtod(optarg, NULL);
            break;

        case 'B':
            d_beam = strtod(optarg, NULL);
            break;

        case 'g':
            g = strtod(optarg, NULL);
            break;

        case 'i':
            x = strtod(optarg, NULL);
            if (x < 0 || x > 90) {
                fprintf(stderr, "Incident angle must be between 0 and 90 degrees\n");
                exit(EXIT_FAILURE);
            }
            else
                mu0 = cos(x * M_PI / 180.0);
            break;

        case 'n':
            n_slab = strtod(optarg, NULL);
            break;

        case 'm':
            machine_output = 1;
            break;

        case 'N':
            n_slide = strtod(optarg, NULL);
            break;

        case 'p':
            n_photons = (long) strtod(optarg, NULL);
            break;

        case 'P':
            dr_port = strtod(optarg, NULL);
            dt_port = dr_port;
            break;

        default:
        case 'h':
        case 'v':
            print_usage();
            break;
        }
    }

    if (machine_output == 0) {
        printf("Albedo                  %10.5f\n",a);
        printf("Optical Depth           %10.5f\n",b);
        printf("Anisotropy              %10.5f\n",g);
        printf("Indices of Refraction\n");
        printf("                  slab  %10.5f\n",n_slab);
        printf("             top slide  %10.5f\n",n_slide);
        printf("          bottom slide  %10.5f\n",n_slide);
        printf("Port and Beam Diameter\n");
        printf("       reflection port  %10.5f mm\n",dr_port);
        printf("     transmission port  %10.5f mm\n",dt_port);
        printf("                  beam  %10.5f mm\n",d_beam);
        printf("Incidence angle         %10.5f\n",acos(mu0)*180.0/M_PI);
        printf("Cos of incidence angle  %10.5f\n",mu0);
        printf("\n");
        printf("  URx    \t   UTx    \t   URU    \t   UTU\n");
    }

    MC_Radial(n_photons, a, b, g, n_slab, n_slide, collimated, mu0, t_sample,
        t_slide, mua_slide, dr_port, dt_port, d_beam, &mc_ur1, &mc_ut1, &mc_ur1_lost, &mc_ut1_lost);

    MC_Radial(n_photons, a, b, g, n_slab, n_slide, diffuse, mu0, t_sample,
        t_slide, mua_slide, dr_port, dt_port, d_beam, &mc_uru, &mc_utu, &mc_uru_lost, &mc_utu_lost);

    ez_RT_Oblique(12, n_slab, n_slide, n_slide, a, b, g, mu0, &URx, &UTx, &URU, &UTU);

    if (machine_output == 0) {
        printf("%9.5f \t%9.5f \t%9.5f \t%9.5f \tMC Calc\n", mc_ur1, mc_ut1, mc_uru, mc_utu);
        printf("%9.5f \t%9.5f \t%9.5f \t%9.5f \tAD Calc\n", URx, UTx, URU, UTU);
        printf("%9.5f \t%9.5f \t%9.5f \t%9.5f \tMC Loss\n", mc_ur1_lost, mc_ut1_lost,
            mc_uru_lost, mc_utu_lost);
//        ez_RT(12, n_slab, n_slide, n_slide, a, b, g, &URx, &UTx, &URU, &UTU);
//        printf("%9.5f \t%9.5f \t%9.5f \t%9.5f \tAD Calc UR1 and UT1\n", URx, UTx, URU, UTU);
        printf
            ("-----------------------------------------------------------------------\n");
        printf("\n");
    } else {
        printf("%9.5f %9.5f %9.5f %9.5f ", mc_ur1, mc_ut1, mc_uru, mc_utu);
        printf("%9.5f %9.5f %9.5f %9.5f ", URx, UTx, URU, UTU);
        printf("%9.5f %9.5f %9.5f %9.5f" , mc_ur1_lost, mc_ut1_lost, mc_uru_lost, mc_utu_lost);
    }

    exit(EXIT_SUCCESS);
}
