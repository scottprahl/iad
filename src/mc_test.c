#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "ad_globl.h"
#include "ad_prime.h"
#include "ad_cone.h"
#include "iad_type.h"
#include "mc_lost.h"

int main(int argc, char **argv)
{
    struct AD_slab_type s;
    double ur1, ut1, uru, utu;
    double ad_ur1, ad_ut1, ad_uru, ad_utu;
    double d = 1;
    double d_slide = 0;
    long n_photons = 1000000;
    int n_quad = 32;
    s.a = 0.0;
    s.b = 0.5;
    s.g = 0.0;
    s.phase_function = HENYEY_GREENSTEIN;
    s.n_slab = 1.0;
    s.n_top_slide = 1.0;
    s.n_bottom_slide = 1.0;
    s.b_top_slide = 0;
    s.b_bottom_slide = 0;
    s.cos_angle = 1;

    printf("   a     b      g      n     n    |   d    dslide |     UR1             UT1             URU             UTU\n");
    printf("                     slab  slide  |   mm     mm   |  AD     MC       AD     MC       AD     MC       AD     MC\n");


    s.a = 0.0;
    s.b = 0.0;
    s.g = 0.0;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.a = 0.0;
    s.b = 0.0;
    s.g = 0.0;
    s.n_slab = 2.5;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.a = 0.0;
    s.b = 0.0;
    s.g = 0.0;
    s.n_slab = 1.5;
    s.n_top_slide = 1.5;
    s.n_bottom_slide = 1.5;
    d_slide = 0.0;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.a = 0.0;
    s.b = 0.0;
    s.g = 0.0;
    s.n_slab = 1.4;
    s.n_top_slide = 1.5;
    s.n_bottom_slide = 1.5;
    d_slide = 0.0;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.a = 0.0;
    s.b = 0.0;
    s.g = 0.0;
    s.n_slab = 1.4;
    s.n_top_slide = 1.5;
    s.n_bottom_slide = 1.5;
    d_slide = 0.15;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);
    printf("\n");

    s.a = 0.0;
    s.b = 0.5;
    s.g = 0.0;
    s.n_slab = 1.0;
    s.n_top_slide = 1.0;
    s.n_bottom_slide = 1.0;
    d_slide = 0;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.a = 0.0;
    s.b = 0.5;
    s.g = 0.0;
    s.n_slab = 1.5;
    s.n_top_slide = 1.0;
    s.n_bottom_slide = 1.0;
    d_slide = 0;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.a = 0.0;
    s.b = 0.5;
    s.g = 0.0;
    s.n_slab = 1.5;
    s.n_top_slide = 1.5;
    s.n_bottom_slide = 1.5;
    d_slide = 0;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.a = 0.0;
    s.b = 0.5;
    s.g = 0.0;
    s.n_slab = 1.4;
    s.n_top_slide = 1.5;
    s.n_bottom_slide = 1.5;
    d_slide = 0.15;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);
    printf("\n");

    s.a = 0.5;
    s.b = 0.5;
    s.g = 0.0;
    s.n_slab = 1.0;
    s.n_top_slide = 1.0;
    s.n_bottom_slide = 1.0;
    d_slide = 0;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.a = 0.5;
    s.b = 0.5;
    s.g = 0.0;
    s.n_slab = 1.5;
    s.n_top_slide = 1.0;
    s.n_bottom_slide = 1.0;
    d_slide = 0;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.a = 0.5;
    s.b = 0.5;
    s.g = 0.0;
    s.n_slab = 1.5;
    s.n_top_slide = 1.5;
    s.n_bottom_slide = 1.5;
    d_slide = 0;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.a = 0.5;
    s.b = 0.5;
    s.g = 0.0;
    s.n_slab = 1.4;
    s.n_top_slide = 1.5;
    s.n_bottom_slide = 1.5;
    d_slide = 0.15;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);
    printf("\n");

    s.b = 100.0;
    s.n_slab = 1.5;
    s.n_top_slide = 1.0;
    s.n_bottom_slide = 1.0;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.n_slab = 1.5;
    s.n_top_slide = 1.5;
    s.n_bottom_slide = 1.5;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.n_slab = 1.3;
    s.n_top_slide = 1.5;
    s.n_bottom_slide = 1.5;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);
    printf("\n");

    s.a = 0.5;
    s.b = 1.0;
    s.n_slab = 1.0;
    s.n_top_slide = 1.0;
    s.n_bottom_slide = 1.0;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.g = 0.5;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.n_slab = 1.5;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);
    printf("\n");

    s.g = 0.9;
    s.b = 4;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    printf("\n");

    d_slide=1;
    s.a = 0.9;
    s.b = 4.0;
    s.n_slab = 1.0;
    s.n_top_slide = 1.0;
    s.n_bottom_slide = 1.0;
    s.g = 0;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.g = 0.9;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.n_slab = 1.5;
    s.n_top_slide = 1.0;
    s.n_bottom_slide = 1.0;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.n_slab = 1.5;
    s.n_top_slide = 1.5;
    s.n_bottom_slide = 1.5;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.n_slab = 1.3;
    s.n_top_slide = 1.5;
    s.n_bottom_slide = 1.5;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    printf("\n");

    s.a = 0.999;
    s.b = 4.0;
    s.n_slab = 1.0;
    s.n_top_slide = 1.0;
    s.n_bottom_slide = 1.0;
    s.g = 0;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.g = 0.9;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.n_slab = 1.5;
    s.n_top_slide = 1.0;
    s.n_bottom_slide = 1.0;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.n_slab = 1.5;
    s.n_top_slide = 1.5;
    s.n_bottom_slide = 1.5;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    s.n_slab = 1.3;
    s.n_top_slide = 1.5;
    s.n_bottom_slide = 1.5;
    MC_RT(s, n_photons, d, d_slide, &ur1,  &ut1,  &uru,  &utu);
    RT(n_quad, &s, &ad_ur1, &ad_ut1, &ad_uru, &ad_utu);
    printf("%5.4f %5.1f %5.4f %5.4f %5.4f | ", s.a, s.b, s.g, s.n_slab, s.n_top_slide);
    printf("%5.4f %5.4f | ", d, d_slide);
    printf("%5.4f %5.4f   %5.4f %5.4f   ", ad_ur1, ur1, ad_ut1, ut1);
    printf("%5.4f %5.4f   %5.4f %5.4f\n", ad_uru, uru, ad_utu, utu);

    printf("\n");
}
