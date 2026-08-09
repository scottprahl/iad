/* Regression tests for the library entry points.

   The command-line program refuses to enter the Monte Carlo loop unless a
   sphere was described, but the library lets a caller ask for lost-light
   runs directly.  Nothing exercised that path, and for a long time it
   answered by deriving the sample ports from sphere parameters that
   measure_OK had never validated.  With no sphere those come out zero wide,
   a port of zero diameter catches nothing, every photon is scored as lost,
   and the caller is handed coefficients built on that fiction.

   These tests link against the objects in this directory rather than an
   installed libiad, so they always describe the tree they were built from.

   Run with no arguments; a nonzero exit status means something regressed. */

#include <stdio.h>
#include <math.h>

#include "ad_globl.h"
#include "iad_type.h"
#include "iad_pub.h"
#include "mc_lost.h"

/* Spheres_Inverse_RT is declared in lib_iad.h, which ctangle produces as a
   side effect of building iad_pub.c.  Depending on that header here would
   make this test re-run ctangle on iad_pub.w, overwriting the tidied
   iad_pub.c and iad_pub.h with raw output every time the suite ran.  The
   prototype is repeated instead; it is the array-based interface meant for
   callers who cannot see the C structures, so it does not change. */

void Spheres_Inverse_RT(double *setup,
    double *analysis, double *sphere_r, double *sphere_t, double *measurements, double *results);

#define TEST_PHOTONS 50000

static int failures = 0;

static void complain(const char *what)
{
    printf("FAIL: %s\n", what);
    failures++;
}

/* Port areas are stored as a fraction of their own sphere, so a diameter has
   to be turned into one before it can be handed to measure_type. */

static double port_fraction(double d_port, double d_sphere)
{
    double ratio = d_port / d_sphere / 2.0;
    return ratio * ratio;
}

/* A thick, strongly scattering slab spreads light far enough sideways that a
   port of ordinary size loses a measurable amount of it. */

static void build_sample(struct measure_type *m, struct invert_type *r)
{
    Initialize_Measure(m);

    m->slab_index = 1.4;
    m->slab_thickness = 5.0;
    m->slab_top_slide_index = 1.0;
    m->slab_top_slide_thickness = 0.0;
    m->slab_bottom_slide_index = 1.0;
    m->slab_bottom_slide_thickness = 0.0;
    m->slab_cos_angle = 1.0;
    m->d_beam = 2.0;
    m->method = SUBSTITUTION;

    m->d_sphere_r = 200.0;
    m->d_sphere_t = 200.0;
    m->as_r = port_fraction(8.0, m->d_sphere_r);
    m->as_t = port_fraction(8.0, m->d_sphere_t);

    Initialize_Result(*m, r, TRUE);
    r->a = 0.995;
    r->b = 30.0;
    r->g = 0.0;
}

/* Without a sphere there is no port, so there is nothing for a photon to
   miss.  Every output must come back untouched at zero -- not merely small,
   exactly zero -- because no simulation was run to produce anything else. */

static void no_sphere_means_no_monte_carlo(void)
{
    struct measure_type m;
    struct invert_type r;
    double ur1 = -1, ut1 = -1, uru = -1, utu = -1;
    double ur1_lost = -1, ut1_lost = -1, uru_lost = -1, utu_lost = -1;

    build_sample(&m, &r);
    m.num_spheres = 0;

    MC_Lost(m, r, TEST_PHOTONS, &ur1, &ut1, &uru, &utu, &ur1_lost, &ut1_lost, &uru_lost, &utu_lost);

    if (ur1_lost != 0.0 || ut1_lost != 0.0 || uru_lost != 0.0 || utu_lost != 0.0)
        complain("no spheres: MC_Lost reported lost light with no port to lose it through");

    if (ur1 != 0.0 || ut1 != 0.0 || uru != 0.0 || utu != 0.0)
        complain("no spheres: MC_Lost left the reflectances holding stale values");

    printf("  no spheres          : losses %.4f %.4f %.4f %.4f\n", ur1_lost, ut1_lost, uru_lost, utu_lost);
}

/* The guard must not fire when a sphere really was described, or the lost
   light correction would quietly vanish for every measurement. */

static void one_sphere_still_runs(void)
{
    struct measure_type m;
    struct invert_type r;
    double ur1, ut1, uru, utu;
    double ur1_lost, ut1_lost, uru_lost, utu_lost;

    build_sample(&m, &r);
    m.num_spheres = 1;

    MC_Lost(m, r, TEST_PHOTONS, &ur1, &ut1, &uru, &utu, &ur1_lost, &ut1_lost, &uru_lost, &utu_lost);

    if (ur1_lost <= 0.0 || ut1_lost <= 0.0)
        complain("one sphere: collimated lost light came back zero");

    if (uru_lost <= 0.0 || utu_lost <= 0.0)
        complain("one sphere: diffuse lost light came back zero");

    printf("  one sphere          : UR1=%.4f UT1=%.4f URU=%.4f UTU=%.4f\n", ur1_lost, ut1_lost, uru_lost, utu_lost);
}

/* One sphere means one sphere at a time, so the reflectance may be measured
   on a different sphere from the transmittance.  Each measurement is scored
   against its own sample port: as_r with d_sphere_r for reflection, as_t
   with d_sphere_t for transmission.  Narrowing one port must move that
   measurement and leave the other alone.  Comparing as_r to as_t directly
   would be a bug, since the fractions belong to different spheres. */

static void each_side_uses_its_own_port(void)
{
    struct measure_type m;
    struct invert_type r;
    double wide_ur1, wide_ut1, wide_uru, wide_utu;
    double wide_ur1_lost, wide_ut1_lost, wide_uru_lost, wide_utu_lost;
    double ur1, ut1, uru, utu;
    double ur1_lost, ut1_lost, uru_lost, utu_lost;

    build_sample(&m, &r);
    m.num_spheres = 1;
    MC_Lost(m, r, TEST_PHOTONS, &wide_ur1, &wide_ut1, &wide_uru, &wide_utu,
        &wide_ur1_lost, &wide_ut1_lost, &wide_uru_lost, &wide_utu_lost);

    /* shrink only the transmission sample port */
    build_sample(&m, &r);
    m.num_spheres = 1;
    m.as_t = port_fraction(4.0, m.d_sphere_t);
    MC_Lost(m, r, TEST_PHOTONS, &ur1, &ut1, &uru, &utu, &ur1_lost, &ut1_lost, &uru_lost, &utu_lost);

    if (ut1_lost <= wide_ut1_lost)
        complain("halving the transmission sample port did not increase UT1 loss");
    if (ur1_lost != wide_ur1_lost)
        complain("the transmission sample port changed the reflection loss");

    printf("  T port  8 -> 4 mm   : UT1 %.4f -> %.4f, UR1 held at %.4f\n", wide_ut1_lost, ut1_lost, ur1_lost);

    /* shrink only the reflection sample port */
    build_sample(&m, &r);
    m.num_spheres = 1;
    m.as_r = port_fraction(4.0, m.d_sphere_r);
    MC_Lost(m, r, TEST_PHOTONS, &ur1, &ut1, &uru, &utu, &ur1_lost, &ut1_lost, &uru_lost, &utu_lost);

    if (ur1_lost <= wide_ur1_lost)
        complain("halving the reflection sample port did not increase UR1 loss");
    if (ut1_lost != wide_ut1_lost)
        complain("the reflection sample port changed the transmission loss");

    printf("  R port  8 -> 4 mm   : UR1 %.4f -> %.4f, UT1 held at %.4f\n", wide_ur1_lost, ur1_lost, ut1_lost);
}

/* The two spheres may differ in size while describing the same physical
   port.  Because the stored fraction is relative to its own sphere, the same
   8 mm port is a different number in each block, and a routine that
   compared the fractions instead of the diameters would see a mismatch that
   is not there.  Halving the sphere and keeping the port must not move the
   answer. */

static void port_size_follows_the_diameter(void)
{
    struct measure_type m;
    struct invert_type r;
    double big_ur1, big_ut1, big_uru, big_utu;
    double big_ur1_lost, big_ut1_lost, big_uru_lost, big_utu_lost;
    double ur1, ut1, uru, utu;
    double ur1_lost, ut1_lost, uru_lost, utu_lost;

    build_sample(&m, &r);
    m.num_spheres = 1;
    MC_Lost(m, r, TEST_PHOTONS, &big_ur1, &big_ut1, &big_uru, &big_utu,
        &big_ur1_lost, &big_ut1_lost, &big_uru_lost, &big_utu_lost);

    build_sample(&m, &r);
    m.num_spheres = 1;
    m.d_sphere_t = 100.0;
    m.as_t = port_fraction(8.0, m.d_sphere_t);
    MC_Lost(m, r, TEST_PHOTONS, &ur1, &ut1, &uru, &utu, &ur1_lost, &ut1_lost, &uru_lost, &utu_lost);

    if (ut1_lost != big_ut1_lost)
        complain("the same port in a smaller sphere changed the transmission loss");

    printf("  T sphere 200->100 mm: same  8 mm port, UT1 held at %.4f\n", ut1_lost);
}

/* Spheres_Inverse_RT takes the number of Monte Carlo runs straight from the
   caller.  With no sphere described that request has to be dropped, so the
   recovered coefficients must not depend on how many runs were asked for.

   The library keeps the adaptive grid in module globals, so repeating an
   identical inversion in one process drifts by around a part in a million.
   The failure this guards against was a factor of five, so the comparison
   is made with a tolerance well below the one and well above the other. */

static void library_ignores_runs_without_a_sphere(void)
{
    double setup[19], analysis[2], sphere_r[8], sphere_t[8], measurements[3];
    double none[4], five[4];
    int i;

    for (i = 0; i < 19; i++)
        setup[i] = 0.0;
    for (i = 0; i < 8; i++) {
        sphere_r[i] = 0.0;
        sphere_t[i] = 0.0;
    }

    setup[0] = 1.4;             /* sample index                       */
    setup[1] = 1.0;             /* top slide index                    */
    setup[2] = 1.0;             /* sample thickness                   */
    setup[4] = 5.0;             /* beam diameter                      */
    setup[5] = 1.0;             /* reflectance standard               */
    setup[6] = 0.0;             /* no spheres                         */
    setup[7] = 200.0;           /* the sphere diameters divide the    */
    setup[12] = 200.0;          /* port fractions and cannot be zero  */
    setup[11] = 0.97;
    setup[16] = 0.97;
    setup[18] = 100000;         /* photons                            */

    analysis[0] = 8;            /* quadrature points                  */
    measurements[0] = 0.4;
    measurements[1] = 0.18;
    measurements[2] = 0.0;

    analysis[1] = 0;
    Spheres_Inverse_RT(setup, analysis, sphere_r, sphere_t, measurements, none);

    analysis[1] = 5;
    Spheres_Inverse_RT(setup, analysis, sphere_r, sphere_t, measurements, five);

    if (none[3] != 0 || five[3] != 0)
        complain("no spheres: the inversion reported an error");

    if (fabs(none[0] - five[0]) > 1e-4 || fabs(none[1] - five[1]) > 1e-4)
        complain("no spheres: asking for Monte Carlo runs changed the coefficients");

    printf("  0 vs 5 MC runs      : mu_a %.6f/%.6f  mu_s' %.6f/%.6f\n", none[0], five[0], none[1], five[1]);
}

int main(void)
{
    printf("== library entry points and lost-light ports\n");

    no_sphere_means_no_monte_carlo();
    one_sphere_still_runs();
    each_side_uses_its_own_port();
    port_size_follows_the_diameter();
    library_ignores_runs_without_a_sphere();

    if (failures) {
        printf("%d check(s) failed\n", failures);
        return 1;
    }

    printf("  all checks passed\n");
    return 0;
}
