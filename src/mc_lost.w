@** Monte Carlo lost light.

This module estimates the fraction of reflected and transmitted light that
misses the integrating-sphere ports because photons leave the finite sample
away from the illuminated area.  The adding-doubling calculation is one
dimensional; this Monte Carlo correction supplies the radial information
needed by |iad| when sample geometry matters.

The sample is represented as an air-glass-sample-glass-air stack.  Photons
are launched from the top of the first slide, transported through Fresnel
boundaries and the scattering slab, and finally scored by radial distance at
the reflection or transmission port.

@(mc_lost.c@>=
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "ad_globl.h"
#include "ad_frsnl.h"
#include "iad_type.h"
#include "iad_util.h"

@<Monte Carlo constants@>@;
@<Monte Carlo global state@>@;

@<Definition for |kiss_rand|@>@;
@<Definition for |kiss_rand_seed|@>@;
@<Definition for |set_photon_seed|@>@;
@<Definition for |MC_Set_Seed|@>@;
@<Definition for |next_photon_seed|@>@;
@<Definition for |rand_zero_one|@>@;
@<Definition for |rand_one_one|@>@;
@<Definition for |fresnel|@>@;
@<Definition for |refract|@>@;
@<Definition for |scatter|@>@;
@<Definition for |launch_point|@>@;
@<Definition for |launch_direction|@>@;
@<Definition for |roulette|@>@;
@<Definition for |move_in_sample|@>@;
@<Definition for |move_at_boundary|@>@;
@<Definition for |move_in_slide|@>@;
@<Definition for |milliseconds|@>@;
@<Definition for |add_to_reflectance_array|@>@;
@<Definition for |add_to_transmittance_array|@>@;
@<Definition for |MC_Radial|@>@;
@<Definition for |MC_Lost|@>@;
@<Definition for |MC_RT|@>@;
@<Definition for |MC_Print_RT_Arrays|@>@;
@<Definition for |MC_Lost_With_Stderr|@>@;

@ The header exposes the Monte Carlo services used by |iad| and by the
standalone |mc_lost| diagnostic program.

@(mc_lost.h@>=
@<Prototype for |MC_Set_Seed|@>;
@<Prototype for |MC_Lost|@>;
@<Prototype for |MC_RT|@>;
@<Prototype for |MC_Radial|@>;
@<Prototype for |MC_Print_RT_Arrays|@>;
@<Prototype for |MC_Lost_With_Stderr|@>;

@*1 Constants and state.

Photons with very small weight are handled by Russian roulette rather than
being followed forever.  The radial arrays are intentionally global so the
standalone executable can request a radial dump after a simulation.

@ @<Monte Carlo constants@>=
#define MIN_WEIGHT 0.0001
#define N_RADIAL_BINS  1001
#define RADIAL_BIN_SIZE 0.02
#define CLOSE(x, y) (fabs((x) - (y)) < 1e-8)

@ The KISS generator state is process-global.  |photon_seed| is advanced once
per photon so repeated simulations can compare nearby optical properties with
the same sequence of photon histories.

@ @<Monte Carlo global state@>=
unsigned long photon_seed = 12345678;

int print_radial_arrays = FALSE;
double R_radial[N_RADIAL_BINS] = { 0 };
double T_radial[N_RADIAL_BINS] = { 0 };

unsigned long kiss_rand_max = ULONG_MAX;
unsigned long kiss_x = 123456789;
unsigned long kiss_y = 362436000;
unsigned long kiss_z = 521288629;
unsigned long kiss_c = 7654321;

@*1 Random numbers.

The random stream uses George Marsaglia's KISS generator: a congruential
generator, a shift-register generator, and two multiply-with-carry generators
combined into one long-period stream.

@ @<Definition for |kiss_rand|@>=
static unsigned long kiss_rand(void)
{
    unsigned long long t, a = 698769069ULL;

    kiss_x = 69069 * kiss_x + 12345;
    kiss_y ^= (kiss_y << 13);
    kiss_y ^= (kiss_y >> 17);
    kiss_y ^= (kiss_y << 5);
    t = a * kiss_z + kiss_c;
    kiss_c = (t >> 32);
    kiss_z = t;
    return kiss_x + kiss_y + kiss_z;
}

@ Knuth's recurrence expands a single seed into the four KISS state variables.

@ @<Definition for |kiss_rand_seed|@>=
static void kiss_rand_seed(unsigned long seed)
{
    static const unsigned long K = 1812433253UL;
    kiss_c = K * (seed ^ (seed >> 30)) + 1;
    kiss_x = K * (kiss_c ^ (kiss_c >> 30)) + 2;
    kiss_y = K * (kiss_x ^ (kiss_x >> 30)) + 3;
    kiss_z = K * (kiss_y ^ (kiss_y >> 30)) + 5;
}

@ @<Definition for |set_photon_seed|@>=
static inline void set_photon_seed(unsigned long new_seed)
{
    photon_seed = new_seed;
}

@ A public seed value of zero requests a non-repeatable seed based on the
current clock time.  Any other seed makes the Monte Carlo sequence repeatable.

@ @<Prototype for |MC_Set_Seed|@>=
void MC_Set_Seed(unsigned long seed)

@ @<Definition for |MC_Set_Seed|@>=
@<Prototype for |MC_Set_Seed|@>
{
    if (seed == 0)
        set_photon_seed((unsigned long) time(NULL));
    else
        set_photon_seed(seed);
}

@ Each photon restarts the KISS sequence from the next photon seed.  This
reduces noise when successive inverse iterations compare nearby optical
properties.

@ @<Definition for |next_photon_seed|@>=
static inline void next_photon_seed(void)
{
    photon_seed = (1812433253UL * photon_seed) & 0xffffffffUL;
    kiss_rand_seed(photon_seed);
    kiss_rand();
    kiss_rand();
    kiss_rand();
}

@ Return a uniform random number in $(0,1]$.  Zero is excluded because it is
used inside logarithms for free-path sampling.

@ @<Definition for |rand_zero_one|@>=
static double rand_zero_one(void)
{
    unsigned long x;
    double xi;

    do {
        x = kiss_rand();
    }
    while (x == 0);

    xi = ((double) x) / ((double) kiss_rand_max);

    return xi;
}

@ @<Definition for |rand_one_one|@>=
static double rand_one_one(void)
{
    return 2.0 * rand_zero_one() - 1.0;
}

@*1 Optical events.

|fresnel| returns the unpolarized Fresnel reflectance for direction cosine
|nu_i| at an interface from index |n_i| to |n_t|.

@ @<Definition for |fresnel|@>=
static double fresnel(double n_i, double n_t, double nu_i)
{
    double nu_t, ratio, temp, temp1;

    if (n_i == n_t)
        return 0.0;

    nu_i = fabs(nu_i);
    if (nu_i == 1.0)
        return sqr((n_i - n_t) / (n_i + n_t));

    ratio = n_i / n_t;
    temp = 1.0 - ratio * ratio * (1.0 - nu_i * nu_i);
    if (temp < 0)
        return 1.0;

    nu_t = sqrt(temp);
    temp = ratio * nu_t;
    temp1 = (nu_i - temp) / (nu_i + temp);
    temp = ratio * nu_i;
    temp = (nu_t - temp) / (nu_t + temp);
    return (temp1 * temp1 + temp * temp) / 2.0;
}

@ Refract the direction cosines across a plane boundary.  The assertions
check Snell's law and preserve unit direction in debug builds.

@ @<Definition for |refract|@>=
static void refract(double n_i, double n_t, double *u, double *v, double *w)
{
    double nu, c;
#ifndef NDEBUG
    double wi = *w;
#endif

    if (n_i == n_t)
        return;

    c = n_i / n_t;
    nu = (*w) * c;

    *u *= c;
    *v *= c;
    if (*w < 0)
        *w = -sqrt(1.0 - sqr(c) + sqr(nu));
    else
        *w = sqrt(1.0 - sqr(c) + sqr(nu));

    assert(CLOSE(n_i * sin(acos(wi)), n_t * sin(acos(*w))));
    assert(((*w) * wi) > 0);
    assert(CLOSE(sqr(*u) + sqr(*v) + sqr(*w), 1.0));
}

@ Henyey-Greenstein scattering changes the direction but not the photon
position.  The isotropic case is sampled directly.

@ @<Definition for |scatter|@>=
static void scatter(double g, double *u, double *v, double *w)
{
    double t1, t2, t3, mu, uu, vv, ww;

    do {
        t1 = rand_one_one();
        t2 = rand_one_one();
        t3 = t1 * t1 + t2 * t2;
    }
    while (t3 > 1);

    if (g == 0) {
        *u = 2.0 * t1 * sqrt(1.0 - t3);
        *v = 2.0 * t2 * sqrt(1.0 - t3);
        *w = 1.0 - 2.0 * t3;
        return;
    }

    mu = (1 - g * g) / (1 - g + 2.0 * g * rand_zero_one());
    mu = (1 + g * g - mu * mu) / 2.0 / g;

    uu = *u;
    vv = *v;
    ww = *w;

    if (fabs(ww) < 0.9) {
        *u = mu * uu + sqrt((1 - mu * mu) / (1 - ww * ww) / t3) * (t1 * uu * ww - t2 * vv);
        *v = mu * vv + sqrt((1 - mu * mu) / (1 - ww * ww) / t3) * (t1 * vv * ww + t2 * uu);
        *w = mu * ww - sqrt((1 - mu * mu) * (1 - ww * ww) / t3) * t1;
    }
    else {
        *u = mu * uu + sqrt((1 - mu * mu) / (1 - vv * vv) / t3) * (t1 * uu * vv + t2 * ww);
        *v = mu * vv - sqrt((1 - mu * mu) * (1 - vv * vv) / t3) * t1;
        *w = mu * ww + sqrt((1 - mu * mu) / (1 - vv * vv) / t3) * (t1 * vv * ww - t2 * uu);
    }
}

@*1 Launching and moving photons.

The initial point is on the top face of the top slide.  A positive beam radius
selects a uniform point on a disk; a zero radius launches every photon on axis.

@ @<Definition for |launch_point|@>=
static void launch_point(double *x, double *y, double *z, double beam_radius, double t_slide)
{
    *x = 0;
    *y = 0;
    *z = -t_slide;

    if (beam_radius > 0) {
        double a, b;
        do {
            a = rand_one_one();
            b = rand_one_one();
        }
        while (a * a + b * b > 1);

        *x = a * beam_radius;
        *y = b * beam_radius;
    }
}

@ A collimated beam uses |mu| as the incident direction cosine.  Diffuse
illumination is sampled over the upper hemisphere.

@ @<Definition for |launch_direction|@>=
static void launch_direction(double *u, double *v, double *w, double cos_cone_angle, double mu)
{
    double phi;
    if (cos_cone_angle == COLLIMATED) {
        *u = sqrt(1 - mu * mu);
        *v = 0;
        *w = mu;
    }
    else {
        *w = sqrt(rand_zero_one());
        phi = 2.0 * M_PI * rand_zero_one();
        *u = cos(phi) * sqrt(1 - sqr(*w));
        *v = sin(phi) * sqrt(1 - sqr(*w));
    }
}

@ Russian roulette terminates very low-weight photons without biasing the
energy balance.  The residual-weight correction is kept so final tallies can
be normalized by the surviving effective launch weight.

@ @<Definition for |roulette|@>=
static void roulette(double *weight, double *residual_weight)
{
    if (*weight > MIN_WEIGHT || *weight == 0.0)
        return;

    *residual_weight += *weight;
    if (rand_zero_one() < 0.1)
        *weight *= 10;
    else
        *weight = 0;
    *residual_weight -= *weight;
}

@ @<Definition for |move_in_sample|@>=
static void move_in_sample(double mu_t, double *x, double *y, double *z, double u, double v, double w)
{
    double step = -log(rand_zero_one()) / mu_t;
    *x += step * u;
    *y += step * v;
    *z += step * w;
}

@ @<Definition for |move_at_boundary|@>=
static void move_at_boundary(double *u, double *v, double *w, double n_i, double n_t)
{
    double r = fresnel(n_i, n_t, *w);
    if (rand_zero_one() < r) {
        *w = -(*w);
    }
    else {
        refract(n_i, n_t, u, v, w);
    }
}

@ |move_in_slide| handles an entire glass layer, including repeated internal
reflections and absorption.  The direction on entry must point from
|z_start| toward |z_end|.

@ @<Definition for |move_in_slide|@>=
static void
move_in_slide(double *x, double *y, double *z, double *u, double *v,
    double *w, double *weight, double z_start, double z_end, double b_slide, double n1, double n2, double n3)
{
    double i_x, i_y, i_z, d_bounce, r1, r2, absorbed_loss_per_bounce;

    if (z_start == z_end) {
        move_at_boundary(u, v, w, n1, n3);
        return;
    }

    assert((z_end - z_start) * (*w) > 0);

    r1 = fresnel(n1, n2, *w);
    if (rand_zero_one() <= r1) {
        *w = -(*w);
        return;
    }

    i_x = *u;
    i_y = *v;
    i_z = *w;
    refract(n1, n2, &i_x, &i_y, &i_z);

    r2 = fresnel(n2, n3, i_z);
    d_bounce = fabs((z_end - z_start) / i_z);
    assert(d_bounce > 0);

    absorbed_loss_per_bounce = exp(-fabs(b_slide / i_z));

    while (1) {
        *x += d_bounce * i_x;
        *y += d_bounce * i_y;
        *z += d_bounce * i_z;
        *weight *= absorbed_loss_per_bounce;
        assert(fabs(*z - z_end) < 1e-8);
        assert((z_end - z_start) * (i_z) > 0);

        if (rand_zero_one() > r2) {
            refract(n1, n3, u, v, w);
            return;
        }

        i_z *= -1;
        *x += d_bounce * i_x;
        *y += d_bounce * i_y;
        *z += d_bounce * i_z;
        *weight *= absorbed_loss_per_bounce;
        assert(fabs(*z - z_start) < 1e-8);
        assert((z_end - z_start) * i_z < 0);

        if (rand_zero_one() > r1) {
            *w = -(*w);
            return;
        }

        i_z *= -1;
    }
}

@ @<Definition for |milliseconds|@>=
static double milliseconds(clock_t start_time)
{
    double t;
    clock_t finish_time = clock();
    t = 1000 * (double) (finish_time - start_time) / CLOCKS_PER_SEC;
    return t;
}

@*1 Radial scoring.

These two helpers add photon weight to the radial arrays and return the
photon radius so the caller can decide whether the photon missed the port.

@ @<Definition for |add_to_reflectance_array|@>=
double add_to_reflectance_array(double x, double y, double z, double w, double weight)
{
    double r = sqrt(x * x + y * y);
    int bin = (int) (r / RADIAL_BIN_SIZE);

    if (bin > N_RADIAL_BINS - 1)
        bin = N_RADIAL_BINS - 1;

    assert(w < 0);
    assert(weight != 0);

    R_radial[bin] += weight;
    return r;
}

@ @<Definition for |add_to_transmittance_array|@>=
double add_to_transmittance_array(double x, double y, double z, double w, double weight)
{
    double r = sqrt(x * x + y * y);
    int bin = (int) (r / RADIAL_BIN_SIZE);

    if (bin > N_RADIAL_BINS - 1)
        bin = N_RADIAL_BINS - 1;

    assert(w > 0);
    assert(weight != 0);

    T_radial[bin] += weight;
    return r;
}

@*1 Transport simulation.

|MC_Radial| is the workhorse.  It returns total reflected and transmitted
fractions and the portions of those totals outside the reflection and
transmission port radii.

@ @<Prototype for |MC_Radial|@>=
void MC_Radial(long photons, double a, double b, double g, double n_sample,
    double n_slide, double cos_cone_angle, double cos_incidence,
    double t_sample, double t_slide, double b_slide,
    double dr_port, double dt_port, double d_beam, double *r_total, double *t_total, double *r_lost, double *t_lost)

@ @<Definition for |MC_Radial|@>=
@<Prototype for |MC_Radial|@>
{
    double x, y, z, u, v, w, weight;
    long i, total_photons;
    double total_weight, distance_remaining, r_beam, mu_t, r, total;

    double r_port_radius = dr_port / 2.0;
    double t_port_radius = dt_port / 2.0;
    double residual_weight = 0.0;
    double total_time = 0;
    clock_t start_time = 0;
    double absorbed = 0;

#ifndef NDEBUG
    double cos_critical = Cos_Critical_Angle(n_sample, 1.0);
    fprintf(stderr, "cos_incidence = %10.5f\n", cos_cone_angle);
    fprintf(stderr, "cos_critical = %10.5f\n", cos_critical);
    fprintf(stderr, "d_beam   = %10.5f\n", d_beam);
    fprintf(stderr, "t_sample = %10.5f\n", t_sample);
    fprintf(stderr, "t_slide  = %10.5f\n", t_slide);
    fprintf(stderr, "n_sample = %10.5f\n", n_sample);
    fprintf(stderr, "n_slide  = %10.5f\n", n_slide);
#endif

    start_time = clock();

    if (cos_cone_angle == COLLIMATED)
        r_beam = d_beam / 2.0;
    else
        r_beam = dr_port / 2.0;

    if (photons >= 0)
        total_photons = photons;
    else {
        total_photons = 1000000;
        total_time = labs(photons);
    }

    *r_lost = 0;
    *t_lost = 0;
    *r_total = 0;
    *t_total = 0;
    total_weight = 0.0;

    if (b < 1e-5)
        b = 1e-5;
    if (b > 1000)
        b = 1000;
    mu_t = b / t_sample;

    for (i = 1; i <= total_photons; i++) {
        next_photon_seed();
        launch_point(&x, &y, &z, r_beam, t_slide);
        launch_direction(&u, &v, &w, cos_cone_angle, cos_incidence);
        assert(w > 0 && z == -t_slide);

        weight = 1;
        total_weight += weight;

        move_in_slide(&x, &y, &z, &u, &v, &w, &weight, z, z + t_slide, b_slide, 1.0, n_slide, n_sample);

        if (w < 0) {
            r = add_to_reflectance_array(x, y, z, w, weight);
            *r_total += weight;
            if (r > r_port_radius)
                *r_lost += weight;
            weight = 0;
            continue;
        }

        assert(w > 0);
        assert(cos_critical < w);
        assert(CLOSE(z, 0));
        assert(weight == 1);
        assert(CLOSE(sqr(u) + sqr(v) + sqr(w), 1));

        while (weight > 0) {
            move_in_sample(mu_t, &x, &y, &z, u, v, w);

            assert(weight <= 1);

            while (z < 0 || z > t_sample) {
                distance_remaining = z / w;
                if (z > t_sample)
                    distance_remaining -= t_sample / w;

                x -= u * distance_remaining;
                y -= v * distance_remaining;
                z -= w * distance_remaining;

                assert(CLOSE(z, 0) || CLOSE(z, t_sample));

                if (w > 0)
                    move_in_slide(&x, &y, &z, &u, &v, &w, &weight, z, z + t_slide, b_slide, n_sample, n_slide, 1.0);
                else
                    move_in_slide(&x, &y, &z, &u, &v, &w, &weight, z, z - t_slide, b_slide, n_sample, n_slide, 1.0);

                if (CLOSE(z, -t_slide) && w < 0) {
                    r = add_to_reflectance_array(x, y, z, w, weight);
                    *r_total += weight;
                    if (r > r_port_radius)
                        *r_lost += weight;
                    weight = 0;
                    break;
                }

                if (CLOSE(z, t_sample + t_slide) && w > 0) {
                    r = add_to_transmittance_array(x, y, z, w, weight);
                    *t_total += weight;
                    if (r > t_port_radius)
                        *t_lost += weight;
                    weight = 0;
                    break;
                }

                assert((CLOSE(z, 0) && w > 0)
                    || (CLOSE(z, t_sample) && w < 0));
                x += u * distance_remaining;
                y += v * distance_remaining;
                z += w * distance_remaining;
            }

            if (weight > 0) {
                assert(z >= 0 && z <= t_sample);
                absorbed += (1 - a) * weight;
                weight *= a;
                roulette(&weight, &residual_weight);
                scatter(g, &u, &v, &w);

                assert(CLOSE(u * u + v * v + w * w, 1.0));
            }
        }

        if (total_time && total_time < milliseconds(start_time))
            break;
    }

    @<Normalize radial Monte Carlo tallies@>@;

    if (print_radial_arrays) {
        @<Print radial Monte Carlo arrays@>@;
    }
}

@ @<Normalize radial Monte Carlo tallies@>=
total = total_weight - residual_weight;
    *r_lost /= total;
    *t_lost /= total;
    *r_total /= total;
    *t_total /= total;

@ @<Print radial Monte Carlo arrays@>=
{
    absorbed /= total;
    fprintf(stderr, "%10d # number of bins\n", N_RADIAL_BINS);
    fprintf(stderr, "%10.5f # bin size [mm]\n", RADIAL_BIN_SIZE);
    fprintf(stderr, "# %10.5f total photons\n", total_weight);
    fprintf(stderr, "# %10.5f residual photons\n", residual_weight);
    fprintf(stderr, "# %10.5f total reflected light\n", *r_total);
    fprintf(stderr, "# %10.5f total transmitted light\n", *t_total);
    fprintf(stderr, "# %10.5f total absorbed light\n", absorbed);
    fprintf(stderr, "# %10.5f total light\n", *r_total + *t_total + absorbed);
    fprintf(stderr, "#    r    \t    R(r)    \t   T(r)\n");
    fprintf(stderr, "#    mm    \t    W/mm2    \t   W/mm2\n");

    for (i = 0; i < N_RADIAL_BINS; i++) {
        double area = M_PI * sqr(RADIAL_BIN_SIZE) * (2 * i + 1);
        fprintf(stderr, "%10.5f\t", i * RADIAL_BIN_SIZE);
        fprintf(stderr, "%10.5f\t", R_radial[i] / total / area);
        fprintf(stderr, "%10.5f\n", T_radial[i] / total / area);
    }
}

@*1 Public entry points.

|MC_Lost| translates the current |iad| measurement and inverse records into
the geometry needed by |MC_Radial|.  It assumes the reflection and
transmission sample ports use the top-slide parameters.

@ @<Prototype for |MC_Lost|@>=
void MC_Lost(struct measure_type m, struct invert_type r, long n_photons,
    double *ur1, double *ut1, double *uru, double *utu,
    double *ur1_lost, double *ut1_lost, double *uru_lost, double *utu_lost)

@ @<Definition for |MC_Lost|@>=
@<Prototype for |MC_Lost|@>
{
    double n_sample = m.slab_index;
    double n_slide = m.slab_top_slide_index;
    double t_sample = m.slab_thickness;
    double t_slide = m.slab_top_slide_thickness;
    double b_slide = m.slab_top_slide_b;

    double dr_port = sqrt(m.as_r) * 2 * m.d_sphere_r;
    double dt_port = sqrt(m.as_t) * 2 * m.d_sphere_t;
    double d_beam = m.d_beam;
    double mu = m.slab_cos_angle;

    if (t_slide == 0.0)
        n_slide = 1.0;
    if (n_slide == 1.0)
        t_slide = 0.0;

    MC_Radial(n_photons / 2, r.a, r.b, r.g, n_sample, n_slide,
        COLLIMATED, mu, t_sample, t_slide, b_slide, dr_port, dt_port, d_beam, ur1, ut1, ur1_lost, ut1_lost);

    *uru_lost = 0;
    *utu_lost = 0;

    if (m.method == SUBSTITUTION)
        MC_Radial(n_photons / 2, r.a, r.b, r.g, n_sample, n_slide,
            DIFFUSE, mu, t_sample, t_slide, b_slide, dr_port, dt_port, d_beam, uru, utu, uru_lost, utu_lost);

    if (*ur1_lost < 0 || *ut1_lost < 0 || *uru_lost < 0 || *utu_lost < 0) {
        exit(EXIT_FAILURE);
    }
}

@ |MC_RT| is a broad-port Monte Carlo check of the adding-doubling
reflection and transmission values.

@ @<Prototype for |MC_RT|@>=
void MC_RT(struct AD_slab_type s, long n_photons, double t_sample,
    double t_slide, double *UR1, double *UT1, double *URU, double *UTU)

@ @<Definition for |MC_RT|@>=
@<Prototype for |MC_RT|@>
{
    double ur1_lost, ut1_lost, uru_lost, utu_lost;
    double dr_port = 1000;
    double dt_port = 1000;
    double d_beam = 0.0;
    double mu = s.cos_angle;

    set_photon_seed(12345);

    MC_Radial(n_photons / 2, s.a, s.b, s.g, s.n_slab, s.n_top_slide,
        COLLIMATED, mu, t_sample, t_slide, s.b_top_slide, dr_port, dt_port, d_beam, UR1, UT1, &ur1_lost, &ut1_lost);

    MC_Radial(n_photons / 2, s.a, s.b, s.g, s.n_slab, s.n_top_slide,
        DIFFUSE, mu, t_sample, t_slide, s.b_top_slide, dr_port, dt_port, d_beam, URU, UTU, &uru_lost, &utu_lost);
}

@ @<Prototype for |MC_Print_RT_Arrays|@>=
void MC_Print_RT_Arrays(int status)

@ @<Definition for |MC_Print_RT_Arrays|@>=
@<Prototype for |MC_Print_RT_Arrays|@>
{
    print_radial_arrays = status;
}

@ Repeating the lost-light calculation provides a simple standard error for
diagnostic and uncertainty-estimation workflows.

@ @<Prototype for |MC_Lost_With_Stderr|@>=
void MC_Lost_With_Stderr(struct measure_type m, struct invert_type r,
    long n_photons, int n_repeats,
    double *ur1, double *ut1, double *uru, double *utu,
    double *mean_ur1_lost, double *mean_ut1_lost,
    double *mean_uru_lost, double *mean_utu_lost,
    double *se_ur1_lost, double *se_ut1_lost, double *se_uru_lost, double *se_utu_lost)

@ @<Definition for |MC_Lost_With_Stderr|@>=
@<Prototype for |MC_Lost_With_Stderr|@>
{
    int i;
    long n_per_run;
    double sum_ur1 = 0, sum_ut1 = 0, sum_uru = 0, sum_utu = 0;
    double sum2_ur1 = 0, sum2_ut1 = 0, sum2_uru = 0, sum2_utu = 0;
    double v_ur1, v_ut1, v_uru, v_utu;
    double cur_ur1_lost, cur_ut1_lost, cur_uru_lost, cur_utu_lost;
    double n = (double) n_repeats;

    if (n_repeats < 1)
        n_repeats = 1;
    n_per_run = n_photons / n_repeats;
    if (n_per_run < 1)
        n_per_run = 1;

    for (i = 0; i < n_repeats; i++) {
        MC_Lost(m, r, n_per_run, ur1, ut1, uru, utu, &cur_ur1_lost, &cur_ut1_lost, &cur_uru_lost, &cur_utu_lost);
        sum_ur1 += cur_ur1_lost;
        sum_ut1 += cur_ut1_lost;
        sum_uru += cur_uru_lost;
        sum_utu += cur_utu_lost;
        sum2_ur1 += cur_ur1_lost * cur_ur1_lost;
        sum2_ut1 += cur_ut1_lost * cur_ut1_lost;
        sum2_uru += cur_uru_lost * cur_uru_lost;
        sum2_utu += cur_utu_lost * cur_utu_lost;
    }

    *mean_ur1_lost = sum_ur1 / n;
    *mean_ut1_lost = sum_ut1 / n;
    *mean_uru_lost = sum_uru / n;
    *mean_utu_lost = sum_utu / n;

    if (n_repeats > 1) {
        v_ur1 = (sum2_ur1 - n * (*mean_ur1_lost) * (*mean_ur1_lost)) / (n - 1.0);
        v_ut1 = (sum2_ut1 - n * (*mean_ut1_lost) * (*mean_ut1_lost)) / (n - 1.0);
        v_uru = (sum2_uru - n * (*mean_uru_lost) * (*mean_uru_lost)) / (n - 1.0);
        v_utu = (sum2_utu - n * (*mean_utu_lost) * (*mean_utu_lost)) / (n - 1.0);
        *se_ur1_lost = sqrt(fmax(0.0, v_ur1) / n);
        *se_ut1_lost = sqrt(fmax(0.0, v_ut1) / n);
        *se_uru_lost = sqrt(fmax(0.0, v_uru) / n);
        *se_utu_lost = sqrt(fmax(0.0, v_utu) / n);
    }
    else {
        *se_ur1_lost = 0.0;
        *se_ut1_lost = 0.0;
        *se_uru_lost = 0.0;
        *se_utu_lost = 0.0;
    }
}
