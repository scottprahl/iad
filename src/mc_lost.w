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

@ The header exposes the Monte Carlo services used by |iad| and by the
standalone |mc_lost| diagnostic program.

@(mc_lost.h@>=
@<Prototype for |MC_Set_Seed|@>;
@<Prototype for |MC_Lost|@>;
@<Prototype for |MC_RT|@>;
@<Prototype for |MC_Radial|@>;
@<Prototype for |MC_Print_RT_Arrays|@>;

@*1 Constants and state.

Photons with very small weight are handled by Russian roulette rather than
being followed forever.  The radial arrays are intentionally global so the
standalone executable can request a radial dump after a simulation.  They
describe one simulation, so |MC_Radial| clears them before it launches any
photons; without that they would accumulate every simulation made during the
run, and a dump taken after the second of two calls would report the sum of
both.

@ @<Monte Carlo constants@>=
#define MIN_WEIGHT 0.0001
#define N_RADIAL_BINS  1001
#define RADIAL_BIN_SIZE 0.02
#define CLOSE(x, y) (fabs((x) - (y)) < 1e-8)
#define DIFFUSE_SEED_OFFSET 0x9e3779b9UL

@ The KISS generator state is process-global.  |photon_seed| is advanced once
per photon so repeated simulations can compare nearby optical properties with
the same sequence of photon histories.

@ @<Monte Carlo global state@>=
unsigned long photon_seed = 12345678;
unsigned long lost_base_seed = 12345678;

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
|lost_base_seed| remembers the choice so |MC_Lost| can return to it.

@ @<Prototype for |MC_Set_Seed|@>=
void MC_Set_Seed(unsigned long seed)

@ @<Definition for |MC_Set_Seed|@>=
@<Prototype for |MC_Set_Seed|@>
{
    if (seed == 0)
        lost_base_seed = (unsigned long) time(NULL);
    else
        lost_base_seed = seed;

    set_photon_seed(lost_base_seed);
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

A photon that has just crossed the top slide and entered the sample has not
scattered yet, but it may already weigh less than one: an absorbing slide
attenuates it on the way through.  The assertion there checks that the weight
is still positive and has not grown, not that it is exactly one, which held
only while every slide was perfectly clear.

@ The two slides are described separately.  The top slide is the one the
light reaches first; the bottom slide is the one on the far side of the
sample.  A slide is absent when its index is one, and then its thickness must
be zero so that the geometry below collapses to a bare surface.

@<Prototype for |MC_Radial|@>=
void MC_Radial(long photons, double a, double b, double g, double n_sample,
    double n_top_slide, double n_bottom_slide,
    double cos_cone_angle, double cos_incidence,
    double t_sample, double t_top_slide, double t_bottom_slide,
    double b_top_slide, double b_bottom_slide,
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
    fprintf(stderr, "illumination  = %s\n",
            (cos_cone_angle == COLLIMATED) ? "collimated" : "diffuse");
    fprintf(stderr, "cos_incidence = %10.5f\n", cos_incidence);
    fprintf(stderr, "cos_critical = %10.5f\n", cos_critical);
    fprintf(stderr, "d_beam   = %10.5f\n", d_beam);
    fprintf(stderr, "t_sample = %10.5f\n", t_sample);
    fprintf(stderr, "t_top    = %10.5f\n", t_top_slide);
    fprintf(stderr, "t_bottom = %10.5f\n", t_bottom_slide);
    fprintf(stderr, "n_sample = %10.5f\n", n_sample);
    fprintf(stderr, "n_top    = %10.5f\n", n_top_slide);
    fprintf(stderr, "n_bottom = %10.5f\n", n_bottom_slide);
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

    for (i = 0; i < N_RADIAL_BINS; i++) {
        R_radial[i] = 0;
        T_radial[i] = 0;
    }

    if (b < 1e-5)
        b = 1e-5;
    if (b > 1000)
        b = 1000;
    mu_t = b / t_sample;

    for (i = 1; i <= total_photons; i++) {
        next_photon_seed();
        launch_point(&x, &y, &z, r_beam, t_top_slide);
        launch_direction(&u, &v, &w, cos_cone_angle, cos_incidence);
        assert(w > 0 && z == -t_top_slide);

        weight = 1;
        total_weight += weight;

        move_in_slide(&x, &y, &z, &u, &v, &w, &weight, z, z + t_top_slide,
            b_top_slide, 1.0, n_top_slide, n_sample);

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
        assert(weight > 0 && weight <= 1);
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
                    move_in_slide(&x, &y, &z, &u, &v, &w, &weight, z, z + t_bottom_slide,
                        b_bottom_slide, n_sample, n_bottom_slide, 1.0);
                else
                    move_in_slide(&x, &y, &z, &u, &v, &w, &weight, z, z - t_top_slide,
                        b_top_slide, n_sample, n_top_slide, 1.0);

                if (CLOSE(z, -t_top_slide) && w < 0) {
                    r = add_to_reflectance_array(x, y, z, w, weight);
                    *r_total += weight;
                    if (r > r_port_radius)
                        *r_lost += weight;
                    weight = 0;
                    break;
                }

                if (CLOSE(z, t_sample + t_bottom_slide) && w > 0) {
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

Each call restarts from |lost_base_seed| so that photon $i$ follows the same
random history every time.  This is the whole reason |next_photon_seed|
reseeds once per photon.  The inverse routine calls |MC_Lost| repeatedly with
optical properties that barely change, and it converges by watching how much
the lost-light estimate moved.  Drawing a fresh slice of the random stream
each call would put the full Monte Carlo standard error into that difference
and hide the real change; sharing the histories cancels most of it.  The two
simulations need independent streams, so the diffuse run is offset from the
collimated one.  Photon counts may differ between calls --- the shorter run
is then a prefix of the longer one, which keeps the estimates correlated.

@ @<Prototype for |MC_Lost|@>=
void MC_Lost(struct measure_type m, struct invert_type r, long n_photons,
    double *ur1, double *ut1, double *uru, double *utu,
    struct lost_type *lost_r, struct lost_type *lost_t, double *utu_lost)

@ @<Definition for |MC_Lost|@>=
@<Prototype for |MC_Lost|@>
{
    double n_sample = m.slab_index;
    double t_sample = m.slab_thickness;

    double n_top = m.slab_top_slide_index;
    double t_top = m.slab_top_slide_thickness;
    double b_top = m.slab_top_slide_b;

    double n_bottom = m.slab_bottom_slide_index;
    double t_bottom = m.slab_bottom_slide_thickness;
    double b_bottom = m.slab_bottom_slide_b;

    double dr_port = sqrt(m.as_r) * 2 * m.d_sphere_r;
    double dt_port = sqrt(m.as_t) * 2 * m.d_sphere_t;
    double d_beam = m.d_beam;
    double mu = m.slab_cos_angle;

    int slides_differ;

    @<Give up when no sphere describes the ports@>@;

    @<Reconcile each slide index with its thickness@>@;

    slides_differ = (n_top != n_bottom) || (t_top != t_bottom) || (b_top != b_bottom);

    set_photon_seed(lost_base_seed);
    MC_Radial(n_photons / 2, r.a, r.b, r.g, n_sample, n_top, n_bottom,
        COLLIMATED, mu, t_sample, t_top, t_bottom, b_top, b_bottom,
        dr_port, dt_port, d_beam, ur1, ut1, &lost_r->direct, &lost_t->direct);

    lost_r->diffuse = 0;
    lost_t->diffuse = 0;
    *utu_lost = 0;

    if (m.method == SUBSTITUTION) {
        set_photon_seed(lost_base_seed + DIFFUSE_SEED_OFFSET);
        MC_Radial(n_photons / 2, r.a, r.b, r.g, n_sample, n_top, n_bottom,
            DIFFUSE, mu, t_sample, t_top, t_bottom, b_top, b_bottom,
            dr_port, dt_port, d_beam, uru, utu, &lost_r->diffuse, utu_lost);

        @<Find the diffuse reflectance loss seen by the transmission sphere@>@;
    }

    if (m.flip_sample && slides_differ)
        @<Take the transmission losses from the flipped sample@>@;

    if (lost_r->direct < 0 || lost_t->direct < 0 ||
        lost_r->diffuse < 0 || *utu_lost < 0) {
        exit(EXIT_FAILURE);
    }
}

@ Port sizes come from the sphere descriptions and from nowhere else, so with
no sphere there is no port and no way to say what missed it.  |measure_OK|
only validates the sphere parameters when |num_spheres| is nonzero, which
means |as_r| and |as_t| may still hold whatever |Initialize_Measure| left
behind.  Deriving |dr_port| and |dt_port| from those gives a port of zero
diameter, and a port of zero diameter catches nothing: every photon is scored
as lost and the inversion is handed a correction that is pure fiction.

The command line already declines to enter the Monte Carlo loop without a
sphere.  The library entry points take the number of runs straight from the
caller, so the same guard belongs here too, where every caller passes.

Nothing is simulated, so the losses are zero.  The reflectances and
transmittances are zeroed as well rather than left holding whatever the
caller had in those variables.

@<Give up when no sphere describes the ports@>=

    if (m.num_spheres <= 0) {
        *ur1 = 0;
        *ut1 = 0;
        *uru = 0;
        *utu = 0;
        lost_r->direct = 0;
        lost_r->diffuse = 0;
        lost_t->direct = 0;
        lost_t->diffuse = 0;
        *utu_lost = 0;
        return;
    }

@ The geometry needs index and thickness to agree about whether a slide is
there at all.  A slide of zero thickness cannot refract, and a slide of index
one is not a slide, so each face is reconciled on its own.  Doing this per
face is the whole point: reconciling only the top and then reusing it for
both is what made \.{-G t} model slides on both faces and \.{-G b} model none.

@<Reconcile each slide index with its thickness@>=

    if (t_top == 0.0)
        n_top = 1.0;
    if (n_top == 1.0)
        t_top = 0.0;

    if (t_bottom == 0.0)
        n_bottom = 1.0;
    if (n_bottom == 1.0)
        t_bottom = 0.0;

@ Both spheres bounce diffuse light off the sample and both lose some of it
out past the rim of their own sample port, so each needs its own figure.  Two
things separate them.

The first is the port.  The run above floods the reflection port and scores
what comes back inside it, which is what |M_R| wants.  The transmission sphere
asks the same question through a port that need not be the same size, and with
one sphere at a time it often is not.  A port half the diameter loses a great
deal more light, so handing the transmission sphere the reflection number was
wrong by roughly the size of the whole correction.

The second is the face.  The beam enters the top of the sample and the
transmission sphere collects what comes out the bottom, so the light rattling
around inside that sphere strikes the sample from {\it below}.  The loss it
suffers is the loss for diffuse light entering the bottom face, which is why
the slides are exchanged below exactly as they are for a flipped sample.  How
much this matters depends entirely on how different the two faces are: with a
slide on top and none underneath the two figures stand at 0.23 and 0.07, a
factor of three apart.  The total |uru| barely moves --- reciprocity holds it
near 0.742 either way --- but its radial spread does not have to be reciprocal,
and it is the spread that decides what clears the rim.

So the second run is skipped only when the geometry makes it redundant: equal
ports, which |measure_OK| guarantees for two spheres, and slides alike, which
makes the two faces interchangeable.  The photon seed is reset to the value
the first diffuse run used, so the two share their random histories and any
difference between them reflects the change in geometry rather than noise.

@<Find the diffuse reflectance loss seen by the transmission sphere@>=
{
    if (dt_port == dr_port && !slides_differ) {
        lost_t->diffuse = lost_r->diffuse;
    } else {
        double flooded_uru, flooded_utu, flooded_utu_lost;

        set_photon_seed(lost_base_seed + DIFFUSE_SEED_OFFSET);
        MC_Radial(n_photons / 2, r.a, r.b, r.g, n_sample, n_bottom, n_top,
            DIFFUSE, mu, t_sample, t_bottom, t_top, b_bottom, b_top,
            dt_port, dt_port, d_beam, &flooded_uru, &flooded_utu,
            &lost_t->diffuse, &flooded_utu_lost);
    }
}

@ \.{-G n} and \.{-G f} describe a sample that is physically turned over
between the two measurements, so the slide stays either against the sphere or
away from it in both.  The sample is stored as the reflection measurement
sees it, which is why the reflectances above need no attention.  The
transmission measurement sees the other face first.

Adding-doubling can ignore this because total transmittance is reciprocal ---
|RT_Flip| swaps the slides, recomputes, and gets the same number back.  Lost
light is not so obliging.  The total that gets through is reciprocal but its
radial spread is not, since it matters a great deal which face carries the
glass the light has to cross.  Measured over two million photons, |ut1_lost|
is 0.0317 with the slide on top and 0.0324 with it on the bottom, while
|ur1_lost| moves from 0.0336 to 0.0307.  Only |utu_lost| comes out the same
either way, diffuse light being symmetric in both directions.

So the transmission losses, and only those, are taken from a second pair of
simulations with the slides exchanged.  This mirrors |RT_Flip| exactly.  The
extra work is skipped unless the slides really do differ, so the ordinary
cases --- no slides, matching slides, every \.{.rxt} file --- cost nothing.

@<Take the transmission losses from the flipped sample@>=
{
    double flipped_ur1, flipped_uru, flipped_ur1_lost, flipped_uru_lost;

    set_photon_seed(lost_base_seed);
    MC_Radial(n_photons / 2, r.a, r.b, r.g, n_sample, n_bottom, n_top,
        COLLIMATED, mu, t_sample, t_bottom, t_top, b_bottom, b_top,
        dr_port, dt_port, d_beam, &flipped_ur1, ut1, &flipped_ur1_lost, &lost_t->direct);

    if (m.method == SUBSTITUTION) {
        set_photon_seed(lost_base_seed + DIFFUSE_SEED_OFFSET);
        MC_Radial(n_photons / 2, r.a, r.b, r.g, n_sample, n_bottom, n_top,
            DIFFUSE, mu, t_sample, t_bottom, t_top, b_bottom, b_top,
            dr_port, dt_port, d_beam, &flipped_uru, utu, &flipped_uru_lost, utu_lost);
    }
}

@ |MC_RT| is a broad-port Monte Carlo check of the adding-doubling
reflection and transmission values.

The slab already carries both slide indices and both slide absorptions; only
the thicknesses have to be supplied, since |AD_slab_type| has no room for
them.

@ @<Prototype for |MC_RT|@>=
void MC_RT(struct AD_slab_type s, long n_photons, double t_sample,
    double t_top_slide, double t_bottom_slide,
    double *UR1, double *UT1, double *URU, double *UTU)

@ @<Definition for |MC_RT|@>=
@<Prototype for |MC_RT|@>
{
    double ur1_lost, ut1_lost, uru_lost, utu_lost;
    double dr_port = 1000;
    double dt_port = 1000;
    double d_beam = 0.0;
    double mu = s.cos_angle;

    set_photon_seed(12345);

    MC_Radial(n_photons / 2, s.a, s.b, s.g, s.n_slab, s.n_top_slide, s.n_bottom_slide,
        COLLIMATED, mu, t_sample, t_top_slide, t_bottom_slide,
        s.b_top_slide, s.b_bottom_slide, dr_port, dt_port, d_beam, UR1, UT1, &ur1_lost, &ut1_lost);

    MC_Radial(n_photons / 2, s.a, s.b, s.g, s.n_slab, s.n_top_slide, s.n_bottom_slide,
        DIFFUSE, mu, t_sample, t_top_slide, t_bottom_slide,
        s.b_top_slide, s.b_bottom_slide, dr_port, dt_port, d_beam, URU, UTU, &uru_lost, &utu_lost);
}

@ @<Prototype for |MC_Print_RT_Arrays|@>=
void MC_Print_RT_Arrays(int status)

@ @<Definition for |MC_Print_RT_Arrays|@>=
@<Prototype for |MC_Print_RT_Arrays|@>
{
    print_radial_arrays = status;
}
