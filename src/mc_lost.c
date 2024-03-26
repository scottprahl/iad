#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "ad_globl.h"
#include "iad_type.h"

#define MIN_WEIGHT 0.001

unsigned long photon_seed = 12345678;

/*
 * `kiss_rand_max` represents the maximum value that `kiss_rand` can return,
 * corresponding to ULONG_MAX, ensuring compatibility across different platforms.
 *
 * Variables kiss_x, kiss_y, kiss_z, and kiss_c are state variables that must be
 * initialized with seed values. The function `kiss_rand` computes and returns
 * the next random number in the sequence, maintaining state between calls.
 */
unsigned long kiss_rand_max = ULONG_MAX;
unsigned long kiss_x = 123456789;
unsigned long kiss_y = 362436000;
unsigned long kiss_z = 521288629;
unsigned long kiss_c = 7654321;

/*
 * KISS (Keep It Simple Stupid) Random Number Generator
 *
 * This RNG combines multiple algorithms for a balanced approach to randomness:
 * - Two Multiply-With-Carry (MWC) generators
 * - A 3-shift register (SHR3)
 * - A Congruential generator (CONG)
 *
 * Operations involve addition, exclusive-or, and bit shifts to produce a sequence
 * with a long period (~2^123), designed by George Marsaglia. This makes the KISS
 * generator suitable for simulations requiring high-quality randomness.
 */
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

/*
 * Initializes the state variables of the KISS random number generator using
 * a seed value based on Donald Knuth's algorithms. This seeding approach
 * leverages a linear congruential generator with multiplier 'K', derived from
 * Knuth's suggestions for generating initial states. The seed influences
 * the initial values of kiss_x, kiss_y, kiss_z, and kiss_c, ensuring that
 * different seeds produce different sequences of random numbers.
 *
 * Each state variable is initialized in a chained manner, where each depends
 * on the previous one, starting with 'kiss_c' based directly on the provided
 * 'seed'. This method helps in dispersing the seed's influence throughout
 * all state variables, enhancing the randomness of initial values.
 *
 * The constant 'K' and the specific transformations applied to 'seed' and the
 * state variables follow recommended practices for initializing random number
 * generators to avoid early patterns or correlations in the generated sequence.
 */
static void kiss_rand_seed(unsigned long seed)
{
    static const unsigned long K = 1812433253UL;
    kiss_c = K * (seed ^ (seed >> 30)) + 1;
    kiss_x = K * (kiss_c ^ (kiss_c >> 30)) + 2;
    kiss_y = K * (kiss_x ^ (kiss_x >> 30)) + 3;
    kiss_z = K * (kiss_y ^ (kiss_y >> 30)) + 5;
}

/* Each photon restarts the KISS random number generator.
   We are interested in subtle changes between successive
   photon paths.  The best way to track these changes is to
   use the same set of random numbers for each photon path.
   new_photon_seed() will use the next random number in a
   sequence.
   We have another random number generator that produces a
   sequence of random numbers for each photon.  This means
   that photon 1001 will reset the KISS generator to the same
   starting number.
*/

/*
 * Configures the KISS random number generator.
 * This allows consistent random number sequences for each simulation and
 * therefore minimizes variation between runs due to small changes in the
 * optical properties arising in successive runs.
 */
static inline void set_photon_seed(unsigned long new_seed)
{
    photon_seed = new_seed;
}

/*
 * The function `set_photon_seed` initializes the seed for photon-specific
 * random number generation, ensuring a distinct starting point for each
 * photon path. `next_photon_seed` updates the seed for the next photon,
 * maintaining the sequence's continuity and resetting the KISS generator
 * with a new seed derived from the previous one. This setup guarantees
 * that each photon, such as photon number 1001, consistently starts with
 * the same initial conditions in the random number sequence, allowing for
 * repeatable and comparable simulations across different photon paths.
 */
static inline void next_photon_seed(void)
{
    photon_seed = (1812433253UL * photon_seed) & 0xffffffffUL;
    kiss_rand_seed(photon_seed);
    kiss_rand();
    kiss_rand();
    kiss_rand();
}

/* random number between 0.0 and 1.0 ... EXCLUDING 0.0 */
static double rand_zero_one(void)
{
    unsigned long x;
    double xi;

    do {
        x = kiss_rand();
    } while (x == 0);

    xi = ((double) x) / ((double) kiss_rand_max);

    return xi;
}

/* random number between -1.0 and 1.0 */
static double rand_one_one(void)
{
    return 2.0 * rand_zero_one() - 1.0;
}

static double sqr(double x)
{
    return x * x;
}

/* Fresnel Reflectance
   incident    medium has index n_i and angle nu_i on
   transmitted medium has index n_t and angle nu_t
*/
static double fresnel(double n_i, double n_t, double nu_i)
{
    double nu_t, ratio, temp, temp1;

    /* matched boundaries? */
    if (n_i == n_t)
        return 0.0;

    /* normal illumination? */
    nu_i = fabs(nu_i);
    if (nu_i == 1.0)
        return sqr((n_i - n_t) / (n_i + n_t));

    /* total internal reflection? */
    ratio = n_i / n_t;
    temp = 1.0 - ratio * ratio * (1.0 - nu_i * nu_i);
    if (temp < 0) return 1.0;

    /* Now calculate Fresnel reflection */
    nu_t = sqrt(temp);
    temp = ratio * nu_t;
    temp1 = (nu_i - temp) / (nu_i + temp);
    temp = ratio * nu_i;
    temp = (nu_t - temp) / (nu_t + temp);
    return (temp1 * temp1 + temp * temp) / 2.0;
}

/* refract a photon into or out of the medium across a z-plane*/
static void refract(double n_i, double n_t, double *u, double *v, double *w)
{
    double nu, c;

    if (n_i == n_t) return;

    c = n_i / n_t;
    nu = (*w) * c;

    *u *= c;
    *v *= c;
    if (*w < 0)
        *w = -sqrt(1.0 - sqr(c) + sqr(nu));
    else
        *w = sqrt(1.0 - sqr(c) + sqr(nu));
}

/* Scatter photon and establish new direction */
static void scatter(double g, double *u, double *v, double *w)
{
    double t1, t2, t3, mu, uu, vv, ww;

    do {
        t1 = rand_one_one();
        t2 = rand_one_one();
        t3 = t1 * t1 + t2 * t2;
    } while (t3 > 1);

    if (g == 0) {               /* isotropic */
        if (t3 == 0) {
            *u = 1;
            *v = 0;
            *w = 0;
        }
        else {
            *u = 2.0 * t3 - 1.0;
            t3 = sqrt((1.0 - sqr(*u)) / t3);
            *v = t1 * t3;
            *w = t2 * t3;
        }
        return;
    }

    mu = (1 - g * g) / (1 - g + 2.0 * g * rand_zero_one());
    mu = (1 + g * g - mu * mu) / 2.0 / g;

    uu = *u;
    vv = *v;
    ww = *w;

    if (fabs(ww) < 0.9) {
        *u = mu * uu + sqrt((1 - mu * mu) / (1 -
                ww * ww) / t3) * (t1 * uu * ww - t2 * vv);
        *v = mu * vv + sqrt((1 - mu * mu) / (1 -
                ww * ww) / t3) * (t1 * vv * ww + t2 * uu);
        *w = mu * ww - sqrt((1 - mu * mu) * (1 - ww * ww) / t3) * t1;
    }
    else {
        *u = mu * uu + sqrt((1 - mu * mu) / (1 -
                vv * vv) / t3) * (t1 * uu * vv + t2 * ww);
        *v = mu * vv - sqrt((1 - mu * mu) * (1 - vv * vv) / t3) * t1;
        *w = mu * ww + sqrt((1 - mu * mu) / (1 -
                vv * vv) / t3) * (t1 * vv * ww - t2 * uu);
    }
}

/* set photon launch point*/
static void launch_point(double *x, double *y, double *z, double beam_radius)
{
    *x = 0;
    *y = 0;
    *z = 0;

    if (beam_radius > 0) {      /* uniform distribution over a disk */
        double a, b;
        do {
            a = rand_one_one();
            b = rand_one_one();
        } while (a * a + b * b > 1);

        *x = a * beam_radius;
        *y = b * beam_radius;
    }
}

/*
 * Initializes the launch direction of a photon.
 * If collimated=True, the photon is launched at an angle mu=cos(theta) from the z-axis.
 * If not collimated, the photon's direction is uniformly distributed over the hemisphere,
 * ensuring a diffuse isotropic emission.
 */
static void launch_direction(double *u, double *v, double *w, int collimated, double mu)
{
    double phi;
    if (collimated) {
        *u = sqrt(1 - mu * mu);
        *v = 0;
        *w = mu;
    } else {
        *w = sqrt(rand_zero_one());
        phi = 2.0 * M_PI * rand_zero_one();
        *u = cos(phi) * sqrt(1 - sqr(*w));
        *v = sin(phi) * sqrt(1 - sqr(*w));
    }
}

/*
 * Implements the Russian roulette technique to unbiasedly terminate a photon's path.
 * If the photon's weight is below MIN_WEIGHT (but not zero), it has a chance to be terminated
 * or to have its weight amplified and continue its path. This decision is made stochastically,
 * giving the photon a 10% chance to survive (weight multiplied by 10) and a 90% chance
 * to be terminated (weight set to 0). This approach ensures unbiased termination by
 * conserving energy: the residual weight is adjusted to account for the absorbed or
 * scattered light, to maintain an exact overall energy balance.
 */
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

/*
 * Advances the photon by a random step within the sample, determined by the
 * medium's total attenuation coefficient (mu_t). The step size follows an
 * exponential distribution, reflecting the natural path lengths in a
 * scattering and absorbing medium. The photon's new position (x, y, z) is
 * updated according to this step, based on its current direction cosines
 * (u, v, w).
 */
static void move_in_sample(double mu_t, double *x, double *y, double *z,
                           double u, double v, double w)
{
    double step = -log(rand_zero_one()) / mu_t;
    *x += step * u;
    *y += step * v;
    *z += step * w;
}

/*
 * On entry the photon is on one boundary of the slide and the direction is
 * into the slide.  On exit, x and y will be the location that the photon
 * leaves the slide and the direction w should be correct.
 * The value of z does not need to be updated because if the ray exits
 * then it does not matter.  If the ray doesn't exit, then it is reflected
 * back into the medium and should have the same value as when it entered
 */
static void move_in_slide(double *x, double *y, double u, double v, double *w,
    double *weight, double t_slide, double mua_slide, double n1, double n2,
    double n3)
{
    double i_x, i_y, i_z, d_bounce, r1, r2, absorbed_light_per_bounce;

    r1 = fresnel(n1, n2, *w);

    if (rand_zero_one() <= r1) {
        /* specularly reflected photon */
        *w = -(*w);
        return;
    }

    /* direction cosines in the slide */
    i_x = u;
    i_y = v;
    i_z = *w;
    refract(n1, n2, &i_x, &i_y, &i_z);

    r2 = fresnel(n2, n3, i_z);
    d_bounce = fabs(t_slide / i_z);
    absorbed_light_per_bounce = 1 - exp(-mua_slide * d_bounce);

    /* photon is in slide and bounces in the slide until it exits */
    while (1) {
        *x += d_bounce * i_x;
        *y += d_bounce * i_y;
        *weight *= 1 - absorbed_light_per_bounce;

        if (rand_zero_one() > r2)
            return;

        *x += d_bounce * i_x;
        *y += d_bounce * i_y;
        *weight *= 1 - absorbed_light_per_bounce;

        if (rand_zero_one() > r1) {
            *w = -(*w);
            return;
        }
    }
}

static double milliseconds(clock_t start_time)
{
    double t;
    clock_t finish_time = clock();
    t = 1000 * (double) (finish_time - start_time) / CLOCKS_PER_SEC;
    return t;
}

void MC_Radial(long photons, double a, double b, double g, double n_sample,
    double n_slide, int collimated, double cos_incidence,
    double t_sample, double t_slide, double mua_slide, double dr_port, double dt_port,
    double d_beam, double *r_total, double *t_total, double *r_lost,
    double *t_lost)
{

    double x, y, z, u, v, w, weight;
    long i, total_photons;
    double residual_weight = 0.0;
    double r_beam;
    double mu_t = b / t_sample;
    double rr_port_squared = dr_port * dr_port / 4.0;
    double rt_port_squared = dt_port * dt_port / 4.0;
    double total_weight, distance_remaining;
    double total_time = 0;
    clock_t start_time = 0;
#ifndef NDEBUG
    double ww;
#endif
    /* high resolution clock ... may be zero at start of process */
    start_time = clock();

    if (collimated)
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

    if (a==0 && collimated)
        return; /* no possible lost light (some diffuse utu losses possible tho) */

    if (dr_port == 0 && dt_port == 0)  /* no port, no loss! */
        return;

    total_weight = 0.0;

    for (i = 1; i <= total_photons; i++) {
        next_photon_seed();
        launch_point(&x, &y, &z, r_beam);
        launch_direction(&u, &v, &w, collimated, cos_incidence);
        weight = 1;
        total_weight += weight;

        /* first interaction from air to slide to sample */
        move_in_slide(&x, &y, u, v, &w, &weight, t_slide, mua_slide, 1.0,
            n_slide, n_sample);

        if (w < 0) {
            /* specularly reflected from surface */
            *r_total += weight;
            if (x * x + y * y > rr_port_squared)
                *r_lost += weight;
            weight = 0;
            continue;
        }

#ifndef NDEBUG
        ww = w;
#endif
        assert(w > 0);
        assert(z == 0);
        assert(weight == 1);

        /* enter the sample */
        refract(1.0, n_sample, &u, &v, &w);

        assert(fabs(u * u + v * v + w * w - 1.0) < 1e-8);
        while (weight > 0) {

            /* move to next scattering or absorption event */
            move_in_sample(mu_t, &x, &y, &z, u, v, w);

            assert(weight <= 1);

            /* deal with photon that may exit the sample */
            while (z < 0 || z > t_sample) {

                /* move back to top or bottom of sample */
                if (z < 0) {
                    distance_remaining = z / w;
                    z = 0.0;
                }
                else {
                    distance_remaining = (z - t_sample) / w;
                    z = t_sample;
                }
                assert(distance_remaining >= 0);

                /* move x,y coordinates back to intersection with boundary */
                x -= u * distance_remaining;
                y -= v * distance_remaining;

                assert(z == 0 || z == t_sample);

                /* w will change sign upon reflection back into sample */
                move_in_slide(&x, &y, u, v, &w, &weight, t_slide, mua_slide,
                    n_sample, n_slide, 1.0);

                /* use direction w to determine if photon is transmitted or reflected */
                if (z == 0) {

                    if (w < 0) {        /* at top surface and moving up */
                        assert(z == 0);
                        *r_total += weight;
                        if (x * x + y * y > rr_port_squared)
                            *r_lost += weight;
                        weight = 0;
                        break;
                    }

                }
                else {

                    if (w > 0) {        /* at bottom and moving down */
                        assert(z == t_sample);
                        *t_total += weight;
                        if (x * x + y * y > rt_port_squared)
                            *t_lost += weight;
                        weight = 0;
                        break;
                    }
                }

                /* move the distance_remaining amount back into sample */
                x += u * distance_remaining;
                y += v * distance_remaining;
                z += w * distance_remaining;
            }

            assert(z >= 0 && z <= t_sample);
            weight *= a;
            roulette(&weight, &residual_weight);
            scatter(g, &u, &v, &w);
        }

        if (total_time && total_time < milliseconds(start_time))
            break;
    }

    *r_lost /= (total_weight + residual_weight);
    *t_lost /= (total_weight + residual_weight);
    *r_total /= (total_weight + residual_weight);
    *t_total /= (total_weight + residual_weight);
}

/* This assumes that the sample port diameter for the reflection and transmission
   measurements are the same. */
void MC_Lost(struct measure_type m, struct invert_type r, long n_photons,
    double *ur1, double *ut1, double *uru, double *utu,
    double *ur1_lost, double *ut1_lost, double *uru_lost, double *utu_lost)
{
    double mua_slide;
    int collimated = 1;
    int diffuse = 0;

    double n_sample = m.slab_index;
    double n_slide = m.slab_top_slide_index;
    double t_sample = m.slab_thickness;
    double t_slide = m.slab_top_slide_thickness;
    double dr_port = sqrt(m.as_r) * 2 * m.d_sphere_r;    /* sample for rsphere */
    double dt_port = sqrt(m.as_t) * 2 * m.d_sphere_t;    /* sample for tsphere */
    double d_beam = m.d_beam;
    double mu = m.slab_cos_angle;

    /* no slide if thickness is zero or the index is 1.0 */
    if (t_slide == 0.0)
        n_slide = 1.0;
    if (n_slide == 1.0)
        t_slide = 0.0;

    if (t_slide == 0)
        mua_slide = 0;
    else
        mua_slide = m.slab_top_slide_b / t_slide;

    // fprintf(stderr, "  reflection port diameter = %8.3f mm\n", dr_port);
    // fprintf(stderr, "transmission port diameter = %8.3f mm\n", dt_port);

    set_photon_seed(12345);
    MC_Radial(n_photons / 2, r.a, r.b, r.g, n_sample, n_slide, collimated, mu,
        t_sample, t_slide, mua_slide, dr_port, dt_port, d_beam, ur1, ut1, ur1_lost,
        ut1_lost);

    *uru_lost = 0;
    *utu_lost = 0;
    if (m.method == SUBSTITUTION)
        MC_Radial(n_photons / 2, r.a, r.b, r.g, n_sample, n_slide, diffuse, mu,
            t_sample, t_slide, mua_slide, dr_port, dt_port, d_beam, uru, utu, uru_lost,
            utu_lost);

    if (*ur1_lost < 0)
        *ur1_lost = 0;
    if (*ut1_lost < 0)
        *ut1_lost = 0;
    if (*uru_lost < 0)
        *uru_lost = 0;
    if (*utu_lost < 0)
        *utu_lost = 0;
}

void MC_RT(struct AD_slab_type s, long n_photons, double *UR1, double *UT1,
    double *URU, double *UTU)
{
    double ur1_lost, ut1_lost, uru_lost, utu_lost;
    int collimated = 1;
    int diffuse = 0;
    double t_sample = 1.0;
    double dr_port = 1000;
    double dt_port = 1000;
    double d_beam = 0.0;
    double t_slide = 0.0;
    double mu = s.cos_angle;
    double mua_slide = 0.0;

    set_photon_seed(12345);
    MC_Radial(n_photons / 2, s.a, s.b, s.g, s.n_slab, s.n_top_slide,
        collimated, mu, t_sample, t_slide, mua_slide, dr_port, dt_port, d_beam, UR1, UT1,
        &ur1_lost, &ut1_lost);
    MC_Radial(n_photons / 2, s.a, s.b, s.g, s.n_slab, s.n_top_slide, diffuse,
        mu, t_sample, t_slide, mua_slide, dr_port, dt_port, d_beam, URU, UTU, &uru_lost,
        &utu_lost);

}
