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

#define MIN_WEIGHT 0.0001
#define N_RADIAL_BINS  1001
#define RADIAL_BIN_SIZE 0.02    /* mm */
#define CLOSE(x, y) (fabs((x) - (y)) < 1e-8)

unsigned long photon_seed = 12345678;

int print_radial_arrays = FALSE;
double R_radial[N_RADIAL_BINS] = { 0 };
double T_radial[N_RADIAL_BINS] = { 0 };

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
    }
    while (x == 0);

    xi = ((double) x) / ((double) kiss_rand_max);

    return xi;
}

/* random number between -1.0 and 1.0 */
static double rand_one_one(void)
{
    return 2.0 * rand_zero_one() - 1.0;
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
    if (temp < 0)
        return 1.0;

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

    // Snell's law
    assert(CLOSE(n_i * sin(acos(wi)), n_t * sin(acos(*w))));

    // Same direction
    assert(((*w) * wi) > 0);

    // direction cosines still have unit magnitude
    assert(CLOSE(sqr(*u) + sqr(*v) + sqr(*w), 1.0));
}

/* Scatter photon and establish new direction */
static void scatter(double g, double *u, double *v, double *w)
{
    double t1, t2, t3, mu, uu, vv, ww;

    do {
        t1 = rand_one_one();
        t2 = rand_one_one();
        t3 = t1 * t1 + t2 * t2;
    }
    while (t3 > 1);

    if (g == 0) {               /* isotropic */
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

/* set photon launch point at the top of the slide */
static void launch_point(double *x, double *y, double *z, double beam_radius, double t_slide)
{
    *x = 0;
    *y = 0;
    *z = -t_slide;

    if (beam_radius > 0) {      /* uniform distribution over a disk */
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

/*
 * Initializes the launch direction of a photon.
 * Right now only oblique collimated angles and diffuse angles are supported
 * If cos_cone_angle=1, the "cone" is a collimated beam and the photon is launched
 * at an angle mu=cos(theta) from the z-axis.
 * If cos_cone_angle=0, the photon's direction is uniformly distributed over the hemisphere,
 * ensuring a diffuse isotropic emission.  The angle mu is irrelevant
 * The cone angle should be added at some point.
 */
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
static void move_in_sample(double mu_t, double *x, double *y, double *z, double u, double v, double w)
{
    double step = -log(rand_zero_one()) / mu_t;
    *x += step * u;
    *y += step * v;
    *z += step * w;
}

static void move_at_boundary(double *u, double *v, double *w, double n_i, double n_t)
{
    // does the photon bounce or get refracted
    double r = fresnel(n_i, n_t, *w);
    if (rand_zero_one() < r) {
        *w = -(*w);
    }
    else {
        refract(n_i, n_t, u, v, w);
    }
}

/*
 * On entry the photon is on one boundary of the slide and the direction is
 * into the slide.  The photon is moving from n1 to n2 to n3.  n2 should always
 * be the index of refraction of the slide (it may equal n1 or n3 in some cases).
 * On exit, x and y will be the location that the photon
 * leaves the slide and the directions u,v,w should be correct.
 * On exit value of z will be on one of the faces of the glass slide.
 *
 */

static void
move_in_slide(double *x, double *y, double *z, double *u, double *v,
    double *w, double *weight, double z_start, double z_end, double b_slide, double n1, double n2, double n3)
{
    double i_x, i_y, i_z, d_bounce, r1, r2, absorbed_loss_per_bounce;

    // case when slide thickness is zero
    if (z_start == z_end) {
        move_at_boundary(u, v, w, n1, n3);
        return;
    }

    // ensure photon is travelling the right direction
    assert((z_end - z_start) * (*w) > 0);

    // reflectivity as we move from n1 to n2
    r1 = fresnel(n1, n2, *w);

    // handle photons that are reflected back at the n1/n2 boundary
    if (rand_zero_one() <= r1) {
        *w = -(*w);
        return;
    }

    // direction cosines within glass
    i_x = *u;
    i_y = *v;
    i_z = *w;
    refract(n1, n2, &i_x, &i_y, &i_z);

    // reflectivity as we move from n2 to n3
    r2 = fresnel(n2, n3, i_z);

    // the distance and absorption for each bounce
    d_bounce = fabs((z_end - z_start) / i_z);
    assert(d_bounce > 0);

    absorbed_loss_per_bounce = exp(-fabs(b_slide / i_z));

    // photon has entered the slide and we bounce it until it exits
    while (1) {
        // move to the n2/n3 boundary
        *x += d_bounce * i_x;
        *y += d_bounce * i_y;
        *z += d_bounce * i_z;
        *weight *= absorbed_loss_per_bounce;
        assert(fabs(*z - z_end) < 1e-8);
        assert((z_end - z_start) * (i_z) > 0);

        // transmitted through the n2/n3 boundary
        if (rand_zero_one() > r2) {
            //fprintf(stderr,"t0 x=%10.5f y=%10.5f z=%10.5f u=%10.5f v=%10.5f w=%10.5f\n", *x, *y, *z, *u, *v, *w);
            //fprintf(stderr, "n_i*sin(theta_i) = %10.5f\n", n1*sin(acos(*w)));
            refract(n1, n3, u, v, w);
            //fprintf(stderr, "n_t*sin(theta_t) = %10.5f\n", n3*sin(acos(*w)));
            //fprintf(stderr,"t1 x=%10.5f y=%10.5f z=%10.5f u=%10.5f v=%10.5f w=%10.5f\n", *x, *y, *z, *u, *v, *w);
            return;
        }

        // bounce back towards the n1/n2 boundary
        i_z *= -1;
        *x += d_bounce * i_x;
        *y += d_bounce * i_y;
        *z += d_bounce * i_z;
        *weight *= absorbed_loss_per_bounce;
        assert(fabs(*z - z_start) < 1e-8);
        assert((z_end - z_start) * i_z < 0);

        // transmitted through (bounced from) the n1/n2 boundary
        if (rand_zero_one() > r1) {
            //fprintf(stderr,"r0 x=%10.5f y=%10.5f z=%10.5f u=%10.5f v=%10.5f w=%10.5f\n", *x, *y, *z, *u, *v, *w);
            *w = -(*w);
            //fprintf(stderr,"r1 x=%10.5f y=%10.5f z=%10.5f u=%10.5f v=%10.5f w=%10.5f\n", *x, *y, *z, *u, *v, *w);
            return;
        }

        // we bounce back towards the n2/n3 boundary
        i_z *= -1;
    }
}


static double milliseconds(clock_t start_time)
{
    double t;
    clock_t finish_time = clock();
    t = 1000 * (double) (finish_time - start_time) / CLOCKS_PER_SEC;
    return t;
}

double add_to_reflectance_array(double x, double y, double z, double w, double weight)
{
    double r = sqrt(x * x + y * y);
    int bin = (int) (r / RADIAL_BIN_SIZE);

    if (bin > N_RADIAL_BINS - 1)
        bin = N_RADIAL_BINS - 1;

    // photon must be going up
    assert(w < 0);

    // weight must be non-zero
    assert(weight != 0);

    R_radial[bin] += weight;
    return r;
}

double add_to_transmittance_array(double x, double y, double z, double w, double weight)
{
    double r = sqrt(x * x + y * y);
    int bin = (int) (r / RADIAL_BIN_SIZE);

    if (bin > N_RADIAL_BINS - 1)
        bin = N_RADIAL_BINS - 1;

    // photon must be going down
    assert(w > 0);

    // photon weight should not be zero
    assert(weight != 0);

    T_radial[bin] += weight;
    return r;
}

/*
 * We have a glass-sample-glass sandwich.
 *    -------------  z = -t_slide
 *        glass
 *    -------------  z = 0
 *        sample
 *    -------------  z = t_sample
 *        glass
 *    -------------  z = t_sample + t_slide
 *
 * Photons start and leave the sample from the top at z=-t_slide
 * the can leave the sample when z= t_sample + t_slide
 */

void
MC_Radial(long photons, double a, double b, double g, double n_sample,
    double n_slide, double cos_cone_angle, double cos_incidence,
    double t_sample, double t_slide, double b_slide, double dr_port,
    double dt_port, double d_beam, double *r_total, double *t_total, double *r_lost, double *t_lost)
{
    double x, y, z, u, v, w, weight;
    long i, total_photons;
    double total_weight, distance_remaining, r_beam, mu_t, r;

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


    /* high resolution clock ... may be zero at start of process */
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

    // reset all accumulators to zero
    *r_lost = 0;
    *t_lost = 0;
    *r_total = 0;
    *t_total = 0;
    total_weight = 0.0;

    // avoid divison by zero or infinity when finding mu_t
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

        // first interaction from air to slide to sample
        //fprintf(stderr," 1 x=%10.5f y=%10.5f z=%10.5f u=%10.5f v=%10.5f w=%10.5f\n", x, y, z, u, v, w);
        move_in_slide(&x, &y, &z, &u, &v, &w, &weight, z, z + t_slide, b_slide, 1.0, n_slide, n_sample);
        //fprintf(stderr," 2 x=%10.5f y=%10.5f z=%10.5f u=%10.5f v=%10.5f w=%10.5f\n", x, y, z, u, v, w);

        // handle the case when the photon is reflected by the slide or surface
        if (w < 0) {
            r = add_to_reflectance_array(x, y, z, w, weight);
            *r_total += weight;
            if (r > r_port_radius)
                *r_lost += weight;
            weight = 0;
            continue;
        }

        // sanity checks for the photon
        assert(w > 0);
        assert(cos_critical < w);
        assert(CLOSE(z, 0));
        assert(weight == 1);

        // ensure direction cosines are valid
        assert(CLOSE(sqr(u) + sqr(v) + sqr(w), 1));

        while (weight > 0) {

            // move to next scattering or absorption event
            move_in_sample(mu_t, &x, &y, &z, u, v, w);

            assert(weight <= 1);

            // handle photons that have left the sample
            while (z < 0 || z > t_sample) {

                // distance of the photon past the top or bottom of sample
                distance_remaining = z / w;
                if (z > t_sample)
                    distance_remaining -= t_sample / w;

                // update the x,y coordinates also
                x -= u * distance_remaining;
                y -= v * distance_remaining;
                z -= w * distance_remaining;

                assert(CLOSE(z, 0) || CLOSE(z, t_sample));

                // move photon for n_sample -> n_slide -> 1 */
                //fprintf(stderr," 3 x=%10.5f y=%10.5f z=%10.5f u=%10.5f v=%10.5f w=%10.5f\n", x, y, z, u, v, w);
                if (w > 0)
                    move_in_slide(&x, &y, &z, &u, &v, &w, &weight, z, z + t_slide, b_slide, n_sample, n_slide, 1.0);
                else
                    move_in_slide(&x, &y, &z, &u, &v, &w, &weight, z, z - t_slide, b_slide, n_sample, n_slide, 1.0);
                //fprintf(stderr," 4 x=%10.5f y=%10.5f z=%10.5f u=%10.5f v=%10.5f w=%10.5f\n", x, y, z, u, v, w);


                // photon exiting from the top
                if (CLOSE(z, -t_slide) && w < 0) {
                    r = add_to_reflectance_array(x, y, z, w, weight);
                    *r_total += weight;
                    if (r > r_port_radius)
                        *r_lost += weight;
                    weight = 0;
                    break;
                }

                // photon exiting from the bottom
                if (CLOSE(z, t_sample + t_slide) && w > 0) {
                    r = add_to_transmittance_array(x, y, z, w, weight);
                    *t_total += weight;
                    if (r > t_port_radius)
                        *t_lost += weight;
                    weight = 0;
                    break;
                }

                // move the distance_remaining amount back into sample
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

                // ensure direction cosines are valid
                assert(CLOSE(u * u + v * v + w * w, 1.0));
            }
        }

        if (total_time && total_time < milliseconds(start_time))
            break;
    }

    double total = total_weight - residual_weight;
    *r_lost /= total;
    *t_lost /= total;
    *r_total /= total;
    *t_total /= total;

    if (print_radial_arrays) {
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
}

/* This assumes that the sample port diameter for the reflection and transmission
   measurements are the same. it also assumes top and bottom slides are the same */
void
MC_Lost(struct measure_type m, struct invert_type r, long n_photons,
    double *ur1, double *ut1, double *uru, double *utu,
    double *ur1_lost, double *ut1_lost, double *uru_lost, double *utu_lost)
{
    double n_sample = m.slab_index;
    double n_slide = m.slab_top_slide_index;
    double t_sample = m.slab_thickness;
    double t_slide = m.slab_top_slide_thickness;
    double b_slide = m.slab_top_slide_b;

    /* retrieve sample port diameters for each sphere */
    double dr_port = sqrt(m.as_r) * 2 * m.d_sphere_r;
    double dt_port = sqrt(m.as_t) * 2 * m.d_sphere_t;
    double d_beam = m.d_beam;
    double mu = m.slab_cos_angle;

    /* no slide if thickness is zero or the index is 1.0 */
    if (t_slide == 0.0)
        n_slide = 1.0;
    if (n_slide == 1.0)
        t_slide = 0.0;

//    set_photon_seed(12345);
    MC_Radial(n_photons / 2, r.a, r.b, r.g, n_sample, n_slide,
        COLLIMATED, mu, t_sample, t_slide, b_slide, dr_port, dt_port, d_beam, ur1, ut1, ur1_lost, ut1_lost);

    *uru_lost = 0;
    *utu_lost = 0;

    /* uru_lost and utu_lost do not apply when dual beam is used */
    if (m.method == SUBSTITUTION)
        MC_Radial(n_photons / 2, r.a, r.b, r.g, n_sample, n_slide,
            DIFFUSE, mu, t_sample, t_slide, b_slide, dr_port, dt_port, d_beam, uru, utu, uru_lost, utu_lost);

    if (*ur1_lost < 0 || *ut1_lost < 0 || *uru_lost < 0 || *utu_lost < 0) {
        //fprintf(stderr, "One or more of MC lost light calculations is negative!\n");
        exit(EXIT_FAILURE);
    }
}

void
MC_RT(struct AD_slab_type s, long n_photons, double t_sample, double t_slide,
    double *UR1, double *UT1, double *URU, double *UTU)
{
    double ur1_lost, ut1_lost, uru_lost, utu_lost;
    double dr_port = 1000;      /* collect all light with in 1 meter port! */
    double dt_port = 1000;
    double d_beam = 0.0;
    double mu = s.cos_angle;

    set_photon_seed(12345);

    MC_Radial(n_photons / 2, s.a, s.b, s.g, s.n_slab, s.n_top_slide,
        COLLIMATED, mu, t_sample, t_slide, s.b_top_slide, dr_port, dt_port, d_beam, UR1, UT1, &ur1_lost, &ut1_lost);

    MC_Radial(n_photons / 2, s.a, s.b, s.g, s.n_slab, s.n_top_slide,
        DIFFUSE, mu, t_sample, t_slide, s.b_top_slide, dr_port, dt_port, d_beam, URU, UTU, &uru_lost, &utu_lost);
}

void MC_Print_RT_Arrays(int status)
{
    print_radial_arrays = status;
}
