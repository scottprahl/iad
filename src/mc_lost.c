#include <assert.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "ad_globl.h"
#include "iad_type.h"

#define MIN_WEIGHT 0.001
#define PI         3.1415926535

/* The KISS generator, (Keep It Simple Stupid), is
   designed to combine the two multiply-with-carry
   generators in MWC with the 3-shift register SHR3 and
   the congruential generator CONG, using addition and
   exclusive-or. Period about 2^123.  George Marsaglia*/

unsigned long kiss_rand_max = ULONG_MAX;
unsigned long kiss_x=123456789;
unsigned long kiss_y=362436000;
unsigned long kiss_z=521288629;
unsigned long kiss_c=7654321;
 
unsigned long kiss_rand()
{
    unsigned long long t, a=698769069ULL;
    
    kiss_x  = 69069*kiss_x+12345;
    kiss_y ^= ( kiss_y << 13 ); 
    kiss_y ^= ( kiss_y >> 17 ); 
    kiss_y ^= ( kiss_y <<  5 );
    t       = a*kiss_z+kiss_c; 
    kiss_c  = (t>>32);
    kiss_z  = t;
    return  kiss_x+kiss_y+kiss_z;
} 

/* use Knuth's algorithms to set initial values */
void kiss_rand_seed(unsigned long seed)
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
unsigned long photon_seed = 12345678;

static inline void set_photon_seed (unsigned long new_seed)
{
    photon_seed = new_seed;
}

static inline void next_photon_seed (void)
{
    photon_seed = (1812433253UL * photon_seed) & 0xffffffffUL;
    kiss_rand_seed(photon_seed);
    kiss_rand();
    kiss_rand();
    kiss_rand();
}

/* more or less portable random number initialization 
   call this with a non-zero seed to start the sequence from a fixed value
   otherwise pass seed=0 so that the program starts with a random number
   
   auto-initialization fails if called faster than once per second.  Since
   this may happen for successive MC_Radial() calls, the random number sequence 
   IS NOT restarted.  The idea is that one seed determines the entire subsequent
   sequence of events.
*/
   
void my_randomize(unsigned long seed)
{
    static int initialized = 0;
    
    if (!initialized) {
        initialized = 1;
        
        if (seed==0) /* time in seconds from start of Unix epoch */
            seed = (unsigned long) time(NULL);
        
        photon_seed = seed;
    }
}


/* random number between 0.0 and 1.0 ... EXCLUDING 0.0 */
static double rand_zero_one(void)
{
    unsigned long x;
    double xi;
    
    do {
        x = kiss_rand();
    } while (x==0);
    
    xi = ((double)x)/((double) kiss_rand_max); 
    
    return xi;
}

/* random number between -1.0 and 1.0 */
static double rand_one_one(void)
{
    return 2.0*rand_zero_one() - 1.0;
}


static double sqr(double x)
{
    return x*x;
}

/* Fresnel Reflectance 
   incident    medium has index n_i and angle nu_i on
   transmitted medium has index n_t and angle nu_t
*/
static double fresnel(double n_i,double n_t,double nu_i)
{
    double nu_t,ratio,temp,temp1;
    
    /* matched boundaries? */
    if (n_i == n_t) return 0.0;
    
    /* normal illumination? */
    nu_i = fabs(nu_i);
    if (nu_i == 1.0) return sqr((n_i-n_t)/(n_i+n_t));
    
    /* total internal reflection? */
    ratio = n_i/n_t;
    temp = 1.0-ratio*ratio*(1.0-nu_i*nu_i);
    if (temp < 0) return 1.0;
    
    /* Now calculate Fresnel reflection */
    nu_t = sqrt(temp);  
    temp = ratio*nu_t;
    temp1= (nu_i-temp)/(nu_i+temp);
    temp = ratio*nu_i;
    temp = (nu_t-temp)/(nu_t+temp);
    return (temp1*temp1+temp*temp)/2.0 ;
}

double cos_critical_angle(double n_i, double n_t)
{
    if (n_t >= n_i)
        return 0.0;
    else 
        return sqrt(1.0 - sqr(n_t/n_i));
}

/* refract a photon into or out of the medium across a z-plane*/
void refract(double n_i, double n_t, double *u, double *v, double *w)
{
    double nu, c;
    
#ifndef NDEBUG
    double wold;
    wold = *w;
#endif

    if (n_i == n_t) return;
    c    = n_i/n_t;
    nu   = *w * c; 
    
    assert(n_i<n_t || fabs(*w) > cos_critical_angle(n_i,n_t));
    *u *= c;
    *v *= c;
    if (*w<0)
        *w  = -sqrt(1.0 - c * c  + nu * nu);
    else
        *w  =  sqrt(1.0 - c * c  + nu * nu);
    
    /* Snell's law valid */
    assert(fabs(n_i*sin(acos(wold))-n_t*sin(acos(*w)))<1e-8);
    
    /* same direction */
    assert(*w*wold>0);
    
    /* still a unit vector */
    assert(fabs(sqr(*u)+sqr(*v)+sqr(*w)-1.0)<1e-8);
}

/* Scatter photon and establish new direction */
void scatter(double g, double *u, double *v, double *w)
{
    double t1, t2, t3, mu, uu, vv, ww;

    do { 
        t1 = rand_one_one(); 
        t2 = rand_one_one(); 
        t3 = t1*t1+t2*t2;
    } while (t3>1);

    if (g==0) {  /* isotropic */
        if (t3 == 0) {
            *u = 1;
            *v = 0;
            *w = 0;
        } else {
            *u = 2.0 * t3 - 1.0;
            t3 = sqrt((1.0 - (*u) * (*u))/t3);
            *v = t1 * t3;
            *w = t2 * t3;
        }
        return;
    } 

    mu = (1-g*g)/(1-g+2.0*g*rand_zero_one());
    mu = (1 + g*g-mu*mu)/2.0/g;

    uu = *u;
    vv = *v;
    ww = *w;

    if ( fabs(ww) < 0.9 ) {  
        *u = mu * uu + sqrt((1-mu*mu)/(1-ww*ww)/t3) * (t1*uu*ww-t2*vv);
        *v = mu * vv + sqrt((1-mu*mu)/(1-ww*ww)/t3) * (t1*vv*ww+t2*uu);
        *w = mu * ww - sqrt((1-mu*mu)*(1-ww*ww)/t3) * t1;
    } else {
        *u = mu * uu + sqrt((1-mu*mu)/(1-vv*vv)/t3) * (t1*uu*vv + t2*ww);
        *v = mu * vv - sqrt((1-mu*mu)*(1-vv*vv)/t3) * t1;
        *w = mu * ww + sqrt((1-mu*mu)/(1-vv*vv)/t3) * (t1*vv*ww - t2*uu);
    }
}

/* set photon launch point*/
static void launch_point(double *x, double *y, double *z, double beam_radius)
{
    *x = 0; 
    *y = 0; 
    *z = 0; 

    if (beam_radius > 0) {  /* uniform distribution over a disk */
        double a, b;
        do {
            a  = rand_one_one(); 
            b  = rand_one_one();
        } while (a*a+b*b>1);
        
        *x = a * beam_radius;
        *y = b * beam_radius; 
    }
}

static void launch_direction(double *u, double *v, double *w, int collimated, double mu)
{
    if (!collimated) {
        do
            scatter(0.0,u,v,w);
        while (*w == 0);
        *w = fabs(*w);
    } else {
        *u = sqrt(1-mu*mu);
        *v = 0;
        *w = mu;
    }
    
}

static void roulette(double *weight, double *residual_weight)
{
    if (*weight>MIN_WEIGHT || *weight==0.0) return;

    *residual_weight += *weight;
    if (rand_zero_one() < 0.1)
        *weight *=10;
    else
        *weight = 0;
    *residual_weight -= *weight;
}

static void move_in_sample(double mu_t, double *x, double *y, double *z, double u, double v, double w)
{
    double step = -log(rand_zero_one()) / mu_t;

    *x += step * u;
    *y += step * v;
    *z += step * w;
}

/*On entry p is on one boundary of the slide and the direction is
  into the slide on entry.  On exit, x and y will be the location
  at exit, but z is not touched.  The direction w should be correct.
  The value of z does not need to be updated because if the ray exits
  then it does not matter.  If the ray doesn't exit, then it is reflected
  back into the medium and should have the same value as when it entered
 */
static void move_in_slide(double *x, double *y, double u, double v, double *w, 
                          double d_slide, double n1, double n2, double n3)
{
    double i_x, i_y, i_z, d_bounce, r1, r2;
    
    r1 = fresnel(n1, n2, *w);

    /* Does the ray makes it into the slide?  The test is <=
       to ensure that when r1==1, the ray is *always* reflected */
    if (rand_zero_one() <= r1)  { 
        *w = -(*w); 
        return;
    }
    
    /* direction cosines in the slide */
    i_x = u;
    i_y = v;
    i_z = *w;
    refract(n1, n2, &i_x, &i_y, &i_z);
    
    r2 = fresnel(n2, n3, i_z);
    d_bounce = fabs(d_slide/i_z);

    /* photon is in slide and bounces in the slide until it exits */
    while (1) {     
        *x += d_bounce * i_x;
        *y += d_bounce * i_y;
        
        if (rand_zero_one() > r2) 
            return;
        
        *x += d_bounce * i_x;
        *y += d_bounce * i_y;

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
    t = 1000*(double)(finish_time-start_time)/CLOCKS_PER_SEC;
    return t;    
}

static void MC_Radial(long photons, double a, double b, double g, double n_sample, 
                      double n_slide, int collimated, double cos_incidence,
                      double d_sample, double d_slide, double d_port, double d_beam,
                      double *r_total, double *t_total, double *r_lost, double *t_lost)
{

    double x, y, z, u, v, w, weight;
    long i,total_photons;
    double residual_weight = 0.0;
    double r_beam;
    double mu_t = b / d_sample;
    double r_port_squared = d_port * d_port / 4.0;
    double total_weight,extra;
    double total_time=0;
    clock_t start_time=0;
#ifndef NDEBUG
	double ww;
#endif        
    /* high resolution clock ... may be zero at start of process */
    start_time = clock();   
        
    if (collimated) 
        r_beam = d_beam / 2.0;
    else
        r_beam = d_port / 2.0;

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
        
    for (i = 1; i <= total_photons; i++) {
        next_photon_seed();
        launch_point(&x,&y,&z,r_beam);
        launch_direction(&u,&v,&w,collimated,cos_incidence);
        weight = w;
        total_weight += weight;
        
        move_in_slide(&x, &y, u, v, &w, d_slide, 1.0, n_slide, n_sample);
        
        if ( w < 0) {
            *r_total += weight;
            if (x*x + y*y > r_port_squared) *r_lost += weight;
            continue;
        } 
        
#ifndef NDEBUG
        ww = w;
#endif
        refract(1.0, n_sample, &u, &v, &w);
        assert(fabs(sin(acos(ww)) - n_sample*sin(acos(w)))<1e-8 || ww == -w);
        assert(fabs(u*u+v*v+w*w-1.0) < 1e-8);
                     
        while (weight>0) {
            
            move_in_sample(mu_t,&x,&y,&z,u,v,w);
            
            while (z < 0 || z > d_sample) {             
                /* move back to top or bottom of sample */
                if (z < 0) {
                    extra = z/w;
                    z  = 0.0;
                } else {
                    extra = (z-d_sample)/w;
                    z  = d_sample;
                } 
                x -= u * extra;
                y -= v * extra;
                    
                /* w changes sign upon reflection back into sample */
                move_in_slide(&x, &y, u, v, &w, d_slide, n_sample, n_slide, 1.0);

                /* use direction w to determine if photon is transmitted out of sample */
                if (z == 0) { 
                    
                    if (w < 0) {  /* still moving up */
                        *r_total += weight;
                        if (x*x + y*y > r_port_squared) *r_lost += weight;
                        weight = 0;
                        break;
                    } 
                    
                } else { 
                
                    if (w > 0) { /* still moving down */
                        *t_total += weight;
                        if (x*x + y*y > r_port_squared) *t_lost += weight;
                        weight = 0;
                        break;
                    } 
                }
                
                /* move the remaining amount back into sample */
                x += u * extra;
                y += v * extra;
                z += w * extra;
            }
            
            weight *= a;
            roulette(&weight, &residual_weight);
            scatter(g, &u, &v, &w);
        }
        
        if (total_time && total_time < milliseconds(start_time))
            break;
        
    }
    
    *r_lost  /= (total_weight + residual_weight);
    *t_lost  /= (total_weight + residual_weight);
    *r_total /= (total_weight + residual_weight);
    *t_total /= (total_weight + residual_weight);
}

void MC_Lost(struct measure_type m, struct invert_type r,  long n_photons,
             double *ur1,      double *ut1,      double *uru,      double *utu, 
             double *ur1_lost, double *ut1_lost, double *uru_lost, double *utu_lost)
{
    int collimated = 1;
    int diffuse = 0;

    double n_sample = m.slab_index;
    double n_slide  = m.slab_top_slide_index;
    double d_sample = m.slab_thickness;
    double d_slide  = m.slab_top_slide_thickness;
    double d_port   = sqrt(m.as_r)*2*m.d_sphere_r;
    double d_beam   = m.d_beam;
    double mu       = m.slab_cos_angle;
    
    /* no slide if thickness is zero or the index is 1.0 */
    if (d_slide == 0.0) n_slide = 1.0;
    if (n_slide == 1.0) d_slide = 0.0;
    
    set_photon_seed(12345);
    MC_Radial(n_photons/2, r.a, r.b, r.g, n_sample, n_slide, collimated, mu, d_sample, d_slide, d_port, d_beam, ur1, ut1, ur1_lost, ut1_lost);
    
    *uru_lost = 0;
    *utu_lost = 0;
    if (m.method==SUBSTITUTION) 
    	MC_Radial(n_photons/2, r.a, r.b, r.g, n_sample, n_slide, diffuse,    mu, d_sample, d_slide, d_port, d_beam, uru, utu, uru_lost, utu_lost);
    
    if (*ur1_lost<0) *ur1_lost = 0;
    if (*ut1_lost<0) *ut1_lost = 0;
    if (*uru_lost<0) *uru_lost = 0;
    if (*utu_lost<0) *utu_lost = 0;
}

void MC_RT(struct AD_slab_type s, long n_photons, double *UR1, double *UT1, double *URU, double *UTU)
{
    double ur1_lost, ut1_lost, uru_lost, utu_lost;
    int collimated = 1;
    int diffuse = 0;
    double d_sample = 1.0;
    double d_port   = 1000;
    double d_beam   = 0.0;
    double d_slide  = 0.0;
    double mu       = s.cos_angle;
    
    set_photon_seed(12345);
    MC_Radial(n_photons/2, s.a, s.b, s.g, s.n_slab, s.n_top_slide, collimated, mu, d_sample, d_slide, d_port, d_beam, UR1, UT1, &ur1_lost, &ut1_lost);
    MC_Radial(n_photons/2, s.a, s.b, s.g, s.n_slab, s.n_top_slide, diffuse,    mu, d_sample, d_slide, d_port, d_beam, URU, UTU, &uru_lost, &utu_lost);
    
}

/*
static void print_usage(void)
{
    fprintf(stderr, "lost1 \n\n");
    fprintf(stderr, "lost finds the reflection and transmission from optical properties\n\n");
    fprintf(stderr, "Usage:  lost [options] input\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h               display help\n");
    fprintf(stderr, "  -a #             albedo (0-1)\n");
    fprintf(stderr, "  -b #             optical thickness (>0)\n");
    fprintf(stderr, "  -g #             scattering anisotropy (-1 to 1)\n");
    fprintf(stderr, "  -n #             specify index of refraction of slab\n");
    fprintf(stderr, "  -s #             specify index of refraction of slide\n");
    fprintf(stderr, "  -N #             number of photons\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "  lost -m data                     data.rt in machine readable format\n");
    fprintf(stderr, "  ad data -o out.txt             out.txt is the \n");
    fprintf(stderr, "  lost -a 0.3                    a=0.3, b=2.0, g=0.0, n=1.0\n");
    exit(0);
}

int main(int argc, char** argv) 
{
    char c;
    double ur1_lost, ut1_lost, uru_lost, utu_lost;
    double ur1, ut1, uru, utu;
    int collimated = 1;
    int diffuse = 0;
    double d_sample = 1.0;
    double d_port   = 10.0;
    double d_beam   = 5.0;
    double d_slide  = 1.0;
    double a = 0.5;
    double b = 2;
    double g = 0.0;
    double n_slab = 1.0;
    double n_top_slide = 1.0;
    int n_photons = 1024*1024;

    while ((c = getopt(argc, argv, "h?va:b:g:n:s:q:")) != EOF) {
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
                n_top_slide = strtod(optarg, NULL);
                break;
                
            default:
            case 'h':
            case '?':
                print_usage();
                break;
        }
    }
    
    MC_Radial(n_photons, a, b, g, n_slab, n_top_slide, collimated, d_sample, d_slide, d_port, d_beam, &ur1, &ut1, &ur1_lost, &ut1_lost);
    MC_Radial(n_photons, a, b, g, n_slab, n_top_slide, diffuse,    d_sample, d_slide, d_port, d_beam, &uru, &utu, &uru_lost, &utu_lost);

    printf("   UR1    \t   UT1    \t   URU    \t   UTU\n");
    printf("%9.5f \t%9.5f \t%9.5f \t%9.5f\n", ur1, ut1, uru, utu);
    printf("%9.5f \t%9.5f \t%9.5f \t%9.5f\n", ur1_lost, ut1_lost, uru_lost, utu_lost);
    return 0;
}

*/
