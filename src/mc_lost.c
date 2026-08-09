/*1:*/
#line 14 "mc_lost.w"

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

/*4:*/
#line 74 "mc_lost.w"

#define MIN_WEIGHT 0.0001
#define N_RADIAL_BINS  1001
#define RADIAL_BIN_SIZE 0.02
#define CLOSE(x, y) (fabs((x) - (y)) < 1e-8)
#define DIFFUSE_SEED_OFFSET 0x9e3779b9UL

/*:4*/
#line 27 "mc_lost.w"

/*6:*/
#line 85 "mc_lost.w"

unsigned long photon_seed= 12345678;
unsigned long lost_base_seed= 12345678;

int print_radial_arrays= FALSE;
double R_radial[N_RADIAL_BINS]= {0};
double T_radial[N_RADIAL_BINS]= {0};

unsigned long kiss_rand_max= ULONG_MAX;
unsigned long kiss_x= 123456789;
unsigned long kiss_y= 362436000;
unsigned long kiss_z= 521288629;
unsigned long kiss_c= 7654321;

/*:6*/
#line 28 "mc_lost.w"


/*8:*/
#line 105 "mc_lost.w"

static unsigned long kiss_rand(void)
{
unsigned long long t,a= 698769069ULL;

kiss_x= 69069*kiss_x+12345;
kiss_y^= (kiss_y<<13);
kiss_y^= (kiss_y>>17);
kiss_y^= (kiss_y<<5);
t= a*kiss_z+kiss_c;
kiss_c= (t>>32);
kiss_z= t;
return kiss_x+kiss_y+kiss_z;
}

/*:8*/
#line 30 "mc_lost.w"

/*10:*/
#line 122 "mc_lost.w"

static void kiss_rand_seed(unsigned long seed)
{
static const unsigned long K= 1812433253UL;
kiss_c= K*(seed^(seed>>30))+1;
kiss_x= K*(kiss_c^(kiss_c>>30))+2;
kiss_y= K*(kiss_x^(kiss_x>>30))+3;
kiss_z= K*(kiss_y^(kiss_y>>30))+5;
}

/*:10*/
#line 31 "mc_lost.w"

/*11:*/
#line 132 "mc_lost.w"

static inline void set_photon_seed(unsigned long new_seed)
{
photon_seed= new_seed;
}

/*:11*/
#line 32 "mc_lost.w"

/*14:*/
#line 145 "mc_lost.w"

/*13:*/
#line 142 "mc_lost.w"

void MC_Set_Seed(unsigned long seed)

/*:13*/
#line 146 "mc_lost.w"

{
if(seed==0)
lost_base_seed= (unsigned long)time(NULL);
else
lost_base_seed= seed;

set_photon_seed(lost_base_seed);
}

/*:14*/
#line 33 "mc_lost.w"

/*16:*/
#line 160 "mc_lost.w"

static inline void next_photon_seed(void)
{
photon_seed= (1812433253UL*photon_seed)&0xffffffffUL;
kiss_rand_seed(photon_seed);
kiss_rand();
kiss_rand();
kiss_rand();
}

/*:16*/
#line 34 "mc_lost.w"

/*18:*/
#line 173 "mc_lost.w"

static double rand_zero_one(void)
{
unsigned long x;
double xi;

do{
x= kiss_rand();
}
while(x==0);

xi= ((double)x)/((double)kiss_rand_max);

return xi;
}

/*:18*/
#line 35 "mc_lost.w"

/*19:*/
#line 189 "mc_lost.w"

static double rand_one_one(void)
{
return 2.0*rand_zero_one()-1.0;
}

/*:19*/
#line 36 "mc_lost.w"

/*21:*/
#line 200 "mc_lost.w"

static double fresnel(double n_i,double n_t,double nu_i)
{
double nu_t,ratio,temp,temp1;

if(n_i==n_t)
return 0.0;

nu_i= fabs(nu_i);
if(nu_i==1.0)
return sqr((n_i-n_t)/(n_i+n_t));

ratio= n_i/n_t;
temp= 1.0-ratio*ratio*(1.0-nu_i*nu_i);
if(temp<0)
return 1.0;

nu_t= sqrt(temp);
temp= ratio*nu_t;
temp1= (nu_i-temp)/(nu_i+temp);
temp= ratio*nu_i;
temp= (nu_t-temp)/(nu_t+temp);
return(temp1*temp1+temp*temp)/2.0;
}

/*:21*/
#line 37 "mc_lost.w"

/*23:*/
#line 228 "mc_lost.w"

static void refract(double n_i,double n_t,double*u,double*v,double*w)
{
double nu,c;
#ifndef NDEBUG
double wi= *w;
#endif
#line 235 "mc_lost.w"

if(n_i==n_t)
return;

c= n_i/n_t;
nu= (*w)*c;

*u*= c;
*v*= c;
if(*w<0)
*w= -sqrt(1.0-sqr(c)+sqr(nu));
else
*w= sqrt(1.0-sqr(c)+sqr(nu));

assert(CLOSE(n_i*sin(acos(wi)),n_t*sin(acos(*w))));
assert(((*w)*wi)> 0);
assert(CLOSE(sqr(*u)+sqr(*v)+sqr(*w),1.0));
}

/*:23*/
#line 38 "mc_lost.w"

/*25:*/
#line 257 "mc_lost.w"

static void scatter(double g,double*u,double*v,double*w)
{
double t1,t2,t3,mu,uu,vv,ww;

do{
t1= rand_one_one();
t2= rand_one_one();
t3= t1*t1+t2*t2;
}
while(t3> 1);

if(g==0){
*u= 2.0*t1*sqrt(1.0-t3);
*v= 2.0*t2*sqrt(1.0-t3);
*w= 1.0-2.0*t3;
return;
}

mu= (1-g*g)/(1-g+2.0*g*rand_zero_one());
mu= (1+g*g-mu*mu)/2.0/g;

uu= *u;
vv= *v;
ww= *w;

if(fabs(ww)<0.9){
*u= mu*uu+sqrt((1-mu*mu)/(1-ww*ww)/t3)*(t1*uu*ww-t2*vv);
*v= mu*vv+sqrt((1-mu*mu)/(1-ww*ww)/t3)*(t1*vv*ww+t2*uu);
*w= mu*ww-sqrt((1-mu*mu)*(1-ww*ww)/t3)*t1;
}
else{
*u= mu*uu+sqrt((1-mu*mu)/(1-vv*vv)/t3)*(t1*uu*vv+t2*ww);
*v= mu*vv-sqrt((1-mu*mu)*(1-vv*vv)/t3)*t1;
*w= mu*ww+sqrt((1-mu*mu)/(1-vv*vv)/t3)*(t1*vv*ww-t2*uu);
}
}

/*:25*/
#line 39 "mc_lost.w"

/*27:*/
#line 300 "mc_lost.w"

static void launch_point(double*x,double*y,double*z,double beam_radius,double t_slide)
{
*x= 0;
*y= 0;
*z= -t_slide;

if(beam_radius> 0){
double a,b;
do{
a= rand_one_one();
b= rand_one_one();
}
while(a*a+b*b> 1);

*x= a*beam_radius;
*y= b*beam_radius;
}
}

/*:27*/
#line 40 "mc_lost.w"

/*29:*/
#line 323 "mc_lost.w"

static void launch_direction(double*u,double*v,double*w,double cos_cone_angle,double mu)
{
double phi;
if(cos_cone_angle==COLLIMATED){
*u= sqrt(1-mu*mu);
*v= 0;
*w= mu;
}
else{
*w= sqrt(rand_zero_one());
phi= 2.0*M_PI*rand_zero_one();
*u= cos(phi)*sqrt(1-sqr(*w));
*v= sin(phi)*sqrt(1-sqr(*w));
}
}

/*:29*/
#line 41 "mc_lost.w"

/*31:*/
#line 344 "mc_lost.w"

static void roulette(double*weight,double*residual_weight)
{
if(*weight> MIN_WEIGHT||*weight==0.0)
return;

*residual_weight+= *weight;
if(rand_zero_one()<0.1)
*weight*= 10;
else
*weight= 0;
*residual_weight-= *weight;
}

/*:31*/
#line 42 "mc_lost.w"

/*32:*/
#line 358 "mc_lost.w"

static void move_in_sample(double mu_t,double*x,double*y,double*z,double u,double v,double w)
{
double step= -log(rand_zero_one())/mu_t;
*x+= step*u;
*y+= step*v;
*z+= step*w;
}

/*:32*/
#line 43 "mc_lost.w"

/*33:*/
#line 367 "mc_lost.w"

static void move_at_boundary(double*u,double*v,double*w,double n_i,double n_t)
{
double r= fresnel(n_i,n_t,*w);
if(rand_zero_one()<r){
*w= -(*w);
}
else{
refract(n_i,n_t,u,v,w);
}
}

/*:33*/
#line 44 "mc_lost.w"

/*35:*/
#line 383 "mc_lost.w"

static void
move_in_slide(double*x,double*y,double*z,double*u,double*v,
double*w,double*weight,double z_start,double z_end,double b_slide,double n1,double n2,double n3)
{
double i_x,i_y,i_z,d_bounce,r1,r2,absorbed_loss_per_bounce;

if(z_start==z_end){
move_at_boundary(u,v,w,n1,n3);
return;
}

assert((z_end-z_start)*(*w)> 0);

r1= fresnel(n1,n2,*w);
if(rand_zero_one()<=r1){
*w= -(*w);
return;
}

i_x= *u;
i_y= *v;
i_z= *w;
refract(n1,n2,&i_x,&i_y,&i_z);

r2= fresnel(n2,n3,i_z);
d_bounce= fabs((z_end-z_start)/i_z);
assert(d_bounce> 0);

absorbed_loss_per_bounce= exp(-fabs(b_slide/i_z));

while(1){
*x+= d_bounce*i_x;
*y+= d_bounce*i_y;
*z+= d_bounce*i_z;
*weight*= absorbed_loss_per_bounce;
assert(fabs(*z-z_end)<1e-8);
assert((z_end-z_start)*(i_z)> 0);

if(rand_zero_one()> r2){
refract(n1,n3,u,v,w);
return;
}

i_z*= -1;
*x+= d_bounce*i_x;
*y+= d_bounce*i_y;
*z+= d_bounce*i_z;
*weight*= absorbed_loss_per_bounce;
assert(fabs(*z-z_start)<1e-8);
assert((z_end-z_start)*i_z<0);

if(rand_zero_one()> r1){
*w= -(*w);
return;
}

i_z*= -1;
}
}

/*:35*/
#line 45 "mc_lost.w"

/*36:*/
#line 444 "mc_lost.w"

static double milliseconds(clock_t start_time)
{
double t;
clock_t finish_time= clock();
t= 1000*(double)(finish_time-start_time)/CLOCKS_PER_SEC;
return t;
}

/*:36*/
#line 46 "mc_lost.w"

/*38:*/
#line 458 "mc_lost.w"

double add_to_reflectance_array(double x,double y,double z,double w,double weight)
{
double r= sqrt(x*x+y*y);
int bin= (int)(r/RADIAL_BIN_SIZE);

if(bin> N_RADIAL_BINS-1)
bin= N_RADIAL_BINS-1;

assert(w<0);
assert(weight!=0);

R_radial[bin]+= weight;
return r;
}

/*:38*/
#line 47 "mc_lost.w"

/*39:*/
#line 474 "mc_lost.w"

double add_to_transmittance_array(double x,double y,double z,double w,double weight)
{
double r= sqrt(x*x+y*y);
int bin= (int)(r/RADIAL_BIN_SIZE);

if(bin> N_RADIAL_BINS-1)
bin= N_RADIAL_BINS-1;

assert(w> 0);
assert(weight!=0);

T_radial[bin]+= weight;
return r;
}

/*:39*/
#line 48 "mc_lost.w"

/*42:*/
#line 515 "mc_lost.w"

/*41:*/
#line 507 "mc_lost.w"

void MC_Radial(long photons,double a,double b,double g,double n_sample,
double n_top_slide,double n_bottom_slide,
double cos_cone_angle,double cos_incidence,
double t_sample,double t_top_slide,double t_bottom_slide,
double b_top_slide,double b_bottom_slide,
double dr_port,double dt_port,double d_beam,double*r_total,double*t_total,double*r_lost,double*t_lost)

/*:41*/
#line 516 "mc_lost.w"

{
double x,y,z,u,v,w,weight;
long i,total_photons;
double total_weight,distance_remaining,r_beam,mu_t,r,total;

double r_port_radius= dr_port/2.0;
double t_port_radius= dt_port/2.0;
double residual_weight= 0.0;
double total_time= 0;
clock_t start_time= 0;
double absorbed= 0;

#ifndef NDEBUG
double cos_critical= Cos_Critical_Angle(n_sample,1.0);
fprintf(stderr,"illumination  = %s\n",
(cos_cone_angle==COLLIMATED)?"collimated":"diffuse");
fprintf(stderr,"cos_incidence = %10.5f\n",cos_incidence);
fprintf(stderr,"cos_critical = %10.5f\n",cos_critical);
fprintf(stderr,"d_beam   = %10.5f\n",d_beam);
fprintf(stderr,"t_sample = %10.5f\n",t_sample);
fprintf(stderr,"t_top    = %10.5f\n",t_top_slide);
fprintf(stderr,"t_bottom = %10.5f\n",t_bottom_slide);
fprintf(stderr,"n_sample = %10.5f\n",n_sample);
fprintf(stderr,"n_top    = %10.5f\n",n_top_slide);
fprintf(stderr,"n_bottom = %10.5f\n",n_bottom_slide);
#endif
#line 543 "mc_lost.w"

start_time= clock();

if(cos_cone_angle==COLLIMATED)
r_beam= d_beam/2.0;
else
r_beam= dr_port/2.0;

if(photons>=0)
total_photons= photons;
else{
total_photons= 1000000;
total_time= labs(photons);
}

*r_lost= 0;
*t_lost= 0;
*r_total= 0;
*t_total= 0;
total_weight= 0.0;

for(i= 0;i<N_RADIAL_BINS;i++){
R_radial[i]= 0;
T_radial[i]= 0;
}

if(b<1e-5)
b= 1e-5;
if(b> 1000)
b= 1000;
mu_t= b/t_sample;

for(i= 1;i<=total_photons;i++){
next_photon_seed();
launch_point(&x,&y,&z,r_beam,t_top_slide);
launch_direction(&u,&v,&w,cos_cone_angle,cos_incidence);
assert(w> 0&&z==-t_top_slide);

weight= 1;
total_weight+= weight;

move_in_slide(&x,&y,&z,&u,&v,&w,&weight,z,z+t_top_slide,
b_top_slide,1.0,n_top_slide,n_sample);

if(w<0){
r= add_to_reflectance_array(x,y,z,w,weight);
*r_total+= weight;
if(r> r_port_radius)
*r_lost+= weight;
weight= 0;
continue;
}

assert(w> 0);
assert(cos_critical<w);
assert(CLOSE(z,0));
assert(weight> 0&&weight<=1);
assert(CLOSE(sqr(u)+sqr(v)+sqr(w),1));

while(weight> 0){
move_in_sample(mu_t,&x,&y,&z,u,v,w);

assert(weight<=1);

while(z<0||z> t_sample){
distance_remaining= z/w;
if(z> t_sample)
distance_remaining-= t_sample/w;

x-= u*distance_remaining;
y-= v*distance_remaining;
z-= w*distance_remaining;

assert(CLOSE(z,0)||CLOSE(z,t_sample));

if(w> 0)
move_in_slide(&x,&y,&z,&u,&v,&w,&weight,z,z+t_bottom_slide,
b_bottom_slide,n_sample,n_bottom_slide,1.0);
else
move_in_slide(&x,&y,&z,&u,&v,&w,&weight,z,z-t_top_slide,
b_top_slide,n_sample,n_top_slide,1.0);

if(CLOSE(z,-t_top_slide)&&w<0){
r= add_to_reflectance_array(x,y,z,w,weight);
*r_total+= weight;
if(r> r_port_radius)
*r_lost+= weight;
weight= 0;
break;
}

if(CLOSE(z,t_sample+t_bottom_slide)&&w> 0){
r= add_to_transmittance_array(x,y,z,w,weight);
*t_total+= weight;
if(r> t_port_radius)
*t_lost+= weight;
weight= 0;
break;
}

assert((CLOSE(z,0)&&w> 0)
||(CLOSE(z,t_sample)&&w<0));
x+= u*distance_remaining;
y+= v*distance_remaining;
z+= w*distance_remaining;
}

if(weight> 0){
assert(z>=0&&z<=t_sample);
absorbed+= (1-a)*weight;
weight*= a;
roulette(&weight,&residual_weight);
scatter(g,&u,&v,&w);

assert(CLOSE(u*u+v*v+w*w,1.0));
}
}

if(total_time&&total_time<milliseconds(start_time))
break;
}

/*43:*/
#line 672 "mc_lost.w"

total= total_weight-residual_weight;
*r_lost/= total;
*t_lost/= total;
*r_total/= total;
*t_total/= total;

/*:43*/
#line 665 "mc_lost.w"


if(print_radial_arrays){
/*44:*/
#line 679 "mc_lost.w"

{
absorbed/= total;
fprintf(stderr,"%10d # number of bins\n",N_RADIAL_BINS);
fprintf(stderr,"%10.5f # bin size [mm]\n",RADIAL_BIN_SIZE);
fprintf(stderr,"# %10.5f total photons\n",total_weight);
fprintf(stderr,"# %10.5f residual photons\n",residual_weight);
fprintf(stderr,"# %10.5f total reflected light\n",*r_total);
fprintf(stderr,"# %10.5f total transmitted light\n",*t_total);
fprintf(stderr,"# %10.5f total absorbed light\n",absorbed);
fprintf(stderr,"# %10.5f total light\n",*r_total+*t_total+absorbed);
fprintf(stderr,"#    r    \t    R(r)    \t   T(r)\n");
fprintf(stderr,"#    mm    \t    W/mm2    \t   W/mm2\n");

for(i= 0;i<N_RADIAL_BINS;i++){
double area= M_PI*sqr(RADIAL_BIN_SIZE)*(2*i+1);
fprintf(stderr,"%10.5f\t",i*RADIAL_BIN_SIZE);
fprintf(stderr,"%10.5f\t",R_radial[i]/total/area);
fprintf(stderr,"%10.5f\n",T_radial[i]/total/area);
}
}

/*:44*/
#line 668 "mc_lost.w"

}
}

/*:42*/
#line 49 "mc_lost.w"

/*47:*/
#line 723 "mc_lost.w"

/*46:*/
#line 718 "mc_lost.w"

void MC_Lost(struct measure_type m,struct invert_type r,long n_photons,
double*ur1,double*ut1,double*uru,double*utu,
double*ur1_lost,double*ut1_lost,double*uru_lost,double*utu_lost)

/*:46*/
#line 724 "mc_lost.w"

{
double n_sample= m.slab_index;
double t_sample= m.slab_thickness;

double n_top= m.slab_top_slide_index;
double t_top= m.slab_top_slide_thickness;
double b_top= m.slab_top_slide_b;

double n_bottom= m.slab_bottom_slide_index;
double t_bottom= m.slab_bottom_slide_thickness;
double b_bottom= m.slab_bottom_slide_b;

double dr_port= sqrt(m.as_r)*2*m.d_sphere_r;
double dt_port= sqrt(m.as_t)*2*m.d_sphere_t;
double d_beam= m.d_beam;
double mu= m.slab_cos_angle;

int slides_differ;

/*48:*/
#line 778 "mc_lost.w"


if(t_top==0.0)
n_top= 1.0;
if(n_top==1.0)
t_top= 0.0;

if(t_bottom==0.0)
n_bottom= 1.0;
if(n_bottom==1.0)
t_bottom= 0.0;

/*:48*/
#line 744 "mc_lost.w"

/*49:*/
#line 806 "mc_lost.w"


if(m.num_spheres<2||dt_port<=0.0)
dt_port= dr_port;

/*:49*/
#line 745 "mc_lost.w"


slides_differ= (n_top!=n_bottom)||(t_top!=t_bottom)||(b_top!=b_bottom);

set_photon_seed(lost_base_seed);
MC_Radial(n_photons/2,r.a,r.b,r.g,n_sample,n_top,n_bottom,
COLLIMATED,mu,t_sample,t_top,t_bottom,b_top,b_bottom,
dr_port,dt_port,d_beam,ur1,ut1,ur1_lost,ut1_lost);

*uru_lost= 0;
*utu_lost= 0;

if(m.method==SUBSTITUTION){
set_photon_seed(lost_base_seed+DIFFUSE_SEED_OFFSET);
MC_Radial(n_photons/2,r.a,r.b,r.g,n_sample,n_top,n_bottom,
DIFFUSE,mu,t_sample,t_top,t_bottom,b_top,b_bottom,
dr_port,dt_port,d_beam,uru,utu,uru_lost,utu_lost);
}

if(m.flip_sample&&slides_differ)
/*50:*/
#line 831 "mc_lost.w"

{
double flipped_ur1,flipped_uru,flipped_ur1_lost,flipped_uru_lost;

set_photon_seed(lost_base_seed);
MC_Radial(n_photons/2,r.a,r.b,r.g,n_sample,n_bottom,n_top,
COLLIMATED,mu,t_sample,t_bottom,t_top,b_bottom,b_top,
dr_port,dt_port,d_beam,&flipped_ur1,ut1,&flipped_ur1_lost,ut1_lost);

if(m.method==SUBSTITUTION){
set_photon_seed(lost_base_seed+DIFFUSE_SEED_OFFSET);
MC_Radial(n_photons/2,r.a,r.b,r.g,n_sample,n_bottom,n_top,
DIFFUSE,mu,t_sample,t_bottom,t_top,b_bottom,b_top,
dr_port,dt_port,d_beam,&flipped_uru,utu,&flipped_uru_lost,utu_lost);
}
}

/*:50*/
#line 765 "mc_lost.w"


if(*ur1_lost<0||*ut1_lost<0||*uru_lost<0||*utu_lost<0){
exit(EXIT_FAILURE);
}
}

/*:47*/
#line 50 "mc_lost.w"

/*53:*/
#line 860 "mc_lost.w"

/*52:*/
#line 855 "mc_lost.w"

void MC_RT(struct AD_slab_type s,long n_photons,double t_sample,
double t_top_slide,double t_bottom_slide,
double*UR1,double*UT1,double*URU,double*UTU)

/*:52*/
#line 861 "mc_lost.w"

{
double ur1_lost,ut1_lost,uru_lost,utu_lost;
double dr_port= 1000;
double dt_port= 1000;
double d_beam= 0.0;
double mu= s.cos_angle;

set_photon_seed(12345);

MC_Radial(n_photons/2,s.a,s.b,s.g,s.n_slab,s.n_top_slide,s.n_bottom_slide,
COLLIMATED,mu,t_sample,t_top_slide,t_bottom_slide,
s.b_top_slide,s.b_bottom_slide,dr_port,dt_port,d_beam,UR1,UT1,&ur1_lost,&ut1_lost);

MC_Radial(n_photons/2,s.a,s.b,s.g,s.n_slab,s.n_top_slide,s.n_bottom_slide,
DIFFUSE,mu,t_sample,t_top_slide,t_bottom_slide,
s.b_top_slide,s.b_bottom_slide,dr_port,dt_port,d_beam,URU,UTU,&uru_lost,&utu_lost);
}

/*:53*/
#line 51 "mc_lost.w"

/*55:*/
#line 883 "mc_lost.w"

/*54:*/
#line 880 "mc_lost.w"

void MC_Print_RT_Arrays(int status)

/*:54*/
#line 884 "mc_lost.w"

{
print_radial_arrays= status;
}/*:55*/
#line 52 "mc_lost.w"


/*:1*/
