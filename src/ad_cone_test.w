@*1 AD Cone Test.
This file provides tests |RT_Cone|.  There is no way to completely
verify this code, but many tests of self-consistency are possible.  
That is what this code does and it should turn up any mistakes in the
implementation.  This data will ultimately be used to verify a layered
Monte Carlo implementation and consequently that will be the ultimate
test.  I just want to have some confidence in the correctness of this
code.

@(ad_cone_test.c@>=
#include <math.h>
#include <float.h>
#include <stdio.h>
#include "nr_util.h"
#include "ad_globl.h"
#include "ad_prime.h"
#include "ad_matrx.h"
#include "ad_frsnl.h"
#include "ad_prime.h"
#include "ad_cone.h"

@<Definition for |PrintTestResults|@>@;
@<Definition for |PrintUnityResults|@>@;
@<Definition for |RT_Cone_Main|@>@;
	
@ A simple utility routine to print the results nicely.
@<Definition for |PrintTestResults|@>=
static void PrintTestResults(int test, int cas, struct AD_slab_type *slab,
            double aUR1, double aUT1, double aURU, double aUTU,
            double bUR1, double bUT1, double bURU, double bUTU)
{
	printf("\nTest:%d.%d\n", test, cas);
	printf("Cone angle           %10.5f\n", acos(slab->cos_angle)*180/3.1415926);
	printf("Cosine of cone angle %10.5f\n", slab->cos_angle);
	printf("Albedo               %10.5f\n", slab->a);
	printf("Optical Depth        %10.5f\n", slab->b);
	printf("Anisotropy           %10.5f\n", slab->g);
	printf("Index for slab       %10.5f\n", slab->n_slab);
	printf("Index for top slide  %10.5f\n", slab->n_top_slide);
	printf("Index for bot slide  %10.5f\n", slab->n_bottom_slide);
	printf("            truth        cone\n");
	printf("UR1     %10.5f    %10.5f\n", aUR1, bUR1);
	printf("UT1     %10.5f    %10.5f\n", aUT1, bUT1);
	printf("URU     %10.5f    %10.5f\n", aURU, bURU);
	printf("UTU     %10.5f    %10.5f\n", aUTU, bUTU);
}

@ A simple utility routine to print the results nicely.
@<Definition for |PrintUnityResults|@>=
static void PrintUnityResults(int test, int cas, struct AD_slab_type *slab,
            double aUR1, double aUT1, double aURU, double aUTU,
            double bUR1, double bUT1, double bURU, double bUTU)
{
	double denom = 1-slab->cos_angle*slab->cos_angle;
	
	printf("\nTest:%d.%d\n", test, cas);
	printf("Cone angle           %10.5f\n", acos(slab->cos_angle)*180/3.1415926);
	printf("Cosine of cone angle %10.5f\n", slab->cos_angle);
	printf("Albedo               %10.5f\n", slab->a);
	printf("Optical Depth        %10.5f\n", slab->b);
	printf("Anisotropy           %10.5f\n", slab->g);
	printf("Index for slab       %10.5f\n", slab->n_slab);
	printf("Index for top slide  %10.5f\n", slab->n_top_slide);
	printf("Index for bot slide  %10.5f\n", slab->n_bottom_slide);
	printf("            truth        cone\n");
	printf("UR1               %10.5f\n", aUR1);
	printf("UT1               %10.5f\n", aUT1);
	printf("UR1+UT1                  %10.5f\n", aUR1+aUT1);
	printf("URU               %10.5f\n", aURU);
	printf("UTU               %10.5f\n", aUTU);
	printf("URU+UTU                  %10.5f\n", aURU+aUTU);
    printf("rc + rd/(1-mu^2) = %10.5f\n", bUR1 - (bUR1-aUR1)/denom);
    printf("tc + td/(1-mu^2) = %10.5f\n", bUT1 - (bUT1-aUT1)/denom);
    printf("           total = %10.5f\n", bUR1 - (bUR1-aUR1)/denom + 
                                      bUT1 - (bUT1-aUT1)/denom);
    printf("rc + rd/(1-mu^2) = %10.5f\n", bURU - (bURU-aURU)/denom);
    printf("tc + td/(1-mu^2) = %10.5f\n", bUTU - (bUTU-aUTU)/denom);
    printf("           total = %10.5f\n", bURU - (bURU-aURU)/denom + 
                                      bUTU - (bUTU-aUTU)/denom);
}

@ @<Definition for |RT_Cone_Main|@>=
int main (int argc, char **argv)
{
double aUR1, aURU, aUT1, aUTU, bUR1, bURU, bUT1, bUTU;
struct AD_slab_type slab;
int N=24;
double mua,musp,mus,d;

	@<Tests with full cone@>@; 
	@<Tests with no scattering@>@; 
	@<Tests with no absorption@>@;
	@<Tests with absorption and scattering@>@;
	@<Tests for Paulo@>@;
	@<Tests that vary g@>@;
	return 0;
}

@ The first set of tests just calls |RT_Cone| with no cone (normal irradiance)
or with a full cone (diffuse irradiance).  We begin simple and progressively
turn on more and more functionality.

@<Tests with full cone@>=


slab.n_slab = 1.0;
slab.n_top_slide = 1.0;
slab.n_bottom_slide = 1.0;
slab.b_top_slide = 0;
slab.b_bottom_slide = 0;
slab.a = 0.0;
slab.b = 0.1;
slab.g = 0.0;
slab.phase_function = HENYEY_GREENSTEIN;

slab.cos_angle = 1;
RT(N, &slab, &aUR1, &aUT1, &aURU, &aUTU);
slab.cos_angle = 0;
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(1,1,&slab,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

slab.a=0.5;
slab.b=0.5;
slab.g=0.875;
slab.cos_angle = 1;
RT(N, &slab, &aUR1, &aUT1, &aURU, &aUTU);
slab.cos_angle = 0;
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(1,2,&slab,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

slab.a=0.0;
slab.b=0.1;
slab.g=0.875;
slab.n_slab = 1.4;
slab.cos_angle = 1;
RT(N, &slab, &aUR1, &aUT1, &aURU, &aUTU);
slab.cos_angle = 0;
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(1,3,&slab,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

slab.a=0.5;
slab.b=0.5;
slab.g=0.875;
slab.n_slab = 1.4;
slab.cos_angle = 1;
RT(N, &slab, &aUR1, &aUT1, &aURU, &aUTU);
slab.cos_angle = 0;
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(1,4,&slab,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

slab.n_top_slide = 1.5;
slab.cos_angle = 1;
RT(N, &slab, &aUR1, &aUT1, &aURU, &aUTU);
slab.cos_angle = 0;
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(1,5,&slab,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

slab.n_bottom_slide = 1.6;
slab.cos_angle = 1;
RT(N, &slab, &aUR1, &aUT1, &aURU, &aUTU);
slab.cos_angle = 0;
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(1,6,&slab,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

@ Now for the first time I let the cosine of the cone angle be something
other than zero or one.  Now figuring out the proper values for the 
reflection and transmission for an absorbing only medium requires a
little bit of math.

Let us assume uniform illumination on an absorbing only slab.  Furthermore,
assume that the index of refraction of the slab is the same as its
environment.  In this case, the transmission through the sample is
given by
 
For diffuse irradiance at the surface, the transmitted light is
$$
T = 2\int_{0}^{\pi\over2} \exp\left[-{\mu_a z\over\cos\theta}\right]
				 \, \sin\theta\cos\theta d\theta
$$
or
$$
T = 2\int_0^1 \exp\left[-{\mu_a d\over\nu'}\right]\nu\, d\nu'
$$
Now if the illumination is only uniform over a cone with an angle
whose cosine is $\nu$, then the
transmission is
$$
T = 2\int_\nu^1 \exp\left[-{\mu_a d\over\nu'}\right]\nu\, d\nu'
$$
or 
$$
T = (\nu'-\mu_a d)\nu'\exp(-\mu_a d/\nu') + (\mu_a d)^2 \Ei(-\mu_a d/\nu')\right|_{\nu'=\nu}^1
$$
one might also be interested that when $d=0$ then 
$$
T= 1 - \nu^2
$$

@<Tests with no scattering@>=
	
slab.n_slab = 1.0;
slab.n_top_slide = 1.0;
slab.n_bottom_slide = 1.0;
slab.b_top_slide = 0;
slab.b_bottom_slide = 0;
slab.a = 0.0;
slab.b = 0.0000001;
slab.g = 0.0;
slab.phase_function = HENYEY_GREENSTEIN;

slab.cos_angle = 0.2;
RT(N, &slab, &aUR1, &aUT1, &aURU, &aUTU);
aURU = 0.0;
aUTU = 1-slab.cos_angle*slab.cos_angle;
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(2,1,&slab,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

slab.b = 1.0;
RT(N, &slab, &aUR1, &aUT1, &aURU, &aUTU);
aURU = 0.0;
aUTU = 0.219314;
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(2,2,&slab,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

@ The tests for no absorption are just to ensure that the total
reflection and transmission is unity for a whole set of cases.

@<Tests with no absorption@>=

slab.n_slab = 1.0;
slab.n_top_slide = 1.0;
slab.n_bottom_slide = 1.0;
slab.b_top_slide = 0;
slab.b_bottom_slide = 0;
slab.a = 1.0;
slab.b = 1.0;
slab.g = 0.0;
slab.phase_function = HENYEY_GREENSTEIN;

printf("*****************************************************************\n");
printf("These tests don't quite add up because they are only approximate \n");
slab.cos_angle = 0.2;
slab.a = 0.0;
RT_Cone(N, &slab, CONE, &aUR1, &aUT1, &aURU, &aUTU);
slab.a = 1.0;
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
PrintUnityResults(3,1,&slab,bUR1,bUT1,bURU,bUTU,aUR1,aUT1,aURU,aUTU);

slab.cos_angle = 0.5;
slab.a = 0.0;
RT_Cone(N, &slab, CONE, &aUR1, &aUT1, &aURU, &aUTU);
slab.a = 1.0;
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
PrintUnityResults(3,1,&slab,bUR1,bUT1,bURU,bUTU,aUR1,aUT1,aURU,aUTU);

slab.cos_angle = 0.8;
slab.a = 0.0;
RT_Cone(N, &slab, CONE, &aUR1, &aUT1, &aURU, &aUTU);
slab.a = 1.0;
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
PrintUnityResults(3,1,&slab,bUR1,bUT1,bURU,bUTU,aUR1,aUT1,aURU,aUTU);

slab.cos_angle = 0.3;
slab.g = 0.875;
slab.a = 0.0;
RT_Cone(N, &slab, CONE, &aUR1, &aUT1, &aURU, &aUTU);
slab.a = 1.0;
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
PrintUnityResults(3,1,&slab,bUR1,bUT1,bURU,bUTU,aUR1,aUT1,aURU,aUTU);

slab.n_slab = 1.4;
slab.a = 0.0;
RT_Cone(N, &slab, CONE, &aUR1, &aUT1, &aURU, &aUTU);
slab.a = 1.0;
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
PrintUnityResults(3,1,&slab,bUR1,bUT1,bURU,bUTU,aUR1,aUT1,aURU,aUTU);

slab.n_top_slide = 1.5;
slab.a = 0.0;
RT_Cone(N, &slab, CONE, &aUR1, &aUT1, &aURU, &aUTU);
slab.a = 1.0;
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
PrintUnityResults(3,1,&slab,bUR1,bUT1,bURU,bUTU,aUR1,aUT1,aURU,aUTU);

slab.n_bottom_slide = 1.6;
slab.a = 0.0;
RT_Cone(N, &slab, CONE, &aUR1, &aUT1, &aURU, &aUTU);
slab.a = 1.0;
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
PrintUnityResults(3,1,&slab,bUR1,bUT1,bURU,bUTU,aUR1,aUT1,aURU,aUTU);

@ Tests with scattering and absorption

@<Tests with absorption and scattering@>=

printf("*****************************************************************\n");
printf("These tests still need some work and will still be approximate \n");
slab.cos_angle = 0.4;
slab.b = 4.0;
slab.a = 0.99;
RT(N, &slab, &aUR1, &aUT1, &aURU, &aUTU);
aUR1 *= (1 - 0.4*0.4);
aUT1 *= (1 - 0.4*0.4);
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(4,1,&slab,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

slab.cos_angle = 0.0;
slab.n_slab = 1.5;
slab.b = 50;
RT(N, &slab, &aUR1, &aUT1, &aURU, &aUTU);
aURU = Diffuse_Glass_R(1.0,1.0,1.5);
aUTU = 0.0;
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(4,2,&slab,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

@ Finally, the motivation for all these tests.  We want Paulo to have
something to test his code with.

@<Tests for Paulo@>=
printf("*****************************************************************\n");
printf("These tests are for Paulo \n");
printf("The truth values are not correct\n");
slab.n_slab = 1.0;
slab.n_top_slide = 1.0;
slab.n_bottom_slide = 1.0;
slab.b_top_slide = 0;
slab.b_bottom_slide = 0;
slab.a = 0.99;
slab.b = 100.0;
slab.g = 0.9;
slab.phase_function = HENYEY_GREENSTEIN;

slab.cos_angle = cos(20*3.1415926/180);
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(5,1,&slab,0.0,0.0,0.0,0.0,bUR1,bUT1,bURU,bUTU);

slab.n_slab = 1.4;
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(5,2,&slab,0.0,0.0,0.0,0.0,bUR1,bUT1,bURU,bUTU);

@ @<Tests that vary g@>=

printf("\nTests that compare fraction of light collected\n");
printf("for the same absorption and reduced scattering but\n");
printf("different sphere sizes with same port size.\n");

slab.n_slab = 1.333;
slab.n_top_slide = 1.5;
slab.n_bottom_slide = 1.5;
slab.cos_angle = 1;
slab.phase_function = HENYEY_GREENSTEIN;
mua=3.9;
musp=11.5;

slab.g=0.1;
d=0.08;
while (slab.g<0.9) {
mus=musp/(1-slab.g);
slab.a=mus/(mua+mus);
slab.b=d*(mua+mus);
RT(N, &slab, &aUR1, &aUT1, &aURU, &aUTU);

slab.cos_angle = cos(atan2(10,150));
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
/*PrintTestResults(9,1,&slab,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);*/
printf("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",d,150.0,slab.g,slab.b,aUT1,bUT1);
slab.cos_angle = cos(atan2(10,70));
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
/*PrintTestResults(9,2,&slab,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);*/
printf("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",d,70.0,slab.g,slab.b,aUT1,bUT1);
slab.g += 0.2;
}
printf("\n");

slab.g=0.1;
d=0.105;
while (slab.g<0.9) {
mus=musp/(1-slab.g);
slab.a=mus/(mua+mus);
slab.b=d*(mua+mus);
RT(N, &slab, &aUR1, &aUT1, &aURU, &aUTU);

slab.cos_angle = cos(atan2(10,150));
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
/*PrintTestResults(9,1,&slab,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);*/
printf("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",d,150.0,slab.g,slab.b,aUT1,bUT1);
slab.cos_angle = cos(atan2(10,70));
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
/*PrintTestResults(9,2,&slab,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);*/
printf("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",d,70.0,slab.g,slab.b,aUT1,bUT1);
slab.g += 0.2;
}
printf("\n");

slab.g=0.1;
d=0.2;
while (slab.g<0.9) {
mus=musp/(1-slab.g);
slab.a=mus/(mua+mus);
slab.b=d*(mua+mus);
RT(N, &slab, &aUR1, &aUT1, &aURU, &aUTU);

slab.cos_angle = cos(atan2(10,150));
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
/*PrintTestResults(9,1,&slab,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);*/
printf("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",d,150.0,slab.g,slab.b,aUT1,bUT1);
slab.cos_angle = cos(atan2(10,70));
RT_Cone(N, &slab, CONE, &bUR1, &bUT1, &bURU, &bUTU);
/*PrintTestResults(9,2,&slab,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);*/
printf("%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",d,70.0,slab.g,slab.b,aUT1,bUT1);
slab.g += 0.2;
}
