@*1 AD Layers Test.
This file provides tests |RT_Layers|.  There is no way to completely
verify this code, but many tests of self-consistency are possible.  
That is what this code does and it should turn up any mistakes in the
implementation.  This data will ultimately be used to verify a layered
Monte Carlo implementation and consequently that will be the ultimate
test.  I just want to have some confidence in the correctness of this
code.

@(ad_layers_test.c@>=
#include <math.h>
#include <float.h>
#include <stdio.h>
#include "nr_util.h"
#include "ad_globl.h"
#include "ad_prime.h"
#include "ad_matrx.h"
#include "ad_prime.h"
#include "ad_layers.h"

@<Definition for |PrintTestResults|@>@;
@<Definition for |RT_Layers_Main|@>@;
	
@ A simple utility routine to print the results nicely.
@<Definition for |PrintTestResults|@>=
static void PrintTestResults(int test, int cas,
            double aUR1, double aUT1, double aURU, double aUTU,
            double bUR1, double bUT1, double bURU, double bUTU)
{
	printf("\nTest:%d.%d\n", test, cas);
	printf("            truth        layers\n");
	printf("UR1     %10.5f    %10.5f\n", aUR1, bUR1);
	printf("UT1     %10.5f    %10.5f\n", aUT1, bUT1);
	printf("URU     %10.5f    %10.5f\n", aURU, bURU);
	printf("UTU     %10.5f    %10.5f\n", aUTU, bUTU);
}

@ @<Definition for |RT_Layers_Main|@>=
int main (int argc, char **argv)
{
double aUR1, aURU, aUT1, aUTU, bUR1, bURU, bUT1, bUTU;
double dUR1,dUT1,dURU,dUTU,cUR1,cUT1,cURU,cUTU;

double a[15], b[15], g[15];
int i;
struct AD_slab_type slab;
int N=32;

	@<Tests with a single layers@>@;
	@<Tests with similar layers@>@;
	@<Tests with clear layers@>@;
	@<Tests with absorbing bounding layers@>@;
	@<Tests with absorbing and clear layers@>@;
	@<Tests with layers reversed@>@;
	@<Tests for Yinchu@>@;
	return 0;
}

@ The first set of tests just call |RT_Layers| with a single slab
instead of the usual call to |RT|.  We start with unscattering samples,
then add boundaries to the unscattering samples.  Finally we add
scattering back in and repeat the same set of tests.

@<Tests with a single layers@>=

a[0]=0.0;
b[0]=0.1;
g[0]=0.875;
ez_RT(N, 1.0, 1.0, 1.0, a[0], b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.0, 1.0, 1.0, 1, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(1,1,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

a[0]=0.0;
b[0]=0.1;
g[0]=0.875;
ez_RT(N, 1.4, 1.0, 1.0, a[0], b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.4, 1.0, 1.0, 1, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(1,2,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

a[0]=0.0;
b[0]=0.1;
g[0]=0.875;
ez_RT(N, 1.4, 1.5, 1.5, a[0], b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.4, 1.5, 1.5, 1, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(1,3,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

a[0]=0.0;
b[0]=0.1;
g[0]=0.875;
ez_RT(N, 1.4, 1.5, 1.6, a[0], b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.4, 1.5, 1.6, 1, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(1,4,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

a[0]=0.5;
b[0]=0.1;
g[0]=0.875;
ez_RT(N, 1.0, 1.0, 1.0, a[0], b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.0, 1.0, 1.0, 1, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(1,5,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

a[0]=0.5;
b[0]=0.1;
g[0]=0.875;
ez_RT(N, 1.4, 1.0, 1.0, a[0], b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.4, 1.0, 1.0, 1, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(1,6,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

a[0]=0.5;
b[0]=0.1;
g[0]=0.875;
ez_RT(N, 1.4, 1.5, 1.5, a[0], b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.4, 1.5, 1.5, 1, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(1,7,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

a[0]=0.5;
b[0]=0.1;
g[0]=0.875;
ez_RT(N, 1.4, 1.5, 1.6, a[0], b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.4, 1.5, 1.6, 1, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(1,8,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

@ The next set of tests just divide a slab with thickness $b$ into $m$ 
individual slabs with thickness $\Delta b = b/m$.  By comparing the
results of |RT_Layers| with the result for a simple homogeneous slab
we can hope to catch obvious bugs.  This is done for  2 and 5 
individual layers.  We then let the slab have an index of refraction
greater than 1.0 and try again.  Then we add two identical slides,
then two non-equal slides.

@<Tests with similar layers@>=

for (i=0; i<2; i++) {
	a[i]=0.5;
	b[i]=0.05;
	g[i]=0.875;
}
ez_RT(N, 1.0, 1.0, 1.0, a[0], 2*b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.0, 1.0, 1.0, 2, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(2,1,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

ez_RT(N, 1.4, 1.0, 1.0, a[0], 2*b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.4, 1.0, 1.0, 2, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(2,2,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

ez_RT(N, 1.4, 1.5, 1.5, a[0], 2*b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.4, 1.5, 1.5, 2, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(2,3,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

ez_RT(N, 1.4, 1.5, 1.6, a[0], 2*b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.4, 1.5, 1.6, 2, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(2,4,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

for (i=0; i<5; i++) {
	a[i]=0.5;
	b[i]=0.02;
	g[i]=0.875;
}
ez_RT(N, 1.0, 1.0, 1.0, a[0], 5*b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.0, 1.0, 1.0, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(2,5,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

ez_RT(N, 1.4, 1.0, 1.0, a[0], 5*b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.4, 1.0, 1.0, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(2,6,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

ez_RT(N, 1.4, 1.5, 1.5, a[0], 5*b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.4, 1.5, 1.5, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(2,7,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

ez_RT(N, 1.4, 1.5, 1.6, a[0], 5*b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.4, 1.5, 1.6, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(2,8,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

@ This set of tests with incorporate perfectly clear layers included
at various depths within the slab.  The clear layers should not 
affect the result.

@<Tests with clear layers@>=

for (i=0; i<5; i++) {
	a[i]=0.5;
	b[i]=0.02;
	g[i]=0.875;
}

b[0]=0.0;
ez_RT(N, 1.0, 1.0, 1.0, a[0], 4*b[1], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.0, 1.0, 1.0, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(3,1,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

ez_RT(N, 1.4, 1.5, 1.6, a[0], 4*b[1], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.4, 1.5, 1.6, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(3,2,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

for (i=0; i<5; i++) {
	a[i]=0.5;
	b[i]=0.02;
	g[i]=0.875;
}

b[4]=0.0;
ez_RT(N, 1.0, 1.0, 1.0, a[0], 4*b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.0, 1.0, 1.0, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(3,3,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

ez_RT(N, 1.4, 1.5, 1.6, a[0], 4*b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.4, 1.5, 1.6, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(3,4,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

for (i=0; i<5; i++) {
	a[i]=0.5;
	b[i]=0.02;
	g[i]=0.875;
}

b[3]=0.0;
ez_RT(N, 1.0, 1.0, 1.0, a[0], 4*b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.0, 1.0, 1.0, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(3,5,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

ez_RT(N, 1.4, 1.5, 1.6, a[0], 4*b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.4, 1.5, 1.6, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(3,6,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

for (i=0; i<5; i++) {
	a[i]=0.5;
	b[i]=0.02;
	g[i]=0.875;
}

b[2]=0.0;
b[4]=0.0;
ez_RT(N, 1.0, 1.0, 1.0, a[0], 3*b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.0, 1.0, 1.0, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(3,7,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

ez_RT(N, 1.4, 1.5, 1.6, a[0], 3*b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.4, 1.5, 1.6, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(3,8,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

@ For this set of tests, I leverage the code added for absorbing
slides.  This way I can add a single absorbing layer above, below, both
above and below.  I cannot really do the case when the index of refraction
of the slide differs from the slab, but I can certainly do the case
when both the slide and the slab have equal non-unity indices of refraction.

@<Tests with absorbing bounding layers@>=

slab.n_slab = 1.0;
slab.n_top_slide = 1.0;
slab.n_bottom_slide = 1.0;
slab.b_top_slide = 0.2;
slab.b_bottom_slide = 0.1;
slab.a = 0.9;
slab.b = 2.0;
slab.g = 0.0;
slab.phase_function = HENYEY_GREENSTEIN;
slab.cos_angle = 1.0;

a[0]=0.0;
b[0]=0.2;
g[0]=0.0;
a[1]=0.9;
b[1]=2.0;
g[1]=0.0;
a[2]=0.0;
b[2]=0.1;
g[2]=0.0;

RT(N, &slab, &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.0, 1.0, 1.0, 3, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(4,1,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

slab.n_slab = 1.4;
slab.n_top_slide = 1.4;
slab.n_bottom_slide = 1.4;
RT(N, &slab, &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.4, 1.4, 1.4, 3, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(4,2,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

@ The fourth set of tests includes both absorbing layers above and below
with various clear layers interspersed within the slab.

@<Tests with absorbing and clear layers@>=

slab.n_slab = 1.0;
slab.n_top_slide = 1.0;
slab.n_bottom_slide = 1.0;
slab.b_top_slide = 0.02;
slab.b_bottom_slide = 0.02;
slab.a = 0.5;
slab.b = 0.06;
slab.g = 0.875;
slab.phase_function = HENYEY_GREENSTEIN;
slab.cos_angle = 1.0;

/* set up abs/scat/clear/scat/clear/scat/abs */
for (i=0; i<7; i++) {
	a[i]=0.5;
	b[i]=0.02;
	g[i]=0.875;
}

a[0]=0.0;	/* top layer only absorbs */
b[2]=0.0;	/* third layer is clear   */
b[4]=0.0;	/* fifth layer is clear   */
a[6]=0.0;   /* seventh layer is abs   */

RT(N, &slab, &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.0, 1.0, 1.0, 7, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(5,1,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

slab.n_slab = 1.4;
slab.n_top_slide = 1.4;
slab.n_bottom_slide = 1.4;
RT(N, &slab, &aUR1, &aUT1, &aURU, &aUTU);
RT_Layers(N, 1.4, 1.4, 1.4, 7, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(5,2,aUR1,aUT1,aURU,aUTU,bUR1,bUT1,bURU,bUTU);

@ The fifth set of tests involves a series of scattering and absorbing layers
in which the reflection and transmission are found for light going downwards.
The order of the layers is then reversed and the downward light in the
first case is compared with the upward light in the second case.

@<Tests with layers reversed@>=
for (i=0; i<5; i++) {
	a[i]=0.5;
	b[i]=0.1;
	g[i]=0.875;
}

a[0]=0.1;
a[1]=0.4;
a[2]=0.5;
a[3]=0.3;
b[4]=3;
a[4]=0.99;

RT_Layers_All(N, 1.4, 1.4, 1.4, 5, a, b, g, \
&aUR1, &aUT1, &aURU, &aUTU, &bUR1, &bUT1, &bURU, &bUTU);

for (i=0; i<5; i++) {
	a[i]=0.5;
	b[i]=0.1;
	g[i]=0.875;
}
a[0]=0.99;
b[0]=3;
a[1]=0.4;
a[2]=0.5;
a[3]=0.3;
a[4]=0.1;

RT_Layers_All(N, 1.4, 1.4, 1.4, 5, a, b, g,\
 &cUR1, &cUT1, &cURU, &cUTU, &dUR1, &dUT1, &dURU, &dUTU);
PrintTestResults(6,1,aUR1,aUT1,aURU,aUTU,dUR1,dUT1,dURU,dUTU);
PrintTestResults(6,2,bUR1,bUT1,bURU,bUTU,cUR1,cUT1,cURU,cUTU);

@ This set of tests is so that Yinchu can check her multi-layered Monte
Carlo code.

@<Tests for Yinchu@>=

for (i=0; i<5; i++) {
	a[i]=0.95;
	b[i]=0.5;
	g[i]=0.875;
}
a[0]=0.1;

RT_Layers(N, 1.0, 1.0, 1.0, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(7,1,bUR1,bUT1,bURU,bUTU,0,0,0,0);

RT_Layers(64, 1.4, 1.0, 1.0, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(7,2,bUR1,bUT1,bURU,bUTU,0,0,0,0);

for (i=0; i<5; i++) {
	a[i]=0.95;
	b[i]=0.5;
	g[i]=0.875;
}
a[1]=0.1;

RT_Layers(N, 1.0, 1.0, 1.0, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(7,3,bUR1,bUT1,bURU,bUTU,0,0,0,0);

RT_Layers(64, 1.4, 1.0, 1.0, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(7,4,bUR1,bUT1,bURU,bUTU,0,0,0,0);

for (i=0; i<5; i++) {
        a[i]=0.95;
        b[i]=0.05;
        g[i]=0.875;
}

RT_Layers(N, 1.0, 1.0, 1.0, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(7,5,bUR1,bUT1,bURU,bUTU,0,0,0,0);

RT_Layers(64, 1.4, 1.0, 1.0, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
PrintTestResults(7,6,bUR1,bUT1,bURU,bUTU,0,0,0,0);

