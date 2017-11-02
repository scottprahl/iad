@*1 AD Global Variables.
Global Routines and Variables.
Changed version to reflect bug fix in the Fresnel routine section.

Revised in May 1995 to allow slides to absorb and various modifications
to improve the way that the file looks.

Revision May 1996 to remove uninitialized tfluence

Revision May 1998 to improve |wrarray|.

@(ad_globl.c@>=
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ad_globl.h"
#include "ad_frsnl.h"

	@<Global variables for adding-doubling@>@;
	@<Definition for |Zero_Layer|@>@;
	@<Definition for |AD_error|@>@;
	@<Definition for |URU_and_UR1|@>@;
	@<Definition for |URU_and_UR1_Cone|@>@;
	@<Definition for |URU_and_URx_Cone|@>@;
	@<Definition for |UFU_and_UF1|@>@;
	@<Definition for |wrmatrix|@>@;
	@<Definition for |wrarray|@>@;

@ @(ad_globl.h@>=
	@h
 	@<Types to export from AD Globals@>@;
	@<External variables to export from AD Globals@>@;
	@<Prototype for |Zero_Layer|@>;
	@<Prototype for |AD_error|@>;
	@<Prototype for |URU_and_UR1|@>;
	@<Prototype for |URU_and_UR1_Cone|@>;
	@<Prototype for |URU_and_URx_Cone|@>;
	@<Prototype for |UFU_and_UF1|@>;
	@<Prototype for |wrmatrix|@>;
	@<Prototype for |wrarray|@>;
	
@*2 Constants.

This is Version 2.0.0 of the adding-doubling code.  (The inverse adding-doubling 
code may have a different version number.)

@ The number of quadrature points determines how accurately the integrals
are performed.  Larger numbers of quadrature points lead to more accurate
solutions.  Fewer points yield much faster computations since the computation
time is proportional to $n^3$ or $n^2\ln n$ because an $n\times n$ matrix 
must be inverted.  

For most practical
purposes four quadrature points is plenty.  However, if you need very accurate reflection and
transmission values, then increase the number of quadrature points.  For example,
if you want to verify a Monte Carlo implementation, then just crank the number
up to 16 or 32 and you are almost certain to get 5 significant digits in your
answer.  

The number of quadrature points does not need to be a power of 2, but it should
be an even number.  If it isn't then somewhere in the bowels of this program
it will get changed.  Finally, if you are unsure of how accurate a solution is, then
increase the number of quadrature points and repeat the algorithm.

There is no intrinsic reason that the maximum number of quadrature points is limited
to 128.  If you have enough memory then this number can be increased.  But if
you have read the stuff above, my feeling is, why bother?  

@d MAX_QUAD_PTS     128
@d DEFAULT_QUAD_PTS  4

@ The two permissible phase functions are isotropic and Henyey-Greenstein.

@d ISOTROPIC 0
@d HENYEY_GREENSTEIN 1

@ The last two constants are related to the details of how the initial
adding-doubling layer is generated.  It is very unlikely that these will
ever be used by anyone.

@d DIAMOND         0
@d INFINITESIMAL_GENERATOR  1

@ This last define is so that intermediate values can be generated
during the calculation of the initial layer matrices.  It is named
after Martin Hammer who requested it.

@d MARTIN_HAMMER 1

@ And finally something for whether the light is conical or oblique

@d CONE      1
@d OBLIQUE   0

@*2 Types.

The fundamental structure for an adding-doubling calculation keeps all the
details of the optical properties of the sample together.  The sample is
bounded by a glass slide above and below.  The glass slides have indicies
of refraction |n_top_slide| and |n_bottom_slide|.  The glass slides may
absorb light, in which case |b_top_slide| or |b_bottom_slide| may be non-zero.

The albedo of the slab is denoted |a|, the optical thickness of the slab by
$b=(\mu_a+\mu_s) d$, and the average cosine of the phase function by |g|.
The phase function of the slab is restricted to just isotropic and Henyey-Greenstein
phase functions at the moment.

@<Types to export from AD Globals@>=

typedef struct AD_slab_type {
  double a; 
  double b; 
  double g; 
  int phase_function; 
  double n_slab;	
  double n_top_slide; 
  double n_bottom_slide; 
  double b_top_slide;
  double b_bottom_slide;
  double cos_angle;
} slab_type;

@ @<Types to export from AD Globals@>=
typedef struct AD_method_type {
  int quad_pts;
  double a_calc, b_calc, g_calc, b_thinnest;
} method_type;

@ The |Martin_Hammer| variable only exists to print internal
results when testing.  Its only a integer and doesn't take
up much space so here it is.

@<Global variables for adding-doubling@>=
#define AD_GLOBAL_SOURCE
double angle[MAX_QUAD_PTS+1];
double weight[MAX_QUAD_PTS+1];
double twoaw[MAX_QUAD_PTS+1];
int Martin_Hammer=0;

@ @<External variables to export from AD Globals@>=
#ifndef AD_GLOBAL_SOURCE
extern double angle[MAX_QUAD_PTS+1];
extern double weight[MAX_QUAD_PTS+1];
extern double twoaw[MAX_QUAD_PTS+1];
extern int Martin_Hammer;

#endif

@*2 Global routines.
My standard error handler
@<Prototype for |AD_error|@>=
void AD_error(char error_text[])

@ @<Definition for |AD_error|@>=
	@<Prototype for |AD_error|@>
{
  fprintf(stderr,"Adding-Doubling error\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

@ @<Prototype for |Zero_Layer|@>=
void Zero_Layer(int n, double **r, double **t)
	
@ @<Definition for |Zero_Layer|@>=
	@<Prototype for |Zero_Layer|@>
{
  int i, j;

  for (i = 1; i <= n; i++) 
    for (j = 1; j <= n; j++) {
 	  t[i][j] = 0.0;
      r[i][j] = 0.0;
    }

  for (i = 1; i <= n; i++) 
	t[i][i] = 1 / twoaw[i];
}

@ Figure out the reflection for collimated irradiance returning within
a cone whose cosine is |mu|.  Note that |mu| is defined on the air side
of the slab and that |mu| is the cosine of the angle that the cone makes
with the normal to the slab,
$$
\hbox{UR1} \equiv \int_\mu^1 R(\nu',1) 2\nu'\,d\nu'
$$
Similarly for irradiance characterized by diffuse light within a cone
one can calculate the amount of reflectance returing within that cone
as
$$
\hbox{URU} \equiv n^2 \int_\mu^1 \int_0^1 R(\nu',\nu'') 2\nu'\,d\nu' 2\nu''\,d\nu''
$$
where, $n^2$ term is to account for the $n^2$ law of radiance.

@<Prototype for |URU_and_UR1_Cone|@>=
void URU_and_UR1_Cone(int n, double n_slab, double mu, double **R, double *URU, double *UR1)

@ @<Definition for |URU_and_UR1_Cone|@>=
	@<Prototype for |URU_and_UR1_Cone|@>
{
  int i, j, last_j;
  double mu_slab;
  double temp = 0.0;

  if (n_slab == 1) 
  	  mu_slab = mu;
  else 
      mu_slab = sqrt(n_slab*n_slab-1+mu*mu)/n_slab;
  
  last_j = 1;
  while (angle[last_j] <= mu_slab)
    last_j++;

  *URU = 0.0;
  for (i = 1; i <= n; i++) {
    temp = 0.0;
    for (j = last_j; j <= n; j++)
      temp += R[i][j] * twoaw[j];
    *URU += temp * twoaw[i];
  }
  *UR1 = temp;
  *URU *= n_slab * n_slab / (1 - mu * mu);
}

@ Figure out the reflection for oblique irradiance returning from a layer
Note that |mu| is the cosine of the angle that the cone makes
with the normal to the slab in air,
$$
\hbox{URx} = \int_\mu^1 R(\nu',\mu) 2\nu'\,d\nu'
$$
For diffuse irradiance over the cone, the total flux back |URU| is somewhat
arbitrarily chosen as the that flux returning in the same cone.  Specifically
as
$$
\hbox{URU} = n^2 \int_\mu^1 \int_\mu^1 R(\nu',\nu'') 2\nu'\,d\nu' 2\nu''\,d\nu''
$$
where, $n^2$ term is to account for the $n^2$ law of radiance. (If you want
the total flux returning within a cone for uniform diffuse illumination then 
use |URU_and UR1_Cone|.)  

@<Prototype for |URU_and_URx_Cone|@>=
void URU_and_URx_Cone(int n, double n_slab, double mu, double **R, double *URU, double *URx)

@ @<Definition for |URU_and_URx_Cone|@>=
	@<Prototype for |URU_and_URx_Cone|@>
{
	int i, j, cone_index;
	double mu_slab, urx, delta, closest_delta;
	double degrees = 180.0/3.1415926535;
	
	mu_slab = sqrt(n_slab*n_slab-1+mu*mu)/n_slab;
	
	closest_delta = 1;
	cone_index = n;
	
	for (i=n; i>=1; i--) {
		delta = fabs(angle[i] - mu_slab);
		if (delta < closest_delta) {
			closest_delta = delta;
			cone_index = i;
		}
	}
	
	if (fabs(angle[cone_index] - mu_slab) > 1e-5) {
		fprintf(stderr, "Something is wrong with the quadrature\n");
		fprintf(stderr, "theta_i = %5.2f degrees or ", acos(mu)*degrees);
		fprintf(stderr, "cos(theta_i) = %8.5f\n", mu);
		fprintf(stderr, "theta_t = %5.2f degrees or ", acos(mu_slab)*degrees);
		fprintf(stderr, "cos(theta_t) = %8.5f\n", mu_slab);
		fprintf(stderr, " index  degrees cosine\n");
		for (i=n; i>=1; i--) {
		    fprintf(stderr, " %5d   %5.2f ", i, acos(angle[i])*degrees);
		    fprintf(stderr, " %8.5f\n", angle[i]);
		}

		fprintf(stderr, "Closest quadrature angle is i=%5d ", cone_index);
		fprintf(stderr, "or cos(theta)=%8.5f\n", angle[cone_index]);
		fprintf(stderr, "Assuming normal incidence\n");
	}
	
    *URU = 0.0;
	for (i = 1; i <= n; i++) {
		urx = 0.0;
		for (j = 1; j <= n; j++)
			urx += R[i][j] * twoaw[j];
		
		*URU += urx * twoaw[i];
		if (i == cone_index) *URx = urx;
	}
	*URU *= n_slab * n_slab;
}

@ Just add up all the angles up to the critical angle.  This is 
a commonly used convenience function to easily calculate |UR1| and
|URU|.  We select the entire range of angles by passing $\cos(\pi/2)= 0$
to the |URU_and_UR1_Cone| routine.

@<Prototype for |URU_and_UR1|@>=
void URU_and_UR1(int n, double n_slab, double **R, double *URU, double *UR1)

@ @<Definition for |URU_and_UR1|@>=
	@<Prototype for |URU_and_UR1|@>
{
  URU_and_UR1_Cone(n, n_slab, 0.0, R, URU, UR1);
}


@ @<Prototype for |UFU_and_UF1|@>=
void UFU_and_UF1(int n, double n_slab, 
                        double **Lup, double **Ldown, double *UFU, double *UF1)

@ @<Definition for |UFU_and_UF1|@>=
	@<Prototype for |UFU_and_UF1|@>
{
  int i, j;
  double temp = 0.0;

  *UFU = 0.0;
  for (j = 1; j <= n; j++) {
    temp = 0.0;
    for (i = 1; i <= n; i++)
      temp += (Lup[i][j] + Ldown[i][j]) * 2 * weight[i];
    *UFU += twoaw[j] * temp;
  }
  *UF1 = temp * n_slab * n_slab;
  *UFU *= n_slab * n_slab / 2;
}


@ @<Prototype for |wrmatrix|@>=
 	void wrmatrix(int n, double **a)

@ @<Definition for |wrmatrix|@>=
	@<Prototype for |wrmatrix|@>
{
  int i, j;
  double tflux, flux;
  
  printf("%9.5f", 0.0);
  for (i = 1; i <= n; i++)
    printf("%9.5f", angle[i]);
    
  printf("     flux\n" );

  tflux = 0.0;
  for (i = 1; i <= n; i++) {
    printf("%9.5f", angle[i]);
    for (j = 1; j <= n; j++)
      if ((a[i][j]>10) || (a[i][j]<-10))
	  	printf("    *****");
	  else
	    printf("%9.5f", a[i][j]);
    flux = 0.0;
    for (j = 1; j <= n; j++)
      if ((a[i][j]<10) && (a[i][j]>-10)) flux += a[i][j] * twoaw[j];
    printf("%9.5f\n", flux);
    tflux += flux * twoaw[i];
  }
 
  printf("%9s", "flux   ");
  for (i = 1; i <= n; i++) {
    flux = 0.0;
    for (j = 1; j <= n; j++)
      if ((a[j][i]<10) && (a[j][i]>-10)) flux += a[j][i] * twoaw[j];
    printf("%9.5f", flux);
  }
  printf("%9.5f\n", tflux);
  for (i = 1; i <= (n + 2); i++)
    printf("*********");
  printf("\n\n");
}


@ @<Prototype for |wrarray|@>=
void wrarray(int n, double *a)

@ @<Definition for |wrarray|@>=
	@<Prototype for |wrarray|@>
{
  int i;
  double sum;

  for (i = 1; i <= n; i++)
    printf("%9.5f", angle[i]);
  printf("%9s\n", " angles");

  sum = 0.0;
  for (i = 1; i <= n; i++) {
      if (a[i]>10 || a[i]<-10)
	  	printf("    *****");
	  else
    	printf("%9.5f", a[i]);
    if (a[i]<10 && a[i]<-10) sum += a[i];
  }
  printf("%9.5f", sum);
  printf("%9s\n", " (natural)");

  sum = 0.0;
  for (i = 1; i <= n; i++) {
      if (a[i]>10 || a[i]<-10)
	  	printf("    *****");
	  else
    	printf("%9.5f", a[i]/twoaw[i]);
    if (a[i]<10 && a[i]<-10) sum += a[i];
  }
  printf("%9.5f", sum);
  printf("%9s\n", "*2aw");
  for (i = 1; i <= (n + 2); i++)
    printf("*********");
  printf("\n\n");
}

@ Just print out an array without mucking @<Prototype for |swrarray|@>=
void swrarray(int n, double *a)

@ @<Definition for |swrarray|@>=
	@<Prototype for |swrarray|@>
{
  int i;
  double sum;

  for (i = 1; i <= n; i++)
    printf("%9.5f", angle[i]);
  printf("%9s\n", "*2aw");
  sum = 0.0;
  for (i = 1; i <= n; i++) {
      if (a[i]>10 || a[i]<-10)
	  	printf("    *****");
	  else
    	printf("%9.5f", a[i]/twoaw[i]);
    if (a[i]<10 && a[i]<-10) sum += a[i];
  }
  printf("%9.5f\n", sum);
  for (i = 1; i <= (n + 2); i++)
    printf("*********");
  printf("\n\n");
}

