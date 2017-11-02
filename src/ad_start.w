@*1 AD Start.
This has the routines for forming the initial matrix to start off an
adding-doubling calculation.  

Added printing of intermediate results for Martin Hammer.

@c
#include <math.h>
#include <float.h>
#include <stdio.h>

#include "ad_frsnl.h"
#include "ad_globl.h"
#include "ad_matrx.h"
#include "ad_phase.h"
#include "ad_radau.h"
#include "ad_start.h"
#include "nr_gaulg.h"
#include "nr_util.h"

	@<Definition for |Get_Start_Depth|@>@;
	@<Definition for |Quadrature|@>@;
	@<Definition for |Choose_Method|@>@;
	@<Definition for |Choose_Cone_Method|@>@;
	@<Definition for |Get_IGI_Layer|@>@;
	@<Definition for |Get_Diamond_Layer|@>@;
	@<Definition for |Init_Layer|@>@;

@ @(ad_start.h@>=
	@<Prototype for |Get_Start_Depth|@>;
	@<Prototype for |Choose_Method|@>;
	@<Prototype for |Choose_Cone_Method|@>;
	@<Prototype for |Init_Layer|@>;
	@<Prototype for |Quadrature|@>;

@*2 Basic routines.

This file contains the three procedures which must be called before any doubling 
may take place.  They should be called in the following order: 

	|Choose_Method|			---  to fill the method record 

	|Quadrature|				---  to calculate the quad angles and weights 

	code to initialize |angle|, |weight|, and |twoaw|
	
	|Init_Layer|					---  to calculate the thin layer |R| and |T| 

	|Double_Until| --- to obtain |R| and |T| for the desired thickness 


@  |Get_Start_Depth| selects the best minimum starting thickness 
to start the doubling process.
The criterion is based on an assessment of the (1) round-off error, 
(2) the angular initialization error, and (3) the thickness 
initialization error.  Wiscombe concluded that an optimal
starting thickness depends on the smallest quadrature angle, and
recommends that when either the infinitesimal generator or diamond
initialization methods are used then the initial thickness is optimal 
when type 2 and 3 errors are comparable, or when 
$$
         d \approx \mu 
$$
Note that round-off is important when the starting thickness is less than 
|1e-4| for diamond initialization and less than
|1e-8| for infinitesimal generator initialization assuming
about 14 significant digits of accuracy.
  
Since the final thickness is determined by repeated doubling, the
starting thickness is found by dividing by 2 until the starting thickness is
less than $\mu$.  Also we make checks for a layer with zero thickness
and one that infinitely thick.

@<Prototype for |Get_Start_Depth|@>=
	double Get_Start_Depth(double mu, double d)

@ @<Definition for |Get_Start_Depth|@>=
	@<Prototype for |Get_Start_Depth|@>

{
	if (d <= 0) 
		return 0.0;

	if (d == HUGE_VAL)
		return (mu / 2.0);

	while (d > mu) d /= 2;
	
	return d;
}

@*2 Quadrature.

@ This returns the quadrature angles using Radau quadrature over the 
interval 0 to 1 if there is no critical angle for total internal reflection
in the slab.  If there is a critical angle whose cosine is $\mu_c$ then
Radau quadrature points are chosen from 0 to $\mu_c$ and Radau
quadrature points over the interval $\mu_c$ to 1.

@<Prototype for |Quadrature|@>=
	 void Quadrature(int n, double n_slab, double *x, double *w)

@ @<Definition for |Quadrature|@>=
	@<Prototype for |Quadrature|@>
{
	int i, nby2;
	double *x1, *w1;
	double mu_c;

	if (n_slab == 1) {
		Radau(0.0, 1.0, x, w, n);
        return;
	}

	mu_c = Cos_Critical_Angle(n_slab,1.0);
	nby2 = n / 2;
	gauleg(0.0, mu_c, x, w, nby2);

	x1 = dvector(1,nby2);
	w1 = dvector(1,nby2);
	Radau(mu_c, 1.0, x1, w1, nby2);
	for (i = 1; i <= nby2; i++) {
		x[nby2 + i] = x1[i];
		w[nby2 + i] = w1[i];
	}
	free_dvector(x1,1,nby2);
	free_dvector(w1,1,nby2);
}

@ |Choose_Method| fills the method structure with correct values for
|a_calc|, |b_calc|, |g_calc|, and |b_thinnest| based on the delta-M
method.  Furthermore, the quadrature angles and weights are also calculated.
Before calling this routines |method.quad_pts| must be set to some 
multiple of 2.  If this routine is not called then it is up to you
to

1. to fill the method record appropriately 

2. call |Quadrature|  

3. fill global arrays |angle|, |weight|, and |twoaw| 

4. determine the thickness of the thinnest layer 


@<Prototype for |Choose_Method|@>=
	void Choose_Method(struct AD_slab_type * slab, struct AD_method_type *method)

@ @<Definition for |Choose_Method|@>=
	@<Prototype for |Choose_Method|@>

{
	double af;
	int i, n;

	if (0<slab->cos_angle && slab->cos_angle < 1) {
		Choose_Cone_Method(slab,method);
		return;
	}
	
	n = method->quad_pts;
	af = pow(slab->g, n)*slab->a;
	method->a_calc = (slab->a - af) / (1 - af);
	method->b_calc = (1 - af) * slab->b;
	method->g_calc = slab->g;
	
	Quadrature(n, slab->n_slab, angle, weight);

	for (i = 1; i <=n; i++)
		twoaw[i] = 2 * angle[i] * weight[i];

	method->b_thinnest = Get_Start_Depth(angle[1],method->b_calc);
}

@ |Choose_Cone_Method| adds the ability to specify a specific quadrature angle
so that accurate estimates of the reflection and transmission might be made
for when the light returning in a particular cone is of interest.  This code
mimicks the usual |Choose_Method| above, and in fact explicitly uses it for
a couple of special cases.  

@<Prototype for |Choose_Cone_Method|@>=
	void Choose_Cone_Method(struct AD_slab_type *slab, 
							struct AD_method_type *method)

@ @<Definition for |Choose_Cone_Method|@>=
	@<Prototype for |Choose_Cone_Method|@>

{
	double af, *angle1, *weight1, cos_crit_angle,mu;
	int i, n, nby2, nby3;

	n = method->quad_pts;
	af = pow(slab->g, n)*slab->a;
	method->a_calc = (slab->a - af) / (1 - af);
	method->b_calc = (1 - af) * slab->b;
	method->g_calc = slab->g;

	@<Special case when cosine is zero@>@;
	@<Special case when no index of refraction change@>@;
	@<Gaussian quadrature from 0 to the critical angle@>@;
	@<Radau quadrature from the critical angle to the cone angle@>@;
	@<Radau quadrature from the cone angle to 1@>@;
}

@ @<print angles@>=

@ @<debug print angles@>=
	{
		printf("****Cone Angle          = %6.2f degrees, Cosine()=%6.4f\n",
		   			acos(slab->cos_angle)*180.0/3.14159,slab->cos_angle);
		double sum=0;
		for (i = 1; i <=n; i++) {
			sum += twoaw[i];
			printf("%02d theta=%6.2f cos(theta)=%6.4f w=%6.4f 2aw=%6.4f\n",
			   i, acos(angle[i])/3.1415926*180.0, angle[i], weight[i],twoaw[i]);
		}
		printf("twoaw sum = %8.4f\n",sum);
	}

@ When the cone angle is zero or ninety degrees then we can just use the
standard method for choosing the quadrature points.

@<Special case when cosine is zero@>=
	if (slab->cos_angle == 0 || slab->cos_angle == 1) {
		Choose_Method(slab, method);
		@<print angles@>@;
		return;
	}
	
@ When there is no index of refraction change, there is no critical angle
to worry about.  Since we want the cone angle to be included as one of our
angles, we use Radau quadrature.  That way both the cone angle and 
perpendicular angles are included.

@<Special case when no index of refraction change@>=
	if (slab->n_slab == 1 && slab->n_top_slide == 1 && slab->n_bottom_slide == 1) {
		nby2 = n / 2;
		Radau(0.0, slab->cos_angle, angle, weight, nby2);
	
		angle1 = dvector(1, nby2);
		weight1 = dvector(1, nby2);
		Radau(slab->cos_angle, 1.0, angle1, weight1, nby2);
		
		for (i = 1; i <= nby2; i++) {
			angle[nby2 + i] = angle1[i];
			weight[nby2 + i] = weight1[i];
		}
		free_dvector(angle1,1,nby2);
		free_dvector(weight1,1,nby2);
	
		for (i = 1; i <=n; i++)
			twoaw[i] = 2 * angle[i] * weight[i];
	
		method->b_thinnest = Get_Start_Depth(angle[1],method->b_calc);
		
		@<print angles@>@;

		return;
	}
	
@ Now we need to include three angles, the critical angle, the cone
angle, and perpendicular.  Now the important angles are the ones in
the slab.  So we calculate the cosine of the critical angle in the 
slab and cosine of the cone angle in the slab.

The critical angle will always be greater than the cone angle in the
slab and therefore the cosine of the critical angle will always be 
less than the cosine of the cone angle.  Thus we will integrate from
zero to the cosine of the critical angle (using Gaussian quadrature
to avoid either endpoint) then from the critical angle to the cone
angle (using Radau quadrature so that the cosine angle will be 
included) and finally from the cone angle to 1 (again using Radau
quadrature so that 1 will be included).

@<Gaussian quadrature from 0 to the critical angle@>=
	cos_crit_angle = Cos_Critical_Angle(slab->n_slab,1.0);
	nby3 = n / 3;
	gauleg(0.0, cos_crit_angle, angle, weight, nby3);

@ @<Radau quadrature from the critical angle to the cone angle@>=

    mu = sqrt(slab->n_slab*slab->n_slab-1+slab->cos_angle*slab->cos_angle)/slab->n_slab;
	angle1 = dvector(1,nby3);
	weight1 = dvector(1,nby3);
	Radau(cos_crit_angle, mu, angle1, weight1, nby3);
	for (i = 1; i <= nby3; i++) {
		angle[nby3 + i] = angle1[i];
		weight[nby3 + i] = weight1[i];
	}

@ @<Radau quadrature from the cone angle to 1@>=
	Radau(mu, 1.0, angle1, weight1, nby3);
	for (i = 1; i <= nby3; i++) {
		angle[nby3 * 2 + i] = angle1[i];
		weight[nby3 * 2 + i] = weight1[i];
	}
	free_dvector(angle1,1,nby3);
	free_dvector(weight1,1,nby3);

	for (i = 1; i <=n; i++)
		twoaw[i] = 2 * angle[i] * weight[i];

	method->b_thinnest = Get_Start_Depth(angle[1],method->b_calc);

	@<print angles@>@;


@*2 Initialization.

The basic idea behind diamond initialization is to rewrite the 
time-independent, one-dimensional, az\-i\-muth\-al\-ly av\-er\-aged, 
radiative transport equation
$$
\nu {\partial L(\tau,\nu)\over\partial \tau} + L(\tau,\nu)
    = {a\over2}\int_{-1}^1 h(\nu,\nu')L(\tau,\nu')\, d\nu'
$$
in a discrete form as
$$
\pm\nu_i{\partial L(\tau,\pm\nu_i)\over\partial \tau} + L(\tau,\pm\nu_i)
= {a\over2}\sum_{j=1}^M w_j \left[ 
   h(\nu_i,\nu_j) L(\tau,\pm\nu_i) + h(\nu_i,-\nu_j) L(\tau,\mp\nu_i)\right]
$$
When this equation is integrated over a thin layer from $\tau^*_0$ to $\tau^*_1$
then get
$$\displaylines{\quad
\pm\nu_i[L(\tau^*_1,\pm\nu_i)-L(\tau^*_0,\pm\nu_i)]+d L_{1/2}(\pm\nu_i)
\hfill\cr
\hfill{} = {a\over2}\sum_{j=1}^M w_j d\left[ 
   h(\nu_i,\nu_j) L_{1/2}(\pm\nu_i) + h(\nu_i,-\nu_j) L_{1/2}(\mp\nu_i)\right]
\quad\cr}
$$
where $d=\tau^*_1-\tau^*_0$.
The integrated radiance $L_{1/2}(\nu)$ is
$$
L_{1/2}(\nu)\equiv{1\over\Delta\tau^*} \int_{\tau^*_0}^{\tau^*_1}L(\tau,\nu)\,d\tau
$$

Exactly how this integral is approximated determines the type of initialization.
Wiscombe evaluated a number of initialization methods and found two that were
useful.  These are the infinitesimal generator and the diamond methods.  The
infinitesmal generator initialization makes the approximation
$$
L_{1/2}(-\nu)=L(\tau^*_1,-\nu)
\qquad\qquad
L_{1/2}(\nu)=L(\tau^*_0,\nu)
$$
and the diamond initialization assumes
$$
L_{1/2}(\nu)={1\over2} [L(\tau^*_0,\nu)+L(\tau^*_1,\nu)]
$$

@*2 Infinitesmal Generator Initialization.

@ |Get_IGI_Layer| generates the starting matrix with the inifinitesimal generator method.
The accuracy is $O(d)$ and assumes that the average irradiance upwards is
equal to that travelling downwards at the top and the average radiance upwards
equals that moving upwards from the bottom.
$$
L_{1/2}(-\nu)=L(\tau^*_1,-\nu)
\qquad\qquad
L_{1/2}(\nu)=L(\tau^*_0,\nu)
$$
After manipulation, Wiscombe obtains these
basic formulas for the infinitesimal generator method,
$$
R = {\hat R} d\qquad\qquad T = I - {\hat T}d
$$
where $d$ is the optical thickness of the layer and $I$ is the
identity matrix.  The values for $\hat R$ and $\hat T$ are given by
$$
\hat R = {a\over2} M^{-1} h^{+-} W \qquad\qquad \hat T = M^{-1}(I - {a\over2}h^{++} W)
$$
where $M$ and $W$ are diagonal matrices composed of the quadrature angles and
their corresponding weights.  Therefore
$$
\hat R_{ij} = {a\over2\mu_i} h^{+-}_{ij} w_j\qquad\qquad 
\hat T_{ij} = {\delta_{ij}\over\mu_i}-{a\over2\mu_i} h^{++}_{ij} w_j
$$
and
$$
R_{ij} = {ad\over2\mu_i} h^{+-}_{ij} w_j\qquad\qquad 
T_{ij} = {ad\over2\mu_i} h^{++}_{ij} + \delta_{ij}\left[1-{d\over\mu_i}\right]
$$

This would be fine, but the way that the reflection and transmission matrices are
set-up requires that each we multiply each matrix on the right by $1/(2\mu_jw_j)$.
Putting things together we get
$$
R_{ij} = {ad\over4\mu_i\mu_j} h^{+-}_{ij}
$$
and
$$
T_{ij} = {ad\over4\mu_i\mu_j} h^{++}_{ij} + 
         {\delta_{ij}\over2\mu_i w_i}\left[1-{d\over\mu_i}\right]
$$

@<Definition for |Get_IGI_Layer|@>=
static void Get_IGI_Layer(struct AD_method_type method, double **h, double **R, double **T)
{
	int i, j, n;
	double a, c, d, temp;

	a = method.a_calc;
	d = method.b_thinnest;
	n = method.quad_pts;

	for (j = 1; j <= n; j++) {
		temp = a * d / 4 / angle[j];
		for (i = 1; i <= n; i++) {
			c = temp/angle[i];
			R[i][j] = c*h[i][-j];
			T[i][j] = c*h[i][j];
		}
		T[j][j] += (1-d/angle[j])/twoaw[j];
	}
}

@*2 Diamond Initialization.

It should be noted up front that the implementation contained herein is
somewhat cryptic.  Much of the complexity comes from using the tricks
in the appendix A of Wiscombe's paper (``On initialization, error and flux
conservation in the doubling method.'')  After spending a whole day tracking
down a small error in the calculation of the reflection matrix, I will spend a
few moments trying to improve the documentation for this whole section.
It should be apparent that this is no substitute for reading the paper.

The advantage of the diamond initialization method is that its accuracy is
of the order of the square of the optical thickness $O(d^2)$.  This means that
much thicker starting layers and retain good starting accuracy.  This reduces
the number of doubling steps that are required.  However, if the layer thickness
is too thin then the accuracy gets much worse because errors in the numerical
precision start to affect the results.

|Get_Diamond_Layer| generates the starting matrix with the diamond method.
This implies that the integral can be replaced by a simple average of the radiances 
at the top and bottom of the layer,
$$
L_{1/2}(\nu)={1\over2} [L(\tau^*_0,\nu)+L(\tau^*_1,\nu)]
$$

@<Definition for |Get_Diamond_Layer|@>=
static void Get_Diamond_Layer(struct AD_method_type method, double **h, double **R, double **T)
{
	@<Local variables and initialization@>@;
	@<Find |r| and |t|@>@;
	@<Find |C=r/(1+t)|@>@;
	@<Find |G=0.5(1+t-C r)|@>@;
    @<print |r|, |t|, and |g| for Martin Hammer@>@;
	@<Calculate |R| and |T|@>@;
	@<Free up memory@>@;	
}


@ This diamond initialization
method uses the same $\hat R$ and $\hat T$ as was used for infinitesimal
generator method.  However, we want to form the |r| and |t|
$$
r = {d\over2}\hat R \qquad\qquad t = {d\over 2} \hat T
$$
Recall that
$$
\hat R_{ij} = {a\over2\mu_i} h^{+-}_{ij} w_j\qquad\qquad 
\hat T_{ij} = {\delta_{ij}\over\mu_i}-{a\over2\mu_i} h^{++}_{ij} w_j
$$
therefore
$$
r_{ij} = {adw_j\over4\mu_i} h^{+-}_{ij}\qquad\qquad 
t_{ij} = \delta_{ij}{d\over2\mu_i}-{adw_j\over4\mu_i} h^{++}_{ij}
$$
If you happen to be wondering why right multiplication by
$1/(2\mu_j w_j)$ is not needed, you would be a thinking sort 
of person.  Division by $1/(2\mu_j w_j)$ is not needed until the
final values for |R| and |T| are formed. 

@<Find |r| and |t|@>=

	for (j = 1; j <= n; j++) {
		temp = a * d * weight[j]/4;
		for (i = 1; i <= n; i++) {
			c = temp/angle[i];
			R[i][j] = c*h[i][-j];
			T[i][j] = -c*h[i][j];
		}
		T[j][j] += d/(2*angle[j]);
	}


@  Wiscombe points out (in Appendix A), that the matrix inversions can be avoided
by noting that if we want $C$ from the combination
$$
C= r(I+t)^{-1}
$$
then one needs only solve the system
$$
(I+t)^T C^T = r^T
$$
for $C$.  This is done in the routine |Left_Inverse_Multiply|.
We just need to create $A=I+T$ and fire it off to
|Left_Inverse_Multiply|.  Actually, Wiscome goes on to suggest
a faster method that takes advantage of the column
oriented structure of storage on the computer.  Since
we are using the Numerical Recipes scheme, I don't think
that his refinement will prove faster because it involves
more multiplications and divisions.  (Actually, that improvement
was exactly what the bug in the program was.  I included
the required multiplications and voil\'a! It  worked.)

@<Find |C=r/(1+t)|@>=

	for (i = 1; i <= n; i++) {
		for (j = 1; j <= n; j++) 
			A[i][j] = T[i][j];
		A[i][i] += 1.0;
	}

	Left_Inverse_Multiply(n, A, R, C);

@ Here the matrix 
$$
G = {1\over2}(I+t-C r)
$$
is formed.  

@<Find |G=0.5(1+t-C r)|@>=
	Matrix_Multiply(n, C, R, G);	
	for (i = 1; i <= n; i++) {
		for (j = 1; j <= n; j++) 
			G[i][j] = (T[i][j] - G[i][j])/2;
		G[i][i] += 0.5;
		}

@ To print intermediate results for Chapter 4 of AJ's book, 
then it is necessary to print things from within |Get_Diamond_Layer|.
Martin Hammer requested that I provide these results.  Since
this is the only time that they are of interest, they are only
printed when both the compiler define |MARTIN_HAMMER| is defined,
and when the variable |Martin_Hammer!=0|.

@<print |r|, |t|, and |g| for Martin Hammer@>=

#ifdef MARTIN_HAMMER
{
double **Ginv, **G2;

if (Martin_Hammer!=0) {
	printf("A from equation 5.55\n");
	wrmatrix(n,T);
	
	printf("B from equation 5.55\n");
	wrmatrix(n,R);
	
	Ginv=dmatrix(1,n,1,n);
	G2=dmatrix(1,n,1,n);
	
	for (i = 1; i <= n; i++) {		
		for (j = 1; j <= n; j++) {		
			G2[i][j] = G[i][j]*2.0;
		}
	}
		
	Matrix_Inverse(n,G2,Ginv);
	
	printf("Inverse of G from equation 5.56\n");
	wrmatrix(n,G2);

	printf("G from equation 5.56\n");
	wrmatrix(n,Ginv);

	free_matrix(Ginv,1,n,1,n);
	free_matrix(G2,1,n,1,n);
	}
}
#endif

@ Now we get the part that I really don't understand.  However, I know
that this works.  There are a couple of confusing transposes and 
bizarre incorporation of |twoaw|, but everything hangs together.
Now since the single layer matrices |R| and |T| are the solutions to
the systems of equations
$$
GR=C\qquad\qquad G(t+I) = I
$$
We do the little shuffle and only find the LU decomposition of
|G| once and use it to find both |R| and |T+1|.

@<Calculate |R| and |T|@>=
	Transpose_Matrix(n, G);
	Decomp(n, G, &condition, ipvt);

	if (condition==1e32)
		AD_error("Singular Matrix ... failed in diamond_init\n");
		
	for (i = 1; i <= n; i++) {		
		@<Solve for row of |R|@>@;
		@<Solve for row of |T|@>@;
	}

#ifdef MARTIN_HAMMER
{
double **T2, **Ginv;

if (Martin_Hammer==5) {

	T2=dmatrix(1,n,1,n);
	Ginv=dmatrix(1,n,1,n);
	
	Copy_Matrix(n,T,T2);
	
	for (i = 1; i <= n; i++) {		
		T2[i][i] += 1/twoaw[i];
	}

	for (i = 1; i <= n; i++) {		
		for (j = 1; j <= n; j++) {		
			T2[i][j] *= twoaw[j]*0.5;
		}
	}

	printf("G=(T-1)/2 from equation 5.55\n");
	wrmatrix(n,T2);

	Matrix_Inverse(n,T2,Ginv);
	
	printf("1/G\n");
	wrmatrix(n,Ginv);
	
	free_matrix(T2,1,n,1,n);
	free_matrix(Ginv,1,n,1,n);
	}
}
#endif
	       
@ We use the decomposed form of |G| to find |R|.  Since |G| is
now the LU decomposition of $G^T$, we must pass rows of the
|C| to |Solve| and get rows back.  Note the finess with
$$
\hbox{work}_j = C_{ji}{a_jw_j\over a_iw_i}
$$
To get everything in the right place.  This is discussed in
Wiscombe's appendix.  Finally, we dutifully
put these values back in |R| and divide by $1/(2\mu_j w_j)$
so that |R| will be symmetric and have the proper form.

@<Solve for row of |R|@>=
		for (j = 1; j <= n; j++)	
			work[j] = C[j][i]*twoaw[j]/twoaw[i];		
		Solve(n, G, work, ipvt);
		for (j = 1; j <= n; j++)	
			R[i][j] = work[j]/twoaw[j];

@ We again use the decomposed form of |G| to find |T|.  This is much
simpler since we only need to pass rows of the identity matrix
back and forth.  We again carefully 
put these values back in |T| and divide by $1/(2\mu_j w_j)$
so that |T| is properly formed.  Oh yes, we can't forget to subtract
the identity matrix!

@<Solve for row of |T|@>=
		for (j = 1; j <= n; j++)	
			work[j] = 0;		
		work[i] = 1.0;
		Solve(n, G, work, ipvt);
		for (j = 1; j <= n; j++)	
			T[i][j] = work[j]/twoaw[j];
		T[i][i] -= 1.0/twoaw[i];			/* Subtract Identity Matrix */

@ Pretty standard stuff here.  Allocate memory and 
print a warning if the thickness is too small.

@<Local variables and initialization@>=

	int i, j, n;
	double **A, **G, **C;
	double a, c, d, temp;
	double *work;
	double condition;
	int	*ipvt;
	
	d = method.b_thinnest;
	a = method.a_calc;
	n = method.quad_pts;

	A=dmatrix(1,n,1,n);
	G=dmatrix(1,n,1,n);
	C=dmatrix(1,n,1,n);
	work = dvector(1, n);
	ipvt = ivector(1, n);

	if (d < 1e-4)
		AD_error("**** Roundoff error is a problem--Use IGI method\n");

@
@<Free up memory@>=	
	free_dvector(work,1,n);
	free_ivector(ipvt,1,n);
	free_dmatrix(A,1,n,1,n);
	free_dmatrix(G,1,n,1,n);
	free_dmatrix(C,1,n,1,n);

@*2 Layer Initialization.

@ |Init_Layer| returns reflection and transmission matrices for a thin layer.
Space must previously been allocated for |R| and |T|.

@<Prototype for |Init_Layer|@>=
	void Init_Layer(struct AD_slab_type slab, struct AD_method_type method, double **R, double **T)

@ @<Definition for |Init_Layer|@>=
	@<Prototype for |Init_Layer|@>@;

{
	double **h;
	int n;

	n=method.quad_pts;

	if (slab.b <= 0) {
		Zero_Layer(n, R, T);
		return;
	}

	h = dmatrix(-n, n, -n, n);
	Get_Phi(n, slab.phase_function, method.g_calc, h);

	if (method.b_thinnest < 1e-4 || method.b_thinnest < 0.09 * angle[1])
		Get_IGI_Layer(method, h, R, T);
	else
		Get_Diamond_Layer(method, h, R, T);

	free_dmatrix(h, -n, n, -n, n);
}
