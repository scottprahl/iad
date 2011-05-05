@*1 AD Radau Quadrature.

This global variable is needed because the degree of the
Legendre Polynomial must be known.  The routine |Radau| stores
the correct value in this.

@d NSLICES     512
@d EPS 		1e-16

@(ad_radau.c@>=
@h
	#include "ad_globl.h"
	#include "ad_radau.h"
	#include "nr_rtsaf.h"
	#include "nr_util.h"
	#include "nr_zbrak.h"

    static int local_n_size;

	@<Prototype for |Pn_and_Pnm1|@>;
	@<Prototype for |Pnd|@>;
	@<Prototype for |phi|@>;
	@<Prototype for |phi_and_phiprime|@>;

	@<Definition for |Pn_and_Pnm1|@>@;
	@<Definition for |Pnd|@>@;
	@<Definition for |phi|@>@;
	@<Definition for |phi_and_phiprime|@>@;
	@<Definition for |Radau|@>@;

@ @(ad_radau.h@>=
	@<Prototype for |Radau|@>;

@*2 Introduction.

The adding-doubling method is based on numerical integration
of functions using quadrature,
$$
\int_0^1 f(\nu,\nu') \, d\nu' = \sum_{k=1}^{N} w_k f(x_k)
$$
The values of the quadrature points $x_k$ and the weights
$w_k$ are chosen in such a way that the integral is evaluated
exactly for a polynomial of order $2N-1$ (or possibly $2N-2$ depending
on the quadrature method).  Using $N$ quadrature points (Gaussian)
is equivalent to the spherical harmonic method of order $P_{N-1}$,
i.e. four quadrature points corresponds to the $P_3$
method.  The specific choice of quadrature methods for samples
with mismatched boundaries is described in the next section.

Total internal reflection causes problems by changing the effective
range of integration. Usually, adding-doubling integrals range from
$0$ to $1$, since the angle varies from ${\pi\over2}$ to $0$ and
therefore the cosine varies from $0$ to $1$. The integrations are
calculated using numerical quadrature, and the quadrature angles are
optimized for this range. If the cosine of the critical angle is
denoted by $\nu_c$ for a boundary layer with total internal
reflection, then the effective range of integration is reduced to
$\nu_c$ to $1$ (because the rest of the integration range is now
zero). To maintain integration accuracy, the integral is broken into
two parts and each is evaluated by quadrature over the specified
subrange,
$$
\int_0^1A(\nu,\nu') B(\nu',\nu'')\,d\nu'  =
          \int_0^{\nu_c}A(\nu,\nu') B(\nu',\nu'')\,d\nu' +
             \int_{\nu_c}^1 A(\nu,\nu') B(\nu',\nu'')\,d\nu' .
$$
Here $A(\nu,\nu')$ and $B(\nu,\nu')$ represent reflection or
transmission functions, and clearly if either is identically
zero for values of $\nu$ less than $\nu_c$, the integration range is
reduced. The calculations in this paper used Gaussian quadrature
for the range from $0$ to $\nu_c$, thereby avoiding
calculations at both endpoints (in particular, the angle $\nu=0$ is
avoided, which may cause division by zero). 
Radau quadrature
is used for the range from $\nu_c$ to $1$, so $\nu=1$ could be 
specified as a quadrature point.  Each part of the
integration range gets half of the quadrature points; when no critical
angle exists, Radau quadrature is used over the entire range.

Radau quadrature requires finding the $n$ roots of the following equation
$$
P_{n-1}(x_i) + {x_i-1\over n} P_{n-1}'(x_i) =0
$$
Here $P_n(x)$ is the $n$th Legendre polynomial of order zero and
$P_{n-1}'(x_i)$ is the first derivative of the $n-1$ Legendre polynomial.
These roots are the required quadrature points for the integration range -1 to
1.   The $n$th integration angle $\nu_n$ corresponds with $x_n=-1$ (normal
incidence).

@*2 Radau.
|Radau| calculates the |n| quadrature points $x_i$ and weights $w_i$.

@<Prototype for |Radau|@>=
void Radau(double x1, double x2, double *x, double *w, int n)

@ @<Definition for |Radau|@>=
	@<Prototype for |Radau|@>
{
	
	x[n] = -1.0;
	w[n] = 2.0 / (n * n);

	switch (n) {
		case  2: @<Values for |n==2|@>@;
		case  4: @<Values for |n==4|@>@;
		case  8: @<Values for |n==8|@>@;
		case 16: @<Values for |n==16|@>@;
		default: @<Values for arbitrary |n|@>@;
	}
	@<Scale values@>@;
}

@ The code to scale values is easy.  Radau quadrature is
defined over the range -1 to 1.  Here we just linearly scale
the width of each interval and weight as appropriate.
To modify for the range $\nu_c$ to $1$ the following relations are
needed to find the necessary integration angles  $\nu_i$ and weights $w_i$ 
$$
\nu_i = {1+\nu_c - (1-\nu_c) x_i \over 2}
$$
and
$$
w_i=   {1-\nu_ c\over (1-x_i) \sqrt{P_{n-1}'(x_i)}}
$$

@<Scale values@>=
{
	double xm, xl;
	int i;
	
	xm = (x2 + x1) * 0.5;
	xl = (x2 - x1) * 0.5;

	for (i = 1; i <= n; i++) {
		x[i] = xm - xl * x[i];
		w[i] = xl * w[i];
	}
	
}


@ Here is the method for finding Radau quadrature points for
non-tabulated values.

@<Values for arbitrary |n|@>=
{
	int i, nb, ndiv;
	double z;
	double *xb1, *xb2;

	@<Allocate memory for Radau@>@;
	@<Bracket roots@>@;
	@<Find roots and weights@>@;
	@<Free memory for Radau@>@;
	break;
	}

@ @<Allocate memory for Radau@>=
	xb1 = dvector(1,NSLICES);
	xb2=dvector(1,NSLICES);
	
@ Bracket |n-1| roots, double |ndiv| if not enough roots are found.
@<Bracket roots@>=
	local_n_size = n;

	if (2*n>NSLICES) 
		ndiv = NSLICES;
	else
		ndiv = 2*n;

	do{
		nb = n - 1;
		zbrak(phi, -1.0, 1.0, ndiv, xb1, xb2, &nb);
		ndiv *= 2;
	}
	while (nb < n - 1 && ndiv <= NSLICES);

	if (nb < n-1) AD_error("Cannot find enough roots for Radau quadrature");

@ Find the roots with an accuracy |EPS| and store them in the array |x|.
Put them in backwards so that |x[n]=-1| is in the correct spot.

@<Find roots and weights@>=
	for (i = 1; i < n; i++) {
		double tmp;
		z = rtsafe(phi_and_phiprime, xb1[i], xb2[i], EPS);
		x[n-i] = z;
		tmp = Pnd(n-1,z);
		w[n-i] = 1 / ((1 - z) * tmp * tmp);
	}

@	@<Free memory for Radau@>=
	free_dvector(xb1,1,NSLICES);
	free_dvector(xb2,1,NSLICES);

@ |Pn_and_Pnm1| returns $P_n(x)$ and $P_{n-1}(x)$

@<Prototype for |Pn_and_Pnm1|@>=
static void Pn_and_Pnm1(int n, double x, double *Pnm1, double *Pn)

@ @<Definition for |Pn_and_Pnm1|@>=
	@<Prototype for |Pn_and_Pnm1|@>
{
	int k;
	double Pk, Pkp1;
	double Pkm1=1.0;

	*Pnm1=1.0;
	*Pn=1.0;
	if (x >= 1.0)
		return;
	
	if (x<=-1.0)
		x=-1;

	Pk =  x;
	
	for (k=1; k<n; k++) {
		Pkp1 = ((2*k+1) * x * Pk - k * Pkm1) / (k+1);
		Pkm1 = Pk;
		Pk = Pkp1;
	}
	
	*Pnm1=Pkm1;
	*Pn=Pk;
}

@ To calculate the weights for the quadrature points we
need to evaluate the first derivative of the Legendre polynomial.
To do this we use a recurrence relation given by  H. H. Michels, 
in ``Abscissas and weigh coefficients for
Lobatto quadrature,'' {\it Math Comp\/}, {\bf 17}, 237-244 (1963).

@<Prototype for |Pnd|@>=
static double Pnd(int n, double x)

@ @<Definition for |Pnd|@>=
	@<Prototype for |Pnd|@>
{
	double p, pminus, pplus;
	int i;

	if (x > 1.0){			
		x=1;
	}
	else if (x<-1.0){
		x=-1;
	}

	pminus = 0;
	p =  1;

	if (n <= 0)
		return pminus;
		
	for (i=1; i<n; i++) {
		pplus = ((2*i + 1) * x * p - (i + 1) * pminus)/i;
		pminus = p;
		p = pplus;
	}
	return p;
}

@ To use Newton's method to find the roots of
$$
\phi_{n-1}(x) = {P_{n-1}(x)+P_n(x) \over 1+x}
$$
we need to find the derivative.  This is 
$$
\phi_{n-1}'(x) = {P_{n-1}'(x)+P_n'(x) \over 1+x}-{P_{n-1}(x)+P_n(x) \over (1+x)^2}
$$
Now we can use our recurrence relation
$$
(1-x^2) P_{n-1}'(x) = nxP_{n-1}(x)-nP_n(x)
$$
To eliminate the derivative terms in the above equation
to get
$$
\phi_{n-1}'= {(nx+x-1)P_{n-1}(x)+(nx+2x-n-1)P_n(x)-(n+1)P_{n+1}(x) \over (1-x)(1+x)^2}
$$
The higher order Legendre Polynomial can be eliminated using
$$
(n+1)P_{n+1}(x) = (2n+1) x P_n(x) - n P_{n-1}(x)
$$
to get
$$
\phi_{n-1}'(x) = {(nx+x+n-1)P_{n-1}(x)+(-nx+x-n-1)P_n(x)\over (1-x)(1+x)^2}
$$
And therefore we just call the routine that will return $P_n(x)$ and
$P_{n-1}(x)$ and multiply by the appropriate factors to obtain both
terms.

The only problem is when $x=1$ or $x=-1$.  Then we get this spurious 
division by zero.  So we special case these and evaluate them elsewhere.

@<Prototype for |phi_and_phiprime|@>=
static void phi_and_phiprime(double x, double *phi, double *phiprime)

@ @<Definition for |phi_and_phiprime|@>=
	@<Prototype for |phi_and_phiprime|@>
{
double Pn, Pnm1;
int n;

	n = local_n_size;
	if (x >= 1.0){
		@<Phi and phiprime at |x=1|@>@;
	}
	else if (x<=-1.0) {
		@<Phi and phiprime at |x=-1|@>@;
	}
	else {
		Pn_and_Pnm1(n, x, &Pnm1, &Pn);
		*phi = (Pn+Pnm1)/(1+x);
		*phiprime = ((n*x-1+x+n)*Pnm1+(-n*x+x-n-1)*Pn)/(1+x)/(1+x)/(1-x);
	}
}

@ To find $\phi(1)$ and $\phi'(1)$ we need to recall a few facts
about Legendre polynomials.  First,
$$
P_n(1)=1
$$
Therefore
$$
\phi(1)=1
$$
The value of the first derivative is somewhat trickier.
Recall that
the Legendre polynomials are solutions to
$$
(1-x^2)P_n''(x) - 2x P_n'(x) + n(n+1) P_n(x) = 0
$$
Now if $x=1$ then the first term on the left hand side will be
zero.  Therefore
$$
P_n'(1) =  {n(n+1)\over 2}
$$
Therefore
$$
\phi_{n-1}'(1) =  {n^2-1\over 2}
$$

@<Phi and phiprime at |x=1|@>=
{
	*phi = 1;
	*phiprime = (n*n-1)/2;
	}
	
@ To evaluate $\phi(-1)$ we must return to the original
definition, i.e. 
So $$
\phi_{n-1}(x) = P_{n-1}(x)+{x-1\over n}P'_{n-1}(x)
$$
To evaluate this
we need to remember some stuff, namely that
$$
P_n(-x) = (-1)^n P_n(x) \qquad \hbox{\rm so}\qquad P_n(-1) = (-1)^n
$$
The value of the first derivative is
again obtained from the differential equation and
$$
P_n'(-1) =  -{n(n+1)\over 2} P_n(-1) = (-1)^{n+1} {n(n+1)\over 2}
$$
Now we just substitute to get
$$
\phi_{n-1}(-1) = (-1)^{n-1}\cdot n
$$
The first derivative is more diffficult. Mathematica says that it is
$$
\phi_{n-1}'(-1) = (-1)^n {n(1-n^2)\over 4}
$$

@<Phi and phiprime at |x=-1|@>=
	*phi = n;
	*phiprime = -n*(1-n*n)/4;
	if (n%2!=1) {
		*phi *= -1;
		*phiprime *= -1;
	}

@ For Radau quadrature, we want to find the $n-1$ roots
of 
$$
\phi_{n-1}(x) = P_{n-1}(x)+{x-1\over n}P'_{n-1}(x)
$$
F. B. Hildebrand notes that by using a recurrence formula
this becomes
$$
\phi_{n-1}(x) = {P_{n-1}(x)+P_n(x) \over 1+x}
$$
This is particularly convenient, because we must find
$P_{n-1}(x)$ before we can find $P_n(x)$ and this is
exactly what |Pn_and_Pnm1| does.

It is noteworthy that this routine uses the recurrence
formula
$$
P_{n+1}(x) = {(2n+1) x P_n(x) - n P_{n-1}(x)\over n+1}
$$
to calculate the Legendre polynomial $P_n(x)$.
This recurrence relation is given
in H. H. Michels, ``Abscissas and weight coefficients for
Lobatto quadrature,'' {\it Math Comp\/}, {\bf 17}, 237-244 (1963).

@<Prototype for |phi|@>=
static double phi(double x)

@ @<Definition for |phi|@>=
	@<Prototype for |phi|@>
{
double Pn, Pnm1;

	if (x<=-1.0) {
		if (local_n_size%2 != 1) 
			return(-local_n_size);
		else 
			return(local_n_size);
	}
	
	Pn_and_Pnm1(local_n_size, x, &Pnm1, &Pn);
	return((Pn+Pnm1)/(1+x));
}

@*2 Radau Tables.

Here is a selection of commonly used number of quadrature points.

@ @<Values for |n==2|@>=
	x[1]= 0.3333333333333334;
	w[1]= 1.5000000000000000;
	break;

@ @<Values for |n==4|@>=
	x[3]=-0.5753189235216942;
	x[2]= 0.1810662711185306;
	x[1]= 0.8228240809745921;
	
	w[3]= 0.6576886399601182;
	w[2]= 0.7763869376863437;
	w[1]= 0.4409244223535367;
	break;

@ @<Values for |n==8|@>= 
	x[7]=-0.8874748789261557;
	x[6]=-0.6395186165262152; 
	x[5]=-0.2947505657736607;
	x[4]= 0.0943072526611108; 
	x[3]= 0.4684203544308211;
	x[2]= 0.7706418936781916;
	x[1]= 0.9550412271225750; 

	w[7]= 0.1853581548029793; 
	w[6]= 0.3041306206467856; 
	w[5]= 0.3765175453891186; 
	w[4]= 0.3915721674524935; 
	w[3]= 0.3470147956345014; 
	w[2]= 0.2496479013298649; 
	w[1]= 0.1145088147442572; 
	break;

@ @<Values for |n==16|@>=
	x[15]=-0.9714610905263484;
	x[14]=-0.9054008198116666;
	x[13]=-0.8045734013587561;
	x[12]=-0.6728619212112202;
	x[11]=-0.5153294780626855;
	x[10]=-0.3380303900599197;
	x[ 9]=-0.1477783218133717;
	x[ 8]= 0.0481153830735303;
	x[ 7]= 0.2421226227060438;
	x[ 6]= 0.4267878274849459;
	x[ 5]= 0.5950144898997919;
	x[ 4]= 0.7403379488928179;
	x[ 3]= 0.8571740937696823;
	x[ 2]= 0.9410354027041150;
	x[ 1]= 0.9887186220549766;

	w[15]= 0.0477022269476863;
	w[14]= 0.0839852814449645;
	w[13]= 0.1170203531038591;
	w[12]= 0.1455555452202026;
	w[11]= 0.1684963978499219;
	w[10]= 0.1849617814886653;
	w[10]= 0.1849617814886653;
	w[ 9]= 0.1943190897115679;
	w[ 8]= 0.1962087882390318;
	w[ 7]= 0.1905582942553547;
	w[ 6]= 0.1775847927527395;
	w[ 5]= 0.1577869218042020;
	w[ 4]= 0.1319256999330681;
	w[ 3]= 0.1009956796217840;
	w[ 2]= 0.0661895086101364;
	w[ 1]= 0.0288971390168143;
	break;
