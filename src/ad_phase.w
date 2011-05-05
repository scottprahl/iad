@*1 AD Phase Function.
This section contains all the routines associated with generating
the necessary matrices for Henyey-Greenstein phase functions.
This is the place to put code to implement other phase functions.

\def\shat{{\hat{\bf s}}}
\def\bfh{{\bf h}}
\def\qandq{\qquad{\rm and}\qquad}

@(ad_phase.c@>=
#include <stdlib.h>
#include <math.h>
#include "nr_util.h"
#include "ad_globl.h"
#include "ad_phase.h"

@<Definition for |Get_Phi|@>@;

@ @(ad_phase.h@>=
	@<Prototype for |Get_Phi|@>;

@*2 Redistribution function.
The single scattering phase function $p(\nu)$ for a tissue determines the
amount of light scattered at an angle $\nu=\cos\theta$ from the direction of
incidence.  The
subtended angle $\nu$ is the dot product
of the unit vectors $\shat_i$ and $\shat_j$
$$
\nu=\shat_i\cdot\shat_j=\nu_i\nu_j+\sqrt{1-\nu_i^2}\sqrt{1-\nu_j^2}\cos\phi
$$
where $\shat_i$ is the incident and $\shat_j$ is the scattered light directions

The redistribution function ${\bf h}_{ij}$ determines the fraction of light
scattered from an incidence cone with angle $\nu_i$ into a cone with angle
$\nu_j$.  The redistribution function is calculated by averaging the phase
function over all possible azimuthal angles for fixed angles $\nu_i$ and
$\nu_j$,
$$
h(\nu_i,\nu_j) = {1\over2\pi}
  \int_0^{2\pi} p(\nu_i\nu_j+\sqrt{1-\nu_i^2}\sqrt{1-\nu_j^2}\cos\phi)\,d\phi
$$
Note that the angles $\nu_i$ and $\nu_j$ may also be negative (light
travelling in the opposite direction).  The full redistribution matrix may be
expressed in terms a $2\times2$ matrix of |n|$\times$|n| matrices
$$
\bfh=\left[\matrix{\bfh^{--}&\bfh^{-+}\cr
                   \bfh^{+-}&\bfh^{++}\cr}
\right] 
$$
The first plus or minus sign indicates the sign in front of the incident 
angle and the second is the sign of the direction of the scattered light. 

When the cosine of the angle of incidence or exitance is unity ($\nu_i=1$ or
$\nu_j=1$), then the redistribution function $h(1,\nu_j)$ is equivalent to the phase
function $p(\nu_j)$.  In the case of isotropic scattering, the
redistribution function is a constant
$$
h(\nu_i,\nu_j) = p(\nu) = {1\over4\pi}.
$$
For Henyey-Greenstein scattering, the redistribution function can be expressed
in terms of the complete elliptic integral of the second kind $E(x)$ 
$$
h(\nu_i,\nu_j) = {2\over\pi}{1-g^2\over (\alpha-\gamma)\sqrt{\alpha+\gamma} }
				  \,E\left(\sqrt{2\gamma\over\alpha+\gamma}\,\right)
$$
where $g$ is the average cosine of the Henyey-Greenstein phase function and
$$
\alpha=1+g^2-2 g \nu_i \nu_j 
\qandq
\gamma=2 g \sqrt{1-\nu_i^2}\sqrt{1-\nu_j^2} 
$$
The function $E(x)$ may be calculated using algorithms found in Press {\it et al.\/}
This method of calculating the phase function is slower than the method
that is used in this program.

Other phase functions require numerical integration of the phase
function.  If the phase function is highly anisotropic, then the
integration over the azimuthal angle is  particularly difficult and care
must be taken to ensure that the integration is accurate.   This is
important because errors in the redistribution function enter directly
into the reflection and transmission matrices for thin layers.  Any
errors will be doubled with each successive addition of layers and small
errors will rapidly increase.

@ An alternate way to calculate the redistribution function is the
$\delta$--$M$ method of Wiscombe.  This method works especially
well for highly anisotropic phase functions.  The number of quadrature
points is specified by $M$.  The $\delta$--$M$ method approximates the
true phase function  by a phase function consisting of a Dirac delta
function and $M-1$ Legendre polynomials
$$
p^*(\nu)= 2 g^M\delta(1-\nu) + (1-g^M) \sum_{k=0}^{M-1} (2k+1)\chi_k^* P_k(\nu)
$$
where 
$$
\chi_k^*={\chi_k-g^M\over 1-g^M}
	\qandq
\chi_k = {1\over2}\int_0^1 p(\nu) P_k(\nu) \,d\nu
$$
When the $\delta$--$M$ method substitutes $p^*(\nu)\rightarrow p(\nu)$, 
then both the albedo and optical thickness must also be changed,
$a^*\rightarrow a$ and $\tau^*\rightarrow\tau$.  This approximation is
analogous to the similarity transformation often used to improve the
diffusion approximation by moving a part ($g^M$) of the scattered light
into the unscattered component.  The new optical
thickness and albedo are
$$
\tau^*=(1-ag^M)\tau  \qandq
a^* = a {1-g^M\over1-ag^M}
$$
This is equivalent transforming the scattering coefficient as
$\mu_s^* = \mu_s(1-g^M)$. The redistribution function can now be written
as
$$
h^*(\nu_i,\nu_j) = \sum_{k=0}^{M-1} (2k+1)\chi_k^* P_k(\nu_i)P_k(\nu_j)
$$
For the special case of a Henyey-Greenstein phase function,
$$
\chi_k^*={g^k-g^M\over1-g^M}.
$$

@ Calculate the renormalization matrix for a Henyey-Greenstein phase
function using the delta-M method.  This version has been optimized 
for isotropic and Henyey-Greenstein phase functions.  

@ @<Prototype for |Get_Phi|@>=
void Get_Phi(int n, int phase_function, double g, double **h)

@ @<Definition for |Get_Phi|@>=
	@<Prototype for |Get_Phi|@>
{
	@<Local variables for |Get_Phi|@>@;
	@<Test for bad calling parameters@>@;
	@<Initialize the phase function matrix@>@;
	@<We're done if phase function is isotropic@>@;
	@<Calculate the quadrature coefficients@>@;
	@<Create Legendre Polynomial matrix@>@;
	@<Calculate the coefficients@>@;
	@<Add the symmetric part of the matrix@>@;
	@<Free |p| and |chi|@>@;
}

@ @<Local variables for |Get_Phi|@>=
    int i, j, k;
    double g2M, gk, x;
    double *chi;
    double **p;

@ @<Test for bad calling parameters@>=
    if (g!=0 && phase_function != HENYEY_GREENSTEIN)
	AD_error("Only the Henyey-Greenstein phase function has been implemented\n");

    if (fabs(g) >= 1)
	AD_error("Get_Phi was called with a bad g_calc value");

@ @<Initialize the phase function matrix@>=
    for (i = -n; i <= n; i++)
	for (j = -n; j <= n; j++)
	    h[i][j] = 1;

    /* zero the zero column and zero row */
    for (i = -n; i <= n; i++) {
	h[i][0] = 0.0;
	h[0][i] = 0.0;
    }

@ @<We're done if phase function is isotropic@>=
    if (g == 0) return;

@ To avoid extra calculation let's define
$$
\hbox{|chi[k]|}\,\equiv (2k+1)\chi_k^*
$$
This will slighly simplify things later on

@<Calculate the quadrature coefficients@>=
	chi = dvector(1, n);
	g2M = pow(g, n);
	gk = 1.0;
	for (k = 1; k < n; k++) {
	    gk *= g;
	    chi[k] = (2 * k + 1) * (gk - g2M) / (1 - g2M);
	}

@ Allocate the matrix for the Legendre values 
this is {\it much\/} more efficient than calculating them
as they are needed.  Since the Legendre polynomial
$P_n(x)$ is generated using recurrence relations, 
all Legendre polynomials $P_k(x)$, where $0\le k\le n$
must also be calculated.  Now the formula
$$
h^*(\nu_i,\nu_j) = \sum_{k=0}^{n-1} (2k+1)\chi_k^* P_k(\nu_i)P_k(\nu_j)
$$
requires all those to be found as well.  There are
$2n+1$ values that must be calculated for $-\mu_n\ldots0\ldots\mu_n$
different arguments.
A simple way is just to put all of the necessary values
in a two-dimensional array and define |p[i][j]==|$\,P_i(\mu_j)$.

@<Create Legendre Polynomial matrix@>=
	@<Allocate the polynomial matrix@>@;
	@<Fill in all the unique values@>@;
	@<Fill in the symmetric values@>@;

@ It is not at all clear that zeroing is needed.
@<Allocate the polynomial matrix@>=
	p = dmatrix(0, n, -n, n);

@ Here I use the recurrence relation
$$
P_{k+1}(\mu_j) = {(2k + 1) x P_k(\mu_j) - k P_{k-1}(\mu_j) \over k + 1}
$$
(which should be stable) to find all the values for all the
positive angles.

@<Fill in all the unique values@>=
	for (j = 1; j <= n; j++) {
	    p[0][j] = 1;
	    x = angle[j];
	    p[1][j] = x;
	    for (k = 1; k < n; k++)
		p[k + 1][j] = ((2 * k + 1) * x * p[k][j] - k * p[k - 1][j]) / (k + 1);
	}

@ I make use of the fact that
$$
P_k(-\nu_j) = (-1)^k P_k(\nu_j)
$$
to fill in all the negative angles in the 
phase function matrix.  This eliminates half
the calculation.  I do two at a time.  This way
there does not need to be a flag.  Since I know that
the dimension of the matrix will be even, this should
not be a problem.  If the matrix is not then you have
problems.

@<Fill in the symmetric values@>=
	for (j = 1; j <= n; j++)
	    for (k = 1; k < n; k++) {
		p[k][-j] = -p[k][j];
		k++;
		p[k][-j] = p[k][j];
	    }

@ Just a straightforward calculation of
$$
h^*(\nu_i,\nu_j) = \sum_{k=0}^{n-1} (2k+1)\chi_k^* P_k(\nu_i)P_k(\nu_j)
$$
and since $\chi_0^*=1$ and $P_0(x)=1$ this is
$$
h^*(\nu_i,\nu_j) = 1+ \sum_{k=1}^{n-1} (2k+1)\chi_k^* P_k(\nu_i)P_k(\nu_j)
$$
Since $h$ has many symmetries, there are only about $n^2/4$ unique
entries.  We only need to calculate those.  Oh yeah, recall that
|chi[k]| includes the factor |2k+1| for speed.

@<Calculate the coefficients@>=
	for (i = 1; i <= n; i++) {
	    for (j = i; j <= n; j++) {
		for (k = 1; k < n; k++) {
		    h[i][j] += chi[k] * p[k][i] * p[k][j];
		    h[-i][j] += chi[k] * p[k][-i] * p[k][j];
			}
	    }
	}

@ Several symmetries in the redistribution matrix are used.
to fill in some entries that begin with a negative angle
$$
h(-\nu_i,\nu_j) = h(\nu_j,-\nu_i)
$$
and secondly
$$
h(-\nu_i,-\nu_j) = h(\nu_j,\nu_i)
$$
Next, some entries along the diagonal are filled in using
$$
h(-\nu_i,-\nu_i) = h(\nu_i,\nu_i)
$$
Finally, the lower triangle is filled in using the values
from the upper half using
$$
h(\nu_i,\nu_j) = h(\nu_j,\nu_i)
$$
This could probably be more elegant, but it hurts my brain
to think about it.  This works and should take advantage of
all the symmetries present.

@<Add the symmetric part of the matrix@>=
	for(i=n;i>=2;i--)
		for (j=1; j<i; j++) {
			h[-i][j] = h[-j][i];
			h[-i][-j] = h[j][i];
			}

	for(i=1;i<=n;i++)
		h[-i][-i]=h[i][i];
		
	for (i = -n; i <= n; i++)
	    for (j = i + 1; j <= n; j++)
		h[j][i] = h[i][j];

@ @<Free |p| and |chi|@>=
	free_dmatrix(p, 0, n, -n, n);
	free_dvector(chi, 1, n);