@*1 AD Double.
This has the routines needed to add layers together in various combinations.

\def\bfr{{\bf R}}
\def\bft{{\bf T}}
\def\bfe{{\bf E}}
\def\bfc{{\bf c}}
\def\bfL{{\bf L}}
\def\bfJ{{\bf J}}

@(ad_doubl.c@>=
#include <math.h>
#include <float.h>
#include "nr_util.h"
#include "ad_matrx.h"
#include "ad_globl.h"
#include "ad_doubl.h"

	@<Definition for |Star_Multiply|@>@;
	@<Definition for |Star_One_Minus|@>@;
	@<Definition for |Basic_Add_Layers|@>@;
	@<Definition for |Basic_Add_Layers_With_Sources|@>@;
	@<Definition for |Add|@>@;
	@<Definition for |Add_With_Sources|@>@;
	@<Definition for |Add_Homogeneous|@>@;
	@<Definition for |Double_Once|@>@;
	@<Definition for |Double_Until|@>@;
	@<Definition for |Double_Until_Infinite|@>@;
	@<Definition for |Between|@>@;

@ @(ad_doubl.h@>=
    @<Prototype for |Add|@>;
	@<Prototype for |Add_With_Sources|@>;
	@<Prototype for |Add_Homogeneous|@>;
	@<Prototype for |Double_Once|@>;
	@<Prototype for |Double_Until|@>;
	@<Prototype for |Double_Until_Infinite|@>;
	@<Prototype for |Between|@>;

@*2 Basic Routine to Add Layers Without Sources.

The basic equations for the adding-doubling method (neglecting sources) are
$$
\bft^{02} = \bft^{12}(\bfe-\bfr^{10}\bfr^{12})^{-1}\bft^{01}
$$
$$
\bfr^{20} = \bft^{12}(\bfe-\bfr^{10}\bfr^{12})^{-1} \bfr^{10}\bft^{21}+\bfr^{21} 
$$
$$
\bft^{20} = \bft^{10}(\bfe-\bfr^{12}\bfr^{10})^{-1}\bft^{21}
$$
$$
\bfr^{02} = \bft^{10}(\bfe-\bfr^{12}\bfr^{10})^{-1} \bfr^{12}\bft^{01}+\bfr^{01} 
$$
Upon examination it is clear that the two sets of equations have
the same form.  Therefore if I implement the first two equations, then
the second set can be obtained by suitable switching of the parameters.
Furthermore, these equations assume some of the multiplications are 
star multiplications. Explicitly,
$$
\bft^{02} = \bft^{12}(\bfe-\bfr^{10}\star\bfr^{12})^{-1}\bft^{01}
$$
and
$$
\bfr^{20} = \bft^{12}(\bfe-\bfr^{10}\star\bfr^{12})^{-1} \bfr^{10}\star\bft^{21}+\bfr^{21} 
$$
where the identity matrix $\bfe$ is then
$$
\bfe^{ij}={1\over2\mu_i w_i} \delta_{ij}
$$
where $\delta_{ij}$ is the usual Kronecker delta.  It is noteworthy that if
say $R^{10}\equiv0$, then $\bfe^{-1}\equiv {\bf c}$ and so
$$
\bft^{02} = \bft^{12}\bfc\bft^{01}=\bft^{12}\star\bft^{01}
$$

One goal of this routine was to make it efficient and easy
to use.  It is possible to call this routine with the same pointer for all
the different reflection matrices and the pointer for
the transmission matrices may be the same also.  (The reflection 
and transmission pointers may need to be distinct.  The temporary memory
pointers |a| and |b| must be distinct from each other and distinct from the reflection and
transmission matrices.)

Note:  it should be possible to eliminate the need for the matrix
|b| if |Inverse_Multiply| could be called with an argument list
like |Inverse_Multiply(n,A,B,A)|.  A quick glance at the code suggests
that this would just force the allocation of the matrix into the 
|Inverse_Multiply| routine and no net gain would result.

@<Definition for |Basic_Add_Layers|@>=

static void Basic_Add_Layers(int n,
	double **R10, double **T01, 
	double **R12, double **R21, double **T12, double **T21,
	double **R20, double **T02,
	double **a, double **b)
{
  Star_Multiply(n, R10, R12, a);/* |a=|$\,\bfr^{10}\star\bfr^{12}$ */
  Star_One_Minus(n, a);/* |a=|$\,\bfe-\bfr^{10}\star\bfr^{12}$ */
  Left_Inverse_Multiply(n, a, T12, b);	
          /* |b=|$\,\bft^{12}(\bfe-\bfr^{10}\bfr^{12})^{-1}$ */
  

  Matrix_Multiply(n, b, R10, a);
          /* |a=|$\,\bft^{12}(\bfe-\bfr^{10}\star\bfr^{12})^{-1}\bfr^{10}$ */
  Star_Multiply(n, a, T21, a);
          /* |a=|$\,\bft^{12}(\bfe-\bfr^{10}\star\bfr^{12})^{-1}\bfr^{10}\star\bft^{21}$ */
  Matrix_Sum(n, R21, a, R20);
  Copy_Matrix(n, T01, a);
  Matrix_Multiply(n, b, a, T02);
}

@*2 Basic Routine to Add Layers With Sources.

The adding-doubling equations including source terms $\bfJ$ are identical to 
those given above for the reflection and transmission.  The only difference is
that the source terms must be kept track of separately according to
$$
\bfJ^{02}_+ = \bfJ^{12}_+ + 
\bft^{12}(\bfe-\bfr^{10}\bfr^{12})^{-1}( \bfJ^{01}_+ +\bfr^{10} \bfJ^{21}_-)
$$
and
$$
\bfJ^{20}_+ = \bfJ^{10}_- + 
\bft^{10}(\bfe-\bfr^{12}\bfr^{10})^{-1}( \bfJ^{21}_- +\bfr^{12} \bfJ^{01}_+)
$$
where the $+$ subscript indicates the downward direction and $-$ indicates the
upward direction.  Note that these subscripts are not needed.  Thus we have
$$
\bfJ^{02} = \bfJ^{12} + 
\bft^{12}(\bfe-\bfr^{10}\bfr^{12})^{-1}( \bfJ^{01} +\bfr^{10} \bfJ^{21})
$$
and
$$
\bfJ^{20} = \bfJ^{10} + 
\bft^{10}(\bfe-\bfr^{12}\bfr^{10})^{-1}( \bfJ^{21} +\bfr^{12} \bfJ^{01})
$$

Again, it is apparent that clever switching of the arguments requires that
only one set of equations needs to be calculated.  
These equations assume some of the multiplications are 
star multiplications. Explicitly,
$$
\bfJ^{02} =  \bfJ^{12} +
\bft^{12}(\bfe-\bfr^{10}\star\bfr^{12})^{-1} (\bfJ^{01}+\bfr^{10}\star\bfJ^{21}) 
$$

@<Definition for |Basic_Add_Layers_With_Sources|@>=

static void Basic_Add_Layers_With_Sources(int n,
	double **R10, double **T01, 
	double **R12, double **R21, double **T12, double **T21,
	double **R20, double **T02, 
	double **J01, double **J12, double **J21, double **J02,
	double **a, double **b)
{
  Star_Multiply(n, R10, R12, a);/* |a=|$\,\bfr^{10}\star\bfr^{12}$ */
  Star_One_Minus(n, a);/* |a=|$\,\bfe-\bfr^{10}\star\bfr^{12}$ */
  Left_Inverse_Multiply(n, a, T12, b);	
          /* |b=|$\,\bft^{12}(\bfe-\bfr^{10}\bfr^{12})^{-1}$ */
  

  Matrix_Multiply(n, b, R10, a);
          /* |a=|$\,\bft^{12}(\bfe-\bfr^{10}\star\bfr^{12})^{-1}\bfr^{10}$ */
  Star_Multiply(n, a, T21, a);
          /* |a=|$\,\bft^{12}(\bfe-\bfr^{10}\star\bfr^{12})^{-1}\bfr^{10}\star\bft^{21}$ */
  Matrix_Sum(n, R21, a, R20);
  Copy_Matrix(n, T01, a);
  Matrix_Multiply(n, b, a, T02);
  
  Star_Multiply(n,R10,J21,a);      /* |a=|$\,\bfr^{10}\star\bfJ^{21}$ */
  Matrix_Sum(n,J01,a,a);              /* |a=|$\,\bfJ^{01}+\bfr^{10}\star\bfJ^{21}$ */
  Matrix_Multiply(n,b,a,J02);      
  /* |J02=|$\,\bft^{12}(\bfe-\bfr^{10}\star\bfr^{12})^{-1}
  (\bfJ^{01}+\bfr^{10}\star\bfJ^{21})$ */
  Matrix_Sum(n,J02,J12,J02);
}

@*2 Higher level routines.

@
@<Prototype for |Add|@>=
void Add(int n, 
	double **R01, double **R10, double **T01, double **T10, 
	double **R12, double **R21, double **T12, double **T21, 
	double **R02, double **R20, double **T02, double **T20)

@ |Add| returns the reflection and transmission matrices
for two different layers added together.  These matrices
do not have to be homogeneous.  The output matrices
|R20|, |R02|, |T20|, and |T02| should be distinct from the
input matrices.

@<Definition for |Add|@>=
	@<Prototype for |Add|@>
{
  @<Allocate memory for |a| and |b|@>@;

  Basic_Add_Layers(n, R10, T01, R12, R21, T12, T21, R20, T02, a, b);
  Basic_Add_Layers(n, R12, T21, R10, R01, T10, T01, R02, T20, a, b);
  
  @<Free Memory for |a| and |b|@>@;
}

@
@<Prototype for |Add_With_Sources|@>=
void Add_With_Sources(int n, 
	double **R01, double **R10, double **T01, double **T10, double **J01, double **J10,
	double **R12, double **R21, double **T12, double **T21, double **J12, double **J21, 
	double **R02, double **R20, double **T02, double **T20, double **J02, double **J20)

@ |Add_With_Sources| returns the reflection and transmission matrices
for two different layers added together.  These matrices
do not have to be homogeneous.  The output matrices
|R20|, |R02|, |T20|, |T02|, |J20|, and |J02| should be distinct from the
input matrices.

@<Definition for |Add_With_Sources|@>=
	@<Prototype for |Add_With_Sources|@>
{
  @<Allocate memory for |a| and |b|@>@;

 Basic_Add_Layers_With_Sources(n, R10, T01, R12, R21, T12, T21, R20, T02, J01,J12,J21,J02,a, b);
 Basic_Add_Layers_With_Sources(n, R12, T21, R10, R01, T10, T01, R02, T20, J21,J10,J01,J20, a, b);
  
  @<Free Memory for |a| and |b|@>@;
}

@
@<Prototype for |Add_Homogeneous|@>=
void Add_Homogeneous(int n, 
	double **R01, double **T01, 
	double **R12, double **T12, 
	double **R02, double **T02)

@
@<Definition for |Add_Homogeneous|@>=
	@<Prototype for |Add_Homogeneous|@>
{
  @<Allocate memory for |a| and |b|@>@;

  Basic_Add_Layers(n, R01, T01, R12, R12, T12, T12, R02, T02, a, b);
  
  @<Free Memory for |a| and |b|@>@;
}

@ This just adds a layer to itself.  Couldn't |Basic_Add_Layers| be
used?  It would mean that there would be no restriction on the 
use of variables --- i.e., |R| could be used as both a factor and
as a result.

@<Prototype for |Double_Once|@>=
void Double_Once(int n, double **R, double **T)

@
@<Definition for |Double_Once|@>=
	@<Prototype for |Double_Once|@>
{
  @<Allocate memory for |a| and |b|@>@;
  Basic_Add_Layers(n, R, T, R, R, T, T, R, T, a, b);
  @<Free Memory for |a| and |b|@>@;
}

@ |Double_Until| and |Double_Until_Infinite| are the only ones that really take advantage of
the external allocation of memory from the routine.  I was kind of careful to 
make sure that this routine terminates if bad |start| and |end| values are given
i.e., |end!=start|$\cdot 2^k$.  Futhermore, it should work correctly
if the target thickness is infinite.  I suppose that I could put some
error warnings in...but right now I don't want to take the time.

@<Prototype for |Double_Until|@>=
void Double_Until(int n, double **r, double **t, double start, double end)

@
@<Definition for |Double_Until|@>=
	@<Prototype for |Double_Until|@>
{
	if (end == HUGE_VAL) {
		Double_Until_Infinite(n, r, t);
		return;
		}

   {
   @<Allocate memory for |a| and |b|@>@;
	while (fabs(end-start) > 0.00001 && end > start) {
		Basic_Add_Layers(n, r, t, r, r, t, t, r, t, a, b);
		start *= 2;
	}
   @<Free Memory for |a| and |b|@>@;
   }
}

@ |Double_Until_Infinite| continues doubling until the thickness
of the slab is essentially infinite.  Originally I had defined
infinite as a diffuse transmission less than $10^{-6}$.  However,
 when the albedo is unity, then this is kind of impractical
and I changed the definition of infinity to be that the diffuse
transmission changes by less than one part in $10^{-6}$
after one doubling step.  The more I think about this, the
less sense it makes....

@<Prototype for |Double_Until_Infinite|@>=
void Double_Until_Infinite(int n, double **r, double **t)

@
@<Definition for |Double_Until_Infinite|@>=
	@<Prototype for |Double_Until_Infinite|@>
{
  double oldutu, UTU, UT1;
  
  @<Allocate memory for |a| and |b|@>@;
  UTU = 0.0;
  do {
    oldutu = UTU;
    Basic_Add_Layers(n, r, t, r, r, t, t, r, t, a, b);
    URU_and_UR1(n, 1.0, t, &UTU, &UT1);
  } while (fabs(UTU - oldutu) >= 0.000001);
  
  @<Free Memory for |a| and |b|@>@;

}

@*2 Internal Radiance.

|Between| finds the radiance between two slabs.
This equation for the upward radiance at the interface between two layers 
is
$$
\bfL_- = (\bfe-\bfr^{12}\star\bfr^{10})^{-1} (\bfr^{12}\star\bft^{01}\star\bfL^0_+ 
+ \bft^{21}\star\bfL^2_-)
$$
where $\bfL^0_+$ is the downward radiance on the top layer and 
$\bfL^2_-$ is the upward radiance on the bottom layer.  The
equation for the downward mid-layer radiance can be obtained similarly using 
$$
\bfL_+ = (\bfe-\bfr^{10}\star\bfr^{12})^{-1} (\bft^{01}\star\bfL^0_+
+ \bfr^{10}\star\bft^{21}\star\bfL^2_-)
$$
Now assume that $\bfL^2_-$ is zero.  Then the matrix
$$
\bfL_- = (\bfe-\bfr^{12}\star\bfr^{10})^{-1} \bfr^{12}\star\bft^{01}
$$
can be used to find the downward fluence by simply star multiplying with
the downward irradiance.  Similarly,
$$
\bfL_+ = (\bfe-\bfr^{10}\star\bfr^{12})^{-1} \bft^{01}
$$

@<Prototype for |Between|@>=
void Between(int n, 
	double **R01, double **R10, double **T01, double **T10, 
	double **R12, double **R21, double **T12, double **T21, 
	double **Lup, double **Ldown)

@ @<Definition for |Between|@>=
	@<Prototype for |Between|@>
{
  @<Allocate memory for |a| and |b|@>@;
  
  Star_Multiply(n, R10, R12, a);
  Star_One_Minus(n, a);
  Right_Inverse_Multiply(n, a, T01, Ldown);
  
  Star_Multiply(n, R12, R10, a);
  Star_One_Minus(n, a);
  Right_Inverse_Multiply(n, a, R12, b);
  Star_Multiply(n, b, T01, Lup);
  
  @<Free Memory for |a| and |b|@>@;
}

@*2 Utility routines.

@ Star matrix multiplication $A \star B$ is defined to directly correspond to an 
integration, i.e.
$$
A\star B = \int_0^1A(\mu_,\mu')B(\mu',\mu'')\, 2\mu d\mu
$$
then
$$
A\star B = \sum_j A^{ij}2\mu_jw_j B^{jk}
$$
where $\mu_j$ is the $j$th quadrature angle and $w_j$ is its corresponding weight.  
It is sometimes useful to consider these matrix ``star 
multiplications'' as normal matrix multiplications which include a diagonal matrix $c$
$$
\bfc_{ij}=2\mu_i w_i \delta_{ij}
$$
Thus a matrix star multiplication may be written
$$
A\star B = A\, \bfc\, B
$$
where the multiplications on the RHS of the above equation are usual matrix 
multiplications.  

Since the routine |Matrix_Multiply| that multiplies the matrices |A| and 
|B| to get |C|, allows |A| and |C| to be coincident.  I first find $C=A c$
and then do $C=C\cdot B$.  This allows us to avoid allocating a temporary matrix. 
|A| may occupy the same memory as |C|, but |B| and |C| must be distinct.

@<Definition for |Star_Multiply|@>=
	static void Star_Multiply(int n, double **A, double **B, double **C)
{
  Right_Diagonal_Multiply(n, A, twoaw, C);
  Matrix_Multiply(n, C, B, C);
}

@ This subtracts the matrix |A| from the unit matrix for star 
multiplication.  

@<Definition for |Star_One_Minus|@>=
	static void Star_One_Minus(int n, double **A)
{
  int i, j;

	for (i=1; i<=n; i++) {
		for (j=1; j<=n; j++) 
			A[i][j] *= -1;
		A[i][i] += 1.0/twoaw[i];
	}
	
}

@ @<Allocate memory for |a| and |b|@>=
	double **a, **b;
	
	a = dmatrix(1,n,1,n);
  	b = dmatrix(1,n,1,n);

@
@<Free Memory for |a| and |b|@>=
  free_dmatrix(a,1,n,1,n);
  free_dmatrix(b,1,n,1,n);
