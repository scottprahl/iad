@*1 AD Matrix.

This is a part of the core suite of files for the adding-doubling
program.  Not surprisingly, this program includes routines to
manipulate matrices.  These routines require that the matrices
be stored using the allocation scheme outline in {\it Numerical Recipes\/}
by Press {\it et al.}  I have spent some time optimizing the matrix
multiplication routine |Matrix_Multiply| because roughly half the time
in any adding-doubling calculation is spent doing matrix multiplication.
Lastly, I should mention that all the routines assume a square matrix
of size |n| by |n|.

\def\bfA{{\bf A}}
\def\bfB{{\bf B}}
\def\bfC{{\bf C}}
\def\bfD{{\bf D}}

@(ad_matrx.c@>=

	#include <stddef.h>
	#include <math.h>
	#include "ad_globl.h"
	#include "ad_matrx.h"
	#include "nr_util.h"

	@<Definition for |Copy_Matrix|@>@;
	@<Definition for |One_Minus|@>@;
	@<Definition for |Transpose_Matrix|@>@;
	@<Definition for |Diagonal_To_Matrix|@>@;
	@<Definition for |Right_Diagonal_Multiply|@>@;
	@<Definition for |Left_Diagonal_Multiply|@>@;
	@<Definition for |Matrix_Multiply|@>@;
	@<Definition for |Matrix_Sum|@>@;
	@<Definition for |Solve|@>@;
	@<Definition for |Decomp|@>@;
	@<Definition for |Matrix_Inverse|@>@;
	@<Definition for |Left_Inverse_Multiply|@>@;
	@<Definition for |Right_Inverse_Multiply|@>@;

@ In this module I collect up information that needs to be written to
the header file \.{ad\_matrx.h} so that other source files that want to
make use of the function defined here will have the necessary declarations
available.

@(ad_matrx.h@>=
	@<Prototype for |Copy_Matrix|@>;
	@<Prototype for |One_Minus|@>;
	@<Prototype for |Transpose_Matrix|@>;
	@<Prototype for |Diagonal_To_Matrix|@>;
	@<Prototype for |Right_Diagonal_Multiply|@>;
	@<Prototype for |Left_Diagonal_Multiply|@>;
	@<Prototype for |Matrix_Multiply|@>;
	@<Prototype for |Matrix_Sum|@>;
	@<Prototype for |Solve|@>;
	@<Prototype for |Decomp|@>;
	@<Prototype for |Matrix_Inverse|@>;
	@<Prototype for |Left_Inverse_Multiply|@>;
	@<Prototype for |Right_Inverse_Multiply|@>;

@*2 Simple Matrix Routines.

@ |Copy_Matrix| replaces the matrix |B| by |A| 

@<Prototype for |Copy_Matrix|@>=
void Copy_Matrix(int n, double **A, double **B)

@ @<Definition for |Copy_Matrix|@>=
		@<Prototype for |Copy_Matrix|@>
{
  double *a_ptr, *b_ptr, *a_last;
  
  a_last = &A[n][n];
  a_ptr = &A[1][1];
  b_ptr = &B[1][1];
  
  while (a_ptr<=a_last)
  	*b_ptr++ = *a_ptr++;
}

@ |One_Minus| replaces the matrix A by 1-A 
@<Prototype for |One_Minus|@>=
void One_Minus(int n, double **A)

@ @<Definition for |One_Minus|@>=
		@<Prototype for |One_Minus|@>
{
  int i, j;

	for (i=1; i<=n; i++) {
		for (j=1; j<=n; j++) 
			A[i][j] *= -1;
		A[i][i] += 1.0;
	}
	
}

@ |Transpose_Matrix| transposes a matrix.
@<Prototype for |Transpose_Matrix|@>=
void Transpose_Matrix(int n, double **a)

@ @<Definition for |Transpose_Matrix|@>=
	@<Prototype for |Transpose_Matrix|@>
{
	int             i, j;
	double          swap;

	for (i = 1; i <= n; i++) {
		for (j = i + 1; j <= n; j++) {
			swap = a[i][j];
			a[i][j] = a[j][i];
			a[j][i] = swap;
		}
	}
}

@ |Diagonal_To_Matrix| converts a diagonal array to a matrix
@<Prototype for |Diagonal_To_Matrix|@>=
void Diagonal_To_Matrix(int n, double *Diag, double **Mat)

@ @<Definition for |Diagonal_To_Matrix|@>=
	@<Prototype for |Diagonal_To_Matrix|@>
{
	int             i, j;

	for (i = 1; i <= n; i++){
		for (j = 1; j <= n; j++)
			Mat[i][j] = 0.0;
		Mat[i][i] = Diag[i];
	}
}

@ |Right_Diagonal_Multiply| multiplies the matrix A by the 
diagonal matrix |B|, puts the
 result in |C|.
 |A| and |C| can be the same matrix
 $$
 C \leftarrow A \cdot B
 $$
 Note that |B| is stored as a vector.
 
@<Prototype for |Right_Diagonal_Multiply|@>=
void Right_Diagonal_Multiply(int n, double **A, double *B, double **C)

@ @<Definition for |Right_Diagonal_Multiply|@>=
	@<Prototype for |Right_Diagonal_Multiply|@>
{
	int             i, j;

	for (i = 1; i <= n; i++) 
		for (j = 1; j <= n; j++) 
			C[i][j] = A[i][j]*B[j];
}

@ |Left_Diagonal_Multiply| multiplies the diagonal matrix a by the matrix B, puts the
  result in C. B and C can be the same matrix
  
@<Prototype for |Left_Diagonal_Multiply|@>=
void Left_Diagonal_Multiply(int n, double *A, double **B, double **C)

@ @<Definition for |Left_Diagonal_Multiply|@>=
	@<Prototype for |Left_Diagonal_Multiply|@>
{
	int             i, j;

	for (i = 1; i <= n; i++) 
		for (j = 1; j <= n; j++)
			C[i][j] = A[i] * B[i][j];
}

@ |Matrix_Sum| adds the two matrices A and B, puts the result in C
  The matrices need not be distinct
  
@<Prototype for |Matrix_Sum|@>=
void Matrix_Sum(int n, double **A, double **B, double **C)

@ @<Definition for |Matrix_Sum|@>=
	@<Prototype for |Matrix_Sum|@>
{
	int             i, j;

	for (i = 1; i <= n; i++)
		for (j = 1; j <= n; j++)
			C[i][j] = A[i][j] + B[i][j];
}

@*2 Matrix Multiplication.
	This is the crux of this whole unit at present.  Most of the time in the 
	adding-doubling algorithm is spent doing matrix multiplication and this
	implementation has been optimized using pointers. 

	|Matrix_Multiply| multiplies the two matrices |A| and |B| and puts the result in |C|.
	The following routine requires that |C| does not occupy the same space as
	|B|, but it can be coincident with |A|.   There is no inherent reason that |A|, |B|,
	and |C| must all be $n\times n$ matrices.  However, all the matrices in
	the adding-doubling method are square and I did not want to pass three 
	separate dimensions to this routine.
	
	The usual way matrix multiplication uses an algorithm something
	similar to:

	@<unused fragment one@>=
	for (i = 1; i <= n; i++) {
		for (j = 1; j <= n; j++) {
			C[i][j] = 0.0;
			for (k = 1; k <= n; k++)
				C[i][j] += A[i][k] * B[k][j];
		}
	}

	@ This has the unfortunate problem that the innermost loop indexes
	successive columns of |A| and successive rows of |B|.  Because 
	indexing successive rows requires something other than a unit
	increment of the matrix pointer, a different algorithm is used.  
	In this case,

		@<unused fragment two@>=
	for (i = 1; i <= n; i++) 
		for (j = 1; j <= n; j++) 
			C[i][j] = 0.0;

	for (i = 1; i <= n; i++) {
		for (k = 1; k <= n; k++) {
			for (j = 1; j <= n; j++)
				C[i][j] += A[i][k] * B[k][j];
		}
	}

	@ This particular form of indexing was chosen to take advantage of the 
	row storage of matrices designated by the Numerical Recipes scheme.
	The innermost loop of the matrix multiplication routine
	now only requires unit increments of the matrix pointers |C| and |B|.
	
	Explictly using pointers to the entries in the salient matrices makes
	this routine roughly 20\% faster than when the above implementation is
	used.  Profiling of the code indicates that roughly 45\% of the time spent in an
 	adding-doubling calculation is spent in this one routine.  Therefore even a 
	modest 20\% increase will translate to a ten percent improvement in performance.

	Finally, the algorithm can be improved to allow the pointers to
	|A| and |C| to be the same.    This is sufficient to allow us to avoid 
	allocating an extra matrix here and there.  It can easily be adapted
	to work with ``star'' multiplication by premultiplying using 
	|Right_Diagonal_Multiply|.  The drawbacks are that a vector |D| must
	be allocated on each call.  It is also necessary to copy the data from
	the vector |D| to the output matrix |C|.  

@ @<Prototype for |Matrix_Multiply|@>=
void Matrix_Multiply(int n, double **A, double **B, double **C)

@ @<Definition for |Matrix_Multiply|@>=
	@<Prototype for |Matrix_Multiply|@>
{
	@<Local variables for |Matrix_Multiply|@>@;
	@<Do awkward cases@>@;
	@<Allocate memory for |D|@>@;
	@<Initialization for |Matrix_Multiply|@>@;
	@<Multiplying |A| and |B|@>@;
	@<Free memory for |D|@>@;
}

	@ @<Local variables for |Matrix_Multiply|@>=
	double          *a_ptr, *a_start;
	double 			*b_start, *b_last;
	double 			*c_start, *c_very_last, *c_ptr;
	double			*D;
	double			*d_start, *d_last;
	register double	t, *d_ptr, *b_ptr;
	ptrdiff_t		row;

	@ @<Do awkward cases@>=
	if (n<=0) {
	 	AD_error("Non-positive dimension passed to Matrix_Multiply");
	}
	else if (n==1) {
		C[1][1]=A[1][1]*B[1][1];
		return;
		}

	@  I need a temporary vector equal to the row length of |C| to hold
	intermediate calculations.  This will allow |A| and |C| to point to
	the same matrix and still yield the correct results.
	
	@<Allocate memory for |D|@>=
	D = dvector(1,n);
	
	@ During the initialization phase, I need to know how far it is from
	one row to the next row.  Because of the peculiar way that Numerical
	Recipes allocates the matrices, this may and probably is not equal to
	|n|.  The number of entries is found explicitly by subtracting a pointer
	to the first entry in row one from the first entry in row two.  This assumes
	that the size of the matrix is at least two.  To make this routine bulletproof,
	this would need to be changed---but I do not think it is really necessary.
	
	@<Initialization for |Matrix_Multiply|@>=
	a_start = &A[1][1];
	b_last  = &B[n][1];
	row    = &A[2][1]-a_start;	
	c_very_last = &C[n][n];
	d_start = &D[1];
	d_last =&D[n];
	
	@ There may be a better way of doing this, but I bet it would depend on specific
	knowlege about how zero is stored in the computer.
	
	@<Zero |D|@>=
	d_ptr=d_start;
	while (d_ptr<=d_last)
		*d_ptr++ = 0.0;

	@ Copy the contents of |D| to |C|.  This could potentially be sped up using
	|memmove()| but I just want it to work for now.
	
	@<Copy |D| into |C|@>=
	d_ptr=d_start;
	c_ptr=c_start;
	while (d_ptr<=d_last)
		*c_ptr++ = *d_ptr++;

	@ Here is the heart of the routine.  The first row of |C| is filled
	completely, then the routine goes on to the second row and so on.
	The inner loop is responsible for multiplying |A[i][k]| (represented
	by |t=*a_ptr|) by every element in row |i| and adding it to the
	appropriate element in row |i| of |C|.  
	
	@<Multiplying |A| and |B|@>=
		
	for (c_start=&C[1][1]; c_start<=c_very_last; c_start += row){
		a_ptr = a_start;
	          @<Zero |D|@>@;
	   
		for (b_start = &B[1][1]; b_start<=b_last; b_start += row) {
			t = *a_ptr++;
			b_ptr=b_start;
			d_ptr=d_start;
			while (d_ptr<=d_last)
				*d_ptr++ += t * (*b_ptr++);
		}
	          @<Copy |D| into |C|@>@;
		a_start+= row;
	}

	@ Dump the memory that was allocated.
	@<Free memory for |D|@>=
	free_dvector(D,1,n);

@*2 Matrix Decomposition.

@ @<Prototype for |Decomp|@>=
void Decomp(int n, double **A, double *condition, int *ipvt)

@ |Decomp| decomposes a double matrix by Gaussian elimination and estimates the condition
 of the matrix.
  
Use solve to compute solutions to linear systems
 
On input |n| is the order of the matrix and |A| is the matrix to be 
triangularized. 
 
On output |A| contains an upper triangular matrix $U$ and a permuted 
version of a lower triangular matrix ${\bf I}-{\bf L}$ so that (permutation 
matrix)*A=L*U.  |condition| is an estimate of the condition of |A|.  For the 
linear system $\bf A X=B$, changes in |A| and |B| may cause changes 
condition times as large in $X$.  If |condition+1.0 = condition|, |A| is 
singular to working precision.  |condition| is set to |1.0e32| if exact 
singularity is detected.  |ipvt| is the pivot vector |ipvt(k)| is the index 
of the kth pivot row |ipvt(n) = (-1)|$^{\rm (number of interchanges)}$

@<Definition for |Decomp|@>=
	@<Prototype for |Decomp|@>
{
	double t, anorm;
	int i, j, k, m;

	@<Do |n==1| case @>@; 
	@<Compute 1-norm of |A|@>@; 
	@<Gaussian elimination with partial pivoting @>@;
	@<Check for singularity@>@;
}

@ This should probably be fixed to compute the inverse of a non-zero
1by 1 matrix.

@<Do |n==1| case @>=
	ipvt[n] = 1;

 	if (n==1) {
 		if (A[1][1]==0){
 			AD_error("1 X 1 Matrix is Singular --- i.e. zero");
 			return;}
 	}
 	
@ @<Compute 1-norm of |A|@>=
	anorm = 0.0;
	for (j = 1; j <= n; j++) {
		t = 0.0;
		for (i = 1; i <= n; i++)
			t += fabs(A[i][j]);
		if (t > anorm)
			anorm = t;
	}
	
@ @<Gaussian elimination with partial pivoting @>=
	for (k = 1; k < n; k++) {		
		@<Find pivot@>@;
		@<Compute multipliers@>@;
		 @<Interchange and eliminate by columns@>@;
	} 

@ @<Find pivot@>=
		m = k;
		for (i = k+1; i <= n; i++) 
			if (fabs(A[i][k]) > fabs(A[m][k]))
				m = i;
		
		ipvt[k] = m;
		if (m != k)
			ipvt[n] *= -1;
		t = A[m][k];
		A[m][k] = A[k][k];
		A[k][k] = t;
		
		 /* skip step if pivot is zero */
		if (t == 0)  continue;
		
@ @<Compute multipliers@>=
		for (i = k+1; i <= n; i++)
			A[i][k] /= -t;
			
@ @<Interchange and eliminate by columns@>=
		for (j = k+1; j <= n; j++) {
			t = A[m][j];
			A[m][j] = A[k][j];
			A[k][j] = t;
			if (t == 0) continue;
			for (i = k+1; i <= n; i++)
				A[i][j] += A[i][k] * t;
		}
		
@ @<Check for singularity@>=

	*condition = 1.0;
	for (k=1;k<=n;k++) {
		if (A[k][k]==0.0){
			*condition=1e32;
			return;
		}
	}

@*2 Solving systems of equations.
@
@<Prototype for |Solve|@>=
void Solve(int n, double **A, double *B, int *ipvt)

@
This procedure finds the solution of the linear system $A X=B$
Don't use if |Decomp| has found a singularity

On input |n| is the order of matrix, |A| is the triangularized matrix 
obtained form |Decomp|.  |B| is the right hand side vector and |ipvt| is
the pivot vector obtained from |Decomp|
  
On output |B| is the solution vector $X$.

@<Definition for |Solve|@>=
		@<Prototype for |Solve|@>
{
	int i, k, m;
	double t;

	@<Forward elimination@>@;
	@<Back substitution@>@;
}

@ @<Forward elimination@>=
	for (k = 1; k<n; k++) {  
		m = ipvt[k];
		t = B[m];
		B[m] = B[k];
		B[k] = t;
		for (i = k+1; i<=n; i++)
			B[i] += A[i][k] * t;
	}

@ @<Back substitution@>=
	for (k = n; k>1; k--) { 
		B[k] /= A[k][k];
		t = -B[k];
		for (i=1; i<k; i++)
			B[i] += A[i][k] * t;
	}
	
	B[1] /= A[1][1];


@ Finds the inverse of the matrix |A| (of order |n|) and stores
 the answer in |Ainv|.
 
@<Prototype for |Matrix_Inverse|@>=
void Matrix_Inverse(int n, double **A, double **Ainv)

@ @<Definition for |Matrix_Inverse|@>=
	@<Prototype for |Matrix_Inverse|@>
{
	int            *ipvt;
	int             i, j;
	double         *work;
	double          condition;

	ipvt = ivector(1, n);
	Decomp(n, A, &condition, ipvt);
	if (condition == (condition + 1) || condition == 1e32) {
		free_ivector(ipvt, 1, n);
		AD_error("Singular Matrix ... failed in Inverse_Multiply\n");
	}
	work = dvector(1, n);
	for (i = 1; i <= n; i++) {
		for (j = 1; j <= n; j++)
			work[j] = 0.0;
		work[i] = 1.0;
		Solve(n, A, work, ipvt);
		for (j = 1; j <= n; j++)
			Ainv[j][i] = work[j];
	}

	free_dvector(work, 1, n);
	free_ivector(ipvt, 1, n);
}


@ @<Prototype for |Left_Inverse_Multiply|@>=
void Left_Inverse_Multiply(int n, double **D, double **C, double **A)

@ |Left_Inverse_Multiply| computes $\bfA=\bfC\cdot\bfD^{-1}$ where |A|, |C| 
and |D| are all |n| by |n| matrices.  This is faster than inverting and 
then multiplying by a factor of six.  Space for |A| should be allocated 
before calling this routine.

@<Definition for |Left_Inverse_Multiply|@>=
	@<Prototype for |Left_Inverse_Multiply|@>
{
	int            *ipvt;
	int             i, j;
	double         *work;
	double          condition;

	Transpose_Matrix(n, D);
	ipvt = ivector(1, n);
	Decomp(n, D, &condition, ipvt);

	/* Check for singular result */
	if (condition == (condition + 1) || condition == 1e32) {
		free_ivector(ipvt, 1, n);
		AD_error("Singular Matrix ... failed in Left_Inverse_Multiply\n");
	}
	
	work = dvector(1, n);
	for (i = 1; i <= n; i++) {		/* Cycle through all the row in |C| */
	
		for (j = 1; j <= n; j++)	/* put a row of |C| into work */
			work[j] = C[i][j];		/* and avoid a Transpose Matrix */
		Solve(n, D, work, ipvt);
		for (j = 1; j <= n; j++)	/* Again avoiding a Transpose Matrix */
			A[i][j] = work[j];		/* stuff the results into a row of |A| */
	}

	free_dvector(work, 1, n);
	free_ivector(ipvt, 1, n);
}

@ @<Prototype for |Right_Inverse_Multiply|@>=
void Right_Inverse_Multiply(int n, double **D, double **C, double **A)

@ |Right_Inverse_Multiply| computes $\bfA=\bfD^{-1}\cdot\bfC$ where |A|, |C| 
and |D| are all |n| by |n| matrices.  This is faster than inverting and 
then multiplying by a factor of six.  Space for |A| should be allocated 
before calling this routine.

@<Definition for |Right_Inverse_Multiply|@>=
	@<Prototype for |Right_Inverse_Multiply|@>
{
	int            *ipvt;
	int             i, j;
	double         *work;
	double          condition;

	ipvt = ivector(1, n);
	Decomp(n, D, &condition, ipvt);

	/* Check for singular result */
	if (condition == (condition + 1) || condition == 1e32) {
		free_ivector(ipvt, 1, n);
		AD_error("Singular Matrix ... failed in Right_Inverse_Multiply\n");
	}
	
	work = dvector(1, n);
	for (i = 1; i <= n; i++) {		/* Cycle through all the rows */
	
		for (j = 1; j <= n; j++)	/* put a column of |C| into work */
			work[j] = C[j][i];	
		Solve(n, D, work, ipvt);
		for (j = 1; j <= n; j++)	/* stuff the results into a column of |A| */
			A[j][i] = work[j];		
	}

	free_dvector(work, 1, n);
	free_ivector(ipvt, 1, n);
}
