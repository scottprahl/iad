@*1 AD Layers.
This file provides routines to obtain reflection and transmission values for
normal illumination of several multiple scattering and absorbing layers.

@(ad_layers.c@>=
#include <math.h>
#include <float.h>
#include "nr_util.h"
#include "ad_globl.h"
#include "ad_bound.h"
#include "ad_doubl.h"
#include "ad_prime.h"
#include "ad_matrx.h"
#include "ad_prime.h"

	@<Definition for |RT_Layers_All|@>@;
	@<Definition for |RT_Layers|@>@;

@ @(ad_layers.h@>=
	@h
	@<Prototype for |RT_Layers|@>;
	@<Prototype for |RT_Layers_All|@>;

@*2 RT Layers.
Sometimes you just need to know the total reflection and transmission from a
target consisting of multiple layers.  This is the routine for you.  It adds
a bunch of scattering and absorbing layers together which have the same index
of refraction together.  The top and bottom are possibly bounded by glass slides.
This is not particularly fast, but it should get the job done.

|nlayers| specifies the number of different layers (not including possible glass
slides above and below the composite sample.  The optical properties are passed in
three zero-based arrays of doubles.  For example |a[1]| is the albedo of the 
second layer.  

@<Prototype for |RT_Layers_All|@>=
void RT_Layers_All(int n, 	
	double nslab, 
	double ntopslide, 
	double nbottomslide, 
	int nlayers,
	double a[],
	double b[],
	double g[],
	double *dUR1, double *dUT1, double *dURU, double *dUTU,
	double *uUR1, double *uUT1, double *uURU, double *uUTU)

@ @<Definition for |RT_Layers_All|@>=
	@<Prototype for |RT_Layers_All|@>
{
	@<Declare variables for |RT_Layers|@>@;
	@<Validate layer properties@>@;
	@<Allocate slab memory@>@;
	@<Initialize slab structure@>@;
	@<Allocate and generate top and bottom boundaries@>@;
	@<Initialize composite layer@>@;
	@<Add all composite layers together@>@;
	@<Add top and bottom boundaries@>@;
	@<Free memory for |RT_Layers|@>@;
}

@ Simple sanity checks to ensure values are reasonable.

@<Validate layer properties@>= 

	if (nlayers < 1) return;
	if (nslab < 0) return;
	if (ntopslide < 0) return;
	if (nbottomslide < 0) return;
	
	for (i=0; i<nlayers; i++) {
		if (a[i]<0 || a[i]>1) return;
		if (b[i]<0 ) return;
		if (g[i]<-1 || g[i]>1) return;
	}
		
@ @<Declare variables for |RT_Layers|@>=
	struct AD_slab_type slab;
	struct AD_method_type method;
	double *R01, *R10, *T01, *T10;
	double *R34, *R43, *T34, *T43;
	double **R12, **R21, **T12, **T21;
	double **R23, **R32, **T23, **T32;
	double **R13, **R31, **T13, **T31;
	double **atemp, **btemp;
	int i;
    *dUR1=-1; 
    *dUT1=-1;
	*dURU=-1;
	*dUTU=-1;
	*uUR1=-1;
	*uUT1=-1;
	*uURU=-1;
	*uUTU=-1;

@ @<Allocate slab memory@>=

	R12 = dmatrix(1, n, 1, n);
	R21 = dmatrix(1, n, 1, n);
	T12 = dmatrix(1, n, 1, n);
	T21 = dmatrix(1, n, 1, n);
	
	R23 = dmatrix(1, n, 1, n);
	R32 = dmatrix(1, n, 1, n);
	T23 = dmatrix(1, n, 1, n);
	T32 = dmatrix(1, n, 1, n);
	
	R13 = dmatrix(1, n, 1, n);
	R31 = dmatrix(1, n, 1, n);
	T13 = dmatrix(1, n, 1, n);
	T31 = dmatrix(1, n, 1, n);
	
	atemp = dmatrix(1, n, 1, n);
	btemp = dmatrix(1, n, 1, n);

@ Create the matrices needed for the top and bottom
@<Allocate and generate top and bottom boundaries@>=
	R01 = dvector(1, n);
	R10 = dvector(1, n);
	T01 = dvector(1, n);
	T10 = dvector(1, n);
	Init_Boundary(slab, n, R01, R10, T01, T10, TOP_BOUNDARY);

	R34 = dvector(1, n);
	R43 = dvector(1, n);
	T34 = dvector(1, n);
	T43 = dvector(1, n);
	Init_Boundary(slab, n, R34, R43, T34, T43, BOTTOM_BOUNDARY);

@ We set this to be a clear layer so that the composite layer will be 
created properly.  The index of refraction of the slab is important so
that the quadrature angles will be chosen correctly.

@<Initialize slab structure@>=
	slab.n_slab = nslab;
	slab.n_top_slide = ntopslide;
	slab.n_bottom_slide = nbottomslide;
	slab.b_top_slide = 0;
	slab.b_bottom_slide = 0;
	slab.a = 0.0;
	slab.b = 0.0;
	slab.g = 0.0;
	slab.phase_function = HENYEY_GREENSTEIN;
	slab.cos_angle = 1.0;

@ The composite layer initially has 0\% reflection and 100\% transmission.
We fob the details on how this layer is created to the |RT_Matrices| which
goes to the trouble to initialize |method| and call |Zero_Layer| for us.
Finally, since this optical problem is not reversible (illumination from
below gives a different answer), we need to initialize the upward matrices
as well.  This simplifies the code when adding successive layers.

@<Initialize composite layer@>=
	RT_Matrices(n, &slab, &method, R23, T23);
	Copy_Matrix(n, R23, R32);
	Copy_Matrix(n, T23, T32);

@ Now add the layers together.  Since the composite layer has been initialized to be
a clear layer, we can just add layers to it.  We start from the bottom.  Find
the transport matrices for this layer.  Add this layer to the top of the composite
layer.  This is repeated for each of the layers.

@<Add all composite layers together@>=
	
	while (nlayers >= 1) {
		nlayers--;
		slab.a = a[nlayers];
		slab.b = b[nlayers];
		slab.g = g[nlayers];
		RT_Matrices(n, &slab, &method, R12, T12);
		Add(n, R12, R12, T12, T12, R23, R32, T23, T32, R13, R31, T13, T31);
		Copy_Matrix(n, R13, R23);
		Copy_Matrix(n, R31, R32);
		Copy_Matrix(n, T13, T23);
		Copy_Matrix(n, T31, T32);
	}		
	
@ The only confusing part about this piece of code is that the layer numbering
gets all messed up.  The composite layer is in the 23 matrices.  This gets added
to the top 01 boundary and should be labeled the 03 matrix.  Instead I use the 
already allocated 13 matrices.  This layer is then added to the bottom 34 
matrices and should result in 04 matrices, but once again I use the 23 matrices.
Finally, the total reflectances and transmittances are calculated, so that all
the remains is to free the allocated memory!  Not so hard after all.

@<Add top and bottom boundaries@>=

	Add_Top(n, R01, R10, T01, T10, R23, R32, T23, T32, R13, R31, T13, T31, atemp, btemp);
	Add_Bottom(n, R13, R31, T13, T31, R34, R43, T34, T43, R23, R32, T23, T32, atemp, btemp);
	URU_and_UR1(n, slab.n_slab, R23, dURU, dUR1);
	URU_and_UR1(n, slab.n_slab, R32, uURU, uUR1);
	Transpose_Matrix(n,T23);
	Transpose_Matrix(n,T32);
	URU_and_UR1(n, slab.n_slab, T23, dUTU, dUT1);
	URU_and_UR1(n, slab.n_slab, T32, uUTU, uUT1);

@ @<Free memory for |RT_Layers|@>=
	free_dvector(R01, 1, n);
	free_dvector(R10, 1, n);
	free_dvector(T01, 1, n);
	free_dvector(T10, 1, n);

	free_dmatrix(R12, 1, n, 1, n);
	free_dmatrix(R21, 1, n, 1, n);
	free_dmatrix(T12, 1, n, 1, n);
	free_dmatrix(T21, 1, n, 1, n);

	free_dmatrix(R23, 1, n, 1, n);
	free_dmatrix(R32, 1, n, 1, n);
	free_dmatrix(T23, 1, n, 1, n);
	free_dmatrix(T32, 1, n, 1, n);

	free_dmatrix(R13, 1, n, 1, n);
	free_dmatrix(R31, 1, n, 1, n);
	free_dmatrix(T13, 1, n, 1, n);
	free_dmatrix(T31, 1, n, 1, n);

	free_dmatrix(atemp, 1, n, 1, n);
	free_dmatrix(btemp, 1, n, 1, n);

	free_dvector(R34, 1, n);
	free_dvector(R43, 1, n);
	free_dvector(T34, 1, n);
	free_dvector(T43, 1, n);

@ This just returns the reflection and transmission for light travelling
downwards.  This is most often what is desired.

@<Prototype for |RT_Layers|@>=
void RT_Layers(int n, 	
	double nslab, 
	double ntopslide, 
	double nbottomslide, 
	int nlayers,
	double a[],
	double b[],
	double g[],
	double *UR1, double *UT1, double *URU, double *UTU)

@ @<Definition for |RT_Layers|@>=
	@<Prototype for |RT_Layers|@>
{
double uUR1, uUT1, uURU, uUTU;

RT_Layers_All(n,nslab,ntopslide,nbottomslide,nlayers,a,b,g,\
				UR1,UT1,URU,UTU, &uUR1, &uUT1, &uURU, &uUTU);
}
