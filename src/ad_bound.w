@*1 AD Boundary.

This section has routines associated with incorporating 
boundary conditions into the adding-doubling algorithm.

@(ad_bound.c@>=

#include <math.h>
#include <stdio.h>
#include "ad_globl.h"
#include "ad_bound.h"
#include "ad_frsnl.h"
#include "ad_matrx.h"
#include "nr_util.h"

    @<Prototype for |A_Add_Slide|@>;
    @<Prototype for |B_Add_Slide|@>;

    @<Definition for |Init_Boundary|@>@;
    @<Definition for |Boundary_RT|@>@;
    @<Definition for |Add_Top|@>@;
    @<Definition for |Add_Bottom|@>@;
    @<Definition for |A_Add_Slide|@>@;
    @<Definition for |B_Add_Slide|@>@;
    @<Definition for |Add_Slides|@>@;
    @<Definition for |Sp_RT|@>@;

@ @(ad_bound.h@>=
    @h
    @<Prototype for |Init_Boundary|@>;
    @<Prototype for |Boundary_RT|@>;
    @<Prototype for |Add_Top|@>;
    @<Prototype for |Add_Bottom|@>;
    @<Prototype for |Add_Slides|@>;
    @<Prototype for |Sp_RT|@>;

@*2 Boundary Initialization.

@ |Init_Boundary| creates reflection and transmission matrices to simulate
a boundary.
 If |boundary==TOP_BOUNDARY| then the arrays returned are for the top surface and the  
 labels are as expected i.e. |T01| is the reflection for light from air passing to the slab. 
 Otherwise the calculations are made for the bottom surface and 
 the labels are backwards i.e. |T01 == T32| and |T10 == T23|, where 0 is the first air slide 
 surface, 1 is the slide/slab surface, 2 is the second slide/slab surface, and 3 is the 
 bottom slide/air surface 

@d TOP_BOUNDARY 0 
@d BOTTOM_BOUNDARY 1

@<Prototype for |Init_Boundary|@>=
void Init_Boundary(struct AD_slab_type slab, int n, @/
double *R01, double *R10, double *T01, double *T10, @/
                   char boundary)

@ @<Definition for |Init_Boundary|@>=
    @<Prototype for |Init_Boundary|@>
{
    if (boundary ==TOP_BOUNDARY) {
        Boundary_RT(1.0, slab.n_top_slide, slab.n_slab, n, slab.b_top_slide, R01, T01);
        Boundary_RT(slab.n_slab, slab.n_top_slide, 1.0, n, slab.b_top_slide, R10, T10);
    }   else {
        Boundary_RT(1.0, slab.n_bottom_slide, slab.n_slab, n, slab.b_bottom_slide, R10, T10);
        Boundary_RT(slab.n_slab, slab.n_bottom_slide, 1.0, n, slab.b_bottom_slide, R01, T01);
    }
}


@ |Boundary_RT| computes the diagonal matrix (represented as an array) 
   that characterizes reflection and transmission at an air (0), absorbing glass (1), 
   slab (2) boundary.  The reflection matrix is the same entering or exiting the slab.  
   The transmission matrices should differ by a factor of ($n_{\rm slab}/n_{\rm outside})^4$, 
   due to $n^2$ law of radiance, but there is some inconsistency in the program and
   if I use this principle then regular calculations for $R$ and $T$ don't work and
   the fluence calculations still don't work.  So punted and took all that code out.
   
   The important point that must be remembered is that all the angles in this
   program assume that the angles are those actually in the sample.  This allows
   angles greater that the critical angle to be used.  Everything is fine as long
   as the index of refraction of the incident medium is 1.0.  If this is not the
   case then the angle inside the medium must be figured out.
   
@ @<Prototype for |Boundary_RT|@>=
void Boundary_RT(double n_i, double n_g, double n_t, int n, double b, @/
                 double *R, double *T)

@ @<Definition for |Boundary_RT|@>=
    @<Prototype for |Boundary_RT|@>
{
    int i;
    double refl, trans;
    double mu;
    
    for (i = 1; i <= n; i++) {
        if (n_i==1.0) 
            mu = Cos_Snell(n_t, angle[i], n_i);
        else 
            mu = angle[i];
        
        Absorbing_Glass_RT(n_i, n_g, n_t, mu, b, &refl, &trans);
        R[i]   = refl * twoaw[i];
        T[i]   = trans;
    }

}

@*2 Boundary incorporation algorithms.

The next two routines |A_Add_Slide| and |B_Add_Slide| are modifications of
the full addition algorithms for dissimilar layers.  They are optimized to 
take advantage of the diagonal nature of the boundary matrices.  There are
two algorithms below to facilitate adding slides below and above the sample.

@ |A_Add_Slide| computes the resulting |R20| and |T02| matrices for a glass slide
on top of an inhomogeneous layer characterized by |R12|, |R21|, |T12|, |T21|.
It is ok if |R21==R12| and |T12==T21|.  But I do not think that it is required
by this routine.
The result matrices |R20| and |T02| should be independent of the input matrices
None of the input matrices are changed
 
The critical quantites are
$$
 T_{02}=T_{12} (E-R_{10}R_{12} )^{-1} T_{01} 
$$
and
$$
R_{20}=T_{12} (E-R_{10}R_{12})^{-1} R_{10} T_{21} + R_{21} 
$$

@<Prototype for |A_Add_Slide|@>=
static void A_Add_Slide(int n, double **R12, double **R21, double **T12, double **T21, @/
double *R10, double *T01, double **R20, double **T02, @/
double **atemp, double **btemp)

@ @<Definition for |A_Add_Slide|@>=
    @<Prototype for |A_Add_Slide|@>
{
  double **ctemp;

  ctemp = R20;
  Left_Diagonal_Multiply(n, R10, R12, atemp);
  One_Minus(n, atemp);
  Left_Inverse_Multiply(n, atemp, T12, ctemp);
  Right_Diagonal_Multiply(n, ctemp, T01, T02);
  Right_Diagonal_Multiply(n, ctemp, R10, btemp);
  Matrix_Multiply(n, btemp, T21, atemp);
  Matrix_Sum(n, R21, atemp, R20);
}

@ |B_Add_Slide| computes the resulting |R02| and |T20| matrices for a glass slide
on top of an inhomogeneous layer characterized by |R12|, |R21|, |T12|, |T21|.
It is ok if |R21==R12| and |T12==T21|.  But I do not think that it is required
by this routine.
The result matrices |R02| and |T20| should be independent of the input matrices
None of the input matrices are changed
 
The critical equations are
$$
 T_{20}=T_{10} (E-R_{12}R_{10} )^{-1} T_{21} 
$$
and
$$
R_{02}=T_{10} (E-R_{12}R_{10})^{-1} R_{12} T_{01} + R_{01} 
$$

@<Prototype for |B_Add_Slide|@>=
static void B_Add_Slide(int n, double **R12, double **T21, @/
double *R01, double *R10, double *T01, double *T10, @/
double **R02, double **T20, @/
double **atemp, double **btemp)


@ @<Definition for |B_Add_Slide|@>=
    @<Prototype for |B_Add_Slide|@>
{
  double **ctemp;
  int i;
  ctemp = R02;
 
  Right_Diagonal_Multiply(n, R12, R10, atemp);
  One_Minus(n, atemp);
  Diagonal_To_Matrix(n, T10, btemp);
  Left_Inverse_Multiply(n, atemp, btemp, ctemp);
  Matrix_Multiply(n, ctemp, T21, T20);
  Matrix_Multiply(n, ctemp, R12, btemp);
  Right_Diagonal_Multiply(n, btemp, T01, R02);
  for (i=1;i<=n;i++) 
    R02[i][i] += R01[i]/twoaw[i]/twoaw[i];
}

@*2 Routines to incorporate slides.

@ |Add_Top| calculates the reflection and transmission matrices for a slab
with a boundary placed on top of it.
$$\vbox{
\settabs\+|R01|, |R10|, |T01|, |T10| \qquad &R, T for slide assuming 0=air and 1=slab\cr
\+|n|                  &size of matrix\cr
\+|R01|, |R10|, |T01|, |T10| & R, T for slide assuming 0=air and 1=slab\cr
\+|R12|, |R21|, |T12|, |T21| & R, T for slab  assuming 1=slide and 2=?\cr
\+|R02|, |R20|, |T02|, |T20| & calc R, T for both  assuming 0=air and 2=?\cr
\+|atemp|, |btemp|       & previously allocated temporary storage matrices\cr
}$$

@<Prototype for |Add_Top|@>=
void Add_Top(int n, double *R01, double *R10, double *T01, double *T10, @/
                       double **R12, double **R21, double **T12, double **T21, @/
                       double **R02, double **R20, double **T02, double **T20, @/
                       double **atemp, double **btemp)

@
@<Definition for |Add_Top|@>=
    @<Prototype for |Add_Top|@>
{
  A_Add_Slide(n, R12, R21, T12, T21, R10, T01, R20, T02, atemp, btemp);
  B_Add_Slide(n, R12, T21, R01, R10, T01, T10, R02, T20, atemp, btemp);
}

@ |Add_Bottom| calculates the reflection and transmission matrices for a 
slab with a boundary placed beneath it
$$\vbox{
\settabs\+|R01|, |R10|, |T01|, |T10| \qquad &R, T for slab assuming 0=? and 1=slab bottom\cr
\+|n|                  &size of matrix\cr
\+|R01|, |R10|, |T01|, |T10|  & R, T for slab assuming 0=slab top and 1=slab bottom\cr
\+|R12|, |R21|, |T12|, |T21| & R, T for slide  assuming 1=slab bottom and 2=slide bottom\cr
\+|R02|, |R20|, |T02|, |T20| & calc R, T for both  assuming 0=slab top and 2=slide bottom\cr
\+|atemp|, |btemp|       & previously allocated temporary storage matrices\cr
}$$

@<Prototype for |Add_Bottom|@>=
void Add_Bottom(int n, double **R01, double **R10, double **T01, double **T10, @/
                          double *R12, double *R21, double *T12, double *T21, @/
                          double **R02, double **R20, double **T02, double **T20, @/
                          double **atemp, double **btemp)

@
@<Definition for |Add_Bottom|@>=
        @<Prototype for |Add_Bottom|@>
{
  A_Add_Slide(n, R10, R01, T10, T01, R12, T21, R02, T20, atemp, btemp);
  B_Add_Slide(n, R10, T01, R21, R12, T21, T12, R20, T02, atemp, btemp);
}

@*2 Including identical slides.
|Add_Slides| is optimized for a slab with equal boundaries on each side.
|Add_Slides| calculates the reflection and transmission matrices for 
a slab with the same boundary placed
above and below it.  It is assumed that the slab is homogeneous. 
in this case the resulting |R| and |T|
matrices are independent of direction.  There are no constraints 
on |R01|, |R10|, |T01|, and |T10|.  The 
handles for |R| and |T|cannot be equal to those for |R_total| and |T_total|.

$$\vbox{
\settabs\+|R01|, |R10|, |T01|, |T10| \qquad &R, T for slide assuming 0=air and 1=slab\cr
\+|n|                                       &size of matrix\cr
\+|R01|, |R10|, |T01|, |T10|                & R, T for slide assuming 0=air and 1=slab\cr
\+|R|, |T|                                  & R, T for homogeneous slab\cr
\+|R_total|, |T_total|                      & R, T for all 3 with top = bottom boundary\cr
\+|atemp|, |btemp|                          & temporary storage matrices\cr
}$$

If equal boundary conditions exist on both sides of the slab then, by 
symmetry, the transmission and reflection operator for light travelling from the top to 
the bottom are equal to those for light propagating from the bottom to the top.  
Consequently only one set need be calculated.  This leads to a faster method for 
calculating the reflection and transmission for a slab with equal boundary conditions 
on each side.  Let the top boundary be layer 01, the medium layer 12, and the bottom 
layer 23.  The boundary conditions on each side are equal:  $R_{01}=R_{32}$, 
$R_{10}=R_{23}$, $T_{01}=T_{32}$, and $T_{10}=T_{23}$.  
For example the light reflected from layer 01 (travelling 
from boundary 0 to boundary 1) will equal the amount of light reflected from layer 
32, since there is no physical difference between the two cases.  The switch in the 
numbering arises from the fact that light passes from the medium to the outside at the 
top surface by going from 1 to 0, and from 2 to 3 on the bottom surface.  The 
reflection and transmission for the slab with boundary conditions are $R_{30}$ and $ T_{03}$ 
respectively.  These are given by
$$
T_{02} = T_{12}(E-R_{10}R_{12})^{-1}T_{01}
$$
and
$$
R_{20} = T_{12}(E-R_{10}R_{12})^{-1} R_{10}T_{21}+R_{21} 
$$
and
$$
T_{03} = T_{10}(E-R_{20}R_{10})^{-1}T_{02}
$$
and
$$
R_{30} = T_{10}(E-R_{20}R_{10})^{-1} R_{20}T_{01}+R_{01} 
$$
Further increases in efficiency may be made by exploiting the diagonal nature of the 
reflection and transmission operators for an interface, since most matrix/matrix 
multiplications above become vector/matrix multiplications.

@<Prototype for |Add_Slides|@>=
void Add_Slides(int n, double *R01, double *R10, double *T01, double *T10, @/
                          double **R, double **T, @/
                          double **R_total, double **T_total, @/
                          double **atemp, double **btemp)
                          
@ 

@<Definition for |Add_Slides|@>=
    @<Prototype for |Add_Slides|@>
{
  int i;
  double **R12, **R21, **T12, **T21;
  double temp;

  R12 = R;
  R21 = R;
  T21 = T;
  T12 = T;
  Left_Diagonal_Multiply(n, R10, R12, atemp);
  One_Minus(n, atemp);
  Left_Inverse_Multiply(n, atemp, T12, T_total);
  Right_Diagonal_Multiply(n, T_total, R10, btemp);
  Matrix_Multiply(n, btemp, T21, R_total);
  Matrix_Sum(n, R_total, R21, R_total);

  Right_Diagonal_Multiply(n, R_total, R10, atemp);
  One_Minus(n, atemp);
  Matrix_Inverse(n, atemp, btemp);
  Left_Diagonal_Multiply(n, T10, btemp, atemp);
  Matrix_Multiply(n, atemp, T_total, btemp);
  Right_Diagonal_Multiply(n, btemp, T01, T_total);
  Matrix_Multiply(n, atemp, R_total, btemp);
  Right_Diagonal_Multiply(n, btemp, T01, R_total);
  for (i = 1; i <= n; i++) {
    temp = twoaw[i];
    R_total[i][i] += R01[i] / (temp * temp);
  }
}

@*2 Specular R and T.

|Sp_RT| calculates the specular reflection and transmission for light incident
on a slide-slab-slide sandwich.   The sample is characterized by the 
record |slab|.  The total unscattered reflection and transmission for oblique irradiance
(|urx| and |utx|) together with their companions |uru| and |utu| for diffuse irradiance.
The cosine of the incident angle is specified by |slab.cos_angle|.

The way that this routine calculates the diffuse unscattered quantities based on the global
quadrature angles previously set-up.  Consequently, these estimates are not exact.  In fact
if |n=4| then only two quadrature points will actually be used to figure out the diffuse
reflection and transmission (assuming mismatched boundaries).  

This algorithm is pretty simple.  Since the quadrature angles are all chosen assuming
points {\bf inside} the medium, I must calculate the corresponding angle for light
entering from the outside.  If the the cosine of this angle is greater than zero then the
angle does not correspond to a direction in which light is totally internally reflected.
For this ray, I find the unscattered that would be reflected or transmitted from the
slab.  I multiply this by the quadrature angle and weight |twoaw[i]| to get the
total diffuse reflectance and transmittance.

Oh, yes.  The mysterious multiplication by a factor of |n_slab*n_slab| is required
to account for the $n^2$-law of radiance. 

@<Prototype for |Sp_RT|@>=
void Sp_RT(int n, struct AD_slab_type slab, double *ur1, double *ut1, double *uru, double *utu)

@ @<Definition for |Sp_RT|@>=
    @<Prototype for |Sp_RT|@>
    {
    double mu_outside,r, t;
    int i;
    
    *uru = 0;
    *utu = 0;
    
    for(i=1; i<=n; i++){
        mu_outside=Cos_Snell(slab.n_slab, angle[i],1.0);
        if (mu_outside!=0) {
            Sp_mu_RT(slab.n_top_slide, slab.n_slab, slab.n_bottom_slide, slab.b_top_slide,
                     slab.b, slab.b_bottom_slide, mu_outside, &r, &t);
            *uru += twoaw[i] * r;
            *utu += twoaw[i] * t;
        }
    }
    
    Sp_mu_RT(slab.n_top_slide, slab.n_slab, slab.n_bottom_slide, slab.b_top_slide,
             slab.b, slab.b_bottom_slide, slab.cos_angle, ur1, ut1);

  *uru *= slab.n_slab * slab.n_slab;
  *utu *= slab.n_slab * slab.n_slab;
}
