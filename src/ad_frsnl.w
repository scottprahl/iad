@*1 AD Fresnel.
This is a part of the core suite of files for the adding-doubling
program.  Not surprisingly, this program includes routines to
calculate Fresnel reflection.  

@(ad_frsnl.c@>=
#include <math.h>
#include <float.h>
#include <stdio.h>
#include "ad_frsnl.h"

    @<Prototype for |Fresnel|@>;
    @<Prototype for |R1|@>;

    @<Definition for |Cos_Critical_Angle|@>@;
    @<Definition for |Cos_Snell|@>@;
    @<Definition for |Fresnel|@>@;
    @<Definition for |Glass|@>@;
    @<Definition for |Absorbing_Glass_RT|@>@;
    @<Definition for |R1|@>@;
    @<Definition for |Sp_mu_RT|@>@;
    @<Definition for |Sp_mu_RT_Flip|@>@;
    @<Definition for |Diffuse_Glass_R|@>@;

@ @(ad_frsnl.h@>=
    @<Prototype for |Cos_Critical_Angle|@>;
    @<Prototype for |Cos_Snell|@>;
    @<Prototype for |Absorbing_Glass_RT|@>;
    @<Prototype for |Sp_mu_RT|@>;
    @<Prototype for |Sp_mu_RT_Flip|@>;
    @<Prototype for |Diffuse_Glass_R|@>;
    @<Prototype for |Glass|@>;

@*2 The critical angle.

 |Cos_Critical_Angle| calculates the cosine of the critical angle.
If there is no critical angle then 0.0 is returned (i.e., $\cos(\pi/2)$).
Note that no trigonmetric functions are required.  
Recalling Snell's law
$$
n_i \sin\theta_i = n_t\sin\theta_t
$$
To find the critical angle, let $\theta_t=\pi/2$ and then
$$
\theta_c = \sin^{-1} {n_t\over n_i}
$$
The cosine of this angle is then
$$
\cos\theta_c = \cos \left(\sin^{-1} {n_t\over n_i}\right) = {\sqrt{n_i^2-n_t^2}\over n_i}
$$
or more simply
$$
\cos\theta_c =\sqrt{1-n^2}
$$
where $n=n_t/n_i$.

@<Prototype for |Cos_Critical_Angle|@>=
    double Cos_Critical_Angle(double ni, double nt)

@ @<Definition for |Cos_Critical_Angle|@>=
    @<Prototype for |Cos_Critical_Angle|@>
{
    double x;
    
    if (nt >= ni)
        return 0.0;
    else {
        x = nt/ni;
        x = sqrt(1.0 - x*x);
        return x;
    }
}

@*2 Snell's Law.

|Cos_Snell| returns the cosine of the angle that the light propagates through 
a medium given the cosine of the angle of incidence and the indices of refraction.  
Let the cosine of the angle of incidence be $\mu_t$, the transmitted cosine as $\mu_t$,
the index of refraction of the incident material $n_i$ and that of the transmitted material
be $n_t$.

Snell's law states
$$
n_i \sin\theta_i= n_t \sin\theta_t
$$
but if the angles are expressed as cosines, $\mu_i=\cos\theta_i$ then
$$
n_i \sin(\cos^{-1}\mu_i)=n_t \sin(\cos^{-1}\mu_t)
$$
Solving for $\mu_t$ yields
$$
\mu_t=\cos\{\sin^{-1}[(n_i/n_t)\sin(\cos^{-1}\mu_i)]\}
$$
which is pretty ugly.  However, note that $\sin(\cos^{-1}\mu)=\sqrt{1-\mu^2}$
and the above becomes
$$
\mu_t = \sqrt{1-(n_i/n_t)^2 (1-\mu_i^2)}
$$
and no trigonmetric calls are necessary.  Hooray!

A few final notes.  I check to make sure that the index of refraction of
changes before calculating a bunch of stuff.  This routine should 
not be passed incident angles greater
than the critical angle, but I shall program defensively and test
to make sure that the argument of the |sqrt| function is non-negative.
If it is, then I return $\mu_t=0$ i.e., $\theta_t=90^\circ$.

I also pretest for the common but trivial case of normal incidence.

@<Prototype for |Cos_Snell|@>=
    double Cos_Snell(double n_i, double mu_i, double n_t)

@ @<Definition for |Cos_Snell|@>=
    @<Prototype for |Cos_Snell|@>
{
        double temp;
        
        if (mu_i==1.0) return 1.0;
        
        if (n_i==n_t)
            return mu_i;

        temp = n_i/n_t;
        temp = 1.0-temp*temp*(1.0 - mu_i*mu_i);
        if (temp<0)
            return 0.0;
        else
            return (sqrt(temp));
}

@*2 Fresnel Reflection.

|Fresnel| calculates the specular reflection for light incident at
an angle $\theta_i$ from the normal (having a cosine equal to $\mu_i$) 
in a medium with index of
refraction |n_i| onto a medium with index of refraction |n_t| .

The usual way to calculate the total reflection for unpolarized light is
to use the Fresnel formula
$$
R = {1\over 2}\left[ {\sin^2(\theta_i-\theta_t)\over \sin^2(\theta_i+\theta_t)}
                     +{\tan^2(\theta_i-\theta_t)\over \tan^2(\theta_i+\theta_t)} \right]
$$
where $\theta_i$ and $\theta_t$ represent the angle (from normal) that light is incident
and the angle at which light is transmitted.  
There are several problems with calculating the reflection using this formula.
First, if the angle of incidence is zero, then the formula results in division by zero.
Furthermore, if the angle of incidence is near zero, then the formula is the ratio
of two small numbers and the results can be inaccurate.
Second, if the angle of incidence exceeds the critical angle, then the calculation of
$\theta_t$ results in an attempt to find the arcsine of a quantity greater than
one.  Third, all calculations in this program are based on the cosine of the angle.
This routine forces the calling routine to find $\theta_i=\cos^{-1} \mu$.  
Fourth, the routine also gives problems when the critical angle is exceeded.

Closer inspection reveals that this is the wrong formulation to use.  The formulas that
should be used for parallel and perpendicular polarization are
$$
R_\parallel =\left[{n_t\cos\theta_i-n_i\cos\theta_t\over
                                        n_t\cos\theta_i+n_i\cos\theta_t}\right]^2,
\qquad\qquad
R_\perp =\left[ {n_i\cos\theta_i-n_t\cos\theta_t\over
                                        n_i\cos\theta_i+n_t\cos\theta_t}\right]^2.
$$
The formula for unpolarized light, written in terms of $\mu_i=\cos\theta_i$ and
$\mu_t=\cos\theta_t$ is
$$
R={1\over 2}\left[{n_t\mu_i-n_i\mu_t\over n_t\mu_i+n_i\mu_t}\right]^2
+{1\over 2}\left[{n_i\mu_i-n_t\mu_t\over n_i\mu_i+n_t\mu_t}\right]^2
$$

This formula has the advantage that no trig routines need to be called and that the
case of normal irradiance does not cause division by zero.  Near normal incidence
remains numerically well-conditioned.  In the routine below, I test for matched 
boundaries and normal incidence to eliminate unnecessary calculations.  I also
test for total internal reflection to avoid possible division by zero.  I also
find the ratio of the indices of refraction to avoid an extra multiplication and
several intermediate variables.
\goodbreak

@ @<Prototype for |Fresnel|@>=
    static double Fresnel(double n_i, double n_t, double mu_i)

@ @<Definition for |Fresnel|@>=
@<Prototype for |Fresnel|@>@;
{
  double mu_t,ratio, temp,temp1;

  if (n_i == n_t)
    return 0.0;

  if (mu_i==1.0) {
    temp = (n_i-n_t)/(n_i+n_t);
    return (temp*temp);
}
    
  if (mu_i == 0.0)
    return 1.0;

    mu_t=Cos_Snell(n_i,mu_i,n_t);
    if (mu_t==0.0)
        return 1.0;
    ratio = n_i/n_t;
    temp=ratio*mu_t;
    temp1 = (mu_i-temp)/(mu_i+temp);
    temp= ratio*mu_i;
    temp = (mu_t-temp)/(mu_t+temp);
    return ( (temp1*temp1+temp*temp)/2);
}

@*2 Reflection from a glass slide.

|Glass| calculates the total specular reflection (i.e., including
multiple internal reflections) based on                  
the indices of refraction of the incident medium |n_i|, the glass |n_g|,     
and medium into which the light is transmitted |n_t| for light incident at
an angle from the normal having cosine |mu_i|.                   

In many tissue optics problems, the sample is constrained by a piece of glass
creating an air-glass-tissue sequence.
The adding-doubling formalism can calculate the effect that the layer of glass will
have on the radiative transport properties by including a layer for the glass-tissue
interface and a layer for the air-glass interface.  However, it is simpler to find net
effect of the glass slide and include only one layer for the glass boundary.  

The first time I implemented this routine, I did not include multiple internal
reflections.  After running test cases, it soon became apparent that  the 
percentage errors were way too
big for media with little absorption and scattering.  It is not hard to find the
result for the reflection from a non-absorbing glass layer (equation A2.21 
in my dissertation) in which multiple reflections are properly accounted for
$$
r_g = {r_1 + r_2 - 2  r_1  r_2 \over 1 - r_1  r_2}
$$
Here $r_1$ is the reflection at the air-glass interface and $r_2$ is the
reflection at the glass-sample interface.

There is one pitfall in calculating $r_g$.  When the angle
of incidence exceeds the critical angle then the formula above causes
division by zero.  If this is the case then $r_1=1$ and can easily
be tested for.

To eliminate unnecessary computation, I check to make sure that 
it really is necessary to call the |Fresnel| routine twice.  
It is noteworthy that the formula for $r_g$ works correctly if the
the first boundary is not totally reflecting but the second one is.
Note that $\mu_g$ gets calculated twice
(once in the first call to |Fresnel| and once directly).

@<Prototype for |Glass|@>=
    double Glass(double n_i, double n_g, double n_t, double mu_i)

@ @<Definition for |Glass|@>=
@<Prototype for |Glass|@>@;
{
  double r1, r2, mu_g,temp;

    if (n_i==n_g) return (Fresnel(n_g,n_t,mu_i));
  
    r1 = Fresnel(n_i, n_g, mu_i);
    if (r1 >= 1.0 || n_g==n_t) return r1;

    mu_g=Cos_Snell(n_i,mu_i,n_g);
        r2 = Fresnel(n_g, n_t, mu_g);
    temp = r1*r2;
    temp = (r1 + r2 - 2 * temp) / (1 - temp);
        return temp;
}

@*2 Reflection from an absorbing slide.

|Absorbing_Glass_RT| calculates the total specular reflection and transmission
(i.e., including
multiple internal reflections) based on                  
the indices of refraction of the incident medium |n_i|, the glass |n_g|,     
and medium into which the light is transmitted |n_t| for light incident at
an angle from the normal having cosine |mu_i|.  The optical thickness of
the glass $b=\mu_a d$ is measured normal to the glass.

This routine was generated to help solve a problem with the inverse adding-doubling
program associated with samples with low absorbances.  A particular situation
arises when the slides have significant absorption relative to the sample
absorption.  Anyway, it is not hard to extend the result for non-absorbing slides
to the absorbing case
$$
r = {r_1 + (1-2r_1)r_2 \exp(-2b/\mu_g)  \over 1 - r_1  r_2 \exp(-2b/\mu_g)}
$$
Here $r_1$ is the reflection at the sample-glass interface and $r_2$ is the
reflection at the glass-air interface and $\mu_g$ is the cosine of the
angle inside the glass.  Note that if $b\ne0$ then the reflection depends
on the order of the indices of refraction, otherwise |n_i| and |n_t|
can be switched and the result should be the same.

The corresponding result for transmission is
$$
t = {(1-r_1) (1-r_2) \exp(-b/\mu_g)  \over 1 - r_1  r_2 \exp(-2b/\mu_g)}
$$

There are two potential pitfalls in the calculation.  The first is
when the angle of incidence exceeds the critical angle then the formula above causes
division by zero.  If this is the case, |Fresnel| will return $r_1=1$ and 
this routine responds appropriately.  The second case is when the optical
thickness of the slide is too large.  

I don't worry too much about optimal coding, because this routine does
not get called all that often and also because |Fresnel| is pretty good
at avoiding unnecessary computations.  At worst this routine just has
a couple of extra function calls and a few extra multiplications.

I also check to make sure that the exponent is not too small.

@<Prototype for |Absorbing_Glass_RT|@>=
    void Absorbing_Glass_RT(double n_i, double n_g, double n_t, double mu_i,double b,
                             double *r, double *t)

@ @<Definition for |Absorbing_Glass_RT|@>=
@<Prototype for |Absorbing_Glass_RT|@>@;
{
  double r1, r2, mu_g, expo, denom;
    *t = 0;
    
    *r = Fresnel(n_i, n_g, mu_i);
    if (*r >= 1.0 || b == HUGE_VAL || mu_i == 0.0)  return;

    mu_g=Cos_Snell(n_i,mu_i,n_g);
    r1 = *r;
    r2 = Fresnel(n_g, n_t, mu_g);
    
    if (b==0.0) {
        *r = (r1 + r2 - 2.0 * r1 *r2) / (1- r1 * r2);
        *t = 1.0 - (*r);
    } else {
        expo = - b/mu_g;
        if (2 * expo  <= DBL_MIN_10_EXP * 2.3025851) return;
        expo = exp(expo);

        denom = 1.0-r1*r2*expo*expo;
        *r = (r1 + (1.0 - 2.0 * r1)*r2*expo*expo) / denom;
        *t = (1.0-r1)*(1.0 -r2)*expo / denom;
    }
}


@*2 Unscattered refl and trans for a sample.

@ |Sp_mu_RT_Flip| finds the reflectance to incorporate flipping of the sample.  This
is needed when the sample is flipped between measurements.  

@<Prototype for |Sp_mu_RT_Flip|@>=
void Sp_mu_RT_Flip(int flip, double n_top, double n_slab, double n_bottom, 
                        double tau_top, double tau_slab, double tau_bottom, double mu, 
                        double *r, double *t)

@ @<Definition for |Sp_mu_RT_Flip|@>=
    @<Prototype for |Sp_mu_RT_Flip|@>
{
    Sp_mu_RT(n_top, n_slab, n_bottom, tau_top, tau_slab, tau_bottom, mu, r, t);
    if (flip && n_top != n_bottom && tau_top != tau_bottom) {
    	double correct_r = *r;
    	Sp_mu_RT(n_bottom, n_slab, n_top, tau_bottom, tau_slab, tau_top, mu, r, t);
    	*r = correct_r;
    }
}

@ |Sp_mu_RT| calculates the unscattered reflection and transmission (i.e., specular) 
through a glass-slab-glass sandwich.  Light is incident at
an angle having a cosine |mu| from air onto a possibly absorbing glass plate with index |n_top|
on a sample with index |n_slab| resting on another possibly absorbing glass plate with index
|n_bottom| and then exiting into air again.

The optical thickness of the slab is |tau_slab|.

@<Prototype for |Sp_mu_RT|@>=
    void Sp_mu_RT(double n_top, double n_slab, double n_bottom, 
                        double tau_top, double tau_slab, double tau_bottom, double mu, 
                        double *r, double *t)

@ @<Definition for |Sp_mu_RT|@>=
    @<Prototype for |Sp_mu_RT|@>
{
    double r_top, r_bottom, t_top, t_bottom, mu_slab, beer, denom, temp, mu_in_slab;
    *r=0;
    *t=0;
    Absorbing_Glass_RT(1.0, n_top, n_slab, mu, tau_top, &r_top, &t_top);
    
	mu_in_slab = Cos_Snell(1.0, mu, n_slab);
    Absorbing_Glass_RT(n_slab, n_bottom, 1.0, mu_in_slab, tau_bottom, &r_bottom, &t_bottom);
    
    @<Calculate |beer|@>@;
    @<Calculate |r| and |t|@>@;
}

@ Nothing tricky here except a check to make sure that the
reflection for the top is not equal to that on the bottom before
calculating it again.  I also drop out of the routine if the top
surface is totally reflecting.

@ I am careful here not to cause an underflow error and to avoid
division by zero.

It turns out that I found a small error in this code fragment.  Basically
I misunderstood what one of the values in \.{float.h} represented.  This
version is now correct

@<Calculate |beer|@>=
    mu_slab=Cos_Snell(1.0, mu, n_slab);

  if (mu_slab == 0)
    beer = 0.0;
  else if (tau_slab == HUGE_VAL)
    beer = 0.0;
  else {
    temp = -tau_slab/mu_slab;
    if (2*temp  <= DBL_MIN_10_EXP * 2.3025851) 
        beer = 0.0;
    else
        beer = exp(temp);
    }
    
@ If $r_{\rm top}$ is the reflection for the top and $r_{\rm bottom}$ is that for the
bottom surface then the total reflection will be
$$
r = r_{\rm top} + {r_{\rm bottom}t_{\rm top}^2 \exp(-2\tau/\mu) \over 
1 - r_{\rm top}  r_{\rm bottom} \exp(-2\tau/\mu)}
$$
and the transmission is
$$
t= {t_{\rm top} t_{\rm bottom} \exp(-\tau/\mu)\over 
1 - r_{\rm top}  r_{\rm bottom} \exp(-2\tau/\mu)}
$$
where $\mu$ is the angle inside the slab and $\tau$ is the optical
thickness of the slab.

I have already calculated the reflections and the exponential attenuation, so
I can just plug into the formula after making sure that it is really necessary.
The denominator cannot be zero since I know |r_top<1| and
that |r_bottom| and |beer| are less than or equal to one.

The bug that was fixed was in the calculated reflection.  I omitted a $r_{\rm bottom}$
in the numerator of the fraction used to calculate the reflection. 

@<Calculate |r| and |t|@>=
    if (beer==0.0){
        *r=r_top;
    } else {
        temp = t_top*beer;
        denom = 1 - r_top * r_bottom * beer * beer;
        *r = r_top + r_bottom*temp * temp  / denom;
        *t = t_bottom*temp/ denom;
  }

@*2 Total diffuse reflection.

|R1| calculates the first moment of the Fresnel reflectance using the analytic 
 solution of Walsh.
The integral of the first moment of the Fresnel reflection ($R_1$) 
has been found analytically by Walsh, [see Ryde 1931]
$$
\eqalign{
R_1 &= {1\over2} + {(m-1)(3m+1)\over 6(m+1)^2} 
        +\left[ {m^2(m^2-1)^2\over(m^2+1)^3}\right]\log\left ( {m-1\over m+1} \right)\cr
        &\qquad- {2m^3 (m^2+2m-1)\over (m^2+1)(m^4-1)} +
        \left[ {8m^4(m^4+1)\over(m^2+1)(m^4-1)^2}\right]\log m
}
$$
where Walsh's parameter $m=n_t/n_i$.    This equation is only valid when 
$n_i<n_t$.  If $n_i>n_t$ then using (see Egan and Hilgeman 1973),
$$
{1-R_1(n_i/n_t)\over n_t^2} = {1-R_1(n_t/n_i)\over n_i^2}
$$
 or
 $$
 R(1/m) = 1-m^2[1-R(m)]
 $$
 
@<Prototype for |R1|@>=
    static double R1(double ni, double nt)

@ @<Definition for |R1|@>=
    @<Prototype for |R1|@>
{
  double m, m2, m4, mm1, mp1, r, temp;

    if (ni==nt) 
        return 0.0;
    
    if (ni<nt)
        m=nt/ni;
    else
        m=ni/nt;

    m2 = m * m;
    m4 = m2 * m2;
    mm1 = m - 1;
    mp1 = m + 1;
    temp = (m2 - 1)/(m2 + 1);

    r = 0.5 + mm1 * (3 * m + 1) / 6 / mp1 / mp1;
    r += m2 * temp*temp / (m2 + 1) * log(mm1 / mp1);
    r -= 2 * m * m2 * (m2 + 2 * m - 1) / (m2 + 1) / (m4 - 1);
    r += 8 * m4 * (m4 + 1) / (m2 + 1) / (m4-1)/(m4-1) * log(m);

    if (ni < nt)
        return r;
    else
        return (1 - (1 - r) / m2);
}

@*2 Diffusion reflection from a glass slide.

|Diffuse_Glass_R| returns the total diffuse specular reflection 
from the air-glass-tissue interface

@<Prototype for |Diffuse_Glass_R|@>=
double Diffuse_Glass_R(double nair, double nslide, double nslab)

@ @<Definition for |Diffuse_Glass_R|@>=
    @<Prototype for |Diffuse_Glass_R|@>
{
  double rairglass, rglasstissue, rtemp;

  rairglass = R1(nair, nslide);
  rglasstissue = R1(nslide, nslab);
  rtemp = rairglass * rglasstissue;
  if (rtemp >=1) 
    return 1.0;
  else
    return ((rairglass + rglasstissue - 2 * rtemp) / (1 - rtemp));
}
