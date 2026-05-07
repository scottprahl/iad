@** IAD Utilities.
\def\sgn{\mathop{\rm sgn}\nolimits}

@(iad_util.c@>=
#include <math.h>
#include <float.h>
#include <stdio.h>
#include "nr_util.h"
#include "ad_globl.h"
#include "ad_frsnl.h"
#include "ad_bound.h"
#include "iad_type.h"
#include "iad_calc.h"
#include "iad_pub.h"
#include "iad_util.h"
unsigned long g_util_debugging = 0;

@h
    @<Definition for |What_Is_B|@>@;
    @<Definition for |Estimate_RT|@>@;
    @<Definition for |a2acalc|@>@;
    @<Definition for |acalc2a|@>@;
    @<Definition for |g2gcalc|@>@;
    @<Definition for |gcalc2g|@>@;
    @<Definition for |b2bcalc|@>@;
    @<Definition for |bcalc2b|@>@;
    @<Definition for |Set_Debugging|@>@;
    @<Definition for |Debug|@>@;
    @<Definition for |sqr|@>@;
    @<Definition for |Calculate_Mua_Musp|@>@;
    @<Definition for |Print_Invert_Type|@>@;
    @<Definition for |Print_Measure_Type|@>@;


@ @(iad_util.h@>=
    @<Prototype for |What_Is_B|@>;
    @<Prototype for |Estimate_RT|@>;
    @<Prototype for |a2acalc|@>;
    @<Prototype for |acalc2a|@>;
    @<Prototype for |g2gcalc|@>;
    @<Prototype for |gcalc2g|@>;
    @<Prototype for |b2bcalc|@>;
    @<Prototype for |bcalc2b|@>;
    @<Prototype for |Set_Debugging|@>;
    @<Prototype for |Calculate_Mua_Musp|@>;
    @<Prototype for |Debug|@>;
    @<Prototype for |sqr|@>;
    @<Prototype for |Print_Invert_Type|@>;
    @<Prototype for |Print_Measure_Type|@>;


@*1 Finding optical thickness.

This routine figures out what the optical thickness of a slab
based on the index of refraction of the slab and the amount
of collimated light that gets through it.

It should be pointed out right here in the front that this
routine does not work for diffuse irradiance, but then the whole
concept of estimating the optical depth for diffuse irradiance
is bogus anyway.


@<Prototype for |What_Is_B|@>=
double What_Is_B(struct AD_slab_type slab, double Tu)

@ @<Definition for |What_Is_B|@>=
    @<Prototype for |What_Is_B|@>

{
    double r1, r2, t1, t2, mu_in_slab;

    @<Calculate specular reflection and transmission@>@;
    @<Check for bad values of |Tu|@>@;
    @<Solve if multiple internal reflections are not present@>@;
    @<Find thickness when multiple internal reflections are present@>@;
}


@ The first thing to do is to find the specular reflection for light
interacting with the top and bottom air-glass-sample interfaces.  I make
a simple check to ensure that the the indices are different before
calculating the bottom reflection.  Most of the time the |r1==r2|,
but there are always those annoying special cases.

@<Calculate specular reflection and transmission@>=
    Absorbing_Glass_RT(1.0, slab.n_top_slide, slab.n_slab,
                       slab.cos_angle, slab.b_top_slide, &r1, &t1);

    mu_in_slab = Cos_Snell(1.0, slab.cos_angle, slab.n_slab);

    Absorbing_Glass_RT(slab.n_slab, slab.n_bottom_slide, 1.0,
                       mu_in_slab, slab.b_bottom_slide, &r2, &t2);


@ Bad values for the unscattered transmission are those that
are non-positive, those greater than one, and those greater than
are possible in a non-absorbing medium, i.e.,
$$
T_c > {t_1 t_2 \over1-r_1r_2}
$$
Since this routine has no way to report errors, I just set the
optical thickness to the natural values in these cases.

@<Check for bad values of |Tu|@>=
    if (Tu <= 0)
        return (HUGE_VAL);

    if (Tu >= t1 * t2 / (1 - r1 * r2))
        return (0.001);

@ If either |r1| or |r2==0| then things are very simple
because the sample does not sustain multiple internal reflections
and the unscattered transmission is
$$
T_c = t_1 t_2 \exp(-b/\nu)
$$
where |b| is the optical thickness and $\nu$ is |slab.cos_angle|.  Clearly,
$$
b = - \nu\ln\left({T_c\over t_1 t_2} \right)
$$

@<Solve if multiple internal reflections are not present@>=
    if (r1 == 0 || r2 == 0)
        return (-slab.cos_angle*log(Tu / t1 / t2));


@ Well I kept putting it off, but now comes the time to solve
the following equation for |b|
$$
T_c = {t_1 t_2\exp(-b)\over 1-r_1r_2 \exp(-2b)}
$$
We note immediately that this is a quadratic equation in
$x=\exp(-b)$.
$$
r_1r_2T_c x^2 + t_1 t_2 x -T_c =0
$$
Sufficient tests
have been made above to ensure that none of the coefficients
are exactly zero. However, it is clear that the leading quadratic term has
a much smaller coefficient than the other two.  Since
$r_1$ and $r_2$ are typically about four percent the product is
roughly $10^{-3}$.  The collimated transmission can be very small
and this makes things even worse.  A further complication is that
we need to choose the only positive root.

Now the roots of $ax^2+bx+c=0$ can be found using the
standard quadratic formula,
$$
x = {-b\pm\sqrt{b^2-4ac}\over 2a}
$$
This is very bad for small values of $a$.  Instead I use
$$
q=-{1\over 2} \left[b+\sgn(b)\sqrt{b^2-4ac}\right]
$$
with the two roots
$$
x={q\over a}\qquad\hbox{and}\qquad x= {c\over q}
$$
Substituting our coefficients
$$
q=-{1\over 2} \left[t_1 t_2+\sqrt{t_1^2 t_2^2+4r_1r_2T_c^2}\right]
$$
With some algebra, this can be shown to be
$$
q = -t_1 t_2 \left[1+{r_1r_2T_c^2\over t_1^2 t_2^2}+\cdots \right]
$$
The only positive root is $x=-T_c/q$.  Therefore
$$
x = { 2 T_c \over t_1 t_2 +\sqrt{t_1^2 t_2^2+4r_1r_2T_c^2}}
$$
(Not very pretty, but straightforward enough.)

@<Find thickness when multiple internal reflections are present@>=

{
    double B;

    B = t1 * t2;
    return (-slab.cos_angle*log(2 * Tu /(B+sqrt(B*B+ 4*Tu * Tu * r1 * r2))));
}

@*1 Estimating R and T.

In several places, it is useful to know an {\it estimate\/} for the values of the
reflection and transmission of the sample based on the measurements.  This
routine provides such an estimate, but it currently ignores anything
corrections that might be made for the integrating spheres.

Good values are only really obtainable when |num_measures==3|, otherwise
we need to make pretty strong assumptions about the reflection and transmission
values.  If |num_measures<3|, then we will assume that no collimated light makes it all
the way through the sample.  The specular reflection is then just that for a
semi-infinite sample and $Tu=0$. If |num_measures==1|, then |Td| is also set
to zero.

{\settabs\+\qquad\qquad&variable&\qquad description of variable&\cr%sample line
\+&|rt|  & total reflection\cr
\+&|rc|  & primary or specular reflection\cr
\+&|rd|  & diffuse or scattered reflection\cr
\+&|tt|  & total transmission\cr
\+&|tp|  & primary or unscattered transmission\cr
\+&|td|  & diffuse or scattered transmission\cr}

@<Prototype for |Estimate_RT|@>=
void Estimate_RT(struct measure_type m, struct invert_type r, double *rt, double *tt,
double *rd, double *rc, double *td, double *tc)

@ @<Definition for |Estimate_RT|@>=
    @<Prototype for |Estimate_RT|@>

{
    @<Calculate the unscattered transmission and reflection@>@;
    @<Estimate the backscattered reflection@>@;
    @<Estimate the scattered transmission@>@;
    @<Debug info for estimate RT@>@;
}

@ If there are three measurements then the specular reflection can
be calculated pretty well.  If there are fewer then
the unscattered transmission is assumed to be zero.  This is not
necessarily the case, but after all, this routine only makes estimates
of the various reflection and transmission quantities.

If there are three measurements, the optical thickness of the sample
is required.  Of course if there are three measurements then the
illumination must be collimated and we can call |What_Is_B| to
find out the optical thickness.  We pass this value to a routine
in the \.{fresnel.h} unit and sit back and wait.

All the above is true if sphere corrections are not needed.
Now, we just fob this off on another function.

@<Calculate the unscattered transmission and reflection@>=

    Calculate_Minimum_MR(m,r,rc,tc);


@ Finding the diffuse reflection is now just a matter of checking
whether V1\% contains the specular reflection from the sample or
not and then just adding or subtracting the specular reflection as
appropriate.

@<Estimate the backscattered reflection@>=
    if (m.fraction_of_ru_in_mr) {
        *rt = m.m_r;
        *rd = *rt - m.fraction_of_ru_in_mr * (*rc);
        if (*rd < 0 ) {
            *rd = 0;
            *rc = *rt;
        }
    }
    else {
        *rd = m.m_r;
        *rt = *rd + *rc;
    }

@ The transmission values follow in much the same way as the
diffuse reflection values --- just subtract the specular
transmission from the total transmission.

@<Estimate the scattered transmission@>=
    if (m.m_u > 0)
        *tc = m.m_u;
    *td = m.m_t - m.fraction_of_tu_in_mt * (*tc);
    *tt = *td + *tc;

@ Collect debugging info here

@<Debug info for estimate RT@>=
    if (0 && Debug(DEBUG_SEARCH)) {
        fprintf(stderr,"SEARCH: r_t = %8.5f ",*rt);
        fprintf(stderr,"r_d = %8.5f ",*rd);
        fprintf(stderr,"r_u = %8.5f\n",*rc);

        fprintf(stderr,"SEARCH: t_t = %8.5f ",*tt);
        fprintf(stderr,"t_d = %8.5f ",*td);
        fprintf(stderr,"t_u = %8.5f\n",*tc);
    }

@*1 Transforming properties.
Routines to convert optical properties to calculation space
and back.

@ |a2acalc| maps albedo $a\in(0,1)$ to an unbounded calculation variable
using the logit (log-odds) transform:
$$
a_{calc} = \ln\!\left({a\over 1-a}\right)
$$
This is the same bijection used internally by scipy's bounded Nelder-Mead
for parameters with both a lower and upper bound.  It is much better
conditioned than the previous $(2a-1)/(a(1-a))$ formula: the logit grows
only logarithmically near the boundaries, so the optimizer sees a
well-scaled landscape across the full range $[0,1]$.

@<Prototype for |a2acalc|@>=
    double a2acalc(double a)

@ @<Definition for |a2acalc|@>=
@<Prototype for |a2acalc|@>
{
    if (a <= 0) return -BIG_A_CALC_VALUE;

    if (a >= 1) return BIG_A_CALC_VALUE;

    return log(a / (1.0 - a));
}

@ |acalc2a| is the inverse of |a2acalc|.
The inverse of the logit is the sigmoid (logistic) function:
$$
a = {1\over 1+e^{-a_{calc}}}
$$

@<Prototype for |acalc2a|@>=
    double acalc2a(double acalc)

@ @<Definition for |acalc2a|@>=
@<Prototype for |acalc2a|@>
{
    if (acalc >= BIG_A_CALC_VALUE)
        return 1.0;

    if (acalc <= -BIG_A_CALC_VALUE)
        return 0.0;

    return 1.0 / (1.0 + exp(-acalc));
}

@ |g2gcalc| maps anisotropy $g\in(-1,1)$ to an unbounded calculation
variable using the inverse hyperbolic tangent (atanh):
$$
g_{calc} = \mathop{\rm atanh}(g) = {1\over2}\ln\!\left({1+g\over 1-g}\right)
$$
Like the logit for albedo, this is the natural bijection for a parameter
bounded symmetrically at $\pm1$ and is far better conditioned near the
boundaries than the previous $g/(1-|g|)$ formula.

@<Prototype for |g2gcalc|@>=
double g2gcalc(double g)

@ @<Definition for |g2gcalc|@>=
@<Prototype for |g2gcalc|@>
{
    double gg = g;
    if (g < -MAX_ABS_G) gg = -MAX_ABS_G;
    if (g >  MAX_ABS_G) gg =  MAX_ABS_G;
    return 0.5 * log((1.0 + gg) / (1.0 - gg));
}

@ |gcalc2g| is the inverse of |g2gcalc|.
The inverse of atanh is tanh:
$$
g = \tanh(g_{calc})
$$
@<Prototype for |gcalc2g|@>=
double gcalc2g(double gcalc)

@ @<Definition for |gcalc2g|@>=
@<Prototype for |gcalc2g|@>
{
    return tanh(gcalc);
}

@ |b2bcalc| is used for the optical depth transformations
it is the inverse of |bcalc2b|.  The relation is
$$
b_{calc} = \ln(b)
$$
The only caveats are to ensure that I don't take the logarithm
of something big or non-positive.

@<Prototype for |b2bcalc|@>=
double b2bcalc(double b)

@ @<Definition for |b2bcalc|@>=
@<Prototype for |b2bcalc|@>
{
    if (b ==  HUGE_VAL) return HUGE_VAL;
    if (b <= 0 ) return 0.0;
    return (log(b));
}

@ |bcalc2b| is used for the anisotropy transformations
it is the inverse of |b2bcalc|.  The relation is
$$
b = \exp(b_{calc})
$$
The only tricky part is to ensure that I don't exponentiate
something big and get an overflow error.  In ANSI \Cee\
the maximum value for $x$ such that $10^x$ is in the range
of representable finite floating point numbers (for doubles)
is given by |DBL_MAX_10_EXP|.  Thus if we want to know if
$$
e^{b_{calc}} > 10^x
$$
or
$$
b_{calc}> x\ln(10) \approx 2.3 x
$$
and this is the criterion that I use.

@<Prototype for |bcalc2b|@>=
double bcalc2b(double bcalc)

@ @<Definition for |bcalc2b|@>=
@<Prototype for |bcalc2b|@>
{
    if (bcalc ==  HUGE_VAL) return HUGE_VAL;
    if (bcalc > 2.3 *  DBL_MAX_10_EXP) return HUGE_VAL;
    return (exp(bcalc));
}

@*1 Some debugging stuff.

@ @<Prototype for |Set_Debugging|@>=
    void Set_Debugging(unsigned long debug_level)

@
@<Definition for |Set_Debugging|@>=
    @<Prototype for |Set_Debugging|@>
{
    g_util_debugging = debug_level;
}

@
@<Prototype for |Debug|@>=
    int Debug(unsigned long mask)

@
@<Definition for |Debug|@>=
    @<Prototype for |Debug|@>
{
    if (g_util_debugging & mask)
        return 1;
    else
        return 0;
}

@
@<Prototype for |sqr|@>=
    double sqr(double x)

@
@<Definition for |sqr|@>=
    @<Prototype for |sqr|@>
{
    return x * x;
}


@ Just figure out the damn scattering and absorption.  This is a nuisance because
|b| may be infinite.

@<Prototype for |Calculate_Mua_Musp|@>=
void Calculate_Mua_Musp(struct measure_type m, struct invert_type r,
                        double *mus, double *musp, double *mua)

@
@<Definition for |Calculate_Mua_Musp|@>=
    @<Prototype for |Calculate_Mua_Musp|@>
{
    if (r.b == HUGE_VAL || isinf(r.b)) {

        if (r.a <= 1e-5) {
            *mus = 0.0;
            *musp = 0.0;
            *mua  = 1.0;
            return;
        }

        if (r.default_mus != UNINITIALIZED) {
            *mus = r.default_mus ;
            *musp = r.default_mus * (1-r.g);
            *mua  = r.default_mus/r.a - r.default_mus;
            return;
        }

        if (r.default_mua != UNINITIALIZED) {
            *mus = r.default_mua / (1-r.a) - r.default_mua;
            *musp = (*mus) * (1-r.g);
            *mua  = r.default_mua;
            return;
        }

        *mus = 1.0;
        *musp = (*mus) * (1-r.g);
        *mua  = (1.0 - r.a) / r.a;
        return;
    }

    *mus  = r.a * r.b / m.slab_thickness;
    *musp = (*mus) * (1 - r.g);
    *mua  = (1 - r.a) * r.b / m.slab_thickness;
}


@
@<Prototype for |Print_Invert_Type|@>=
    void Print_Invert_Type(struct invert_type r)

@
@<Definition for |Print_Invert_Type|@>=
    @<Prototype for |Print_Invert_Type|@>
{
    fprintf(stderr, "\n");
    fprintf(stderr,"default  a=%10.5f   b=%10.5f    g=%10.5f\n",
        r.default_a, r.default_b,r.default_g);
    fprintf(stderr,"slab     a=%10.5f   b=%10.5f    g=%10.5f\n",
        r.slab.a, r.slab.b,r.slab.g);
    fprintf(stderr,"n      top=%10.5f mid=%10.5f  bot=%10.5f\n",
        r.slab.n_top_slide, r.slab.n_slab,r.slab.n_bottom_slide);
    fprintf(stderr,"thick  top=%10.5f cos=%10.5f  bot=%10.5f\n",
        r.slab.b_top_slide, r.slab.cos_angle,r.slab.b_bottom_slide);
    fprintf(stderr,"search = %d quadrature points = %d\n", r.search,r.method.quad_pts );
    fprintf(stderr,"default_a = %10.5f\n", r.default_a );
    fprintf(stderr,"default_b = %10.5f\n", r.default_b );
    fprintf(stderr,"default_g = %10.5f\n", r.default_g );
    fprintf(stderr,"default_mua = %10.5f\n", r.default_mua );
    fprintf(stderr,"default_mus = %10.5f\n", r.default_mus );
}

@
@<Prototype for |Print_Measure_Type|@>=
    void Print_Measure_Type(struct measure_type m)

@
@<Definition for |Print_Measure_Type|@>=
    @<Prototype for |Print_Measure_Type|@>
{
    fprintf(stderr, "\n");
    fprintf(stderr,"#                        Beam diameter = %7.1f mm\n", m.d_beam);
    fprintf(stderr,"#                     Sample thickness = %7.1f mm\n",
            m.slab_thickness );
    fprintf(stderr,"#                  Top slide thickness = %7.1f mm\n",
                    m.slab_top_slide_thickness );
    fprintf(stderr,"#               Bottom slide thickness = %7.1f mm\n",
                    m.slab_bottom_slide_thickness );
    fprintf(stderr,"#           Sample index of refraction = %7.3f\n",
            m.slab_index );
    fprintf(stderr,"#        Top slide index of refraction = %7.3f\n",
            m.slab_top_slide_index );
    fprintf(stderr,"#     Bottom slide index of refraction = %7.3f\n",
            m.slab_bottom_slide_index );
    fprintf(stderr,"#    Fraction unscattered light in M_R = %7.1f %%\n",
    m.fraction_of_ru_in_mr*100);
    fprintf(stderr,"#    Fraction unscattered light in M_T = %7.1f %%\n",
    m.fraction_of_tu_in_mt*100);
    fprintf(stderr,"# \n");
    fprintf(stderr,"# Reflection sphere\n");
    fprintf(stderr,"#                      sphere diameter = %7.1f mm\n",
    m.d_sphere_r );
    fprintf(stderr,"#                 sample port diameter = %7.1f mm\n",
    2*m.d_sphere_r*sqrt(m.as_r) );
    fprintf(stderr,"#               entrance port diameter = %7.1f mm\n",
    2*m.d_sphere_r*sqrt(m.at_r) );
    fprintf(stderr,"#               detector port diameter = %7.1f mm\n",
    2*m.d_sphere_r*sqrt(m.ad_r) );
    fprintf(stderr,"#                     wall reflectance = %7.1f %%\n", m.rw_r*100 );
    fprintf(stderr,"#                 standard reflectance = %7.1f %%\n", m.rstd_r*100 );
    fprintf(stderr,"#                 detector reflectance = %7.1f %%\n", m.rd_r*100 );
    fprintf(stderr,"#                              spheres = %7d\n", m.num_spheres );
    fprintf(stderr,"#                             measures = %7d\n", m.num_measures );
    fprintf(stderr,"#                               method = %7d\n", m.method );
    fprintf(stderr,"area_r as=%10.5f  ad=%10.5f    ae=%10.5f  aw=%10.5f\n",
        m.as_r, m.ad_r, m.at_r, m.aw_r);
    fprintf(stderr,"refls  rd=%10.5f  rw=%10.5f  rstd=%10.5f   f=%10.5f\n",
        m.rd_r, m.rw_r, m.rstd_r, m.f_r);
    fprintf(stderr,"area_t as=%10.5f  ad=%10.5f    ae=%10.5f  aw=%10.5f\n",
        m.as_t, m.ad_t, m.at_t, m.aw_t);
    fprintf(stderr,"refls  rd=%10.5f  rw=%10.5f  rstd=%10.5f\n",
        m.rd_t, m.rw_t, m.rstd_t);
    fprintf(stderr,"lost  ur1=%10.5f ut1=%10.5f   uru=%10.5f  utu=%10.5f\n",
        m.ur1_lost, m.ut1_lost, m.utu_lost, m.utu_lost);
}
