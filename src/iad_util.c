/*1:*/
#line 6 "iad_util.w"

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
unsigned long g_util_debugging= 0;

#define BIG_A_VALUE 999999.0
#define SMALL_A_VALUE 0.000001 \


#line 20 "iad_util.w"

/*4:*/
#line 80 "iad_util.w"

/*3:*/
#line 77 "iad_util.w"

double What_Is_B(struct AD_slab_type slab,double Tc)

/*:3*/
#line 81 "iad_util.w"


{
double r1,r2,t1,t2,mu_in_slab;

/*5:*/
#line 99 "iad_util.w"

Absorbing_Glass_RT(1.0,slab.n_top_slide,slab.n_slab,
slab.cos_angle,slab.b_top_slide,&r1,&t1);

mu_in_slab= Cos_Snell(1.0,slab.cos_angle,slab.n_slab);

Absorbing_Glass_RT(slab.n_slab,slab.n_bottom_slide,1.0,
mu_in_slab,slab.b_bottom_slide,&r2,&t2);


/*:5*/
#line 86 "iad_util.w"

/*6:*/
#line 118 "iad_util.w"

if(Tc<=0)
return(HUGE_VAL);

if(Tc>=t1*t2/(1-r1*r2))
return(0.001);

/*:6*/
#line 87 "iad_util.w"

/*7:*/
#line 136 "iad_util.w"

if(r1==0||r2==0)
return(-slab.cos_angle*log(Tc/t1/t2));


/*:7*/
#line 88 "iad_util.w"

/*8:*/
#line 187 "iad_util.w"


{
double B;

B= t1*t2;
return(-slab.cos_angle*log(2*Tc/(B+sqrt(B*B+4*Tc*Tc*r1*r2))));
}

/*:8*/
#line 89 "iad_util.w"

}


/*:4*/
#line 21 "iad_util.w"

/*10:*/
#line 222 "iad_util.w"

/*9:*/
#line 218 "iad_util.w"

void Estimate_RT(struct measure_type m,struct invert_type r,double*rt,double*tt,
double*rd,double*rc,double*td,double*tc)

/*:9*/
#line 223 "iad_util.w"


{
/*11:*/
#line 246 "iad_util.w"


Calculate_Minimum_MR(m,r,rc,tc);


/*:11*/
#line 226 "iad_util.w"

/*12:*/
#line 256 "iad_util.w"

if(m.fraction_of_rc_in_mr){
*rt= m.m_r;
*rd= *rt-m.fraction_of_rc_in_mr*(*rc);
if(*rd<0){
*rd= 0;
*rc= *rt;
}
}
else{
*rd= m.m_r;
*rt= *rd+*rc;
}
if(Debug(DEBUG_SEARCH)){
fprintf(stderr,"        rt = %.5f\n",*rt);
fprintf(stderr,"    est rd = %.5f\n",*rd);
}

/*:12*/
#line 227 "iad_util.w"

/*13:*/
#line 278 "iad_util.w"

if(m.num_measures==1){
*tt= 0.0;
*td= 0.0;
}
else if(m.fraction_of_tc_in_mt){
*tt= m.m_t;
*td= *tt-*tc;
if(*td<0){
*tc= *tt;
*td= 0;
}
}
else{
*td= m.m_t;
*tt= *td+*tc;
}
if(Debug(DEBUG_SEARCH)){
fprintf(stderr,"        tt = %.5f\n",*tt);
fprintf(stderr,"    est td = %.5f\n",*td);
}

/*:13*/
#line 228 "iad_util.w"

}

/*:10*/
#line 22 "iad_util.w"

/*16:*/
#line 316 "iad_util.w"

/*15:*/
#line 313 "iad_util.w"

double a2acalc(double a)

/*:15*/
#line 317 "iad_util.w"

{
if(a<=0)return-BIG_A_VALUE;

if(a>=1)return BIG_A_VALUE;
return((2*a-1)/a/(1-a));
}

/*:16*/
#line 23 "iad_util.w"

/*18:*/
#line 345 "iad_util.w"

/*17:*/
#line 342 "iad_util.w"

double acalc2a(double acalc)

/*:17*/
#line 346 "iad_util.w"

{
if(acalc==BIG_A_VALUE)
return 1.0;
else
if(acalc==-BIG_A_VALUE)
return 0.0;
else
if(fabs(acalc)<SMALL_A_VALUE)
return 0.5;
else
return((-2+acalc+sqrt(acalc*acalc+4))/(2*acalc));
}

/*:18*/
#line 24 "iad_util.w"

/*20:*/
#line 370 "iad_util.w"

/*19:*/
#line 367 "iad_util.w"

double g2gcalc(double g)

/*:19*/
#line 371 "iad_util.w"

{
if(g<=-1)return(-HUGE_VAL);

if(g>=1)return(HUGE_VAL);

return(g/(1-fabs(g)));
}

/*:20*/
#line 25 "iad_util.w"

/*22:*/
#line 388 "iad_util.w"

/*21:*/
#line 385 "iad_util.w"

double gcalc2g(double gcalc)

/*:21*/
#line 389 "iad_util.w"

{
if(gcalc==-HUGE_VAL)return-1.0;
if(gcalc==HUGE_VAL)return 1.0;
return(gcalc/(1+fabs(gcalc)));
}

/*:22*/
#line 26 "iad_util.w"

/*24:*/
#line 407 "iad_util.w"

/*23:*/
#line 404 "iad_util.w"

double b2bcalc(double b)

/*:23*/
#line 408 "iad_util.w"

{
if(b==HUGE_VAL)return HUGE_VAL;
if(b<=0)return 0.0;
return(log(b));
}

/*:24*/
#line 27 "iad_util.w"

/*26:*/
#line 437 "iad_util.w"

/*25:*/
#line 434 "iad_util.w"

double bcalc2b(double bcalc)

/*:25*/
#line 438 "iad_util.w"

{
if(bcalc==HUGE_VAL)return HUGE_VAL;
if(bcalc> 2.3*DBL_MAX_10_EXP)return HUGE_VAL;
return(exp(bcalc));
}

/*:26*/
#line 28 "iad_util.w"

/*28:*/
#line 452 "iad_util.w"

/*27:*/
#line 449 "iad_util.w"

void twoprime(double a,double b,double g,double*ap,double*bp)

/*:27*/
#line 453 "iad_util.w"

{
if(a==1&&g==1)
*ap= 0.0;
else
*ap= (1-g)*a/(1-a*g);

if(b==HUGE_VAL)
*bp= HUGE_VAL;
else
*bp= (1-a*g)*b;
}

/*:28*/
#line 29 "iad_util.w"

/*30:*/
#line 472 "iad_util.w"

/*29:*/
#line 469 "iad_util.w"

void twounprime(double ap,double bp,double g,double*a,double*b)

/*:29*/
#line 473 "iad_util.w"

{
*a= ap/(1-g+ap*g);
if(bp==HUGE_VAL)
*b= HUGE_VAL;
else
*b= (1+ap*g/(1-g))*bp;
}

/*:30*/
#line 30 "iad_util.w"

/*32:*/
#line 491 "iad_util.w"

/*31:*/
#line 488 "iad_util.w"

void abgg2ab(double a1,double b1,double g1,double g2,double*a2,double*b2)

/*:31*/
#line 492 "iad_util.w"

{
double a,b;

twoprime(a1,b1,g1,&a,&b);
twounprime(a,b,g2,a2,b2);
}

/*:32*/
#line 31 "iad_util.w"

/*34:*/
#line 512 "iad_util.w"

/*33:*/
#line 509 "iad_util.w"

void abgb2ag(double a1,double b1,double b2,double*a2,double*g2)

/*:33*/
#line 513 "iad_util.w"

{
if(b1==0||b2==0){
*a2= a1;
*g2= 0;
}

if(b2<b1)
b2= b1;

if(a1==0)*a2= 0.0;
else{
if(a1==1)
*a2= 1.0;
else{
if(b1==0||b2==HUGE_VAL)
*a2= a1;
else
*a2= 1+b1/b2*(a1-1);
}
}
if(*a2==0||b2==0||b2==HUGE_VAL)
*g2= 0.5;
else
*g2= (1-b1/b2)/(*a2);
}

/*:34*/
#line 32 "iad_util.w"

/*41:*/
#line 613 "iad_util.w"

/*40:*/
#line 610 "iad_util.w"

void quick_guess(struct measure_type m,struct invert_type r,double*a,double*b,double*g)

/*:40*/
#line 614 "iad_util.w"

{
double UR1,UT1,rd,td,tc,rc,bprime,aprime,alpha,beta,logr;

Estimate_RT(m,r,&UR1,&UT1,&rd,&rc,&td,&tc);
/*42:*/
#line 636 "iad_util.w"

if(UT1==1)
aprime= 1.0;
else if(rd/(1-UT1)>=0.1)
{
double tmp= (1-rd-UT1)/(1-UT1);
aprime= 1-4.0/9.0*tmp*tmp;
}
else if(rd<0.05&&UT1<0.4)
aprime= 1-(1-10*rd)*(1-10*rd);
else if(rd<0.1&&UT1<0.4)
aprime= 0.5+(rd-0.05)*4;
else
{
double tmp= (1-4*rd-UT1)/(1-UT1);
aprime= 1-tmp*tmp;
}

/*:42*/
#line 619 "iad_util.w"


switch(m.num_measures){
case 1:
/*44:*/
#line 670 "iad_util.w"

*g= r.default_g;
*a= aprime/(1-*g+aprime*(*g));
*b= HUGE_VAL;

/*:44*/
#line 623 "iad_util.w"

break;
case 2:
/*45:*/
#line 675 "iad_util.w"

/*43:*/
#line 654 "iad_util.w"

if(rd<0.01){
bprime= What_Is_B(r.slab,UT1);
fprintf(stderr,"low rd<0.01! ut1=%f aprime=%f bprime=%f\n",UT1,aprime,bprime);
}else if(UT1<=0)
bprime= HUGE_VAL;
else if(UT1> 0.1)
bprime= 2*exp(5*(rd-UT1)*log(2.0));
else{
alpha= 1/log(0.05/1.0);
beta= log(1.0)/log(0.05/1.0);
logr= log(UR1);
bprime= log(UT1)-beta*log(0.05)+beta*logr;
bprime/= alpha*log(0.05)-alpha*logr-1;
}

/*:43*/
#line 676 "iad_util.w"


*g= r.default_g;
*a= aprime/(1-*g+aprime**g);
*b= bprime/(1-*a**g);

/*:45*/
#line 626 "iad_util.w"

break;
case 3:
/*46:*/
#line 682 "iad_util.w"


switch(r.search){
case FIND_A:
/*47:*/
#line 699 "iad_util.w"


*g= r.default_g;
*a= aprime/(1-*g+aprime**g);
*b= What_Is_B(r.slab,m.m_u);

/*:47*/
#line 686 "iad_util.w"

break;
case FIND_B:
/*48:*/
#line 705 "iad_util.w"


*g= r.default_g;
*a= 0.0;
*b= What_Is_B(r.slab,m.m_u);

/*:48*/
#line 689 "iad_util.w"

break;
case FIND_AB:
/*49:*/
#line 711 "iad_util.w"


*g= r.default_g;

if(*g==1)
*a= 0.0;
else
*a= aprime/(1-*g+aprime**g);

/*43:*/
#line 654 "iad_util.w"

if(rd<0.01){
bprime= What_Is_B(r.slab,UT1);
fprintf(stderr,"low rd<0.01! ut1=%f aprime=%f bprime=%f\n",UT1,aprime,bprime);
}else if(UT1<=0)
bprime= HUGE_VAL;
else if(UT1> 0.1)
bprime= 2*exp(5*(rd-UT1)*log(2.0));
else{
alpha= 1/log(0.05/1.0);
beta= log(1.0)/log(0.05/1.0);
logr= log(UR1);
bprime= log(UT1)-beta*log(0.05)+beta*logr;
bprime/= alpha*log(0.05)-alpha*logr-1;
}

/*:43*/
#line 720 "iad_util.w"

if(bprime==HUGE_VAL||*a**g==1)
*b= HUGE_VAL;
else
*b= bprime/(1-*a**g);

/*:49*/
#line 692 "iad_util.w"

break;
case FIND_AG:
/*50:*/
#line 726 "iad_util.w"

*b= What_Is_B(r.slab,m.m_u);
if(*b==HUGE_VAL||*b==0){
*a= aprime;
*g= r.default_g;
}else{
/*43:*/
#line 654 "iad_util.w"

if(rd<0.01){
bprime= What_Is_B(r.slab,UT1);
fprintf(stderr,"low rd<0.01! ut1=%f aprime=%f bprime=%f\n",UT1,aprime,bprime);
}else if(UT1<=0)
bprime= HUGE_VAL;
else if(UT1> 0.1)
bprime= 2*exp(5*(rd-UT1)*log(2.0));
else{
alpha= 1/log(0.05/1.0);
beta= log(1.0)/log(0.05/1.0);
logr= log(UR1);
bprime= log(UT1)-beta*log(0.05)+beta*logr;
bprime/= alpha*log(0.05)-alpha*logr-1;
}

/*:43*/
#line 732 "iad_util.w"

*a= 1+bprime*(aprime-1)/(*b);
if(*a<0.1)
*g= 0.0;
else
*g= (1-bprime/(*b))/(*a);
}

/*:50*/
#line 695 "iad_util.w"

break;
}

/*:46*/
#line 629 "iad_util.w"

break;
}

/*51:*/
#line 740 "iad_util.w"

if(*a<0)
*a= 0.0;
if(*g<0)
*g= 0.0;
else if(*g>=1)
*g= 0.5;

/*:51*/
#line 633 "iad_util.w"

}

/*:41*/
#line 33 "iad_util.w"

/*54:*/
#line 754 "iad_util.w"

/*53:*/
#line 750 "iad_util.w"

void Set_Debugging(unsigned long debug_level)

/*:53*/
#line 755 "iad_util.w"

{
g_util_debugging= debug_level;
}

/*:54*/
#line 34 "iad_util.w"

/*56:*/
#line 765 "iad_util.w"

/*55:*/
#line 761 "iad_util.w"

int Debug(unsigned long mask)

/*:55*/
#line 766 "iad_util.w"

{
if(g_util_debugging&mask)
return 1;
else
return 0;
}

/*:56*/
#line 35 "iad_util.w"

/*58:*/
#line 779 "iad_util.w"

/*57:*/
#line 775 "iad_util.w"

void Print_Invert_Type(struct invert_type r)

/*:57*/
#line 780 "iad_util.w"

{
fprintf(stderr,"\n");
fprintf(stderr,"default  a=%10.5f   b=%10.5f    g=%10.5f\n",
r.default_a,r.default_b,r.default_g);
fprintf(stderr,"slab     a=%10.5f   b=%10.5f    g=%10.5f\n",
r.slab.a,r.slab.b,r.slab.g);
fprintf(stderr,"n      top=%10.5f mid=%10.5f  bot=%10.5f\n",
r.slab.n_top_slide,r.slab.n_slab,r.slab.n_bottom_slide);
fprintf(stderr,"thick  top=%10.5f cos=%10.5f  bot=%10.5f\n",
r.slab.b_top_slide,r.slab.cos_angle,r.slab.b_bottom_slide);
fprintf(stderr,"search = %d quadrature points = %d\n",r.search,r.method.quad_pts);
}

/*:58*/
#line 36 "iad_util.w"

/*60:*/
#line 799 "iad_util.w"

/*59:*/
#line 795 "iad_util.w"

void Print_Measure_Type(struct measure_type m)

/*:59*/
#line 800 "iad_util.w"

{
fprintf(stderr,"\n");
fprintf(stderr,"#                        Beam diameter = %7.1f mm\n",m.d_beam);
fprintf(stderr,"#                     Sample thickness = %7.1f mm\n",
m.slab_thickness);
fprintf(stderr,"#                  Top slide thickness = %7.1f mm\n",
m.slab_top_slide_thickness);
fprintf(stderr,"#               Bottom slide thickness = %7.1f mm\n",
m.slab_bottom_slide_thickness);
fprintf(stderr,"#           Sample index of refraction = %7.3f\n",
m.slab_index);
fprintf(stderr,"#        Top slide index of refraction = %7.3f\n",
m.slab_top_slide_index);
fprintf(stderr,"#     Bottom slide index of refraction = %7.3f\n",
m.slab_bottom_slide_index);
fprintf(stderr,"#    Fraction unscattered light in M_R = %7.1f %%\n",
m.fraction_of_rc_in_mr*100);
fprintf(stderr,"#    Fraction unscattered light in M_T = %7.1f %%\n",
m.fraction_of_tc_in_mt*100);
fprintf(stderr,"# \n");
fprintf(stderr,"# Reflection sphere\n");
fprintf(stderr,"#                      sphere diameter = %7.1f mm\n",
m.d_sphere_r);
fprintf(stderr,"#                 sample port diameter = %7.1f mm\n",
2*m.d_sphere_r*sqrt(m.as_r));
fprintf(stderr,"#               entrance port diameter = %7.1f mm\n",
2*m.d_sphere_r*sqrt(m.ae_r));
fprintf(stderr,"#               detector port diameter = %7.1f mm\n",
2*m.d_sphere_r*sqrt(m.ad_r));
fprintf(stderr,"#                     wall reflectance = %7.1f %%\n",m.rw_r*100);
fprintf(stderr,"#                 standard reflectance = %7.1f %%\n",m.rstd_r*100);
fprintf(stderr,"#                 detector reflectance = %7.1f %%\n",m.rd_r*100);
fprintf(stderr,"#                              spheres = %7d\n",m.num_spheres);
fprintf(stderr,"#                             measures = %7d\n",m.num_measures);
fprintf(stderr,"#                               method = %7d\n",m.method);
fprintf(stderr,"area_r as=%10.5f  ad=%10.5f    ae=%10.5f  aw=%10.5f\n",
m.as_r,m.ad_r,m.ae_r,m.aw_r);
fprintf(stderr,"refls  rd=%10.5f  rw=%10.5f  rstd=%10.5f   f=%10.5f\n",
m.rd_r,m.rw_r,m.rstd_r,m.f_r);
fprintf(stderr,"area_t as=%10.5f  ad=%10.5f    ae=%10.5f  aw=%10.5f\n",
m.as_t,m.ad_t,m.ae_t,m.aw_t);
fprintf(stderr,"refls  rd=%10.5f  rw=%10.5f  rstd=%10.5f   f=%10.5f\n",
m.rd_t,m.rw_t,m.rstd_t,m.f_t);
fprintf(stderr,"lost  ur1=%10.5f ut1=%10.5f   uru=%10.5f  utu=%10.5f\n",
m.ur1_lost,m.ut1_lost,m.utu_lost,m.utu_lost);
}/*:60*/
#line 37 "iad_util.w"



/*:1*/
