/*1:*/
#line 7 "ad_prime.w"

#include <math.h> 
#include <float.h> 
#include <stdio.h> 
#include "nr_util.h"
#include "ad_globl.h"
#include "ad_bound.h"
#include "ad_start.h"
#include "ad_doubl.h"
#include "ad_prime.h"
#include "ad_matrx.h"
#include "ad_cone.h"

/*5:*/
#line 51 "ad_prime.w"

/*4:*/
#line 47 "ad_prime.w"

void RT_Matrices(int n,struct AD_slab_type*slab,struct AD_method_type*method,
double**R,double**T)

/*:4*/
#line 52 "ad_prime.w"

{
double d;

if(n<3)
method->quad_pts= DEFAULT_QUAD_PTS;
else
if(n> MAX_QUAD_PTS)
method->quad_pts= MAX_QUAD_PTS;
else
if((n&1)==1)
method->quad_pts= n/2*2;
else
method->quad_pts= n;

Choose_Method(slab,method);

if(slab->b<=0){
Zero_Layer(n,R,T);
return;
}

n= method->quad_pts;
Init_Layer(*slab,*method,R,T);

if(slab->b==HUGE_VAL)
d= 1.0;
else
d= method->b_thinnest*slab->b/method->b_calc;

Double_Until(n,R,T,d,slab->b);
}


/*:5*/
#line 20 "ad_prime.w"

/*7:*/
#line 108 "ad_prime.w"

/*6:*/
#line 105 "ad_prime.w"

void RT(int n,struct AD_slab_type*slab,double*UR1,double*UT1,double*URU,double*UTU)

/*:6*/
#line 109 "ad_prime.w"

{
/*8:*/
#line 146 "ad_prime.w"

double**R,**T,**R2,**T2;
double*R01,*R10,*T01,*T10;
double*R23,*R32,*T23,*T32;
double**R02,**R20,**T02,**T20;
double**R03,**R30,**T03,**T30;
double**atemp,**btemp;
struct AD_method_type method;
*UR1= -1;
*URU= -1;
*UT1= -1;
*UTU= -1;

/*:8*/
#line 111 "ad_prime.w"


if(slab->cos_angle!=1.0){
RT_Cone(n,slab,OBLIQUE,UR1,UT1,URU,UTU);
return;
}

/*9:*/
#line 159 "ad_prime.w"

if(slab->n_slab<0)return;
if(slab->n_top_slide<0)return;
if(slab->n_bottom_slide<0)return;
if(slab->a<0||slab->a> 1)return;
if(slab->g<-1||slab->g> 1)return;
if(slab->b<0)return;

/*:9*/
#line 118 "ad_prime.w"


/*10:*/
#line 168 "ad_prime.w"

R= dmatrix(1,n,1,n);
T= dmatrix(1,n,1,n);
RT_Matrices(n,slab,&method,R,T);

/*:10*/
#line 120 "ad_prime.w"


if(slab->b==0){
Sp_RT(n,*slab,UR1,UT1,URU,UTU);

}else if(slab->n_slab==1&&slab->n_top_slide==1&&slab->n_bottom_slide==1
&&slab->b_top_slide==0&&slab->b_bottom_slide==0){
/*11:*/
#line 173 "ad_prime.w"

URU_and_UR1(n,slab->n_slab,R,URU,UR1);
URU_and_UR1(n,slab->n_slab,T,UTU,UT1);

/*:11*/
#line 127 "ad_prime.w"

}else if(slab->n_top_slide==slab->n_bottom_slide&&
slab->b_top_slide==0&&slab->b_bottom_slide==0){
/*12:*/
#line 177 "ad_prime.w"

R01= dvector(1,n);
R10= dvector(1,n);
T01= dvector(1,n);
T10= dvector(1,n);
Init_Boundary(*slab,method.quad_pts,R01,R10,T01,T10,TOP_BOUNDARY);

/*:12*/
#line 130 "ad_prime.w"

/*13:*/
#line 184 "ad_prime.w"

atemp= dmatrix(1,n,1,n);
btemp= dmatrix(1,n,1,n);
R2= dmatrix(1,n,1,n);
T2= dmatrix(1,n,1,n);
Add_Slides(n,R01,R10,T01,T10,R,T,R2,T2,atemp,btemp);
URU_and_UR1(n,slab->n_slab,R2,URU,UR1);
URU_and_UR1(n,slab->n_slab,T2,UTU,UT1);
free_dmatrix(atemp,1,n,1,n);
free_dmatrix(btemp,1,n,1,n);
free_dmatrix(R2,1,n,1,n);
free_dmatrix(T2,1,n,1,n);

/*:13*/
#line 131 "ad_prime.w"

/*14:*/
#line 197 "ad_prime.w"

free_dvector(R01,1,n);
free_dvector(R10,1,n);
free_dvector(T01,1,n);
free_dvector(T10,1,n);

/*:14*/
#line 132 "ad_prime.w"

}else{
/*12:*/
#line 177 "ad_prime.w"

R01= dvector(1,n);
R10= dvector(1,n);
T01= dvector(1,n);
T10= dvector(1,n);
Init_Boundary(*slab,method.quad_pts,R01,R10,T01,T10,TOP_BOUNDARY);

/*:12*/
#line 134 "ad_prime.w"

/*15:*/
#line 203 "ad_prime.w"

R23= dvector(1,n);
R32= dvector(1,n);
T23= dvector(1,n);
T32= dvector(1,n);
Init_Boundary(*slab,method.quad_pts,R23,R32,T23,T32,BOTTOM_BOUNDARY);

/*:15*/
#line 135 "ad_prime.w"

/*16:*/
#line 210 "ad_prime.w"


R02= dmatrix(1,n,1,n);
R20= dmatrix(1,n,1,n);
T02= dmatrix(1,n,1,n);
T20= dmatrix(1,n,1,n);
R03= dmatrix(1,n,1,n);
R30= dmatrix(1,n,1,n);
T03= dmatrix(1,n,1,n);
T30= dmatrix(1,n,1,n);
atemp= dmatrix(1,n,1,n);
btemp= dmatrix(1,n,1,n);

/*:16*/
#line 136 "ad_prime.w"

/*17:*/
#line 223 "ad_prime.w"

Add_Top(n,R01,R10,T01,T10,R,R,T,T,R02,R20,T02,T20,atemp,btemp);
Add_Bottom(n,R02,R20,T02,T20,R23,R32,T23,T32,R03,R30,T03,T30,atemp,btemp);
URU_and_UR1(n,slab->n_slab,R03,URU,UR1);
Transpose_Matrix(n,T03);
URU_and_UR1(n,slab->n_slab,T03,UTU,UT1);

/*:17*/
#line 137 "ad_prime.w"

/*18:*/
#line 230 "ad_prime.w"

free_dmatrix(R02,1,n,1,n);
free_dmatrix(R20,1,n,1,n);
free_dmatrix(T02,1,n,1,n);
free_dmatrix(T20,1,n,1,n);

free_dmatrix(R03,1,n,1,n);
free_dmatrix(R30,1,n,1,n);
free_dmatrix(T03,1,n,1,n);
free_dmatrix(T30,1,n,1,n);

free_dmatrix(atemp,1,n,1,n);
free_dmatrix(btemp,1,n,1,n);

/*:18*/
#line 138 "ad_prime.w"

/*19:*/
#line 244 "ad_prime.w"

free_dvector(R23,1,n);
free_dvector(R32,1,n);
free_dvector(T23,1,n);
free_dvector(T32,1,n);


/*:19*/
#line 139 "ad_prime.w"

/*14:*/
#line 197 "ad_prime.w"

free_dvector(R01,1,n);
free_dvector(R10,1,n);
free_dvector(T01,1,n);
free_dvector(T10,1,n);

/*:14*/
#line 140 "ad_prime.w"

}

/*20:*/
#line 251 "ad_prime.w"

free_dmatrix(R,1,n,1,n);
free_dmatrix(T,1,n,1,n);

/*:20*/
#line 143 "ad_prime.w"

}

/*:7*/
#line 21 "ad_prime.w"

/*22:*/
#line 271 "ad_prime.w"

/*21:*/
#line 262 "ad_prime.w"

void ez_RT(int n,double nslab,
double ntopslide,
double nbottomslide,
double a,
double b,
double g,
double*UR1,double*UT1,double*URU,double*UTU)

/*:21*/
#line 272 "ad_prime.w"

{
struct AD_slab_type slab;

slab.n_slab= nslab;
slab.n_top_slide= ntopslide;
slab.n_bottom_slide= nbottomslide;
slab.b_top_slide= 0;
slab.b_bottom_slide= 0;
slab.a= a;
slab.b= b;
slab.g= g;
slab.phase_function= HENYEY_GREENSTEIN;
slab.cos_angle= 1.0;
RT(n,&slab,UR1,UT1,URU,UTU);
}

/*:22*/
#line 22 "ad_prime.w"

/*26:*/
#line 339 "ad_prime.w"

/*25:*/
#line 336 "ad_prime.w"

void RTabs(int n,struct AD_slab_type*slab,double*UR1,double*UT1,double*URU,double*UTU)

/*:25*/
#line 340 "ad_prime.w"

{
/*27:*/
#line 364 "ad_prime.w"

double**R,**T;
double*R01,*R10,*T01,*T10;
double*R23,*R32,*T23,*T32;
double**R02,**R20,**T02,**T20;
double**R03,**R30,**T03,**T30;
double**atemp,**btemp;
struct AD_method_type method;

/*:27*/
#line 342 "ad_prime.w"

double**Rtop,**Ttop,**Rbottom,**Tbottom;
struct AD_slab_type slab1;
double btop,bbottom;

/*10:*/
#line 168 "ad_prime.w"

R= dmatrix(1,n,1,n);
T= dmatrix(1,n,1,n);
RT_Matrices(n,slab,&method,R,T);

/*:10*/
#line 347 "ad_prime.w"

/*28:*/
#line 373 "ad_prime.w"


slab1.b= slab->b_top_slide;
slab1.cos_angle= slab->cos_angle;
slab1.a= 0;
slab1.g= 0;
slab1.phase_function= HENYEY_GREENSTEIN;
slab1.n_slab= slab->n_slab;
slab1.n_top_slide= 1.0;
slab1.n_bottom_slide= 1.0;
slab1.b_top_slide= 0.0;
slab1.b_bottom_slide= 0.0;

Rtop= dmatrix(1,n,1,n);
Ttop= dmatrix(1,n,1,n);

RT_Matrices(n,&slab1,&method,Rtop,Ttop);

/*:28*/
#line 348 "ad_prime.w"

/*29:*/
#line 391 "ad_prime.w"

slab1.b= slab->b_bottom_slide;
slab1.cos_angle= slab->cos_angle;

Rbottom= dmatrix(1,n,1,n);
Tbottom= dmatrix(1,n,1,n);
RT_Matrices(n,&slab1,&method,Rbottom,Tbottom);

/*:29*/
#line 349 "ad_prime.w"

/*16:*/
#line 210 "ad_prime.w"


R02= dmatrix(1,n,1,n);
R20= dmatrix(1,n,1,n);
T02= dmatrix(1,n,1,n);
T20= dmatrix(1,n,1,n);
R03= dmatrix(1,n,1,n);
R30= dmatrix(1,n,1,n);
T03= dmatrix(1,n,1,n);
T30= dmatrix(1,n,1,n);
atemp= dmatrix(1,n,1,n);
btemp= dmatrix(1,n,1,n);

/*:16*/
#line 350 "ad_prime.w"


/*30:*/
#line 399 "ad_prime.w"

btop= slab->b_top_slide;
slab->b_top_slide= 0;
/*12:*/
#line 177 "ad_prime.w"

R01= dvector(1,n);
R10= dvector(1,n);
T01= dvector(1,n);
T10= dvector(1,n);
Init_Boundary(*slab,method.quad_pts,R01,R10,T01,T10,TOP_BOUNDARY);

/*:12*/
#line 402 "ad_prime.w"

slab->b_top_slide= btop;

/*:30*/
#line 352 "ad_prime.w"

/*31:*/
#line 405 "ad_prime.w"

bbottom= slab->b_bottom_slide;
slab->b_bottom_slide= 0;
/*15:*/
#line 203 "ad_prime.w"

R23= dvector(1,n);
R32= dvector(1,n);
T23= dvector(1,n);
T32= dvector(1,n);
Init_Boundary(*slab,method.quad_pts,R23,R32,T23,T32,BOTTOM_BOUNDARY);

/*:15*/
#line 408 "ad_prime.w"

slab->b_bottom_slide= bbottom;

/*:31*/
#line 353 "ad_prime.w"


/*32:*/
#line 411 "ad_prime.w"

Add(n,Rtop,Rtop,Ttop,Ttop,R,R,T,T,R02,R20,T02,T20);
Add(n,R02,R20,T02,T20,Rbottom,Rbottom,Tbottom,Tbottom,R03,R30,T03,T30);
Add_Top(n,R01,R10,T01,T10,R03,R30,T03,T30,R02,R20,T02,T20,atemp,btemp);
Add_Bottom(n,R02,R20,T02,T20,R23,R32,T23,T32,R03,R30,T03,T30,atemp,btemp);
URU_and_UR1(n,slab->n_slab,R03,URU,UR1);
Transpose_Matrix(n,T03);
URU_and_UR1(n,slab->n_slab,T03,UTU,UT1);

/*:32*/
#line 355 "ad_prime.w"


/*18:*/
#line 230 "ad_prime.w"

free_dmatrix(R02,1,n,1,n);
free_dmatrix(R20,1,n,1,n);
free_dmatrix(T02,1,n,1,n);
free_dmatrix(T20,1,n,1,n);

free_dmatrix(R03,1,n,1,n);
free_dmatrix(R30,1,n,1,n);
free_dmatrix(T03,1,n,1,n);
free_dmatrix(T30,1,n,1,n);

free_dmatrix(atemp,1,n,1,n);
free_dmatrix(btemp,1,n,1,n);

/*:18*/
#line 357 "ad_prime.w"

/*19:*/
#line 244 "ad_prime.w"

free_dvector(R23,1,n);
free_dvector(R32,1,n);
free_dvector(T23,1,n);
free_dvector(T32,1,n);


/*:19*/
#line 358 "ad_prime.w"

/*14:*/
#line 197 "ad_prime.w"

free_dvector(R01,1,n);
free_dvector(R10,1,n);
free_dvector(T01,1,n);
free_dvector(T10,1,n);

/*:14*/
#line 359 "ad_prime.w"

/*20:*/
#line 251 "ad_prime.w"

free_dmatrix(R,1,n,1,n);
free_dmatrix(T,1,n,1,n);

/*:20*/
#line 360 "ad_prime.w"

/*33:*/
#line 420 "ad_prime.w"

free_dmatrix(Rtop,1,n,1,n);
free_dmatrix(Ttop,1,n,1,n);
free_dmatrix(Rbottom,1,n,1,n);
free_dmatrix(Tbottom,1,n,1,n);

/*:33*/
#line 361 "ad_prime.w"

}

/*:26*/
#line 23 "ad_prime.w"

/*36:*/
#line 447 "ad_prime.w"

/*35:*/
#line 443 "ad_prime.w"

void Flux_Fluence(int n,struct AD_slab_type*slab,double zmin,double zmax,int intervals,
double*UF1_array,double*UFU_array,double*flux_up,double*flux_down)

/*:35*/
#line 448 "ad_prime.w"

{
/*37:*/
#line 469 "ad_prime.w"


double*R01,*R10,*T01,*T10;
double*R56,*R65,*T56,*T65;
double**R12,**T12;
double**R23,**T23;
double**R34,**T34;
double**R45,**T45;
double**R02,**R20,**T02,**T20;
double**R46,**R64,**T46,**T64;
double**R03,**R30,**T03,**T30;
double**R36,**R63,**T36,**T63;
double**Lup,**Ldown;
double**a,**b;
double flx_down,flx_up,UFU,UF1;
double slab_thickness;

struct AD_method_type method;
int i,j;

/*:37*/
#line 450 "ad_prime.w"


if(intervals> MAX_FLUENCE_INTERVALS)
AD_error("too many intervals requested.  increase the const max_fluence_intervals\n");

/*38:*/
#line 489 "ad_prime.w"

slab_thickness= slab->b;
slab->b= zmin;
R12= dmatrix(1,n,1,n);
T12= dmatrix(1,n,1,n);
RT_Matrices(n,slab,&method,R12,T12);

R01= dvector(1,n);
R10= dvector(1,n);
T01= dvector(1,n);
T10= dvector(1,n);
Init_Boundary(*slab,method.quad_pts,R01,R10,T01,T10,TOP_BOUNDARY);

R20= dmatrix(1,n,1,n);
T20= dmatrix(1,n,1,n);
R02= dmatrix(1,n,1,n);
T02= dmatrix(1,n,1,n);
a= dmatrix(1,n,1,n);
b= dmatrix(1,n,1,n);
Add_Top(n,R01,R10,T01,T10,R12,R12,T12,T12,R02,R20,T02,T20,a,b);

free_dmatrix(R12,1,n,1,n);
free_dmatrix(T12,1,n,1,n);
free_dvector(R01,1,n);
free_dvector(R10,1,n);
free_dvector(T01,1,n);
free_dvector(T10,1,n);

/*:38*/
#line 455 "ad_prime.w"

/*39:*/
#line 517 "ad_prime.w"

slab->b= slab_thickness-zmax;
R45= dmatrix(1,n,1,n);
T45= dmatrix(1,n,1,n);
RT_Matrices(n,slab,&method,R45,T45);
R56= dvector(1,n);
R65= dvector(1,n);
T56= dvector(1,n);
T65= dvector(1,n);
Init_Boundary(*slab,method.quad_pts,R56,R65,T56,T65,BOTTOM_BOUNDARY);
R46= dmatrix(1,n,1,n);
T46= dmatrix(1,n,1,n);
R64= dmatrix(1,n,1,n);
T64= dmatrix(1,n,1,n);
Add_Bottom(n,R45,R45,T45,T45,R56,R65,T56,T65,R46,R64,T46,T64,a,b);
free_dmatrix(R45,1,n,1,n);
free_dmatrix(T45,1,n,1,n);
free_dvector(R56,1,n);
free_dvector(R65,1,n);
free_dvector(T56,1,n);
free_dvector(T65,1,n);
free_dmatrix(a,1,n,1,n);
free_dmatrix(b,1,n,1,n);

/*:39*/
#line 456 "ad_prime.w"

/*40:*/
#line 541 "ad_prime.w"

R23= dmatrix(1,n,1,n);
T23= dmatrix(1,n,1,n);
R03= dmatrix(1,n,1,n);
T03= dmatrix(1,n,1,n);
R30= dmatrix(1,n,1,n);
T30= dmatrix(1,n,1,n);

R34= dmatrix(1,n,1,n);
T34= dmatrix(1,n,1,n);
R63= dmatrix(1,n,1,n);
T63= dmatrix(1,n,1,n);
R36= dmatrix(1,n,1,n);
T36= dmatrix(1,n,1,n);

Lup= dmatrix(1,n,1,n);
Ldown= dmatrix(1,n,1,n);

/*:40*/
#line 457 "ad_prime.w"


for(i= 0;i<=intervals;i++){

/*41:*/
#line 559 "ad_prime.w"

slab->b= (zmax-zmin)/intervals*i;
RT_Matrices(n,slab,&method,R23,T23);
Add(n,R02,R20,T02,T20,R23,R23,T23,T23,R03,R30,T03,T30);

slab->b= (zmax-zmin)-slab->b;
RT_Matrices(n,slab,&method,R34,T34);
Add(n,R34,R34,T34,T34,R46,R64,T46,T64,R36,R63,T36,T63);

Between(n,R03,R30,T03,T30,R36,R63,T36,T63,Lup,Ldown);

/*:41*/
#line 461 "ad_prime.w"


/*42:*/
#line 570 "ad_prime.w"

UFU_and_UF1(n,slab->n_slab,Lup,Ldown,&UFU,&UF1);
UF1_array[i]= UF1;
UFU_array[i]= UFU;

flx_down= 0.0;
flx_up= 0.0;
for(j= 1;j<=n;j++){
flx_down+= twoaw[j]*Ldown[j][n];
flx_up+= twoaw[j]*Lup[j][n];
}
flux_down[i]= flx_down*slab->n_slab*slab->n_slab;
flux_up[i]= flx_up*slab->n_slab*slab->n_slab;

/*:42*/
#line 463 "ad_prime.w"

}

/*43:*/
#line 584 "ad_prime.w"

free_dmatrix(R02,1,n,1,n);
free_dmatrix(T02,1,n,1,n);
free_dmatrix(R20,1,n,1,n);
free_dmatrix(T20,1,n,1,n);

free_dmatrix(R23,1,n,1,n);
free_dmatrix(T23,1,n,1,n);

free_dmatrix(R03,1,n,1,n);
free_dmatrix(T03,1,n,1,n);
free_dmatrix(R30,1,n,1,n);
free_dmatrix(T30,1,n,1,n);

free_dmatrix(R34,1,n,1,n);
free_dmatrix(T34,1,n,1,n);

free_dmatrix(R63,1,n,1,n);
free_dmatrix(T63,1,n,1,n);
free_dmatrix(R36,1,n,1,n);
free_dmatrix(T36,1,n,1,n);

free_dmatrix(R64,1,n,1,n);
free_dmatrix(T64,1,n,1,n);
free_dmatrix(R46,1,n,1,n);
free_dmatrix(T46,1,n,1,n);

free_dmatrix(Lup,1,n,1,n);
free_dmatrix(Ldown,1,n,1,n);
/*:43*/
#line 466 "ad_prime.w"

}

/*:36*/
#line 24 "ad_prime.w"

/*24:*/
#line 306 "ad_prime.w"

/*23:*/
#line 296 "ad_prime.w"

void ez_RT_unscattered(int n,
double nslab,
double ntopslide,
double nbottomslide,
double a,
double b,
double g,
double*UR1,double*UT1,double*URU,double*UTU)

/*:23*/
#line 307 "ad_prime.w"

{
struct AD_slab_type slab;

slab.n_slab= nslab;
slab.n_top_slide= ntopslide;
slab.n_bottom_slide= nbottomslide;
slab.b_top_slide= 0;
slab.b_bottom_slide= 0;
slab.a= a;
slab.b= b;
slab.g= g;
slab.phase_function= HENYEY_GREENSTEIN;
slab.cos_angle= 1.0;
Sp_RT(n,slab,UR1,UT1,URU,UTU);
}

/*:24*/
#line 25 "ad_prime.w"


/*:1*/
