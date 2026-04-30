/*1:*/
#line 3 "iad_find.w"

#include <math.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include "ad_globl.h"
#include "nr_util.h"
#include "nr_mnbrk.h"
#include "nr_brent.h"
#include "nr_amoeb.h"
#include "iad_type.h"
#include "iad_util.h"
#include "iad_calc.h"
#include "iad_agrid.h"
#include "iad_find.h"

#define NUMBER_OF_GUESSES 10

guess_type guess[NUMBER_OF_GUESSES];

static double scipy_bounded_initial_vertex(double value,double lower,double upper)
{
if(value> upper)
value= 2.0*upper-value;
if(value<lower)
value= lower;
if(value> upper)
value= upper;
return value;
}

int compare_guesses(const void*p1,const void*p2)
{
guess_type*g1= (guess_type*)p1;
guess_type*g2= (guess_type*)p2;

if(g1->distance<g2->distance)
return-1;
else if(g1->distance==g2->distance)
return 0;
else
return 1;
}

/*17:*/
#line 366 "iad_find.w"

/*16:*/
#line 363 "iad_find.w"

void U_Find_Ba(struct measure_type m,struct invert_type*r)

/*:16*/
#line 367 "iad_find.w"

{
double ax,bx,cx,fa,fb,fc,ba;

if(Debug(DEBUG_SEARCH)){
fprintf(stderr,"SEARCH: Using U_Find_Bs()");
fprintf(stderr," (mu=%6.4f)",r->slab.cos_angle);
if(r->default_bs!=UNINITIALIZED)
fprintf(stderr,"  default_bs = %8.5f",r->default_bs);
if(r->default_g!=UNINITIALIZED)
fprintf(stderr,"  default_g = %8.5f",r->default_g);
fprintf(stderr,"\n");
}

r->slab.a= 0;
r->slab.g= (r->default_g==UNINITIALIZED)?0:r->default_g;
r->slab.b= (r->default_bs==UNINITIALIZED)?HUGE_VAL:r->default_bs;

if(m.m_t==0){
r->slab.b= HUGE_VAL;
U_Find_A(m,r);
return;
}

Set_Calc_State(m,*r);

ax= b2bcalc(0.1);
bx= b2bcalc(1.0);
mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,Find_Ba_fn,&r->AD_iterations);
r->final_distance= brent(ax,bx,cx,Find_Ba_fn,r->tolerance,&ba,&r->AD_iterations);


r->slab.a= (r->slab.b)/(bcalc2b(ba)+r->slab.b);
r->slab.b= bcalc2b(ba)+r->slab.b;

/*12:*/
#line 280 "iad_find.w"

r->a= r->slab.a;
r->b= r->slab.b;
r->g= r->slab.g;
r->found= (r->tolerance> r->final_distance);
Set_Calc_State(m,*r);

/*:12*/
#line 402 "iad_find.w"

}

/*:17*/
#line 46 "iad_find.w"

/*15:*/
#line 315 "iad_find.w"

/*14:*/
#line 312 "iad_find.w"

void U_Find_Bs(struct measure_type m,struct invert_type*r)

/*:14*/
#line 316 "iad_find.w"

{
double ax,bx,cx,fa,fb,fc,bs;

if(Debug(DEBUG_SEARCH)){
fprintf(stderr,"SEARCH: Using U_Find_Bs()");
fprintf(stderr," (mu=%6.4f)",r->slab.cos_angle);
if(r->default_ba!=UNINITIALIZED)
fprintf(stderr,"  default_ba = %8.5f",r->default_ba);
if(r->default_g!=UNINITIALIZED)
fprintf(stderr,"  default_g = %8.5f",r->default_g);
fprintf(stderr,"\n");
}

if(m.m_t==0){
r->slab.b= HUGE_VAL;
U_Find_A(m,r);
return;
}

r->slab.a= 0;
r->slab.g= (r->default_g==UNINITIALIZED)?0:r->default_g;
r->slab.b= (r->default_ba==UNINITIALIZED)?HUGE_VAL:r->default_ba;

Set_Calc_State(m,*r);

ax= b2bcalc(0.1);
bx= b2bcalc(1.0);
mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,Find_Bs_fn,&r->AD_iterations);
r->final_distance= brent(ax,bx,cx,Find_Bs_fn,r->tolerance,&bs,&r->AD_iterations);


r->slab.a= bcalc2b(bs)/(bcalc2b(bs)+r->slab.b);
r->slab.b= bcalc2b(bs)+r->slab.b;

/*12:*/
#line 280 "iad_find.w"

r->a= r->slab.a;
r->b= r->slab.b;
r->g= r->slab.g;
r->found= (r->tolerance> r->final_distance);
Set_Calc_State(m,*r);

/*:12*/
#line 351 "iad_find.w"

}

/*:15*/
#line 47 "iad_find.w"

/*19:*/
#line 417 "iad_find.w"

/*18:*/
#line 414 "iad_find.w"

void U_Find_A(struct measure_type m,struct invert_type*r)

/*:18*/
#line 418 "iad_find.w"
{
double Rt,Tt,Rd,Ru,Td,Tu;

if(Debug(DEBUG_SEARCH)){
fprintf(stderr,"SEARCH: Using U_Find_A()");
fprintf(stderr," (mu=%6.4f)",r->slab.cos_angle);
if(r->default_b!=UNINITIALIZED)
fprintf(stderr,"  default_b = %8.5f",r->default_b);
if(r->default_g!=UNINITIALIZED)
fprintf(stderr,"  default_g = %8.5f",r->default_g);
fprintf(stderr,"\n");
}

Estimate_RT(m,*r,&Rt,&Tt,&Rd,&Ru,&Td,&Tu);

r->slab.g= (r->default_g==UNINITIALIZED)?0:r->default_g;
r->slab.b= (r->default_b==UNINITIALIZED)?HUGE_VAL:r->default_b;
r->slab.a= 0.0;
r->final_distance= 0.0;
Set_Calc_State(m,*r);

if(Rt> 0.99999){
r->final_distance= Find_A_fn(a2acalc(1.0));
r->slab.a= 1.0;
}else{
double x,ax,bx,cx,fa,fb,fc;

ax= a2acalc(0.3);
bx= a2acalc(0.5);

mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,Find_A_fn,&r->AD_iterations);
r->final_distance= brent(ax,bx,cx,Find_A_fn,r->tolerance,&x,&r->AD_iterations);
r->slab.a= acalc2a(x);
}

/*12:*/
#line 280 "iad_find.w"

r->a= r->slab.a;
r->b= r->slab.b;
r->g= r->slab.g;
r->found= (r->tolerance> r->final_distance);
Set_Calc_State(m,*r);

/*:12*/
#line 453 "iad_find.w"

}

/*:19*/
#line 48 "iad_find.w"

/*23:*/
#line 513 "iad_find.w"

/*22:*/
#line 510 "iad_find.w"

void U_Find_B(struct measure_type m,struct invert_type*r)

/*:22*/
#line 514 "iad_find.w"
{
double Rt,Tt,Rd,Ru,Td,Tu;

if(Debug(DEBUG_SEARCH)){
fprintf(stderr,"SEARCH: Using U_Find_B()");
fprintf(stderr," (mu=%6.4f)",r->slab.cos_angle);
if(r->default_a!=UNINITIALIZED)
fprintf(stderr,"  default_a = %8.5f",r->default_a);
if(r->default_g!=UNINITIALIZED)
fprintf(stderr,"  default_g = %8.5f",r->default_g);
fprintf(stderr,"\n");
}

Estimate_RT(m,*r,&Rt,&Tt,&Rd,&Ru,&Td,&Tu);

r->slab.g= (r->default_g==UNINITIALIZED)?0:r->default_g;
r->slab.a= (r->default_a==UNINITIALIZED)?0:r->default_a;
r->slab.b= 0.5;
r->final_distance= 0.0;
Set_Calc_State(m,*r);

/*24:*/
#line 542 "iad_find.w"

{
double x,ax,bx,cx,fa,fb,fc;

ax= b2bcalc(0.1);
bx= b2bcalc(10);

mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,Find_B_fn,&r->AD_iterations);
r->final_distance= brent(ax,bx,cx,Find_B_fn,r->tolerance,&x,&r->AD_iterations);
r->slab.b= bcalc2b(x);
Set_Calc_State(m,*r);
}


/*:24*/
#line 535 "iad_find.w"


/*12:*/
#line 280 "iad_find.w"

r->a= r->slab.a;
r->b= r->slab.b;
r->g= r->slab.g;
r->found= (r->tolerance> r->final_distance);
Set_Calc_State(m,*r);

/*:12*/
#line 537 "iad_find.w"

}

/*:23*/
#line 49 "iad_find.w"

/*21:*/
#line 461 "iad_find.w"

/*20:*/
#line 458 "iad_find.w"

void U_Find_G(struct measure_type m,struct invert_type*r)

/*:20*/
#line 462 "iad_find.w"
{
double Rt,Tt,Rd,Ru,Td,Tu;
double x,ax,bx,cx,fa,fb,fc;

if(Debug(DEBUG_SEARCH)){
fprintf(stderr,"SEARCH: Using U_Find_G()");
fprintf(stderr," (mu=%6.4f)",r->slab.cos_angle);
if(r->default_a!=UNINITIALIZED)
fprintf(stderr,"  default_a = %8.5f",r->default_a);
if(r->default_b!=UNINITIALIZED)
fprintf(stderr,"  default_b = %8.5f",r->default_b);
fprintf(stderr,"\n");
}

Estimate_RT(m,*r,&Rt,&Tt,&Rd,&Ru,&Td,&Tu);

r->slab.a= (r->default_a==UNINITIALIZED)?0.5:r->default_a;
if(r->default_b!=UNINITIALIZED)
r->slab.b= r->default_b;
else if(m.m_u> 0)
r->slab.b= What_Is_B(r->slab,m.m_u);
else
r->slab.b= HUGE_VAL;

r->slab.g= 0.0;
r->final_distance= 0.0;
Set_Calc_State(m,*r);

ax= g2gcalc(-0.99);
bx= g2gcalc(0.99);

mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,Find_G_fn,&r->AD_iterations);
r->final_distance= brent(ax,bx,cx,Find_G_fn,r->tolerance,&x,&r->AD_iterations);

r->slab.g= gcalc2g(x);
Set_Calc_State(m,*r);


/*12:*/
#line 280 "iad_find.w"

r->a= r->slab.a;
r->b= r->slab.b;
r->g= r->slab.g;
r->found= (r->tolerance> r->final_distance);
Set_Calc_State(m,*r);

/*:12*/
#line 500 "iad_find.w"

}


/*:21*/
#line 50 "iad_find.w"

/*26:*/
#line 572 "iad_find.w"

/*25:*/
#line 569 "iad_find.w"

void U_Find_AG(struct measure_type m,struct invert_type*r)

/*:25*/
#line 573 "iad_find.w"

{
/*5:*/
#line 116 "iad_find.w"

int i,i_best,j_best;
double*x,*y,**p;

x= dvector(1,2);
y= dvector(1,3);
p= dmatrix(1,3,1,2);

/*:5*/
#line 575 "iad_find.w"


if(Debug(DEBUG_SEARCH)){
fprintf(stderr,"SEARCH: Using U_Find_AG()");
fprintf(stderr," mu=%4.2f, ",r->slab.cos_angle);
if(m.num_measures==3)
fprintf(stderr," b= %6.3f  (M_U)",What_Is_B(r->slab,m.m_u));
else if(r->default_b!=UNINITIALIZED)
fprintf(stderr," b = %6.3f (constrained)",r->default_b);
else
fprintf(stderr," b = %6.3f (default)",1.0);
fprintf(stderr,"\n");
}

if(m.num_measures==3)
r->slab.b= What_Is_B(r->slab,m.m_u);
else if(r->default_b==UNINITIALIZED)
r->slab.b= 1;
else
r->slab.b= r->default_b;

Set_Calc_State(m,*r);
/*46:*/
#line 939 "iad_find.w"

{
int n_agrid;
int k;
size_t count;


abg_distance(r->slab.a,r->slab.b,r->slab.g,&(guess[0]));

if(r->MC_iterations> 0){
for(k= 1;k<NUMBER_OF_GUESSES;k++)
guess[k]= guess[0];
count= 1;
}else{

if(!AGrid_Valid(m,*r))AGrid_Build(m,*r);


n_agrid= AGrid_Fill_Guesses(m.m_r,m.m_t,
guess+1,NUMBER_OF_GUESSES-1);
count= (size_t)(1+n_agrid);
qsort((void*)guess,count,sizeof(guess_type),compare_guesses);
}

if(Debug(DEBUG_BEST_GUESS)){
fprintf(stderr,"BEST: AGRID + PREV GUESSES\n");
fprintf(stderr,"BEST:  k      albedo          b          g   distance\n");
for(k= 0;k<(int)count&&k<7;k++){
fprintf(stderr,"BEST:%3d  ",k);
fprintf(stderr,"%10.5f ",guess[k].a);
fprintf(stderr,"%10.5f ",guess[k].b);
fprintf(stderr,"%10.5f ",guess[k].g);
fprintf(stderr,"%10.5f\n",guess[k].distance);
}
}
}/*:46*/
#line 597 "iad_find.w"

/*27:*/
#line 613 "iad_find.w"

{
int k,kk;
double a0= scipy_bounded_initial_vertex(guess[0].a,0.0,1.0);
double g0= scipy_bounded_initial_vertex(guess[0].g,-MAX_ABS_G,MAX_ABS_G);


p[1][1]= scipy_bounded_initial_vertex(a0,0.0,1.0);
p[1][2]= scipy_bounded_initial_vertex(g0,-MAX_ABS_G,MAX_ABS_G);



{
double a_step= (r->mc_simplex_a_step> 0.0)?r->mc_simplex_a_step
:(a0!=0.0?0.05*a0:0.00025);
double g_step= (r->mc_simplex_g_step> 0.0)?r->mc_simplex_g_step
:(g0!=0.0?0.05*g0:0.00025);
p[2][1]= scipy_bounded_initial_vertex(a0+a_step,0.0,1.0);
p[2][2]= scipy_bounded_initial_vertex(g0,-MAX_ABS_G,MAX_ABS_G);
p[3][1]= scipy_bounded_initial_vertex(a0,0.0,1.0);
p[3][2]= scipy_bounded_initial_vertex(g0+g_step,-MAX_ABS_G,MAX_ABS_G);
}


for(k= 1;k<7;k++){
if(guess[0].a!=guess[k].a)
break;
}
for(kk= 1;kk<7;kk++){
if(kk==k)
continue;
if(guess[0].g!=guess[kk].g||guess[k].g!=guess[kk].g)
break;
}

if(Debug(DEBUG_BEST_GUESS)){
fprintf(stderr,"-----------------------------------------------------\n");
fprintf(stderr,"BEST: <1> ");
fprintf(stderr,"%10.5f ",guess[0].a);
fprintf(stderr,"%10.5f ",guess[0].b);
fprintf(stderr,"%10.5f ",guess[0].g);
fprintf(stderr,"%10.5f\n",guess[0].distance);
fprintf(stderr,"BEST: <2> (physical simplex, a+step)\n");
fprintf(stderr,"BEST: <3> (physical simplex, g+step)\n");
fprintf(stderr,"\n");
}
}

/*:27*/
#line 598 "iad_find.w"

/*28:*/
#line 661 "iad_find.w"


for(i= 1;i<=3;i++){
x[1]= p[i][1];
x[2]= p[i][2];
y[i]= Find_AG_fn(x);
}


/*:28*/
#line 599 "iad_find.w"

{
double lo[3],hi[3];
lo[1]= 0.0;hi[1]= 1.0;
lo[2]= -MAX_ABS_G;hi[2]= MAX_ABS_G;
amoeba(p,y,2,r->tolerance,Find_AG_fn,&r->AD_iterations,lo,hi);
}
/*29:*/
#line 674 "iad_find.w"

r->final_distance= 10;
for(i= 1;i<=3;i++){
if(y[i]<r->final_distance){
r->slab.a= p[i][1];
r->slab.g= p[i][2];
r->final_distance= y[i];
}
}


/*:29*/
#line 606 "iad_find.w"

/*11:*/
#line 270 "iad_find.w"

if(r->slab.a> 1.0-1e-3){
guess_type boundary_guess;
abg_distance(1.0,r->slab.b,r->slab.g,&boundary_guess);
if(boundary_guess.distance<r->final_distance){
r->slab.a= 1.0;
r->final_distance= boundary_guess.distance;
}
}

/*:11*/
#line 607 "iad_find.w"

/*13:*/
#line 288 "iad_find.w"

free_dvector(x,1,2);
free_dvector(y,1,3);
free_dmatrix(p,1,3,1,2);

/*:13*/
#line 608 "iad_find.w"


/*12:*/
#line 280 "iad_find.w"

r->a= r->slab.a;
r->b= r->slab.b;
r->g= r->slab.g;
r->found= (r->tolerance> r->final_distance);
Set_Calc_State(m,*r);

/*:12*/
#line 610 "iad_find.w"

}

/*:26*/
#line 51 "iad_find.w"

/*4:*/
#line 81 "iad_find.w"

/*3:*/
#line 78 "iad_find.w"

void U_Find_AB(struct measure_type m,struct invert_type*r)

/*:3*/
#line 82 "iad_find.w"

{
/*5:*/
#line 116 "iad_find.w"

int i,i_best,j_best;
double*x,*y,**p;

x= dvector(1,2);
y= dvector(1,3);
p= dmatrix(1,3,1,2);

/*:5*/
#line 84 "iad_find.w"


if(Debug(DEBUG_SEARCH)){
fprintf(stderr,"SEARCH: Using U_Find_AB()");
fprintf(stderr," mu=%4.2f, g=",r->slab.cos_angle);
if(r->default_g!=UNINITIALIZED)
fprintf(stderr,"  %7.3f (constrained g)",r->default_g);
else
fprintf(stderr,"  %7.3f (default)",0.0);
fprintf(stderr,"\n");
}

r->slab.g= (r->default_g==UNINITIALIZED)?0:r->default_g;
Set_Calc_State(m,*r);

/*46:*/
#line 939 "iad_find.w"

{
int n_agrid;
int k;
size_t count;


abg_distance(r->slab.a,r->slab.b,r->slab.g,&(guess[0]));

if(r->MC_iterations> 0){
for(k= 1;k<NUMBER_OF_GUESSES;k++)
guess[k]= guess[0];
count= 1;
}else{

if(!AGrid_Valid(m,*r))AGrid_Build(m,*r);


n_agrid= AGrid_Fill_Guesses(m.m_r,m.m_t,
guess+1,NUMBER_OF_GUESSES-1);
count= (size_t)(1+n_agrid);
qsort((void*)guess,count,sizeof(guess_type),compare_guesses);
}

if(Debug(DEBUG_BEST_GUESS)){
fprintf(stderr,"BEST: AGRID + PREV GUESSES\n");
fprintf(stderr,"BEST:  k      albedo          b          g   distance\n");
for(k= 0;k<(int)count&&k<7;k++){
fprintf(stderr,"BEST:%3d  ",k);
fprintf(stderr,"%10.5f ",guess[k].a);
fprintf(stderr,"%10.5f ",guess[k].b);
fprintf(stderr,"%10.5f ",guess[k].g);
fprintf(stderr,"%10.5f\n",guess[k].distance);
}
}
}/*:46*/
#line 99 "iad_find.w"

/*7:*/
#line 170 "iad_find.w"

{
int k,kk;
double a0= scipy_bounded_initial_vertex(guess[0].a,0.0,1.0);
double b0= scipy_bounded_initial_vertex((guess[0].b<1e8)?guess[0].b:1.0,0.0,1e10);


p[1][1]= scipy_bounded_initial_vertex(a0,0.0,1.0);
p[1][2]= scipy_bounded_initial_vertex(b0,0.0,1e10);



{
double a_step= (r->mc_simplex_a_step> 0.0)?r->mc_simplex_a_step
:(a0!=0.0?0.05*a0:0.00025);
double b_step= (r->mc_simplex_b_step> 0.0)?r->mc_simplex_b_step
:(b0!=0.0?0.05*b0:0.00025);
p[2][1]= scipy_bounded_initial_vertex(a0+a_step,0.0,1.0);
p[2][2]= scipy_bounded_initial_vertex(b0,0.0,1e10);
p[3][1]= scipy_bounded_initial_vertex(a0,0.0,1.0);
p[3][2]= scipy_bounded_initial_vertex(b0+b_step,0.0,1e10);
}


for(k= 1;k<7;k++){
if(guess[0].a!=guess[k].a)
break;
}
for(kk= 1;kk<7;kk++){
if(k==kk)
continue;
if(guess[0].b!=guess[kk].b||guess[k].b!=guess[kk].b)
break;
}

if(Debug(DEBUG_BEST_GUESS)){
fprintf(stderr,"-----------------------------------------------------\n");
fprintf(stderr,"BEST: <1> ");
fprintf(stderr,"%10.5f ",guess[0].a);
fprintf(stderr,"%10.5f ",guess[0].b);
fprintf(stderr,"%10.5f ",guess[0].g);
fprintf(stderr,"%10.5f\n",guess[0].distance);
fprintf(stderr,"BEST: <2> (physical simplex, a+step)\n");
fprintf(stderr,"BEST: <3> (physical simplex, b+step)\n");
fprintf(stderr,"\n");
}
if(Debug(DEBUG_MC)){
m.ur1_lost= guess[kk].ur1_lost;
m.ut1_lost= guess[kk].ut1_lost;
m.uru_lost= guess[kk].uru_lost;
m.utu_lost= guess[kk].utu_lost;
}
}

/*:7*/
#line 100 "iad_find.w"

/*8:*/
#line 224 "iad_find.w"


for(i= 1;i<=3;i++){
x[1]= p[i][1];
x[2]= p[i][2];
y[i]= Find_AB_fn(x);
}

/*:8*/
#line 101 "iad_find.w"

{
double lo[3],hi[3];
lo[1]= 0.0;hi[1]= 1.0;
lo[2]= 0.0;hi[2]= 1e10;
amoeba(p,y,2,r->tolerance,Find_AB_fn,&r->AD_iterations,lo,hi);
}
/*9:*/
#line 232 "iad_find.w"

r->final_distance= 10;
for(i= 1;i<=3;i++){
if(y[i]<r->final_distance){
r->slab.a= p[i][1];
r->slab.b= p[i][2];
r->final_distance= y[i];
}
}

/*:9*/
#line 108 "iad_find.w"

/*10:*/
#line 247 "iad_find.w"

if(r->slab.a> 1.0-1e-3){
double saved_fd= r->final_distance;
double saved_a= r->slab.a;
double saved_b= r->slab.b;
int saved_iter= r->AD_iterations;
double saved_default_a= r->default_a;
r->default_a= 1.0;
U_Find_B(m,r);
if(r->final_distance>=saved_fd){

r->slab.a= saved_a;
r->slab.b= saved_b;
r->final_distance= saved_fd;
}
r->default_a= saved_default_a;
r->AD_iterations+= saved_iter;
}

/*:10*/
#line 109 "iad_find.w"


/*13:*/
#line 288 "iad_find.w"

free_dvector(x,1,2);
free_dvector(y,1,3);
free_dmatrix(p,1,3,1,2);

/*:13*/
#line 111 "iad_find.w"

/*12:*/
#line 280 "iad_find.w"

r->a= r->slab.a;
r->b= r->slab.b;
r->g= r->slab.g;
r->found= (r->tolerance> r->final_distance);
Set_Calc_State(m,*r);

/*:12*/
#line 112 "iad_find.w"

}

/*:4*/
#line 52 "iad_find.w"

/*31:*/
#line 692 "iad_find.w"

/*30:*/
#line 689 "iad_find.w"

void U_Find_BG(struct measure_type m,struct invert_type*r)

/*:30*/
#line 693 "iad_find.w"

{
/*5:*/
#line 116 "iad_find.w"

int i,i_best,j_best;
double*x,*y,**p;

x= dvector(1,2);
y= dvector(1,3);
p= dmatrix(1,3,1,2);

/*:5*/
#line 695 "iad_find.w"


if(Debug(DEBUG_SEARCH)){
fprintf(stderr,"SEARCH: Using U_Find_BG()");
fprintf(stderr," (mu=%6.4f)",r->slab.cos_angle);
if(r->default_a!=UNINITIALIZED)
fprintf(stderr,"  default_a = %8.5f",r->default_a);
fprintf(stderr,"\n");
}

r->slab.a= (r->default_a==UNINITIALIZED)?0:r->default_a;
Set_Calc_State(m,*r);

/*46:*/
#line 939 "iad_find.w"

{
int n_agrid;
int k;
size_t count;


abg_distance(r->slab.a,r->slab.b,r->slab.g,&(guess[0]));

if(r->MC_iterations> 0){
for(k= 1;k<NUMBER_OF_GUESSES;k++)
guess[k]= guess[0];
count= 1;
}else{

if(!AGrid_Valid(m,*r))AGrid_Build(m,*r);


n_agrid= AGrid_Fill_Guesses(m.m_r,m.m_t,
guess+1,NUMBER_OF_GUESSES-1);
count= (size_t)(1+n_agrid);
qsort((void*)guess,count,sizeof(guess_type),compare_guesses);
}

if(Debug(DEBUG_BEST_GUESS)){
fprintf(stderr,"BEST: AGRID + PREV GUESSES\n");
fprintf(stderr,"BEST:  k      albedo          b          g   distance\n");
for(k= 0;k<(int)count&&k<7;k++){
fprintf(stderr,"BEST:%3d  ",k);
fprintf(stderr,"%10.5f ",guess[k].a);
fprintf(stderr,"%10.5f ",guess[k].b);
fprintf(stderr,"%10.5f ",guess[k].g);
fprintf(stderr,"%10.5f\n",guess[k].distance);
}
}
}/*:46*/
#line 708 "iad_find.w"

/*33:*/
#line 727 "iad_find.w"

{
int k,kk;
double b0= scipy_bounded_initial_vertex((guess[0].b<1e8)?guess[0].b:1.0,0.0,1e10);
double g0= scipy_bounded_initial_vertex(guess[0].g,-MAX_ABS_G,MAX_ABS_G);


p[1][1]= scipy_bounded_initial_vertex(b0,0.0,1e10);
p[1][2]= scipy_bounded_initial_vertex(g0,-MAX_ABS_G,MAX_ABS_G);



{
double b_step= (r->mc_simplex_b_step> 0.0)?r->mc_simplex_b_step
:(b0!=0.0?0.05*b0:0.00025);
double g_step= (r->mc_simplex_g_step> 0.0)?r->mc_simplex_g_step
:(g0!=0.0?0.05*g0:0.00025);
p[2][1]= scipy_bounded_initial_vertex(b0+b_step,0.0,1e10);
p[2][2]= scipy_bounded_initial_vertex(g0,-MAX_ABS_G,MAX_ABS_G);
p[3][1]= scipy_bounded_initial_vertex(b0,0.0,1e10);
p[3][2]= scipy_bounded_initial_vertex(g0+g_step,-MAX_ABS_G,MAX_ABS_G);
}


for(k= 1;k<7;k++){
if(guess[0].b!=guess[k].b)
break;
}
for(kk= 1;kk<7;kk++){
if(kk==k)
continue;
if(guess[0].g!=guess[kk].g||guess[k].g!=guess[kk].g)
break;
}

if(Debug(DEBUG_BEST_GUESS)){
fprintf(stderr,"-----------------------------------------------------\n");
fprintf(stderr,"BEST: <1> ");
fprintf(stderr,"%10.5f ",guess[0].a);
fprintf(stderr,"%10.5f ",guess[0].b);
fprintf(stderr,"%10.5f ",guess[0].g);
fprintf(stderr,"%10.5f\n",guess[0].distance);
fprintf(stderr,"BEST: <2> (physical simplex, b+step)\n");
fprintf(stderr,"BEST: <3> (physical simplex, g+step)\n");
fprintf(stderr,"\n");
}
}

/*:33*/
#line 709 "iad_find.w"

/*34:*/
#line 775 "iad_find.w"


for(i= 1;i<=3;i++){
x[1]= p[i][1];
x[2]= p[i][2];
y[i]= Find_BG_fn(x);
}

/*:34*/
#line 710 "iad_find.w"

{
double lo[3],hi[3];
lo[1]= 0.0;hi[1]= 1e10;
lo[2]= -MAX_ABS_G;hi[2]= MAX_ABS_G;
amoeba(p,y,2,r->tolerance,Find_BG_fn,&r->AD_iterations,lo,hi);
}
/*35:*/
#line 787 "iad_find.w"

r->final_distance= 10;
for(i= 1;i<=3;i++){
if(y[i]<r->final_distance){
r->slab.b= p[i][1];
r->slab.g= p[i][2];
r->final_distance= y[i];
}
}

/*:35*/
#line 717 "iad_find.w"


/*13:*/
#line 288 "iad_find.w"

free_dvector(x,1,2);
free_dvector(y,1,3);
free_dmatrix(p,1,3,1,2);

/*:13*/
#line 719 "iad_find.w"

/*12:*/
#line 280 "iad_find.w"

r->a= r->slab.a;
r->b= r->slab.b;
r->g= r->slab.g;
r->found= (r->tolerance> r->final_distance);
Set_Calc_State(m,*r);

/*:12*/
#line 720 "iad_find.w"

}

/*:31*/
#line 53 "iad_find.w"

/*37:*/
#line 810 "iad_find.w"

/*36:*/
#line 807 "iad_find.w"

void U_Find_BaG(struct measure_type m,struct invert_type*r)

/*:36*/
#line 811 "iad_find.w"

{
/*5:*/
#line 116 "iad_find.w"

int i,i_best,j_best;
double*x,*y,**p;

x= dvector(1,2);
y= dvector(1,3);
p= dmatrix(1,3,1,2);

/*:5*/
#line 813 "iad_find.w"

Set_Calc_State(m,*r);
/*6:*/
#line 133 "iad_find.w"

{

size_t count= NUMBER_OF_GUESSES;

abg_distance(r->slab.a,r->slab.b,r->slab.g,&(guess[0]));

if(!Valid_Grid(m,*r))Fill_Grid(m,*r,1);


Near_Grid_Points(m.m_r,m.m_t,r->search,&i_best,&j_best);
Grid_ABG(i_best,j_best,&(guess[1]));
Grid_ABG(i_best+1,j_best,&(guess[2]));
Grid_ABG(i_best-1,j_best,&(guess[3]));
Grid_ABG(i_best,j_best+1,&(guess[4]));
Grid_ABG(i_best,j_best-1,&(guess[5]));
Grid_ABG(i_best+1,j_best+1,&(guess[6]));
Grid_ABG(i_best-1,j_best-1,&(guess[7]));
Grid_ABG(i_best+1,j_best-1,&(guess[8]));
Grid_ABG(i_best-1,j_best+1,&(guess[9]));

qsort((void*)guess,count,sizeof(guess_type),compare_guesses);

if(Debug(DEBUG_BEST_GUESS)){
int k;
fprintf(stderr,"BEST: GRID GUESSES\n");
fprintf(stderr,"BEST:  k      albedo          b          g   distance\n");
for(k= 0;k<=6;k++){
fprintf(stderr,"BEST:%3d  ",k);
fprintf(stderr,"%10.5f ",guess[k].a);
fprintf(stderr,"%10.5f ",guess[k].b);
fprintf(stderr,"%10.5f ",guess[k].g);
fprintf(stderr,"%10.5f\n",guess[k].distance);
}
}
}

/*:6*/
#line 815 "iad_find.w"

/*38:*/
#line 825 "iad_find.w"


if(guess[0].b> r->default_bs){
p[1][1]= b2bcalc(guess[0].b-r->default_bs);
p[2][1]= b2bcalc(2*(guess[0].b-r->default_bs));
p[3][1]= p[1][1];
}else{
p[1][1]= b2bcalc(0.0001);
p[2][1]= b2bcalc(0.001);
p[3][1]= p[1][1];
}

p[1][2]= g2gcalc(guess[0].g);
p[2][2]= p[1][2];
p[3][2]= g2gcalc(0.9*guess[0].g+0.05);

/*:38*/
#line 816 "iad_find.w"

/*39:*/
#line 841 "iad_find.w"


for(i= 1;i<=3;i++){
x[1]= p[i][1];
x[2]= p[i][2];
y[i]= Find_BaG_fn(x);
}

/*:39*/
#line 817 "iad_find.w"

amoeba(p,y,2,r->tolerance,Find_BaG_fn,&r->AD_iterations,NULL,NULL);
/*40:*/
#line 853 "iad_find.w"

r->final_distance= 10;
for(i= 1;i<=3;i++){
if(y[i]<r->final_distance){
r->slab.b= bcalc2b(p[i][1])+r->default_bs;
r->slab.a= r->default_bs/r->slab.b;
r->slab.g= gcalc2g(p[i][2]);
r->final_distance= y[i];
}
}

/*:40*/
#line 819 "iad_find.w"


/*13:*/
#line 288 "iad_find.w"

free_dvector(x,1,2);
free_dvector(y,1,3);
free_dmatrix(p,1,3,1,2);

/*:13*/
#line 821 "iad_find.w"

/*12:*/
#line 280 "iad_find.w"

r->a= r->slab.a;
r->b= r->slab.b;
r->g= r->slab.g;
r->found= (r->tolerance> r->final_distance);
Set_Calc_State(m,*r);

/*:12*/
#line 822 "iad_find.w"

}

/*:37*/
#line 54 "iad_find.w"

/*42:*/
#line 877 "iad_find.w"

/*41:*/
#line 874 "iad_find.w"

void U_Find_BsG(struct measure_type m,struct invert_type*r)

/*:41*/
#line 878 "iad_find.w"

{
/*5:*/
#line 116 "iad_find.w"

int i,i_best,j_best;
double*x,*y,**p;

x= dvector(1,2);
y= dvector(1,3);
p= dmatrix(1,3,1,2);

/*:5*/
#line 880 "iad_find.w"


if(Debug(DEBUG_SEARCH)){
fprintf(stderr,"SEARCH: Using U_Find_BsG()");
fprintf(stderr," (mu=%6.4f)",r->slab.cos_angle);
if(r->default_ba!=UNINITIALIZED)
fprintf(stderr,"  default_ba = %8.5f",r->default_ba);
fprintf(stderr,"\n");
}

Set_Calc_State(m,*r);
/*6:*/
#line 133 "iad_find.w"

{

size_t count= NUMBER_OF_GUESSES;

abg_distance(r->slab.a,r->slab.b,r->slab.g,&(guess[0]));

if(!Valid_Grid(m,*r))Fill_Grid(m,*r,1);


Near_Grid_Points(m.m_r,m.m_t,r->search,&i_best,&j_best);
Grid_ABG(i_best,j_best,&(guess[1]));
Grid_ABG(i_best+1,j_best,&(guess[2]));
Grid_ABG(i_best-1,j_best,&(guess[3]));
Grid_ABG(i_best,j_best+1,&(guess[4]));
Grid_ABG(i_best,j_best-1,&(guess[5]));
Grid_ABG(i_best+1,j_best+1,&(guess[6]));
Grid_ABG(i_best-1,j_best-1,&(guess[7]));
Grid_ABG(i_best+1,j_best-1,&(guess[8]));
Grid_ABG(i_best-1,j_best+1,&(guess[9]));

qsort((void*)guess,count,sizeof(guess_type),compare_guesses);

if(Debug(DEBUG_BEST_GUESS)){
int k;
fprintf(stderr,"BEST: GRID GUESSES\n");
fprintf(stderr,"BEST:  k      albedo          b          g   distance\n");
for(k= 0;k<=6;k++){
fprintf(stderr,"BEST:%3d  ",k);
fprintf(stderr,"%10.5f ",guess[k].a);
fprintf(stderr,"%10.5f ",guess[k].b);
fprintf(stderr,"%10.5f ",guess[k].g);
fprintf(stderr,"%10.5f\n",guess[k].distance);
}
}
}

/*:6*/
#line 891 "iad_find.w"

/*43:*/
#line 901 "iad_find.w"


p[1][1]= b2bcalc(guess[0].b-r->default_ba);
p[1][2]= g2gcalc(guess[0].g);

p[2][1]= b2bcalc(2*guess[0].b-2*r->default_ba);
p[2][2]= p[1][2];

p[3][1]= p[1][1];
p[3][2]= g2gcalc(0.9*guess[0].g+0.05);


/*:43*/
#line 892 "iad_find.w"

/*44:*/
#line 913 "iad_find.w"


for(i= 1;i<=3;i++){
x[1]= p[i][1];
x[2]= p[i][2];
y[i]= Find_BsG_fn(x);
}

/*:44*/
#line 893 "iad_find.w"

amoeba(p,y,2,r->tolerance,Find_BsG_fn,&r->AD_iterations,NULL,NULL);
/*45:*/
#line 921 "iad_find.w"

r->final_distance= 10;
for(i= 1;i<=3;i++){
if(y[i]<r->final_distance){
r->slab.b= bcalc2b(p[i][1])+r->default_ba;
r->slab.a= 1-r->default_ba/r->slab.b;
r->slab.g= gcalc2g(p[i][2]);
r->final_distance= y[i];
}
}

/*:45*/
#line 895 "iad_find.w"


/*13:*/
#line 288 "iad_find.w"

free_dvector(x,1,2);
free_dvector(y,1,3);
free_dmatrix(p,1,3,1,2);

/*:13*/
#line 897 "iad_find.w"

/*12:*/
#line 280 "iad_find.w"

r->a= r->slab.a;
r->b= r->slab.b;
r->g= r->slab.g;
r->found= (r->tolerance> r->final_distance);
Set_Calc_State(m,*r);

/*:12*/
#line 898 "iad_find.w"

}

/*:42*/
#line 55 "iad_find.w"


/*:1*/
