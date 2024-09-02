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
#include "iad_find.h"

#define NUMBER_OF_GUESSES 10

guess_type guess[NUMBER_OF_GUESSES];

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

/*15:*/
#line 308 "iad_find.w"

/*14:*/
#line 305 "iad_find.w"

void U_Find_Ba(struct measure_type m,struct invert_type*r)

/*:14*/
#line 309 "iad_find.w"

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

/*10:*/
#line 222 "iad_find.w"

r->a= r->slab.a;
r->b= r->slab.b;
r->g= r->slab.g;
r->found= (r->tolerance> r->final_distance);
Set_Calc_State(m,*r);

/*:10*/
#line 344 "iad_find.w"

}

/*:15*/
#line 34 "iad_find.w"

/*13:*/
#line 257 "iad_find.w"

/*12:*/
#line 254 "iad_find.w"

void U_Find_Bs(struct measure_type m,struct invert_type*r)

/*:12*/
#line 258 "iad_find.w"

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

/*10:*/
#line 222 "iad_find.w"

r->a= r->slab.a;
r->b= r->slab.b;
r->g= r->slab.g;
r->found= (r->tolerance> r->final_distance);
Set_Calc_State(m,*r);

/*:10*/
#line 293 "iad_find.w"

}

/*:13*/
#line 35 "iad_find.w"

/*17:*/
#line 359 "iad_find.w"

/*16:*/
#line 356 "iad_find.w"

void U_Find_A(struct measure_type m,struct invert_type*r)

/*:16*/
#line 360 "iad_find.w"
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

/*10:*/
#line 222 "iad_find.w"

r->a= r->slab.a;
r->b= r->slab.b;
r->g= r->slab.g;
r->found= (r->tolerance> r->final_distance);
Set_Calc_State(m,*r);

/*:10*/
#line 395 "iad_find.w"

}

/*:17*/
#line 36 "iad_find.w"

/*21:*/
#line 455 "iad_find.w"

/*20:*/
#line 452 "iad_find.w"

void U_Find_B(struct measure_type m,struct invert_type*r)

/*:20*/
#line 456 "iad_find.w"
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

/*22:*/
#line 484 "iad_find.w"

{
double x,ax,bx,cx,fa,fb,fc;

ax= b2bcalc(0.1);
bx= b2bcalc(10);

mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,Find_B_fn,&r->AD_iterations);
r->final_distance= brent(ax,bx,cx,Find_B_fn,r->tolerance,&x,&r->AD_iterations);
r->slab.b= bcalc2b(x);
Set_Calc_State(m,*r);
}


/*:22*/
#line 477 "iad_find.w"


/*10:*/
#line 222 "iad_find.w"

r->a= r->slab.a;
r->b= r->slab.b;
r->g= r->slab.g;
r->found= (r->tolerance> r->final_distance);
Set_Calc_State(m,*r);

/*:10*/
#line 479 "iad_find.w"

}

/*:21*/
#line 37 "iad_find.w"

/*19:*/
#line 403 "iad_find.w"

/*18:*/
#line 400 "iad_find.w"

void U_Find_G(struct measure_type m,struct invert_type*r)

/*:18*/
#line 404 "iad_find.w"
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


/*10:*/
#line 222 "iad_find.w"

r->a= r->slab.a;
r->b= r->slab.b;
r->g= r->slab.g;
r->found= (r->tolerance> r->final_distance);
Set_Calc_State(m,*r);

/*:10*/
#line 442 "iad_find.w"

}


/*:19*/
#line 38 "iad_find.w"

/*24:*/
#line 514 "iad_find.w"

/*23:*/
#line 511 "iad_find.w"

void U_Find_AG(struct measure_type m,struct invert_type*r)

/*:23*/
#line 515 "iad_find.w"

{
/*5:*/
#line 98 "iad_find.w"

int i,i_best,j_best;
double*x,*y,**p;

x= dvector(1,2);
y= dvector(1,3);
p= dmatrix(1,3,1,2);

/*:5*/
#line 517 "iad_find.w"


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
/*6:*/
#line 115 "iad_find.w"

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
#line 539 "iad_find.w"

/*25:*/
#line 549 "iad_find.w"

{
int k,kk;

p[1][1]= a2acalc(guess[0].a);
p[1][2]= g2gcalc(guess[0].g);

for(k= 1;k<7;k++){
if(guess[0].a!=guess[k].a)
break;
}

p[2][1]= a2acalc(guess[k].a);
p[2][2]= g2gcalc(guess[k].g);

for(kk= 1;kk<7;kk++){
if(kk==k)
continue;
if(guess[0].g!=guess[kk].g||guess[k].g!=guess[kk].g)
break;
}
p[3][1]= a2acalc(guess[kk].a);
p[3][2]= g2gcalc(guess[kk].g);

if(Debug(DEBUG_BEST_GUESS)){
fprintf(stderr,"-----------------------------------------------------\n");
fprintf(stderr,"BEST: <1> ");
fprintf(stderr,"%10.5f ",guess[0].a);
fprintf(stderr,"%10.5f ",guess[0].b);
fprintf(stderr,"%10.5f ",guess[0].g);
fprintf(stderr,"%10.5f\n",guess[0].distance);
fprintf(stderr,"BEST: <2> ");
fprintf(stderr,"%10.5f ",guess[k].a);
fprintf(stderr,"%10.5f ",guess[k].b);
fprintf(stderr,"%10.5f ",guess[k].g);
fprintf(stderr,"%10.5f\n",guess[k].distance);
fprintf(stderr,"BEST: <3> ");
fprintf(stderr,"%10.5f ",guess[kk].a);
fprintf(stderr,"%10.5f ",guess[kk].b);
fprintf(stderr,"%10.5f ",guess[kk].g);
fprintf(stderr,"%10.5f\n",guess[kk].distance);
fprintf(stderr,"\n");
}
}

/*:25*/
#line 540 "iad_find.w"

/*26:*/
#line 594 "iad_find.w"


for(i= 1;i<=3;i++){
x[1]= p[i][1];
x[2]= p[i][2];
y[i]= Find_AG_fn(x);
}


/*:26*/
#line 541 "iad_find.w"

amoeba(p,y,2,r->tolerance,Find_AG_fn,&r->AD_iterations);
/*27:*/
#line 607 "iad_find.w"

r->final_distance= 10;
for(i= 1;i<=3;i++){
if(y[i]<r->final_distance){
r->slab.a= acalc2a(p[i][1]);
r->slab.g= gcalc2g(p[i][2]);
r->final_distance= y[i];
}
}


/*:27*/
#line 543 "iad_find.w"

/*11:*/
#line 230 "iad_find.w"

free_dvector(x,1,2);
free_dvector(y,1,3);
free_dmatrix(p,1,3,1,2);

/*:11*/
#line 544 "iad_find.w"


/*10:*/
#line 222 "iad_find.w"

r->a= r->slab.a;
r->b= r->slab.b;
r->g= r->slab.g;
r->found= (r->tolerance> r->final_distance);
Set_Calc_State(m,*r);

/*:10*/
#line 546 "iad_find.w"

}

/*:24*/
#line 39 "iad_find.w"

/*4:*/
#line 69 "iad_find.w"

/*3:*/
#line 66 "iad_find.w"

void U_Find_AB(struct measure_type m,struct invert_type*r)

/*:3*/
#line 70 "iad_find.w"

{
/*5:*/
#line 98 "iad_find.w"

int i,i_best,j_best;
double*x,*y,**p;

x= dvector(1,2);
y= dvector(1,3);
p= dmatrix(1,3,1,2);

/*:5*/
#line 72 "iad_find.w"


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

/*6:*/
#line 115 "iad_find.w"

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
#line 87 "iad_find.w"

/*7:*/
#line 152 "iad_find.w"

{
int k,kk;

p[1][1]= a2acalc(guess[0].a);
p[1][2]= b2bcalc(guess[0].b);

for(k= 1;k<7;k++){
if(guess[0].a!=guess[k].a)
break;
}

p[2][1]= a2acalc(guess[k].a);
p[2][2]= b2bcalc(guess[k].b);

for(kk= 1;kk<7;kk++){
if(k==kk)
continue;
if(guess[0].b!=guess[kk].b||guess[k].b!=guess[kk].b)
break;
}

p[3][1]= a2acalc(guess[kk].a);
p[3][2]= b2bcalc(guess[kk].b);

if(Debug(DEBUG_BEST_GUESS)){
fprintf(stderr,"-----------------------------------------------------\n");
fprintf(stderr,"BEST: <1> ");
fprintf(stderr,"%10.5f ",guess[0].a);
fprintf(stderr,"%10.5f ",guess[0].b);
fprintf(stderr,"%10.5f ",guess[0].g);
fprintf(stderr,"%10.5f\n",guess[0].distance);
fprintf(stderr,"BEST: <2> ");
fprintf(stderr,"%10.5f ",guess[k].a);
fprintf(stderr,"%10.5f ",guess[k].b);
fprintf(stderr,"%10.5f ",guess[k].g);
fprintf(stderr,"%10.5f\n",guess[k].distance);
fprintf(stderr,"BEST: <3> ");
fprintf(stderr,"%10.5f ",guess[kk].a);
fprintf(stderr,"%10.5f ",guess[kk].b);
fprintf(stderr,"%10.5f ",guess[kk].g);
fprintf(stderr,"%10.5f\n",guess[kk].distance);
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
#line 88 "iad_find.w"

/*8:*/
#line 204 "iad_find.w"


for(i= 1;i<=3;i++){
x[1]= p[i][1];
x[2]= p[i][2];
y[i]= Find_AB_fn(x);
}

/*:8*/
#line 89 "iad_find.w"

amoeba(p,y,2,r->tolerance,Find_AB_fn,&r->AD_iterations);
/*9:*/
#line 212 "iad_find.w"

r->final_distance= 10;
for(i= 1;i<=3;i++){
if(y[i]<r->final_distance){
r->slab.a= acalc2a(p[i][1]);
r->slab.b= bcalc2b(p[i][2]);
r->final_distance= y[i];
}
}

/*:9*/
#line 91 "iad_find.w"


/*11:*/
#line 230 "iad_find.w"

free_dvector(x,1,2);
free_dvector(y,1,3);
free_dmatrix(p,1,3,1,2);

/*:11*/
#line 93 "iad_find.w"

/*10:*/
#line 222 "iad_find.w"

r->a= r->slab.a;
r->b= r->slab.b;
r->g= r->slab.g;
r->found= (r->tolerance> r->final_distance);
Set_Calc_State(m,*r);

/*:10*/
#line 94 "iad_find.w"

}

/*:4*/
#line 40 "iad_find.w"

/*29:*/
#line 625 "iad_find.w"

/*28:*/
#line 622 "iad_find.w"

void U_Find_BG(struct measure_type m,struct invert_type*r)

/*:28*/
#line 626 "iad_find.w"

{
/*5:*/
#line 98 "iad_find.w"

int i,i_best,j_best;
double*x,*y,**p;

x= dvector(1,2);
y= dvector(1,3);
p= dmatrix(1,3,1,2);

/*:5*/
#line 628 "iad_find.w"


if(Debug(DEBUG_SEARCH)){
fprintf(stderr,"SEARCH: Using U_Find_BG()");
fprintf(stderr," (mu=%6.4f)",r->slab.cos_angle);
if(r->default_a!=UNINITIALIZED)
fprintf(stderr,"  default_a = %8.5f",r->default_a);
fprintf(stderr,"\n");
}

r->slab.a= (r->default_a==UNINITIALIZED)?0:r->default_a;
Set_Calc_State(m,*r);

/*6:*/
#line 115 "iad_find.w"

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
#line 641 "iad_find.w"

/*31:*/
#line 655 "iad_find.w"

{
int k,kk;

p[1][1]= b2bcalc(guess[0].b);
p[1][2]= g2gcalc(guess[0].g);

for(k= 1;k<7;k++){
if(guess[0].b!=guess[k].b)
break;
}

p[2][1]= b2bcalc(guess[k].b);
p[2][2]= g2gcalc(guess[k].g);

for(kk= 1;kk<7;kk++){
if(kk==k)
continue;
if(guess[0].g!=guess[kk].g||guess[k].g!=guess[kk].g)
break;
}
p[3][1]= b2bcalc(guess[kk].b);
p[3][2]= g2gcalc(guess[kk].g);

if(Debug(DEBUG_BEST_GUESS)){
fprintf(stderr,"-----------------------------------------------------\n");
fprintf(stderr,"BEST: <1> ");
fprintf(stderr,"%10.5f ",guess[0].a);
fprintf(stderr,"%10.5f ",guess[0].b);
fprintf(stderr,"%10.5f ",guess[0].g);
fprintf(stderr,"%10.5f\n",guess[0].distance);
fprintf(stderr,"BEST: <2> ");
fprintf(stderr,"%10.5f ",guess[k].a);
fprintf(stderr,"%10.5f ",guess[k].b);
fprintf(stderr,"%10.5f ",guess[k].g);
fprintf(stderr,"%10.5f\n",guess[k].distance);
fprintf(stderr,"BEST: <3> ");
fprintf(stderr,"%10.5f ",guess[kk].a);
fprintf(stderr,"%10.5f ",guess[kk].b);
fprintf(stderr,"%10.5f ",guess[kk].g);
fprintf(stderr,"%10.5f\n",guess[kk].distance);
fprintf(stderr,"\n");
}
}

/*:31*/
#line 642 "iad_find.w"

/*32:*/
#line 700 "iad_find.w"


for(i= 1;i<=3;i++){
x[1]= p[i][1];
x[2]= p[i][2];
y[i]= Find_BG_fn(x);
}

/*:32*/
#line 643 "iad_find.w"

amoeba(p,y,2,r->tolerance,Find_BG_fn,&r->AD_iterations);
/*33:*/
#line 712 "iad_find.w"

r->final_distance= 10;
for(i= 1;i<=3;i++){
if(y[i]<r->final_distance){
r->slab.b= bcalc2b(p[i][1]);
r->slab.g= gcalc2g(p[i][2]);
r->final_distance= y[i];
}
}

/*:33*/
#line 645 "iad_find.w"


/*11:*/
#line 230 "iad_find.w"

free_dvector(x,1,2);
free_dvector(y,1,3);
free_dmatrix(p,1,3,1,2);

/*:11*/
#line 647 "iad_find.w"

/*10:*/
#line 222 "iad_find.w"

r->a= r->slab.a;
r->b= r->slab.b;
r->g= r->slab.g;
r->found= (r->tolerance> r->final_distance);
Set_Calc_State(m,*r);

/*:10*/
#line 648 "iad_find.w"

}

/*:29*/
#line 41 "iad_find.w"

/*35:*/
#line 735 "iad_find.w"

/*34:*/
#line 732 "iad_find.w"

void U_Find_BaG(struct measure_type m,struct invert_type*r)

/*:34*/
#line 736 "iad_find.w"

{
/*5:*/
#line 98 "iad_find.w"

int i,i_best,j_best;
double*x,*y,**p;

x= dvector(1,2);
y= dvector(1,3);
p= dmatrix(1,3,1,2);

/*:5*/
#line 738 "iad_find.w"

Set_Calc_State(m,*r);
/*6:*/
#line 115 "iad_find.w"

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
#line 740 "iad_find.w"

/*36:*/
#line 750 "iad_find.w"


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

/*:36*/
#line 741 "iad_find.w"

/*37:*/
#line 766 "iad_find.w"


for(i= 1;i<=3;i++){
x[1]= p[i][1];
x[2]= p[i][2];
y[i]= Find_BaG_fn(x);
}

/*:37*/
#line 742 "iad_find.w"

amoeba(p,y,2,r->tolerance,Find_BaG_fn,&r->AD_iterations);
/*38:*/
#line 778 "iad_find.w"

r->final_distance= 10;
for(i= 1;i<=3;i++){
if(y[i]<r->final_distance){
r->slab.b= bcalc2b(p[i][1])+r->default_bs;
r->slab.a= r->default_bs/r->slab.b;
r->slab.g= gcalc2g(p[i][2]);
r->final_distance= y[i];
}
}

/*:38*/
#line 744 "iad_find.w"


/*11:*/
#line 230 "iad_find.w"

free_dvector(x,1,2);
free_dvector(y,1,3);
free_dmatrix(p,1,3,1,2);

/*:11*/
#line 746 "iad_find.w"

/*10:*/
#line 222 "iad_find.w"

r->a= r->slab.a;
r->b= r->slab.b;
r->g= r->slab.g;
r->found= (r->tolerance> r->final_distance);
Set_Calc_State(m,*r);

/*:10*/
#line 747 "iad_find.w"

}

/*:35*/
#line 42 "iad_find.w"

/*40:*/
#line 802 "iad_find.w"

/*39:*/
#line 799 "iad_find.w"

void U_Find_BsG(struct measure_type m,struct invert_type*r)

/*:39*/
#line 803 "iad_find.w"

{
/*5:*/
#line 98 "iad_find.w"

int i,i_best,j_best;
double*x,*y,**p;

x= dvector(1,2);
y= dvector(1,3);
p= dmatrix(1,3,1,2);

/*:5*/
#line 805 "iad_find.w"


if(Debug(DEBUG_SEARCH)){
fprintf(stderr,"SEARCH: Using U_Find_BsG()");
fprintf(stderr," (mu=%6.4f)",r->slab.cos_angle);
if(r->default_ba!=UNINITIALIZED)
fprintf(stderr,"  default_ba = %8.5f",r->default_ba);
fprintf(stderr,"\n");
}

Set_Calc_State(m,*r);
/*6:*/
#line 115 "iad_find.w"

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
#line 816 "iad_find.w"

/*41:*/
#line 826 "iad_find.w"


p[1][1]= b2bcalc(guess[0].b-r->default_ba);
p[1][2]= g2gcalc(guess[0].g);

p[2][1]= b2bcalc(2*guess[0].b-2*r->default_ba);
p[2][2]= p[1][2];

p[3][1]= p[1][1];
p[3][2]= g2gcalc(0.9*guess[0].g+0.05);


/*:41*/
#line 817 "iad_find.w"

/*42:*/
#line 838 "iad_find.w"


for(i= 1;i<=3;i++){
x[1]= p[i][1];
x[2]= p[i][2];
y[i]= Find_BsG_fn(x);
}

/*:42*/
#line 818 "iad_find.w"

amoeba(p,y,2,r->tolerance,Find_BsG_fn,&r->AD_iterations);
/*43:*/
#line 846 "iad_find.w"

r->final_distance= 10;
for(i= 1;i<=3;i++){
if(y[i]<r->final_distance){
r->slab.b= bcalc2b(p[i][1])+r->default_ba;
r->slab.a= 1-r->default_ba/r->slab.b;
r->slab.g= gcalc2g(p[i][2]);
r->final_distance= y[i];
}
}/*:43*/
#line 820 "iad_find.w"


/*11:*/
#line 230 "iad_find.w"

free_dvector(x,1,2);
free_dvector(y,1,3);
free_dmatrix(p,1,3,1,2);

/*:11*/
#line 822 "iad_find.w"

/*10:*/
#line 222 "iad_find.w"

r->a= r->slab.a;
r->b= r->slab.b;
r->g= r->slab.g;
r->found= (r->tolerance> r->final_distance);
Set_Calc_State(m,*r);

/*:10*/
#line 823 "iad_find.w"

}

/*:40*/
#line 43 "iad_find.w"


/*:1*/
