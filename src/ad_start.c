/*1:*/
#line 7 "./ad_start.w"

#include <math.h> 
#include <float.h> 
#include <stdio.h> 

#include "ad_frsnl.h"
#include "ad_globl.h"
#include "ad_matrx.h"
#include "ad_phase.h"
#include "ad_radau.h"
#include "ad_start.h"
#include "nr_gaulg.h"
#include "nr_util.h"

/*5:*/
#line 77 "./ad_start.w"

/*4:*/
#line 74 "./ad_start.w"

double Get_Start_Depth(double mu,double d)

/*:4*/
#line 78 "./ad_start.w"


{
if(d<=0)
return 0.0;

if(d==HUGE_VAL)
return(mu/2.0);

while(d> mu)d/= 2;

return d;
}

/*:5*/
#line 21 "./ad_start.w"

/*8:*/
#line 103 "./ad_start.w"

/*7:*/
#line 100 "./ad_start.w"

void Quadrature(int n,double n_slab,double*x,double*w)

/*:7*/
#line 104 "./ad_start.w"

{
int i,nby2;
double*x1,*w1;
double mu_c;

if(n_slab==1){
Radau(0.0,1.0,x,w,n);
return;
}

mu_c= Cos_Critical_Angle(n_slab,1.0);
nby2= n/2;
gauleg(0.0,mu_c,x,w,nby2);

x1= dvector(1,nby2);
w1= dvector(1,nby2);
Radau(mu_c,1.0,x1,w1,nby2);
for(i= 1;i<=nby2;i++){
x[nby2+i]= x1[i];
w[nby2+i]= w1[i];
}
free_dvector(x1,1,nby2);
free_dvector(w1,1,nby2);
}

/*:8*/
#line 22 "./ad_start.w"

/*10:*/
#line 149 "./ad_start.w"

/*9:*/
#line 146 "./ad_start.w"

void Choose_Method(struct AD_slab_type*slab,struct AD_method_type*method)

/*:9*/
#line 150 "./ad_start.w"


{
double af;
int i,n;

if(0<slab->cos_angle&&slab->cos_angle<1){
Choose_Cone_Method(slab,method);
return;
}

n= method->quad_pts;
af= pow(slab->g,n)*slab->a;
method->a_calc= (slab->a-af)/(1-af);
method->b_calc= (1-af)*slab->b;
method->g_calc= slab->g;

Quadrature(n,slab->n_slab,angle,weight);

for(i= 1;i<=n;i++)
twoaw[i]= 2*angle[i]*weight[i];

method->b_thinnest= Get_Start_Depth(angle[1],method->b_calc);
}

/*:10*/
#line 23 "./ad_start.w"

/*12:*/
#line 185 "./ad_start.w"

/*11:*/
#line 181 "./ad_start.w"

void Choose_Cone_Method(struct AD_slab_type*slab,
struct AD_method_type*method)

/*:11*/
#line 186 "./ad_start.w"


{
double af,*angle1,*weight1,cos_crit_angle,mu;
int i,n,nby2,nby3;

n= method->quad_pts;
af= pow(slab->g,n)*slab->a;
method->a_calc= (slab->a-af)/(1-af);
method->b_calc= (1-af)*slab->b;
method->g_calc= slab->g;

/*15:*/
#line 223 "./ad_start.w"

if(slab->cos_angle==0||slab->cos_angle==1){
Choose_Method(slab,method);
/*13:*/
#line 205 "./ad_start.w"


/*:13*/
#line 226 "./ad_start.w"

return;
}

/*:15*/
#line 198 "./ad_start.w"

/*16:*/
#line 235 "./ad_start.w"

if(slab->n_slab==1&&slab->n_top_slide==1&&slab->n_bottom_slide==1){
nby2= n/2;
Radau(0.0,slab->cos_angle,angle,weight,nby2);

angle1= dvector(1,nby2);
weight1= dvector(1,nby2);
Radau(slab->cos_angle,1.0,angle1,weight1,nby2);

for(i= 1;i<=nby2;i++){
angle[nby2+i]= angle1[i];
weight[nby2+i]= weight1[i];
}
free_dvector(angle1,1,nby2);
free_dvector(weight1,1,nby2);

for(i= 1;i<=n;i++)
twoaw[i]= 2*angle[i]*weight[i];

method->b_thinnest= Get_Start_Depth(angle[1],method->b_calc);

/*13:*/
#line 205 "./ad_start.w"


/*:13*/
#line 256 "./ad_start.w"


return;
}

/*:16*/
#line 199 "./ad_start.w"

/*17:*/
#line 275 "./ad_start.w"

cos_crit_angle= Cos_Critical_Angle(slab->n_slab,1.0);
nby3= n/3;
gauleg(0.0,cos_crit_angle,angle,weight,nby3);

/*:17*/
#line 200 "./ad_start.w"

/*18:*/
#line 280 "./ad_start.w"


mu= sqrt(slab->n_slab*slab->n_slab-1+slab->cos_angle*slab->cos_angle)/slab->n_slab;
angle1= dvector(1,nby3);
weight1= dvector(1,nby3);
Radau(cos_crit_angle,mu,angle1,weight1,nby3);
for(i= 1;i<=nby3;i++){
angle[nby3+i]= angle1[i];
weight[nby3+i]= weight1[i];
}

/*:18*/
#line 201 "./ad_start.w"

/*19:*/
#line 291 "./ad_start.w"

Radau(mu,1.0,angle1,weight1,nby3);
for(i= 1;i<=nby3;i++){
angle[nby3*2+i]= angle1[i];
weight[nby3*2+i]= weight1[i];
}
free_dvector(angle1,1,nby3);
free_dvector(weight1,1,nby3);

for(i= 1;i<=n;i++)
twoaw[i]= 2*angle[i]*weight[i];

method->b_thinnest= Get_Start_Depth(angle[1],method->b_calc);

/*13:*/
#line 205 "./ad_start.w"


/*:13*/
#line 305 "./ad_start.w"



/*:19*/
#line 202 "./ad_start.w"

}

/*:12*/
#line 24 "./ad_start.w"

/*22:*/
#line 397 "./ad_start.w"

static void Get_IGI_Layer(struct AD_method_type method,double**h,double**R,double**T)
{
int i,j,n;
double a,c,d,temp;

a= method.a_calc;
d= method.b_thinnest;
n= method.quad_pts;

for(j= 1;j<=n;j++){
temp= a*d/4/angle[j];
for(i= 1;i<=n;i++){
c= temp/angle[i];
R[i][j]= c*h[i][-j];
T[i][j]= c*h[i][j];
}
T[j][j]+= (1-d/angle[j])/twoaw[j];
}
}

/*:22*/
#line 25 "./ad_start.w"

/*23:*/
#line 442 "./ad_start.w"

static void Get_Diamond_Layer(struct AD_method_type method,double**h,double**R,double**T)
{
/*31:*/
#line 671 "./ad_start.w"


int i,j,n;
double**A,**G,**C;
double a,c,d,temp;
double*work;
double condition;
int*ipvt;

d= method.b_thinnest;
a= method.a_calc;
n= method.quad_pts;

A= dmatrix(1,n,1,n);
G= dmatrix(1,n,1,n);
C= dmatrix(1,n,1,n);
work= dvector(1,n);
ipvt= ivector(1,n);

if(d<1e-4)
AD_error("**** Roundoff error is a problem--Use IGI method\n");

/*:31*/
#line 445 "./ad_start.w"

/*24:*/
#line 476 "./ad_start.w"


for(j= 1;j<=n;j++){
temp= a*d*weight[j]/4;
for(i= 1;i<=n;i++){
c= temp/angle[i];
R[i][j]= c*h[i][-j];
T[i][j]= -c*h[i][j];
}
T[j][j]+= d/(2*angle[j]);
}


/*:24*/
#line 446 "./ad_start.w"

/*25:*/
#line 509 "./ad_start.w"


for(i= 1;i<=n;i++){
for(j= 1;j<=n;j++)
A[i][j]= T[i][j];
A[i][i]+= 1.0;
}

Left_Inverse_Multiply(n,A,R,C);

/*:25*/
#line 447 "./ad_start.w"

/*26:*/
#line 525 "./ad_start.w"

Matrix_Multiply(n,C,R,G);
for(i= 1;i<=n;i++){
for(j= 1;j<=n;j++)
G[i][j]= (T[i][j]-G[i][j])/2;
G[i][i]+= 0.5;
}

/*:26*/
#line 448 "./ad_start.w"

/*27:*/
#line 540 "./ad_start.w"


#ifdef MARTIN_HAMMER
{
double**Ginv,**G2;

if(Martin_Hammer!=0){
printf("A from equation 5.55\n");
wrmatrix(n,T);

printf("B from equation 5.55\n");
wrmatrix(n,R);

Ginv= dmatrix(1,n,1,n);
G2= dmatrix(1,n,1,n);

for(i= 1;i<=n;i++){
for(j= 1;j<=n;j++){
G2[i][j]= G[i][j]*2.0;
}
}

Matrix_Inverse(n,G2,Ginv);

printf("Inverse of G from equation 5.56\n");
wrmatrix(n,G2);

printf("G from equation 5.56\n");
wrmatrix(n,Ginv);

free_matrix(Ginv,1,n,1,n);
free_matrix(G2,1,n,1,n);
}
}
#endif

/*:27*/
#line 449 "./ad_start.w"

/*28:*/
#line 587 "./ad_start.w"

Transpose_Matrix(n,G);
Decomp(n,G,&condition,ipvt);

if(condition==1e32)
AD_error("Singular Matrix ... failed in diamond_init\n");

for(i= 1;i<=n;i++){
/*29:*/
#line 645 "./ad_start.w"

for(j= 1;j<=n;j++)
work[j]= C[j][i]*twoaw[j]/twoaw[i];
Solve(n,G,work,ipvt);
for(j= 1;j<=n;j++)
R[i][j]= work[j]/twoaw[j];

/*:29*/
#line 595 "./ad_start.w"

/*30:*/
#line 659 "./ad_start.w"

for(j= 1;j<=n;j++)
work[j]= 0;
work[i]= 1.0;
Solve(n,G,work,ipvt);
for(j= 1;j<=n;j++)
T[i][j]= work[j]/twoaw[j];
T[i][i]-= 1.0/twoaw[i];

/*:30*/
#line 596 "./ad_start.w"

}

#ifdef MARTIN_HAMMER
{
double**T2,**Ginv;

if(Martin_Hammer==5){

T2= dmatrix(1,n,1,n);
Ginv= dmatrix(1,n,1,n);

Copy_Matrix(n,T,T2);

for(i= 1;i<=n;i++){
T2[i][i]+= 1/twoaw[i];
}

for(i= 1;i<=n;i++){
for(j= 1;j<=n;j++){
T2[i][j]*= twoaw[j]*0.5;
}
}

printf("G=(T-1)/2 from equation 5.55\n");
wrmatrix(n,T2);

Matrix_Inverse(n,T2,Ginv);

printf("1/G\n");
wrmatrix(n,Ginv);

free_matrix(T2,1,n,1,n);
free_matrix(Ginv,1,n,1,n);
}
}
#endif

/*:28*/
#line 450 "./ad_start.w"

/*32:*/
#line 694 "./ad_start.w"

free_dvector(work,1,n);
free_ivector(ipvt,1,n);
free_dmatrix(A,1,n,1,n);
free_dmatrix(G,1,n,1,n);
free_dmatrix(C,1,n,1,n);

/*:32*/
#line 451 "./ad_start.w"

}


/*:23*/
#line 26 "./ad_start.w"

/*35:*/
#line 709 "./ad_start.w"

/*34:*/
#line 706 "./ad_start.w"

void Init_Layer(struct AD_slab_type slab,struct AD_method_type method,double**R,double**T)

/*:34*/
#line 710 "./ad_start.w"


{
double**h;
int n;

n= method.quad_pts;

if(slab.b<=0){
Zero_Layer(n,R,T);
return;
}

h= dmatrix(-n,n,-n,n);
Get_Phi(n,slab.phase_function,method.g_calc,h);

if(method.b_thinnest<1e-4||method.b_thinnest<0.09*angle[1])
Get_IGI_Layer(method,h,R,T);
else
Get_Diamond_Layer(method,h,R,T);

free_dmatrix(h,-n,n,-n,n);
}/*:35*/
#line 27 "./ad_start.w"


/*:1*/
