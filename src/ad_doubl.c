/*1:*/
#line 11 "ad_doubl.w"

#include <math.h> 
#include <float.h> 
#include <stdio.h> 
#include "nr_util.h"
#include "ad_matrx.h"
#include "ad_globl.h"
#include "ad_doubl.h"

/*21:*/
#line 404 "ad_doubl.w"

static void Star_Multiply(int n,double**A,double**B,double**C)
{
Right_Diagonal_Multiply(n,A,twoaw,C);
Matrix_Multiply(n,C,B,C);
}

/*:21*/
#line 20 "ad_doubl.w"

/*22:*/
#line 414 "ad_doubl.w"

static void Star_One_Minus(int n,double**A)
{
int i,j;

for(i= 1;i<=n;i++){
for(j= 1;j<=n;j++)
A[i][j]*= -1;
A[i][i]+= 1.0/twoaw[i];
}

}

/*:22*/
#line 21 "ad_doubl.w"

/*3:*/
#line 92 "ad_doubl.w"


static void Basic_Add_Layers(int n,
double**R10,double**T01,
double**R12,double**R21,double**T12,double**T21,
double**R20,double**T02,
double**a,double**b)
{
Star_Multiply(n,R10,R12,a);
Star_One_Minus(n,a);
Left_Inverse_Multiply(n,a,T12,b);



Matrix_Multiply(n,b,R10,a);

Star_Multiply(n,a,T21,a);

Matrix_Sum(n,R21,a,R20);
Copy_Matrix(n,T01,a);
Matrix_Multiply(n,b,a,T02);
}

/*:3*/
#line 22 "ad_doubl.w"

/*4:*/
#line 150 "ad_doubl.w"


static void Basic_Add_Layers_With_Sources(int n,
double**R10,double**T01,
double**R12,double**R21,double**T12,double**T21,
double**R20,double**T02,
double**J01,double**J12,double**J21,double**J02,
double**a,double**b)
{
Star_Multiply(n,R10,R12,a);
Star_One_Minus(n,a);
Left_Inverse_Multiply(n,a,T12,b);



Matrix_Multiply(n,b,R10,a);

Star_Multiply(n,a,T21,a);

Matrix_Sum(n,R21,a,R20);
Copy_Matrix(n,T01,a);
Matrix_Multiply(n,b,a,T02);

Star_Multiply(n,R10,J21,a);
Matrix_Sum(n,J01,a,a);
Matrix_Multiply(n,b,a,J02);


Matrix_Sum(n,J02,J12,J02);
}

/*:4*/
#line 23 "ad_doubl.w"

/*7:*/
#line 196 "ad_doubl.w"

/*6:*/
#line 184 "ad_doubl.w"

void Add(int n,
double**R01,double**R10,double**T01,double**T10,
double**R12,double**R21,double**T12,double**T21,
double**R02,double**R20,double**T02,double**T20)

/*:6*/
#line 197 "ad_doubl.w"

{
/*23:*/
#line 427 "ad_doubl.w"

double**a,**b;

a= dmatrix(1,n,1,n);
b= dmatrix(1,n,1,n);

/*:23*/
#line 199 "ad_doubl.w"


Basic_Add_Layers(n,R10,T01,R12,R21,T12,T21,R20,T02,a,b);
Basic_Add_Layers(n,R12,T21,R10,R01,T10,T01,R02,T20,a,b);

/*24:*/
#line 434 "ad_doubl.w"

free_dmatrix(a,1,n,1,n);
free_dmatrix(b,1,n,1,n);/*:24*/
#line 204 "ad_doubl.w"

}

/*:7*/
#line 24 "ad_doubl.w"

/*9:*/
#line 220 "ad_doubl.w"

/*8:*/
#line 208 "ad_doubl.w"

void Add_With_Sources(int n,
double**R01,double**R10,double**T01,double**T10,double**J01,double**J10,
double**R12,double**R21,double**T12,double**T21,double**J12,double**J21,
double**R02,double**R20,double**T02,double**T20,double**J02,double**J20)

/*:8*/
#line 221 "ad_doubl.w"

{
/*23:*/
#line 427 "ad_doubl.w"

double**a,**b;

a= dmatrix(1,n,1,n);
b= dmatrix(1,n,1,n);

/*:23*/
#line 223 "ad_doubl.w"


Basic_Add_Layers_With_Sources(n,R10,T01,R12,R21,T12,T21,R20,T02,J01,J12,J21,J02,a,b);
Basic_Add_Layers_With_Sources(n,R12,T21,R10,R01,T10,T01,R02,T20,J21,J10,J01,J20,a,b);

/*24:*/
#line 434 "ad_doubl.w"

free_dmatrix(a,1,n,1,n);
free_dmatrix(b,1,n,1,n);/*:24*/
#line 228 "ad_doubl.w"

}

/*:9*/
#line 25 "ad_doubl.w"

/*11:*/
#line 239 "ad_doubl.w"

/*10:*/
#line 232 "ad_doubl.w"

void Add_Homogeneous(int n,
double**R01,double**T01,
double**R12,double**T12,
double**R02,double**T02)

/*:10*/
#line 240 "ad_doubl.w"

{
/*23:*/
#line 427 "ad_doubl.w"

double**a,**b;

a= dmatrix(1,n,1,n);
b= dmatrix(1,n,1,n);

/*:23*/
#line 242 "ad_doubl.w"


Basic_Add_Layers(n,R01,T01,R12,R12,T12,T12,R02,T02,a,b);

/*24:*/
#line 434 "ad_doubl.w"

free_dmatrix(a,1,n,1,n);
free_dmatrix(b,1,n,1,n);/*:24*/
#line 246 "ad_doubl.w"

}

/*:11*/
#line 26 "ad_doubl.w"

/*13:*/
#line 258 "ad_doubl.w"

/*12:*/
#line 254 "ad_doubl.w"

void Double_Once(int n,double**R,double**T)

/*:12*/
#line 259 "ad_doubl.w"

{
/*23:*/
#line 427 "ad_doubl.w"

double**a,**b;

a= dmatrix(1,n,1,n);
b= dmatrix(1,n,1,n);

/*:23*/
#line 261 "ad_doubl.w"

Basic_Add_Layers(n,R,T,R,R,T,T,R,T,a,b);
/*24:*/
#line 434 "ad_doubl.w"

free_dmatrix(a,1,n,1,n);
free_dmatrix(b,1,n,1,n);/*:24*/
#line 263 "ad_doubl.w"

}

/*:13*/
#line 27 "ad_doubl.w"

/*15:*/
#line 277 "ad_doubl.w"

/*14:*/
#line 273 "ad_doubl.w"

void Double_Until(int n,double**r,double**t,double start,double end)

/*:14*/
#line 278 "ad_doubl.w"

{
if(end==HUGE_VAL){
Double_Until_Infinite(n,r,t);
return;
}

{
/*23:*/
#line 427 "ad_doubl.w"

double**a,**b;

a= dmatrix(1,n,1,n);
b= dmatrix(1,n,1,n);

/*:23*/
#line 286 "ad_doubl.w"

while(fabs(end-start)> 0.00001&&end> start){
Basic_Add_Layers(n,r,t,r,r,t,t,r,t,a,b);
start*= 2;
}

/*24:*/
#line 434 "ad_doubl.w"

free_dmatrix(a,1,n,1,n);
free_dmatrix(b,1,n,1,n);/*:24*/
#line 292 "ad_doubl.w"

}
}

/*:15*/
#line 28 "ad_doubl.w"

/*17:*/
#line 309 "ad_doubl.w"

/*16:*/
#line 305 "ad_doubl.w"

void Double_Until_Infinite(int n,double**r,double**t)

/*:16*/
#line 310 "ad_doubl.w"

{
double oldutu,UTU,UT1;

/*23:*/
#line 427 "ad_doubl.w"

double**a,**b;

a= dmatrix(1,n,1,n);
b= dmatrix(1,n,1,n);

/*:23*/
#line 314 "ad_doubl.w"

UTU= 0.0;
do{
oldutu= UTU;
Basic_Add_Layers(n,r,t,r,r,t,t,r,t,a,b);
URU_and_UR1(n,1.0,t,&UTU,&UT1);
}while(fabs(UTU-oldutu)>=0.000001);

/*24:*/
#line 434 "ad_doubl.w"

free_dmatrix(a,1,n,1,n);
free_dmatrix(b,1,n,1,n);/*:24*/
#line 322 "ad_doubl.w"


}

/*:17*/
#line 29 "ad_doubl.w"

/*19:*/
#line 358 "ad_doubl.w"

/*18:*/
#line 352 "ad_doubl.w"

void Between(int n,
double**R01,double**R10,double**T01,double**T10,
double**R12,double**R21,double**T12,double**T21,
double**Lup,double**Ldown)

/*:18*/
#line 359 "ad_doubl.w"

{
/*23:*/
#line 427 "ad_doubl.w"

double**a,**b;

a= dmatrix(1,n,1,n);
b= dmatrix(1,n,1,n);

/*:23*/
#line 361 "ad_doubl.w"


Star_Multiply(n,R10,R12,a);
Star_One_Minus(n,a);
Right_Inverse_Multiply(n,a,T01,Ldown);

Star_Multiply(n,R12,R10,a);
Star_One_Minus(n,a);
Right_Inverse_Multiply(n,a,R12,b);
Star_Multiply(n,b,T01,Lup);

/*24:*/
#line 434 "ad_doubl.w"

free_dmatrix(a,1,n,1,n);
free_dmatrix(b,1,n,1,n);/*:24*/
#line 372 "ad_doubl.w"

}

/*:19*/
#line 30 "ad_doubl.w"


/*:1*/
