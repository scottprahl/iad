/*1:*/
#line 18 "ad_matrx.w"


#include<stddef.h> 
#include<math.h> 
#include"ad_globl.h"
#include"ad_matrx.h"
#include"nr_util.h"

/*5:*/
#line 67 "ad_matrx.w"

/*4:*/
#line 64 "ad_matrx.w"

void Copy_Matrix(int n,double**A,double**B)

/*:4*/
#line 68 "ad_matrx.w"

{
double*a_ptr,*b_ptr,*a_last;

a_last= &A[n][n];
a_ptr= &A[1][1];
b_ptr= &B[1][1];

while(a_ptr<=a_last)
*b_ptr++= *a_ptr++;
}

/*:5*/
#line 26 "ad_matrx.w"

/*7:*/
#line 84 "ad_matrx.w"

/*6:*/
#line 81 "ad_matrx.w"

void One_Minus(int n,double**A)

/*:6*/
#line 85 "ad_matrx.w"

{
int i,j;

for(i= 1;i<=n;i++){
for(j= 1;j<=n;j++)
A[i][j]*= -1;
A[i][i]+= 1.0;
}

}

/*:7*/
#line 27 "ad_matrx.w"

/*9:*/
#line 101 "ad_matrx.w"

/*8:*/
#line 98 "ad_matrx.w"

void Transpose_Matrix(int n,double**a)

/*:8*/
#line 102 "ad_matrx.w"

{
int i,j;
double swap;

for(i= 1;i<=n;i++){
for(j= i+1;j<=n;j++){
swap= a[i][j];
a[i][j]= a[j][i];
a[j][i]= swap;
}
}
}

/*:9*/
#line 28 "ad_matrx.w"

/*11:*/
#line 120 "ad_matrx.w"

/*10:*/
#line 117 "ad_matrx.w"

void Diagonal_To_Matrix(int n,double*Diag,double**Mat)

/*:10*/
#line 121 "ad_matrx.w"

{
int i,j;

for(i= 1;i<=n;i++){
for(j= 1;j<=n;j++)
Mat[i][j]= 0.0;
Mat[i][i]= Diag[i];
}
}

/*:11*/
#line 29 "ad_matrx.w"

/*13:*/
#line 144 "ad_matrx.w"

/*12:*/
#line 141 "ad_matrx.w"

void Right_Diagonal_Multiply(int n,double**A,double*B,double**C)

/*:12*/
#line 145 "ad_matrx.w"

{
int i,j;

for(i= 1;i<=n;i++)
for(j= 1;j<=n;j++)
C[i][j]= A[i][j]*B[j];
}

/*:13*/
#line 30 "ad_matrx.w"

/*15:*/
#line 160 "ad_matrx.w"

/*14:*/
#line 157 "ad_matrx.w"

void Left_Diagonal_Multiply(int n,double*A,double**B,double**C)

/*:14*/
#line 161 "ad_matrx.w"

{
int i,j;

for(i= 1;i<=n;i++)
for(j= 1;j<=n;j++)
C[i][j]= A[i]*B[i][j];
}

/*:15*/
#line 31 "ad_matrx.w"

/*22:*/
#line 250 "ad_matrx.w"

/*21:*/
#line 247 "ad_matrx.w"

void Matrix_Multiply(int n,double**A,double**B,double**C)

/*:21*/
#line 251 "ad_matrx.w"

{
/*23:*/
#line 261 "ad_matrx.w"

double*a_ptr,*a_start;
double*b_start,*b_last;
double*c_start,*c_very_last,*c_ptr;
double*D;
double*d_start,*d_last;
register double t,*d_ptr,*b_ptr;
ptrdiff_t row;

/*:23*/
#line 253 "ad_matrx.w"

/*24:*/
#line 270 "ad_matrx.w"

if(n<=0){
AD_error("Non-positive dimension passed to Matrix_Multiply");
}
else if(n==1){
C[1][1]= A[1][1]*B[1][1];
return;
}

/*:24*/
#line 254 "ad_matrx.w"

/*25:*/
#line 283 "ad_matrx.w"

D= dvector(1,n);

/*:25*/
#line 255 "ad_matrx.w"

/*26:*/
#line 294 "ad_matrx.w"

a_start= &A[1][1];
b_last= &B[n][1];
row= &A[2][1]-a_start;
c_very_last= &C[n][n];
d_start= &D[1];
d_last= &D[n];

/*:26*/
#line 256 "ad_matrx.w"

/*29:*/
#line 325 "ad_matrx.w"


for(c_start= &C[1][1];c_start<=c_very_last;c_start+= row){
a_ptr= a_start;
/*27:*/
#line 305 "ad_matrx.w"

d_ptr= d_start;
while(d_ptr<=d_last)
*d_ptr++= 0.0;

/*:27*/
#line 329 "ad_matrx.w"


for(b_start= &B[1][1];b_start<=b_last;b_start+= row){
t= *a_ptr++;
b_ptr= b_start;
d_ptr= d_start;
while(d_ptr<=d_last)
*d_ptr+++= t*(*b_ptr++);
}
/*28:*/
#line 313 "ad_matrx.w"

d_ptr= d_start;
c_ptr= c_start;
while(d_ptr<=d_last)
*c_ptr++= *d_ptr++;

/*:28*/
#line 338 "ad_matrx.w"

a_start+= row;
}

/*:29*/
#line 257 "ad_matrx.w"

/*30:*/
#line 343 "ad_matrx.w"

free_dvector(D,1,n);

/*:30*/
#line 258 "ad_matrx.w"

}

/*:22*/
#line 32 "ad_matrx.w"

/*17:*/
#line 176 "ad_matrx.w"

/*16:*/
#line 173 "ad_matrx.w"

void Matrix_Sum(int n,double**A,double**B,double**C)

/*:16*/
#line 177 "ad_matrx.w"

{
int i,j;

for(i= 1;i<=n;i++)
for(j= 1;j<=n;j++)
C[i][j]= A[i][j]+B[i][j];
}

/*:17*/
#line 33 "ad_matrx.w"

/*43:*/
#line 464 "ad_matrx.w"

/*42:*/
#line 451 "ad_matrx.w"

void Solve(int n,double**A,double*B,int*ipvt)

/*:42*/
#line 465 "ad_matrx.w"

{
int i,k,m;
double t;

/*44:*/
#line 474 "ad_matrx.w"

for(k= 1;k<n;k++){
m= ipvt[k];
t= B[m];
B[m]= B[k];
B[k]= t;
for(i= k+1;i<=n;i++)
B[i]+= A[i][k]*t;
}

/*:44*/
#line 470 "ad_matrx.w"

/*45:*/
#line 484 "ad_matrx.w"

for(k= n;k> 1;k--){
B[k]/= A[k][k];
t= -B[k];
for(i= 1;i<k;i++)
B[i]+= A[i][k]*t;
}

B[1]/= A[1][1];


/*:45*/
#line 471 "ad_matrx.w"

}

/*:43*/
#line 34 "ad_matrx.w"

/*33:*/
#line 368 "ad_matrx.w"

/*32:*/
#line 348 "ad_matrx.w"

void Decomp(int n,double**A,double*condition,int*ipvt)

/*:32*/
#line 369 "ad_matrx.w"

{
double t,anorm;
int i,j,k,m;

/*34:*/
#line 383 "ad_matrx.w"

ipvt[n]= 1;

if(n==1){
if(A[1][1]==0){
AD_error("1 X 1 Matrix is Singular --- i.e. zero");
return;}
}

/*:34*/
#line 374 "ad_matrx.w"

/*35:*/
#line 392 "ad_matrx.w"

anorm= 0.0;
for(j= 1;j<=n;j++){
t= 0.0;
for(i= 1;i<=n;i++)
t+= fabs(A[i][j]);
if(t> anorm)
anorm= t;
}

/*:35*/
#line 375 "ad_matrx.w"

/*36:*/
#line 402 "ad_matrx.w"

for(k= 1;k<n;k++){
/*37:*/
#line 409 "ad_matrx.w"

m= k;
for(i= k+1;i<=n;i++)
if(fabs(A[i][k])> fabs(A[m][k]))
m= i;

ipvt[k]= m;
if(m!=k)
ipvt[n]*= -1;
t= A[m][k];
A[m][k]= A[k][k];
A[k][k]= t;


if(t==0)continue;

/*:37*/
#line 404 "ad_matrx.w"

/*38:*/
#line 425 "ad_matrx.w"

for(i= k+1;i<=n;i++)
A[i][k]/= -t;

/*:38*/
#line 405 "ad_matrx.w"

/*39:*/
#line 429 "ad_matrx.w"

for(j= k+1;j<=n;j++){
t= A[m][j];
A[m][j]= A[k][j];
A[k][j]= t;
if(t==0)continue;
for(i= k+1;i<=n;i++)
A[i][j]+= A[i][k]*t;
}

/*:39*/
#line 406 "ad_matrx.w"

}

/*:36*/
#line 376 "ad_matrx.w"

/*40:*/
#line 439 "ad_matrx.w"


*condition= 1.0;
for(k= 1;k<=n;k++){
if(A[k][k]==0.0){
*condition= 1e32;
return;
}
}

/*:40*/
#line 377 "ad_matrx.w"

}

/*:33*/
#line 35 "ad_matrx.w"

/*47:*/
#line 501 "ad_matrx.w"

/*46:*/
#line 498 "ad_matrx.w"

void Matrix_Inverse(int n,double**A,double**Ainv)

/*:46*/
#line 502 "ad_matrx.w"

{
int*ipvt;
int i,j;
double*work;
double condition;

ipvt= ivector(1,n);
Decomp(n,A,&condition,ipvt);
if(condition==(condition+1)||condition==1e32){
free_ivector(ipvt,1,n);
AD_error("Singular Matrix ... failed in Inverse_Multiply\n");
}
work= dvector(1,n);
for(i= 1;i<=n;i++){
for(j= 1;j<=n;j++)
work[j]= 0.0;
work[i]= 1.0;
Solve(n,A,work,ipvt);
for(j= 1;j<=n;j++)
Ainv[j][i]= work[j];
}

free_dvector(work,1,n);
free_ivector(ipvt,1,n);
}


/*:47*/
#line 36 "ad_matrx.w"

/*49:*/
#line 538 "ad_matrx.w"

/*48:*/
#line 530 "ad_matrx.w"

void Left_Inverse_Multiply(int n,double**D,double**C,double**A)

/*:48*/
#line 539 "ad_matrx.w"

{
int*ipvt;
int i,j;
double*work;
double condition;

Transpose_Matrix(n,D);
ipvt= ivector(1,n);
Decomp(n,D,&condition,ipvt);


if(condition==(condition+1)||condition==1e32){
free_ivector(ipvt,1,n);
AD_error("Singular Matrix ... failed in Left_Inverse_Multiply\n");
}

work= dvector(1,n);
for(i= 1;i<=n;i++){

for(j= 1;j<=n;j++)
work[j]= C[i][j];
Solve(n,D,work,ipvt);
for(j= 1;j<=n;j++)
A[i][j]= work[j];
}

free_dvector(work,1,n);
free_ivector(ipvt,1,n);
}

/*:49*/
#line 37 "ad_matrx.w"

/*51:*/
#line 578 "ad_matrx.w"

/*50:*/
#line 570 "ad_matrx.w"

void Right_Inverse_Multiply(int n,double**D,double**C,double**A)

/*:50*/
#line 579 "ad_matrx.w"

{
int*ipvt;
int i,j;
double*work;
double condition;

ipvt= ivector(1,n);
Decomp(n,D,&condition,ipvt);


if(condition==(condition+1)||condition==1e32){
free_ivector(ipvt,1,n);
AD_error("Singular Matrix ... failed in Right_Inverse_Multiply\n");
}

work= dvector(1,n);
for(i= 1;i<=n;i++){

for(j= 1;j<=n;j++)
work[j]= C[j][i];
Solve(n,D,work,ipvt);
for(j= 1;j<=n;j++)
A[j][i]= work[j];
}

free_dvector(work,1,n);
free_ivector(ipvt,1,n);
}/*:51*/
#line 38 "ad_matrx.w"


/*:1*/
