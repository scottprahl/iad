/*1:*/
#line 10 "ad_radau.w"

#define NSLICES 512
#define EPS 1e-16 \


#line 11 "ad_radau.w"

#include"ad_globl.h"
#include"ad_radau.h"
#include"nr_rtsaf.h"
#include"nr_util.h"
#include"nr_zbrak.h"

static int local_n_size;

/*12:*/
#line 198 "ad_radau.w"

static void Pn_and_Pnm1(int n,double x,double*Pnm1,double*Pn)

/*:12*/
#line 20 "ad_radau.w"
;
/*14:*/
#line 234 "ad_radau.w"

static double Pnd(int n,double x)

/*:14*/
#line 21 "ad_radau.w"
;
/*20:*/
#line 406 "ad_radau.w"

static double phi(double x)

/*:20*/
#line 22 "ad_radau.w"
;
/*16:*/
#line 296 "ad_radau.w"

static void phi_and_phiprime(double x,double*phi,double*phiprime)

/*:16*/
#line 23 "ad_radau.w"
;

/*13:*/
#line 201 "ad_radau.w"

/*12:*/
#line 198 "ad_radau.w"

static void Pn_and_Pnm1(int n,double x,double*Pnm1,double*Pn)

/*:12*/
#line 202 "ad_radau.w"

{
int k;
double Pk,Pkp1;
double Pkm1= 1.0;

*Pnm1= 1.0;
*Pn= 1.0;
if(x>=1.0)
return;

if(x<=-1.0)
x= -1;

Pk= x;

for(k= 1;k<n;k++){
Pkp1= ((2*k+1)*x*Pk-k*Pkm1)/(k+1);
Pkm1= Pk;
Pk= Pkp1;
}

*Pnm1= Pkm1;
*Pn= Pk;
}

/*:13*/
#line 25 "ad_radau.w"

/*15:*/
#line 237 "ad_radau.w"

/*14:*/
#line 234 "ad_radau.w"

static double Pnd(int n,double x)

/*:14*/
#line 238 "ad_radau.w"

{
double p,pminus,pplus;
int i;

if(x> 1.0){
x= 1;
}
else if(x<-1.0){
x= -1;
}

pminus= 0;
p= 1;

if(n<=0)
return pminus;

for(i= 1;i<n;i++){
pplus= ((2*i+1)*x*p-(i+1)*pminus)/i;
pminus= p;
p= pplus;
}
return p;
}

/*:15*/
#line 26 "ad_radau.w"

/*21:*/
#line 409 "ad_radau.w"

/*20:*/
#line 406 "ad_radau.w"

static double phi(double x)

/*:20*/
#line 410 "ad_radau.w"

{
double Pn,Pnm1;

if(x<=-1.0){
if(local_n_size%2!=1)
return(-local_n_size);
else
return(local_n_size);
}

Pn_and_Pnm1(local_n_size,x,&Pnm1,&Pn);
return((Pn+Pnm1)/(1+x));
}

/*:21*/
#line 27 "ad_radau.w"

/*17:*/
#line 299 "ad_radau.w"

/*16:*/
#line 296 "ad_radau.w"

static void phi_and_phiprime(double x,double*phi,double*phiprime)

/*:16*/
#line 300 "ad_radau.w"

{
double Pn,Pnm1;
int n;

n= local_n_size;
if(x>=1.0){
/*18:*/
#line 344 "ad_radau.w"

{
*phi= 1;
*phiprime= (n*n-1)/2;
}

/*:18*/
#line 307 "ad_radau.w"

}
else if(x<=-1.0){
/*19:*/
#line 374 "ad_radau.w"

*phi= n;
*phiprime= -n*(1-n*n)/4;
if(n%2!=1){
*phi*= -1;
*phiprime*= -1;
}

/*:19*/
#line 310 "ad_radau.w"

}
else{
Pn_and_Pnm1(n,x,&Pnm1,&Pn);
*phi= (Pn+Pnm1)/(1+x);
*phiprime= ((n*x-1+x+n)*Pnm1+(-n*x+x-n-1)*Pn)/(1+x)/(1+x)/(1-x);
}
}

/*:17*/
#line 28 "ad_radau.w"

/*5:*/
#line 96 "ad_radau.w"

/*4:*/
#line 93 "ad_radau.w"

void Radau(double x1,double x2,double*x,double*w,int n)

/*:4*/
#line 97 "ad_radau.w"

{

x[n]= -1.0;
w[n]= 2.0/(n*n);

switch(n){
case 2:/*23:*/
#line 429 "ad_radau.w"

x[1]= 0.3333333333333334;
w[1]= 1.5000000000000000;
break;

/*:23*/
#line 104 "ad_radau.w"

case 4:/*24:*/
#line 434 "ad_radau.w"

x[3]= -0.5753189235216942;
x[2]= 0.1810662711185306;
x[1]= 0.8228240809745921;

w[3]= 0.6576886399601182;
w[2]= 0.7763869376863437;
w[1]= 0.4409244223535367;
break;

/*:24*/
#line 105 "ad_radau.w"

case 8:/*25:*/
#line 444 "ad_radau.w"

x[7]= -0.8874748789261557;
x[6]= -0.6395186165262152;
x[5]= -0.2947505657736607;
x[4]= 0.0943072526611108;
x[3]= 0.4684203544308211;
x[2]= 0.7706418936781916;
x[1]= 0.9550412271225750;

w[7]= 0.1853581548029793;
w[6]= 0.3041306206467856;
w[5]= 0.3765175453891186;
w[4]= 0.3915721674524935;
w[3]= 0.3470147956345014;
w[2]= 0.2496479013298649;
w[1]= 0.1145088147442572;
break;

/*:25*/
#line 106 "ad_radau.w"

case 16:/*26:*/
#line 462 "ad_radau.w"

x[15]= -0.9714610905263484;
x[14]= -0.9054008198116666;
x[13]= -0.8045734013587561;
x[12]= -0.6728619212112202;
x[11]= -0.5153294780626855;
x[10]= -0.3380303900599197;
x[9]= -0.1477783218133717;
x[8]= 0.0481153830735303;
x[7]= 0.2421226227060438;
x[6]= 0.4267878274849459;
x[5]= 0.5950144898997919;
x[4]= 0.7403379488928179;
x[3]= 0.8571740937696823;
x[2]= 0.9410354027041150;
x[1]= 0.9887186220549766;

w[15]= 0.0477022269476863;
w[14]= 0.0839852814449645;
w[13]= 0.1170203531038591;
w[12]= 0.1455555452202026;
w[11]= 0.1684963978499219;
w[10]= 0.1849617814886653;
w[10]= 0.1849617814886653;
w[9]= 0.1943190897115679;
w[8]= 0.1962087882390318;
w[7]= 0.1905582942553547;
w[6]= 0.1775847927527395;
w[5]= 0.1577869218042020;
w[4]= 0.1319256999330681;
w[3]= 0.1009956796217840;
w[2]= 0.0661895086101364;
w[1]= 0.0288971390168143;
break;/*:26*/
#line 107 "ad_radau.w"

default:/*7:*/
#line 145 "ad_radau.w"

{
int i,nb,ndiv;
double z;
double*xb1,*xb2;

/*8:*/
#line 158 "ad_radau.w"

xb1= dvector(1,NSLICES);
xb2= dvector(1,NSLICES);

/*:8*/
#line 151 "ad_radau.w"

/*9:*/
#line 163 "ad_radau.w"

local_n_size= n;

if(2*n> NSLICES)
ndiv= NSLICES;
else
ndiv= 2*n;

do{
nb= n-1;
zbrak(phi,-1.0,1.0,ndiv,xb1,xb2,&nb);
ndiv*= 2;
}
while(nb<n-1&&ndiv<=NSLICES);

if(nb<n-1)AD_error("Cannot find enough roots for Radau quadrature");

/*:9*/
#line 152 "ad_radau.w"

/*10:*/
#line 183 "ad_radau.w"

for(i= 1;i<n;i++){
double tmp;
z= rtsafe(phi_and_phiprime,xb1[i],xb2[i],EPS);
x[n-i]= z;
tmp= Pnd(n-1,z);
w[n-i]= 1/((1-z)*tmp*tmp);
}

/*:10*/
#line 153 "ad_radau.w"

/*11:*/
#line 192 "ad_radau.w"

free_dvector(xb1,1,NSLICES);
free_dvector(xb2,1,NSLICES);

/*:11*/
#line 154 "ad_radau.w"

break;
}

/*:7*/
#line 108 "ad_radau.w"

}
/*6:*/
#line 126 "ad_radau.w"

{
double xm,xl;
int i;

xm= (x2+x1)*0.5;
xl= (x2-x1)*0.5;

for(i= 1;i<=n;i++){
x[i]= xm-xl*x[i];
w[i]= xl*w[i];
}

}


/*:6*/
#line 110 "ad_radau.w"

}

/*:5*/
#line 29 "ad_radau.w"


/*:1*/
