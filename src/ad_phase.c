/*1:*/
#line 10 "ad_phase.w"

#include <stdlib.h> 
#include <math.h> 
#include "nr_util.h"
#include "ad_globl.h"
#include "ad_phase.h"

/*7:*/
#line 129 "ad_phase.w"

/*6:*/
#line 126 "ad_phase.w"

void Get_Phi(int n,int phase_function,double g,double**h)

/*:6*/
#line 130 "ad_phase.w"

{
/*8:*/
#line 143 "ad_phase.w"

int i,j,k;
double g2M,gk,x;
double*chi;
double**p;

/*:8*/
#line 132 "ad_phase.w"

/*9:*/
#line 149 "ad_phase.w"

if(g!=0&&phase_function!=HENYEY_GREENSTEIN)
AD_error("Only the Henyey-Greenstein phase function has been implemented\n");

if(fabs(g)>=1)
AD_error("Get_Phi was called with a bad g_calc value");

/*:9*/
#line 133 "ad_phase.w"

/*10:*/
#line 156 "ad_phase.w"

for(i= -n;i<=n;i++)
for(j= -n;j<=n;j++)
h[i][j]= 1;


for(i= -n;i<=n;i++){
h[i][0]= 0.0;
h[0][i]= 0.0;
}

/*:10*/
#line 134 "ad_phase.w"

/*11:*/
#line 167 "ad_phase.w"

if(g==0)return;

/*:11*/
#line 135 "ad_phase.w"

/*12:*/
#line 176 "ad_phase.w"

chi= dvector(1,n);
g2M= pow(g,n);
gk= 1.0;
for(k= 1;k<n;k++){
gk*= g;
chi[k]= (2*k+1)*(gk-g2M)/(1-g2M);
}

/*:12*/
#line 136 "ad_phase.w"

/*13:*/
#line 200 "ad_phase.w"

/*14:*/
#line 206 "ad_phase.w"

p= dmatrix(0,n,-n,n);

/*:14*/
#line 201 "ad_phase.w"

/*15:*/
#line 216 "ad_phase.w"

for(j= 1;j<=n;j++){
p[0][j]= 1;
x= angle[j];
p[1][j]= x;
for(k= 1;k<n;k++)
p[k+1][j]= ((2*k+1)*x*p[k][j]-k*p[k-1][j])/(k+1);
}

/*:15*/
#line 202 "ad_phase.w"

/*16:*/
#line 237 "ad_phase.w"

for(j= 1;j<=n;j++)
for(k= 1;k<n;k++){
p[k][-j]= -p[k][j];
k++;
p[k][-j]= p[k][j];
}

/*:16*/
#line 203 "ad_phase.w"


/*:13*/
#line 137 "ad_phase.w"

/*17:*/
#line 257 "ad_phase.w"

for(i= 1;i<=n;i++){
for(j= i;j<=n;j++){
for(k= 1;k<n;k++){
h[i][j]+= chi[k]*p[k][i]*p[k][j];
h[-i][j]+= chi[k]*p[k][-i]*p[k][j];
}
}
}

/*:17*/
#line 138 "ad_phase.w"

/*18:*/
#line 289 "ad_phase.w"

for(i= n;i>=2;i--)
for(j= 1;j<i;j++){
h[-i][j]= h[-j][i];
h[-i][-j]= h[j][i];
}

for(i= 1;i<=n;i++)
h[-i][-i]= h[i][i];

for(i= -n;i<=n;i++)
for(j= i+1;j<=n;j++)
h[j][i]= h[i][j];

/*:18*/
#line 139 "ad_phase.w"

/*19:*/
#line 303 "ad_phase.w"

free_dmatrix(p,0,n,-n,n);
free_dvector(chi,1,n);/*:19*/
#line 140 "ad_phase.w"

}

/*:7*/
#line 17 "ad_phase.w"


/*:1*/
