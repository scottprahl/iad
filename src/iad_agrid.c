/*2:*/
#line 25 "iad_agrid.w"


#include <math.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include "ad_globl.h"
#include "iad_type.h"
#include "iad_util.h"
#include "iad_calc.h"
#include "iad_agrid.h"

/*4:*/
#line 63 "iad_agrid.w"

#define AGRID_TOL           0.05    
#define AGRID_MIN_DEPTH     2       
#define AGRID_MAX_DEPTH     4       
#define AGRID_INIT_CAP      512     
#define AGRID_CLOSURE_N     4       
#define AGRID_MIN_LOG_B     (-8.0)
#define AGRID_MAX_LOG_B_AB  8.0     
#define AGRID_MAX_LOG_B_BG  10.0    

/*:4*/
#line 37 "iad_agrid.w"

/*5:*/
#line 81 "iad_agrid.w"


typedef struct{
double a,b,g;
double ur1,ut1,uru,utu;
}agrid_entry_t;

static agrid_entry_t*AGrid_Cache= NULL;
static int AGrid_N= 0;
static int AGrid_Cap= 0;
static int AGrid_Search= -1;

/*:5*/
#line 38 "iad_agrid.w"

/*6:*/
#line 99 "iad_agrid.w"


static int AGrid_Initialized= 0;
static double AGrid_slab_n= 0;
static double AGrid_slab_n_top= 0;
static double AGrid_slab_n_bottom= 0;
static double AGrid_slab_b_top= 0;
static double AGrid_slab_b_bottom= 0;
static double AGrid_slab_cos_angle= 0;
static int AGrid_num_spheres= -1;
static int AGrid_num_measures= -1;
static double AGrid_m_u= 0;

static double AGrid_fixed_g= 0;
static double AGrid_fixed_b= 0;
static double AGrid_fixed_a= 0;

/*:6*/
#line 39 "iad_agrid.w"


/*7:*/
#line 122 "iad_agrid.w"

static double agrid_nonlinear_a(double u)
{

return 1.0-(1.0-u)*(1.0-u)*(1.0+2.0*u);
}

/*:7*/
#line 41 "iad_agrid.w"

/*8:*/
#line 129 "iad_agrid.w"

static double agrid_nonlinear_g(double v)
{

double vv= (1.0-v)*(1.0-v)*(1.0+2.0*v);
return(1.0-2.0*vv)*MAX_ABS_G;
}

/*:8*/
#line 42 "iad_agrid.w"

/*9:*/
#line 141 "iad_agrid.w"

static void agrid_make_abg(double u,double v,int search,
double*a,double*b,double*g)
{
switch(search){
case FIND_AB:

*a= agrid_nonlinear_a(u);
*b= exp(v);
*g= AGrid_fixed_g;
break;
case FIND_AG:

*a= agrid_nonlinear_a(u);
*b= AGrid_fixed_b;
*g= agrid_nonlinear_g(v);
break;
case FIND_BG:

*a= AGrid_fixed_a;
*b= exp(u);
*g= agrid_nonlinear_g(v);
break;
default:
*a= 0.5;*b= 1.0;*g= 0.0;
break;
}


if(*b<1e-8)*b= 1e-8;
if(*g<-MAX_ABS_G)*g= -MAX_ABS_G;
if(*g> MAX_ABS_G)*g= MAX_ABS_G;
}

/*:9*/
#line 43 "iad_agrid.w"

/*10:*/
#line 177 "iad_agrid.w"

static void agrid_ensure_capacity(void)
{
if(AGrid_N<AGrid_Cap)return;
AGrid_Cap= (AGrid_Cap==0)?AGRID_INIT_CAP:AGrid_Cap*2;
AGrid_Cache= (agrid_entry_t*)realloc(AGrid_Cache,
(size_t)AGrid_Cap*sizeof(agrid_entry_t));
if(AGrid_Cache==NULL){
fprintf(stderr,"AGrid: out of memory\n");
exit(EXIT_FAILURE);
}
}

/*:10*/
#line 44 "iad_agrid.w"

/*11:*/
#line 195 "iad_agrid.w"

static int agrid_find_entry(double a,double b,double g)
{
int i;
for(i= 0;i<AGrid_N;i++){
if(AGrid_Cache[i].a==a&&
AGrid_Cache[i].b==b&&
AGrid_Cache[i].g==g)
return i;
}
return-1;
}

/*:11*/
#line 45 "iad_agrid.w"

/*12:*/
#line 212 "iad_agrid.w"

static int agrid_add_or_get(double a,double b,double g)
{
int idx= agrid_find_entry(a,b,g);
if(idx>=0)return idx;

agrid_ensure_capacity();
idx= AGrid_N;
AGrid_Cache[idx].a= a;
AGrid_Cache[idx].b= b;
AGrid_Cache[idx].g= g;
abg_eval(a,b,g,
&AGrid_Cache[idx].ur1,
&AGrid_Cache[idx].ut1,
&AGrid_Cache[idx].uru,
&AGrid_Cache[idx].utu);
AGrid_N++;
return idx;
}

/*:12*/
#line 46 "iad_agrid.w"

/*13:*/
#line 244 "iad_agrid.w"

static double agrid_raw_interp_error(int c00,int c10,int c01,int c11,int cc)
{
double mr00,mt00,mr10,mt10,mr01,mt01,mr11,mt11,mrcc,mtcc;
double interp_mr,interp_mt;
agrid_entry_t*e;

e= &AGrid_Cache[c00];
abg_sphere_mr_mt(e->a,e->b,e->g,e->ur1,e->ut1,e->uru,e->utu,&mr00,&mt00);
e= &AGrid_Cache[c10];
abg_sphere_mr_mt(e->a,e->b,e->g,e->ur1,e->ut1,e->uru,e->utu,&mr10,&mt10);
e= &AGrid_Cache[c01];
abg_sphere_mr_mt(e->a,e->b,e->g,e->ur1,e->ut1,e->uru,e->utu,&mr01,&mt01);
e= &AGrid_Cache[c11];
abg_sphere_mr_mt(e->a,e->b,e->g,e->ur1,e->ut1,e->uru,e->utu,&mr11,&mt11);
e= &AGrid_Cache[cc];
abg_sphere_mr_mt(e->a,e->b,e->g,e->ur1,e->ut1,e->uru,e->utu,&mrcc,&mtcc);

interp_mr= 0.25*(mr00+mr10+mr01+mr11);
interp_mt= 0.25*(mt00+mt10+mt01+mt11);
return fabs(mrcc-interp_mr)+fabs(mtcc-interp_mt);
}

/*:13*/
#line 47 "iad_agrid.w"

/*14:*/
#line 269 "iad_agrid.w"

static void agrid_subdivide(double u0,double u1,double v0,double v1,
int depth,int search)
{
double um= 0.5*(u0+u1);
double vm= 0.5*(v0+v1);
double a00,b00,g00;
double a10,b10,g10;
double a01,b01,g01;
double a11,b11,g11;
double amm,bmm,gmm;
int c00,c10,c01,c11,cc;
double err;
int need_split;


agrid_make_abg(u0,v0,search,&a00,&b00,&g00);
agrid_make_abg(u1,v0,search,&a10,&b10,&g10);
agrid_make_abg(u0,v1,search,&a01,&b01,&g01);
agrid_make_abg(u1,v1,search,&a11,&b11,&g11);
agrid_make_abg(um,vm,search,&amm,&bmm,&gmm);

c00= agrid_add_or_get(a00,b00,g00);
c10= agrid_add_or_get(a10,b10,g10);
c01= agrid_add_or_get(a01,b01,g01);
c11= agrid_add_or_get(a11,b11,g11);
cc= agrid_add_or_get(amm,bmm,gmm);

err= agrid_raw_interp_error(c00,c10,c01,c11,cc);
need_split= (depth<AGRID_MIN_DEPTH)||
(err> AGRID_TOL&&depth<AGRID_MAX_DEPTH);

if(need_split){
agrid_subdivide(u0,um,v0,vm,depth+1,search);
agrid_subdivide(um,u1,v0,vm,depth+1,search);
agrid_subdivide(u0,um,vm,v1,depth+1,search);
agrid_subdivide(um,u1,vm,v1,depth+1,search);
}

}

/*:14*/
#line 48 "iad_agrid.w"

/*15:*/
#line 312 "iad_agrid.w"

static void agrid_save_context(struct measure_type m,struct invert_type r)
{
AGrid_slab_n= m.slab_index;
AGrid_slab_n_top= m.slab_top_slide_index;
AGrid_slab_n_bottom= m.slab_bottom_slide_index;
AGrid_slab_b_top= r.slab.b_top_slide;
AGrid_slab_b_bottom= r.slab.b_bottom_slide;
AGrid_slab_cos_angle= m.slab_cos_angle;
AGrid_num_spheres= m.num_spheres;
AGrid_num_measures= m.num_measures;
AGrid_m_u= m.m_u;



AGrid_fixed_g= r.slab.g;
AGrid_fixed_b= r.slab.b;
AGrid_fixed_a= (r.default_a==UNINITIALIZED)?0.0:r.default_a;

AGrid_Search= r.search;
AGrid_Initialized= 1;
}

/*:15*/
#line 49 "iad_agrid.w"

/*17:*/
#line 344 "iad_agrid.w"

/*16:*/
#line 341 "iad_agrid.w"

int AGrid_Valid(struct measure_type m,struct invert_type r)

/*:16*/
#line 345 "iad_agrid.w"

{
if(!AGrid_Initialized||AGrid_N==0)return 0;
if(AGrid_Search!=r.search)return 0;

if(m.slab_index!=AGrid_slab_n)return 0;
if(m.slab_cos_angle!=AGrid_slab_cos_angle)return 0;
if(m.slab_top_slide_index!=AGrid_slab_n_top)return 0;
if(m.slab_bottom_slide_index!=AGrid_slab_n_bottom)return 0;
if(r.slab.b_top_slide!=AGrid_slab_b_top)return 0;
if(r.slab.b_bottom_slide!=AGrid_slab_b_bottom)return 0;

if(m.num_spheres!=AGrid_num_spheres)return 0;
if(m.num_measures!=AGrid_num_measures)return 0;
if(m.num_measures==3&&m.m_u!=AGrid_m_u)return 0;






if(r.search==FIND_AB&&r.slab.g!=AGrid_fixed_g)return 0;
if(r.search==FIND_AG&&r.slab.b!=AGrid_fixed_b)return 0;
if(r.search==FIND_BG){
double fa= (r.default_a==UNINITIALIZED)?0.0:r.default_a;
if(fa!=AGrid_fixed_a)return 0;
}

return 1;
}

/*:17*/
#line 50 "iad_agrid.w"

/*19:*/
#line 382 "iad_agrid.w"

/*18:*/
#line 379 "iad_agrid.w"

void AGrid_Build(struct measure_type m,struct invert_type r)

/*:18*/
#line 383 "iad_agrid.w"

{
double u0,u1,v0,v1;
int search= r.search;


AGrid_N= 0;


if(search==FIND_AB)AGrid_fixed_g= r.slab.g;
if(search==FIND_AG)AGrid_fixed_b= r.slab.b;
if(search==FIND_BG)
AGrid_fixed_a= (r.default_a==UNINITIALIZED)?0.0:r.default_a;


Set_Calc_State(m,r);

if(Debug(DEBUG_GRID)){
switch(search){
case FIND_AB:
fprintf(stderr,"AGRID: Filling AB grid (g=%.5f)\n",AGrid_fixed_g);
break;
case FIND_AG:
fprintf(stderr,"AGRID: Filling AG grid (b=%.5f)\n",AGrid_fixed_b);
break;
case FIND_BG:
fprintf(stderr,"AGRID: Filling BG grid (a=%.5f)\n",AGrid_fixed_a);
break;
default:
fprintf(stderr,"AGRID: Filling grid for search=%d\n",search);
break;
}
}


switch(search){
case FIND_AB:
u0= 0.0;u1= 1.0;
v0= AGRID_MIN_LOG_B;v1= AGRID_MAX_LOG_B_AB;
break;
case FIND_AG:
u0= 0.0;u1= 1.0;
v0= 0.0;v1= 1.0;
break;
case FIND_BG:
u0= AGRID_MIN_LOG_B;u1= AGRID_MAX_LOG_B_BG;
v0= 0.0;v1= 1.0;
break;
default:
return;
}

agrid_subdivide(u0,u1,v0,v1,0,search);

agrid_save_context(m,r);

if(Debug(DEBUG_GRID))
fprintf(stderr,"AGRID: Built %d entries for search=%d\n",AGrid_N,search);
if(Debug(DEBUG_LOST_LIGHT))
fprintf(stderr,"GRID: %d AD evaluations to fill the grid\n",AGrid_N);
}

/*:19*/
#line 51 "iad_agrid.w"

/*21:*/
#line 459 "iad_agrid.w"

/*20:*/
#line 455 "iad_agrid.w"

int AGrid_Fill_Guesses(double m_r,double m_t,
guess_type*guesses,int max_n)

/*:20*/
#line 460 "iad_agrid.w"

{
int i,j,n_out;
double*all_dist;
double best_dist[AGRID_CLOSURE_N];
int best_idx[AGRID_CLOSURE_N];
int n_best;

if(AGrid_N==0)return 0;


all_dist= (double*)malloc((size_t)AGrid_N*sizeof(double));
if(all_dist==NULL){
fprintf(stderr,"AGrid: out of memory in AGrid_Fill_Guesses\n");
exit(EXIT_FAILURE);
}

n_best= 0;
for(i= 0;i<AGrid_N;i++){
double d= abg_stored_distance(AGrid_Cache[i].a,
AGrid_Cache[i].b,
AGrid_Cache[i].g,
AGrid_Cache[i].ur1,
AGrid_Cache[i].ut1,
AGrid_Cache[i].uru,
AGrid_Cache[i].utu);
all_dist[i]= d;
if(Debug(DEBUG_GRID)&&d<0.15){
fprintf(stderr,"AGRID_ALL: a=%8.5f b=%9.5f g=%7.4f ur1=%8.5f ut1=%8.5f d=%9.6f\n",
AGrid_Cache[i].a,AGrid_Cache[i].b,AGrid_Cache[i].g,
AGrid_Cache[i].ur1,AGrid_Cache[i].ut1,d);
}
if(n_best<AGRID_CLOSURE_N){
best_dist[n_best]= d;
best_idx[n_best]= i;
n_best++;
}else{

int worst= 0;
for(j= 1;j<n_best;j++)
if(best_dist[j]> best_dist[worst])worst= j;
if(d<best_dist[worst]){
best_dist[worst]= d;
best_idx[worst]= i;
}
}
}










{
double axis0[AGRID_CLOSURE_N*2];
double axis1[AGRID_CLOSURE_N*2];
int n0= 0,n1= 0;
int k,l,found;
int use_geo0,use_geo1;


use_geo0= (AGrid_Search==FIND_BG);

use_geo1= (AGrid_Search==FIND_AB);



while(n0<AGRID_CLOSURE_N){
double best_d= 1e30;
int best_i= -1;
double vc;
for(i= 0;i<AGrid_N;i++){
switch(AGrid_Search){
case FIND_AB:vc= AGrid_Cache[i].a;break;
case FIND_AG:vc= AGrid_Cache[i].a;break;
case FIND_BG:vc= AGrid_Cache[i].b;break;
default:vc= AGrid_Cache[i].a;break;
}
found= 0;
for(l= 0;l<n0;l++)if(axis0[l]==vc){found= 1;break;}
if(!found&&all_dist[i]<best_d){
best_d= all_dist[i];
best_i= i;
}
}
if(best_i<0)break;
switch(AGrid_Search){
case FIND_AB:axis0[n0++]= AGrid_Cache[best_i].a;break;
case FIND_AG:axis0[n0++]= AGrid_Cache[best_i].a;break;
case FIND_BG:axis0[n0++]= AGrid_Cache[best_i].b;break;
default:axis0[n0++]= AGrid_Cache[best_i].a;break;
}
}


while(n1<AGRID_CLOSURE_N){
double best_d= 1e30;
int best_i= -1;
double vc;
for(i= 0;i<AGrid_N;i++){
switch(AGrid_Search){
case FIND_AB:vc= AGrid_Cache[i].b;break;
case FIND_AG:vc= AGrid_Cache[i].g;break;
case FIND_BG:vc= AGrid_Cache[i].g;break;
default:vc= AGrid_Cache[i].b;break;
}
found= 0;
for(l= 0;l<n1;l++)if(axis1[l]==vc){found= 1;break;}
if(!found&&all_dist[i]<best_d){
best_d= all_dist[i];
best_i= i;
}
}
if(best_i<0)break;
switch(AGrid_Search){
case FIND_AB:axis1[n1++]= AGrid_Cache[best_i].b;break;
case FIND_AG:axis1[n1++]= AGrid_Cache[best_i].g;break;
case FIND_BG:axis1[n1++]= AGrid_Cache[best_i].g;break;
default:axis1[n1++]= AGrid_Cache[best_i].b;break;
}
}





if(n0>=2){
double sorted0[AGRID_CLOSURE_N];
int n_sorted= n0;
for(k= 0;k<n0;k++)sorted0[k]= axis0[k];
for(k= 1;k<n_sorted;k++){
double tmp= sorted0[k];
l= k-1;
while(l>=0&&sorted0[l]> tmp){sorted0[l+1]= sorted0[l];l--;}
sorted0[l+1]= tmp;
}
for(k= 0;k<n_sorted-1&&n0<AGRID_CLOSURE_N*2;k++){
double lo= sorted0[k],hi= sorted0[k+1];
double mid= (use_geo0&&lo> 0&&hi> 0)?sqrt(lo*hi):0.5*(lo+hi);
found= 0;
for(l= 0;l<n0;l++)if(axis0[l]==mid){found= 1;break;}
if(!found)axis0[n0++]= mid;
}
}


if(n1>=2){
double sorted1[AGRID_CLOSURE_N];
int n_sorted= n1;
for(k= 0;k<n1;k++)sorted1[k]= axis1[k];
for(k= 1;k<n_sorted;k++){
double tmp= sorted1[k];
l= k-1;
while(l>=0&&sorted1[l]> tmp){sorted1[l+1]= sorted1[l];l--;}
sorted1[l+1]= tmp;
}
for(k= 0;k<n_sorted-1&&n1<AGRID_CLOSURE_N*2;k++){
double lo= sorted1[k],hi= sorted1[k+1];
double mid= (use_geo1&&lo> 0&&hi> 0)?sqrt(lo*hi):0.5*(lo+hi);
found= 0;
for(l= 0;l<n1;l++)if(axis1[l]==mid){found= 1;break;}
if(!found)axis1[n1++]= mid;
}
}


for(k= 0;k<n0;k++){
for(l= 0;l<n1;l++){
double a,b,g;
double d;
int idx;

switch(AGrid_Search){
case FIND_AB:a= axis0[k];b= axis1[l];g= AGrid_fixed_g;break;
case FIND_AG:a= axis0[k];b= AGrid_fixed_b;g= axis1[l];break;
case FIND_BG:a= AGrid_fixed_a;b= axis0[k];g= axis1[l];break;
default:a= axis0[k];b= axis1[l];g= AGrid_fixed_g;break;
}

if(b<1e-8)b= 1e-8;
if(g<-MAX_ABS_G)g= -MAX_ABS_G;
if(g> MAX_ABS_G)g= MAX_ABS_G;

idx= agrid_add_or_get(a,b,g);
d= abg_stored_distance(AGrid_Cache[idx].a,
AGrid_Cache[idx].b,
AGrid_Cache[idx].g,
AGrid_Cache[idx].ur1,
AGrid_Cache[idx].ut1,
AGrid_Cache[idx].uru,
AGrid_Cache[idx].utu);


if(n_best<AGRID_CLOSURE_N){
best_dist[n_best]= d;
best_idx[n_best]= idx;
n_best++;
}else{
int worst= 0;
for(j= 1;j<n_best;j++)
if(best_dist[j]> best_dist[worst])worst= j;
if(d<best_dist[worst]){
best_dist[worst]= d;
best_idx[worst]= idx;
}
}
}
}
}

free(all_dist);


for(i= 1;i<n_best;i++){
int ti= best_idx[i];
double td= best_dist[i];
j= i-1;
while(j>=0&&best_dist[j]> td){
best_dist[j+1]= best_dist[j];
best_idx[j+1]= best_idx[j];
j--;
}
best_dist[j+1]= td;
best_idx[j+1]= ti;
}


n_out= (n_best<max_n)?n_best:max_n;
for(i= 0;i<n_out;i++){
agrid_entry_t*e= &AGrid_Cache[best_idx[i]];
guesses[i].a= e->a;
guesses[i].b= e->b;
guesses[i].g= e->g;
guesses[i].distance= best_dist[i];
guesses[i].ur1_lost= 0.0;
guesses[i].ut1_lost= 0.0;
guesses[i].uru_lost= 0.0;
guesses[i].utu_lost= 0.0;
}

if(Debug(DEBUG_BEST_GUESS)){
fprintf(stderr,"BEST: AGRID GUESSES\n");
fprintf(stderr,"BEST:  k      albedo          b          g   distance\n");
for(i= 0;i<n_out&&i<7;i++){
fprintf(stderr,"BEST:%3d  ",i);
fprintf(stderr,"%10.5f ",guesses[i].a);
fprintf(stderr,"%10.5f ",guesses[i].b);
fprintf(stderr,"%10.5f ",guesses[i].g);
fprintf(stderr,"%10.5f\n",guesses[i].distance);
}
}

return n_out;
}/*:21*/
#line 52 "iad_agrid.w"


/*:2*/
