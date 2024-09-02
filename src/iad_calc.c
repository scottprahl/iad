/*1:*/
#line 17 "iad_calc.w"


#include <math.h> 
#include <string.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include "nr_util.h"
#include "nr_zbrent.h"
#include "ad_globl.h"
#include "ad_frsnl.h"
#include "ad_prime.h"
#include "iad_type.h"
#include "iad_util.h"
#include "iad_calc.h"
#include "mc_lost.h"

#define ABIT 1e-6
#define A_COLUMN 1
#define B_COLUMN 2
#define G_COLUMN 3
#define URU_COLUMN 4
#define UTU_COLUMN 5
#define UR1_COLUMN 6
#define UT1_COLUMN 7
#define URU_LOST_COLUMN 8
#define UTU_LOST_COLUMN 9
#define UR1_LOST_COLUMN 10
#define UT1_LOST_COLUMN 11

#define REFLECTION_SPHERE 1
#define TRANSMISSION_SPHERE 0
#define T_TRUST_FACTOR 1

#define SWAP(a,b) {double swap= (a);(a)= (b);(b)= swap;}

static int CALCULATING_GRID= 0;
static struct measure_type MM;
static struct invert_type RR;
static struct measure_type MGRID;
static struct invert_type RGRID;
static double**The_Grid= NULL;
static double GG_a;
static double GG_b;
static double GG_g;
static double GG_bs;
static double GG_ba;
static boolean_type The_Grid_Initialized= FALSE;
static boolean_type The_Grid_Search= -1;

/*17:*/
#line 373 "iad_calc.w"

/*16:*/
#line 370 "iad_calc.w"

void Set_Calc_State(struct measure_type m,struct invert_type r)

/*:16*/
#line 374 "iad_calc.w"

{
memcpy(&MM,&m,sizeof(struct measure_type));
memcpy(&RR,&r,sizeof(struct invert_type));
}

/*:17*/
#line 66 "iad_calc.w"

/*19:*/
#line 386 "iad_calc.w"

/*18:*/
#line 383 "iad_calc.w"

void Get_Calc_State(struct measure_type*m,struct invert_type*r)

/*:18*/
#line 387 "iad_calc.w"

{
memcpy(m,&MM,sizeof(struct measure_type));
memcpy(r,&RR,sizeof(struct invert_type));
}

/*:19*/
#line 67 "iad_calc.w"

/*21:*/
#line 399 "iad_calc.w"

/*20:*/
#line 396 "iad_calc.w"

boolean_type Same_Calc_State(struct measure_type m,struct invert_type r)

/*:20*/
#line 400 "iad_calc.w"

{
if(The_Grid==NULL)return FALSE;
if(!The_Grid_Initialized)return FALSE;

if(r.search!=RR.search)return FALSE;
if(r.method.quad_pts!=RR.method.quad_pts)return FALSE;
if(r.slab.a!=RR.slab.a)return FALSE;
if(r.slab.b!=RR.slab.b)return FALSE;
if(r.slab.g!=RR.slab.g)return FALSE;
if(r.slab.phase_function!=RR.slab.phase_function)return FALSE;
if(r.slab.n_slab!=RR.slab.n_slab)return FALSE;
if(r.slab.n_top_slide!=RR.slab.n_top_slide)return FALSE;
if(r.slab.n_bottom_slide!=RR.slab.n_bottom_slide)return FALSE;
if(r.slab.b_top_slide!=RR.slab.b_top_slide)return FALSE;
if(r.slab.b_bottom_slide!=RR.slab.b_bottom_slide)return FALSE;
if(r.slab.cos_angle!=RR.slab.cos_angle)return FALSE;
if((m.num_measures==3)&&(m.m_u!=MGRID.m_u))return(FALSE);
return TRUE;
}

/*:21*/
#line 68 "iad_calc.w"


/*39:*/
#line 762 "iad_calc.w"

void Fill_AB_Grid(struct measure_type m,struct invert_type r)

/*:39*/
#line 70 "iad_calc.w"
;
/*45:*/
#line 848 "iad_calc.w"

void Fill_AG_Grid(struct measure_type m,struct invert_type r)

/*:45*/
#line 71 "iad_calc.w"
;

/*37:*/
#line 657 "iad_calc.w"

/*36:*/
#line 653 "iad_calc.w"

void RT_Flip(int flip,int n,struct AD_slab_type*slab,double*UR1,double*UT1,
double*URU,double*UTU)

/*:36*/
#line 658 "iad_calc.w"

{
double correct_UR1,correct_URU;
RT(n,slab,UR1,UT1,URU,UTU);
if(flip){
correct_UR1= *UR1;
correct_URU= *URU;
SWAP(slab->n_top_slide,slab->n_bottom_slide)
SWAP(slab->b_top_slide,slab->b_bottom_slide)
RT(n,slab,UR1,UT1,URU,UTU);
SWAP(slab->n_top_slide,slab->n_bottom_slide)
SWAP(slab->b_top_slide,slab->b_bottom_slide)
*UR1= correct_UR1;
*URU= correct_URU;
}
}

/*:37*/
#line 73 "iad_calc.w"

/*23:*/
#line 424 "iad_calc.w"

/*22:*/
#line 421 "iad_calc.w"

void Allocate_Grid(search_type s)

/*:22*/
#line 425 "iad_calc.w"

{
(void)s;
The_Grid= dmatrix(0,GRID_SIZE*GRID_SIZE,1,11);
if(The_Grid==NULL)AD_error("unable to allocate the grid matrix");
The_Grid_Initialized= FALSE;
}

/*:23*/
#line 74 "iad_calc.w"

/*27:*/
#line 476 "iad_calc.w"

/*26:*/
#line 473 "iad_calc.w"

boolean_type Valid_Grid(struct measure_type m,struct invert_type r)

/*:26*/
#line 477 "iad_calc.w"

{
int s= r.search;
/*28:*/
#line 487 "iad_calc.w"

if(The_Grid==NULL){
if(Debug(DEBUG_GRID))
fprintf(stderr,"GRID: Fill because NULL\n");
return(FALSE);
}
if(!The_Grid_Initialized){
if(Debug(DEBUG_GRID))
fprintf(stderr,"GRID: Fill because not initialized\n");
return(FALSE);
}

/*:28*//*29:*/
#line 501 "iad_calc.w"

if(The_Grid_Search!=s){
if(Debug(DEBUG_GRID))
fprintf(stderr,"GRID: Fill because search type changed\n");
return(FALSE);
}

/*:29*//*30:*/
#line 510 "iad_calc.w"


if((m.num_measures==3)&&(m.m_u!=MGRID.m_u)){
if(Debug(DEBUG_GRID))
fprintf(stderr,"GRID: Fill because unscattered light changed\n");
return(FALSE);
}

/*:30*//*31:*/
#line 520 "iad_calc.w"

if(m.slab_index!=MGRID.slab_index){
if(Debug(DEBUG_GRID))
fprintf(stderr,"GRID: Fill because slab refractive index changed\n");
return(FALSE);
}
if(m.slab_cos_angle!=MGRID.slab_cos_angle){
if(Debug(DEBUG_GRID))
fprintf(stderr,"GRID: Fill because light angle changed\n");
return(FALSE);
}

if(m.slab_top_slide_index!=MGRID.slab_top_slide_index){
if(Debug(DEBUG_GRID))
fprintf(stderr,"GRID: Fill because top slide index changed\n");
return(FALSE);
}

if(m.slab_bottom_slide_index!=MGRID.slab_bottom_slide_index){
if(Debug(DEBUG_GRID))
fprintf(stderr,"GRID: Fill because bottom slide index changed\n");
return(FALSE);
}

if(s==FIND_AB&&r.slab.g!=RGRID.slab.g){
if(Debug(DEBUG_GRID))
fprintf(stderr,"GRID: Fill because anisotropy changed\n");
return(FALSE);
}

if(s==FIND_AG&&r.slab.b!=RGRID.slab.b){
if(Debug(DEBUG_GRID))
fprintf(stderr,"GRID: Fill because optical depth changed\n");
return(FALSE);
}

if(s==FIND_BsG&&r.default_ba!=RGRID.default_ba){
if(Debug(DEBUG_GRID))
fprintf(stderr,"GRID: Fill because mu_a changed\n");
return(FALSE);
}

if(s==FIND_BaG&&r.default_bs!=RGRID.default_bs){
if(Debug(DEBUG_GRID))
fprintf(stderr,"GRID: Fill because mu_s changed\n");
return(FALSE);
}


/*:31*/
#line 480 "iad_calc.w"


return(TRUE);
}

/*:27*/
#line 75 "iad_calc.w"

/*38:*/
#line 679 "iad_calc.w"

static void fill_grid_entry(int i,int j)
{
double ur1,ut1,uru,utu;
double ur1_lost,ut1_lost,uru_lost,utu_lost;

if(RR.slab.b<=1e-6)RR.slab.b= 1e-6;

if(Debug(DEBUG_GRID_CALC)&&i==0&&j==0){
fprintf(stderr,"+   i   j ");
fprintf(stderr,"      a         b          g     |");
fprintf(stderr,"     M_R        grid  |");
fprintf(stderr,"     M_T        grid\n");
}

if(Debug(DEBUG_EVERY_CALC)){
if(!CALCULATING_GRID)
fprintf(stderr,"a=%8.5f b=%10.5f g=%8.5f ",RR.slab.a,RR.slab.b,RR.slab.g);
else{
if(j==0)fprintf(stderr,".");
if(i+1==GRID_SIZE&&j==0)fprintf(stderr,"\n");
}
}

RT_Flip(MM.flip_sample,RR.method.quad_pts,&RR.slab,&ur1,&ut1,&uru,&utu);

if(Debug(DEBUG_EVERY_CALC)&&!CALCULATING_GRID)
fprintf(stderr,"ur1=%8.5f ut1=%8.5f\n",ur1,ut1);

The_Grid[GRID_SIZE*i+j][A_COLUMN]= RR.slab.a;
The_Grid[GRID_SIZE*i+j][B_COLUMN]= RR.slab.b;
The_Grid[GRID_SIZE*i+j][G_COLUMN]= RR.slab.g;
The_Grid[GRID_SIZE*i+j][UR1_COLUMN]= ur1;
The_Grid[GRID_SIZE*i+j][UT1_COLUMN]= ut1;
The_Grid[GRID_SIZE*i+j][URU_COLUMN]= uru;
The_Grid[GRID_SIZE*i+j][UTU_COLUMN]= utu;

if(Debug(DEBUG_MC)){
if(i==0&&j==0)
fprintf(stderr,"Filling Grid\n%2d ",0);




RR.a= RR.slab.a;
RR.b= RR.slab.b;
RR.g= RR.slab.g;
MC_Lost(MM,RR,100000,&ur1,&ut1,&uru,&utu,&ur1_lost,&ut1_lost,&uru_lost,&utu_lost);

The_Grid[GRID_SIZE*i+j][UR1_LOST_COLUMN]= ur1_lost;
The_Grid[GRID_SIZE*i+j][UT1_LOST_COLUMN]= ut1_lost;
The_Grid[GRID_SIZE*i+j][URU_LOST_COLUMN]= uru_lost;
The_Grid[GRID_SIZE*i+j][UTU_LOST_COLUMN]= utu_lost;





fprintf(stderr,".");
if(j==GRID_SIZE-1){
if(i!=GRID_SIZE-1)
fprintf(stderr,"\n%2d ",i+1);
else
fprintf(stderr,"\n");
}

if(Debug(DEBUG_GRID_CALC)){
fprintf(stderr,"+ %3d %3d ",i,j);
fprintf(stderr,"%10.5f %10.5f %10.5f |",RR.slab.a,RR.slab.b,RR.slab.g);
fprintf(stderr,"%10.5f %10.5f |",MM.m_r,uru);
fprintf(stderr,"%10.5f %10.5f \n",MM.m_t,utu);
}
}
}

/*:38*/
#line 76 "iad_calc.w"

/*55:*/
#line 1056 "iad_calc.w"

/*54:*/
#line 1053 "iad_calc.w"

void Fill_Grid(struct measure_type m,struct invert_type r,int force_new)

/*:54*/
#line 1057 "iad_calc.w"

{
if(force_new||!Same_Calc_State(m,r)){
switch(r.search){
case FIND_AB:
Fill_AB_Grid(m,r);
break;
case FIND_AG:
Fill_AG_Grid(m,r);
break;
case FIND_BG:
Fill_BG_Grid(m,r);
break;
case FIND_BaG:
Fill_BaG_Grid(m,r);
break;
case FIND_BsG:
Fill_BsG_Grid(m,r);
break;
default:
AD_error("Attempt to fill grid for unknown search case.");
}
}

Get_Calc_State(&MGRID,&RGRID);
}

/*:55*/
#line 77 "iad_calc.w"

/*35:*/
#line 607 "iad_calc.w"

/*34:*/
#line 604 "iad_calc.w"

void Near_Grid_Points(double r,double t,search_type s,int*i_min,int*j_min)

/*:34*/
#line 608 "iad_calc.w"

{
int i,j;
double fval;
double smallest= 10.0;
struct measure_type old_mm;
struct invert_type old_rr;
(void)r;
(void)t;
(void)s;

if(Debug(DEBUG_GRID))
fprintf(stderr,"GRID: Finding best grid points\n");
Get_Calc_State(&old_mm,&old_rr);

*i_min= 0;
*j_min= 0;
for(i= 0;i<GRID_SIZE;i++){
for(j= 0;j<GRID_SIZE;j++){

CALCULATING_GRID= 1;
fval= Calculate_Grid_Distance(i,j);
CALCULATING_GRID= 0;

if(fval<smallest){
*i_min= i;
*j_min= j;
smallest= fval;
}
}
}

Set_Calc_State(old_mm,old_rr);
}

/*:35*/
#line 78 "iad_calc.w"

/*40:*/
#line 765 "iad_calc.w"

/*39:*/
#line 762 "iad_calc.w"

void Fill_AB_Grid(struct measure_type m,struct invert_type r)

/*:39*/
#line 766 "iad_calc.w"

{
int i,j;
double min_log_b= -8;
double max_log_b= +8;

if(Debug(DEBUG_GRID))
fprintf(stderr,"GRID: Filling AB grid (g=%.5f)\n",RR.slab.g);

if(The_Grid==NULL)Allocate_Grid(r.search);
/*47:*/
#line 886 "iad_calc.w"

GG_a= 0.0;
GG_b= 0.0;
GG_g= 0.0;
GG_bs= 0.0;
GG_ba= 0.0;


/*:47*/
#line 776 "iad_calc.w"


Set_Calc_State(m,r);

GG_g= RR.slab.g;
for(i= 0;i<GRID_SIZE;i++){
double x= (double)i/(GRID_SIZE-1.0);
RR.slab.b= exp(min_log_b+(max_log_b-min_log_b)*x);
for(j= 0;j<GRID_SIZE;j++){
/*42:*/
#line 813 "iad_calc.w"

{
double x= (double)j/(GRID_SIZE-1.0);
RR.slab.a= 1.0-(1.0-x)*(1.0-x)*(1.0+2.0*x);
}

/*:42*/
#line 785 "iad_calc.w"

fill_grid_entry(i,j);
}
}

The_Grid_Initialized= TRUE;
The_Grid_Search= FIND_AB;
}

/*:40*/
#line 79 "iad_calc.w"

/*46:*/
#line 851 "iad_calc.w"

/*45:*/
#line 848 "iad_calc.w"

void Fill_AG_Grid(struct measure_type m,struct invert_type r)

/*:45*/
#line 852 "iad_calc.w"

{
int i,j;
double max_a= -10;
double min_a= 10;

if(Debug(DEBUG_GRID))
fprintf(stderr,"GRID: Filling AG grid\n");

if(The_Grid==NULL)Allocate_Grid(r.search);
/*47:*/
#line 886 "iad_calc.w"

GG_a= 0.0;
GG_b= 0.0;
GG_g= 0.0;
GG_bs= 0.0;
GG_ba= 0.0;


/*:47*/
#line 862 "iad_calc.w"


Set_Calc_State(m,r);
GG_b= r.slab.b;
for(i= 0;i<GRID_SIZE;i++){
/*43:*/
#line 822 "iad_calc.w"

{
double x= (double)i/(GRID_SIZE-1.0);
double xx= (1.0-x)*(1.0-x)*(1.0+2.0*x);
RR.slab.g= (1.0-2.0*xx)*MAX_ABS_G;
}

/*:43*/
#line 867 "iad_calc.w"

for(j= 0;j<GRID_SIZE;j++){
/*42:*/
#line 813 "iad_calc.w"

{
double x= (double)j/(GRID_SIZE-1.0);
RR.slab.a= 1.0-(1.0-x)*(1.0-x)*(1.0+2.0*x);
}

/*:42*/
#line 869 "iad_calc.w"

fill_grid_entry(i,j);
if(RR.slab.a<min_a)min_a= RR.slab.a;
if(RR.slab.a> max_a)max_a= RR.slab.a;
}
}

if(Debug(DEBUG_GRID)){
fprintf(stderr,"GRID: a        = %9.7f to %9.7f \n",min_a,max_a);
fprintf(stderr,"GRID: b        = %9.5f \n",r.slab.b);
fprintf(stderr,"GRID: g  range = %9.6f to %9.6f \n",-MAX_ABS_G,MAX_ABS_G);
}

The_Grid_Initialized= TRUE;
The_Grid_Search= FIND_AG;
}

/*:46*/
#line 80 "iad_calc.w"

/*49:*/
#line 904 "iad_calc.w"

/*48:*/
#line 901 "iad_calc.w"

void Fill_BG_Grid(struct measure_type m,struct invert_type r)

/*:48*/
#line 905 "iad_calc.w"

{
int i,j;
double min_log_b= -8;
double max_log_b= +10;

if(The_Grid==NULL)Allocate_Grid(r.search);
/*47:*/
#line 886 "iad_calc.w"

GG_a= 0.0;
GG_b= 0.0;
GG_g= 0.0;
GG_bs= 0.0;
GG_ba= 0.0;


/*:47*/
#line 912 "iad_calc.w"


if(Debug(DEBUG_GRID))
fprintf(stderr,"GRID: Filling BG grid\n");

Set_Calc_State(m,r);
RR.slab.a= RR.default_a;
GG_a= RR.slab.a;

for(i= 0;i<GRID_SIZE;i++){
double x= (double)i/(GRID_SIZE-1.0);
RR.slab.b= exp(min_log_b+(max_log_b-min_log_b)*x);
for(j= 0;j<GRID_SIZE;j++){
/*44:*/
#line 829 "iad_calc.w"

{
double x= (double)j/(GRID_SIZE-1.0);
double xx= (1.0-x)*(1.0-x)*(1.0+2.0*x);
RR.slab.g= (1.0-2.0*xx)*MAX_ABS_G;
}

/*:44*/
#line 925 "iad_calc.w"

fill_grid_entry(i,j);
}
}

if(Debug(DEBUG_GRID)){
fprintf(stderr,"GRID: a        = %9.7f \n",RR.default_a);
fprintf(stderr,"GRID: b  range = %9.5f to %9.3f \n",exp(min_log_b),exp(max_log_b));
fprintf(stderr,"GRID: g  range = %9.6f to %9.6f \n",-MAX_ABS_G,MAX_ABS_G);
}

The_Grid_Initialized= TRUE;
The_Grid_Search= FIND_BG;
}

/*:49*/
#line 81 "iad_calc.w"

/*51:*/
#line 948 "iad_calc.w"

/*50:*/
#line 945 "iad_calc.w"

void Fill_BaG_Grid(struct measure_type m,struct invert_type r)

/*:50*/
#line 949 "iad_calc.w"

{
int i,j;
double max_a= -10;
double min_a= 10;
double bs= r.default_bs;
double min_log_ba= -8;
double max_log_ba= +10;

if(The_Grid==NULL)Allocate_Grid(r.search);
/*47:*/
#line 886 "iad_calc.w"

GG_a= 0.0;
GG_b= 0.0;
GG_g= 0.0;
GG_bs= 0.0;
GG_ba= 0.0;


/*:47*/
#line 959 "iad_calc.w"


if(Debug(DEBUG_GRID)){
fprintf(stderr,"GRID: Filling BaG grid\n");
fprintf(stderr,"GRID:       bs = %9.5f\n",bs);
fprintf(stderr,"GRID: ba range = %9.6f to %9.3f \n",exp(min_log_ba),exp(max_log_ba));
}

Set_Calc_State(m,r);
GG_bs= bs;
for(i= 0;i<GRID_SIZE;i++){
double x= (double)i/(GRID_SIZE-1.0);
double ba= exp(min_log_ba+(max_log_ba-min_log_ba)*x);
RR.slab.b= ba+bs;
if(RR.slab.b> 0)
RR.slab.a= bs/RR.slab.b;
else
RR.slab.a= 0;
if(RR.slab.a<0)RR.slab.a= 0;
if(RR.slab.a<min_a)min_a= RR.slab.a;
if(RR.slab.a> max_a)max_a= RR.slab.a;

for(j= 0;j<GRID_SIZE;j++){
/*44:*/
#line 829 "iad_calc.w"

{
double x= (double)j/(GRID_SIZE-1.0);
double xx= (1.0-x)*(1.0-x)*(1.0+2.0*x);
RR.slab.g= (1.0-2.0*xx)*MAX_ABS_G;
}

/*:44*/
#line 982 "iad_calc.w"

fill_grid_entry(i,j);
}
}

if(Debug(DEBUG_GRID)){
fprintf(stderr,"GRID: a        = %9.7f to %9.7f \n",min_a,max_a);
fprintf(stderr,"GRID: b  range = %9.5f to %9.3f \n",exp(min_log_ba)+bs,exp(max_log_ba)+bs);
fprintf(stderr,"GRID: g  range = %9.6f to %9.6f \n",-MAX_ABS_G,MAX_ABS_G);
}

The_Grid_Initialized= TRUE;
The_Grid_Search= FIND_BaG;
}

/*:51*/
#line 82 "iad_calc.w"

/*53:*/
#line 1004 "iad_calc.w"

/*52:*/
#line 1001 "iad_calc.w"

void Fill_BsG_Grid(struct measure_type m,struct invert_type r)

/*:52*/
#line 1005 "iad_calc.w"

{
int i,j;
double max_a= -10;
double min_a= 10;
double ba= r.default_ba;
double min_log_bs= -8;
double max_log_bs= +10;

if(The_Grid==NULL)Allocate_Grid(r.search);
/*47:*/
#line 886 "iad_calc.w"

GG_a= 0.0;
GG_b= 0.0;
GG_g= 0.0;
GG_bs= 0.0;
GG_ba= 0.0;


/*:47*/
#line 1015 "iad_calc.w"


if(Debug(DEBUG_GRID)){
fprintf(stderr,"GRID: Filling BsG grid\n");
fprintf(stderr,"GRID:       ba = %9.5f\n",ba);
fprintf(stderr,"GRID: bs range = %9.6f to %9.3f \n",exp(min_log_bs),exp(max_log_bs));
}

Set_Calc_State(m,r);
GG_ba= RR.default_ba;
for(i= 0;i<GRID_SIZE;i++){
double x= (double)i/(GRID_SIZE-1.0);
double bs= exp(min_log_bs+(max_log_bs-min_log_bs)*x);
RR.slab.b= ba+bs;
if(RR.slab.b> 0)
RR.slab.a= 1-RR.default_ba/RR.slab.b;
else
RR.slab.a= 0;
if(RR.slab.a<0)RR.slab.a= 0;
if(RR.slab.a<min_a)min_a= RR.slab.a;
if(RR.slab.a> max_a)max_a= RR.slab.a;

for(j= 0;j<GRID_SIZE;j++){
/*44:*/
#line 829 "iad_calc.w"

{
double x= (double)j/(GRID_SIZE-1.0);
double xx= (1.0-x)*(1.0-x)*(1.0+2.0*x);
RR.slab.g= (1.0-2.0*xx)*MAX_ABS_G;
}

/*:44*/
#line 1038 "iad_calc.w"

fill_grid_entry(i,j);
}
}

if(Debug(DEBUG_GRID)){
fprintf(stderr,"GRID: a  range = %9.7f to %9.7f \n",min_a,max_a);
fprintf(stderr,"GRID: b  range = %9.5f to %9.3f \n",exp(min_log_bs)+ba,exp(max_log_bs)+ba);
fprintf(stderr,"GRID: g  range = %9.6f to %9.6f \n",-MAX_ABS_G,MAX_ABS_G);
}

The_Grid_Initialized= TRUE;
The_Grid_Search= FIND_BsG;
}

/*:53*/
#line 83 "iad_calc.w"

/*25:*/
#line 439 "iad_calc.w"

/*24:*/
#line 436 "iad_calc.w"

void Grid_ABG(int i,int j,guess_type*guess)

/*:24*/
#line 440 "iad_calc.w"

{
if(0<=i&&i<GRID_SIZE&&0<=j&&j<GRID_SIZE){
guess->a= The_Grid[GRID_SIZE*i+j][A_COLUMN];
guess->b= The_Grid[GRID_SIZE*i+j][B_COLUMN];
guess->g= The_Grid[GRID_SIZE*i+j][G_COLUMN];
guess->ur1_lost= The_Grid[GRID_SIZE*i+j][UR1_LOST_COLUMN];
guess->ut1_lost= The_Grid[GRID_SIZE*i+j][UT1_LOST_COLUMN];
guess->uru_lost= The_Grid[GRID_SIZE*i+j][URU_LOST_COLUMN];
guess->utu_lost= The_Grid[GRID_SIZE*i+j][UTU_LOST_COLUMN];
guess->distance= Calculate_Grid_Distance(i,j);
}else{
guess->a= 0.5;
guess->b= 0.5;
guess->g= 0.5;
guess->ur1_lost= 0;
guess->ut1_lost= 0;
guess->uru_lost= 0;
guess->utu_lost= 0;
guess->distance= 999;
}
}

/*:25*/
#line 84 "iad_calc.w"

/*6:*/
#line 184 "iad_calc.w"

/*5:*/
#line 181 "iad_calc.w"

double Gain(int sphere,struct measure_type m,double uru_sample,double uru_third)

/*:5*/
#line 185 "iad_calc.w"

{
double inv_gain;

if(sphere==REFLECTION_SPHERE){
if(m.baffle_r){
inv_gain= m.rw_r+(m.at_r/m.aw_r)*uru_third;
inv_gain*= m.aw_r+(1-m.at_r)*(m.ad_r*m.rd_r+m.as_r*uru_sample);
inv_gain= 1.0-inv_gain;
}else{
inv_gain= 1.0-m.aw_r*m.rw_r-m.ad_r*m.rd_r-m.as_r*uru_sample-m.at_r*uru_third;
}
}else
if(m.baffle_t){
inv_gain= m.rw_t+(m.at_t/m.aw_t)*uru_third;
inv_gain*= m.aw_t+(1-m.at_t)*(m.ad_t*m.rd_t+m.as_t*uru_sample);
inv_gain= 1.0-inv_gain;
}else{
inv_gain= 1.0-m.aw_t*m.rw_t-m.ad_t*m.rd_t-m.as_t*uru_sample-m.at_t*uru_third;
}
return 1.0/inv_gain;
}

/*:6*/
#line 85 "iad_calc.w"

/*8:*/
#line 224 "iad_calc.w"

/*7:*/
#line 221 "iad_calc.w"

double Gain_11(struct measure_type m,double URU,double tdiffuse)

/*:7*/
#line 225 "iad_calc.w"

{
double G,GP,G11;

G= Gain(REFLECTION_SPHERE,m,URU,0);
GP= Gain(TRANSMISSION_SPHERE,m,URU,0);

G11= G/(1-m.as_r*m.as_t*m.aw_r*m.aw_t*(1-m.at_r)*(1-m.at_t)
*G*GP*tdiffuse*tdiffuse);

return G11;
}

/*:8*/
#line 86 "iad_calc.w"

/*10:*/
#line 249 "iad_calc.w"

/*9:*/
#line 246 "iad_calc.w"

double Gain_22(struct measure_type m,double URU,double tdiffuse)

/*:9*/
#line 250 "iad_calc.w"

{
double G,GP,G22;

G= Gain(REFLECTION_SPHERE,m,URU,0);
GP= Gain(TRANSMISSION_SPHERE,m,URU,0);

G22= GP/(1-m.as_r*m.as_t*m.aw_r*m.aw_t*(1-m.at_r)*(1-m.at_t)
*G*GP*tdiffuse*tdiffuse);

return G22;
}

/*:10*/
#line 87 "iad_calc.w"

/*12:*/
#line 294 "iad_calc.w"

/*11:*/
#line 290 "iad_calc.w"

double Two_Sphere_R(struct measure_type m,
double UR1,double URU,double UT1,double UTU)

/*:11*/
#line 295 "iad_calc.w"

{
double x,GP;
GP= Gain(TRANSMISSION_SPHERE,m,URU,0);

x= m.ad_r*(1-m.at_r)*m.rw_r*Gain_11(m,URU,UTU);
x*= (1-m.f_r)*UR1+m.rw_r*m.f_r+(1-m.f_r)*m.as_t*(1-m.at_t)*m.rw_t*UT1*UTU*GP;
return x;
}

/*:12*/
#line 88 "iad_calc.w"

/*14:*/
#line 331 "iad_calc.w"

/*13:*/
#line 327 "iad_calc.w"

double Two_Sphere_T(struct measure_type m,
double UR1,double URU,double UT1,double UTU)

/*:13*/
#line 332 "iad_calc.w"

{
double x,G;
G= Gain(REFLECTION_SPHERE,m,URU,0);
x= m.ad_t*(1-m.at_t)*m.rw_t*Gain_22(m,URU,UTU);
x*= (1-m.f_r)*UT1+(1-m.at_r)*m.rw_r*m.as_r*UTU*(m.f_r*m.rw_r+(1-m.f_r)*UR1)*G;
return x;
}

/*:14*/
#line 89 "iad_calc.w"


/*61:*/
#line 1184 "iad_calc.w"

/*60:*/
#line 1177 "iad_calc.w"

void Calculate_Distance_With_Corrections(
double UR1,double UT1,
double Ru,double Tu,
double URU,double UTU,
double*M_R,double*M_T,double*dev)

/*:60*/
#line 1185 "iad_calc.w"

{
/*62:*/
#line 1221 "iad_calc.w"

double UR1_calc,UT1_calc,URU_calc,UTU_calc;

URU_calc= URU-MM.uru_lost;
if(URU_calc<0)URU_calc= 0;

UTU_calc= UTU-MM.utu_lost;
if(UTU_calc<0)UTU_calc= 0;

/*:62*//*63:*/
#line 1236 "iad_calc.w"

UR1_calc= UR1-(1.0-MM.fraction_of_ru_in_mr)*Ru-MM.ur1_lost;
if(UR1_calc<0)UR1_calc= 0;

UT1_calc= UT1-(1.0-MM.fraction_of_tu_in_mt)*Tu-MM.ut1_lost;
if(UT1_calc<0)UT1_calc= 0;

/*:63*/
#line 1187 "iad_calc.w"


switch(MM.num_spheres){
case 0:
/*64:*/
#line 1248 "iad_calc.w"

{
*M_R= UR1_calc;
*M_T= UT1_calc;
}

/*:64*/
#line 1191 "iad_calc.w"

break;

case 1:
if(MM.method==COMPARISON){
/*77:*/
#line 1590 "iad_calc.w"

{
*M_R= UR1_calc;
*M_T= UT1_calc;
}

/*:77*/
#line 1196 "iad_calc.w"

}else{
/*71:*/
#line 1371 "iad_calc.w"

double P_std,P,P_0,G,G_0,G_std,r_first,P_ss,P_su;
r_first= 1;

if(MM.baffle_r)r_first= MM.rw_r*(1-MM.at_r);

UR1_calc= UR1-Ru-MM.ur1_lost;
if(UR1_calc<0)UR1_calc= 0;

G_0= Gain(REFLECTION_SPHERE,MM,0.0,0.0);
G= Gain(REFLECTION_SPHERE,MM,URU_calc,0.0);
G_std= Gain(REFLECTION_SPHERE,MM,MM.rstd_r,0.0);

P_std= G_std*(MM.rstd_r*(1-MM.f_r)+MM.f_r*MM.rw_r);
P_0= G_0*(MM.f_r*MM.rw_r);
P_ss= r_first*(UR1_calc*(1-MM.f_r)+MM.f_r*MM.rw_r);
P_su= MM.rw_r*(1-MM.f_r)*MM.fraction_of_ru_in_mr*Ru;
P= G*(P_ss+P_su);

*M_R= MM.rstd_r*(P-P_0)/(P_std-P_0);

if(Debug(DEBUG_SPHERE_GAIN)&&!CALCULATING_GRID){
double G_none,P_none,M_none;
G_none= Gain(REFLECTION_SPHERE,MM,URU,0.0);
P_none= G_none*(P_ss+P_su);
M_none= MM.rstd_r*(P_none-P_0)/(P_std-P_0);
fprintf(stderr,"SPHERE: REFLECTION\n");
fprintf(stderr,"SPHERE:       baffle = %d\n",MM.baffle_r);
fprintf(stderr,"SPHERE:       R_u collected = %5.1f%%\n",MM.fraction_of_ru_in_mr*100);
fprintf(stderr,"SPHERE:       hits sphere wall first = %5.1f%%\n",MM.f_r*100);
fprintf(stderr,"SPHERE:       UR1 = %7.3f   UR1_calc = %7.3f\n",UR1,UR1_calc);
fprintf(stderr,"SPHERE:       URU = %7.3f   URU_calc = %7.3f\n",URU,URU_calc);
fprintf(stderr,"SPHERE:       R_u = %7.3f\n",Ru);
fprintf(stderr,"SPHERE:       G_0 = %7.3f        P_0 = %7.3f\n",G_0,P_0);
fprintf(stderr,"SPHERE:         G = %7.3f          P = %7.3f\n",G,P);
if(MM.uru_lost> 0)
fprintf(stderr,"SPHERE: G_no_lost = %7.3f  P_no_lost = %7.3f\n",G_none,P_none);
fprintf(stderr,"SPHERE:     G_cal = %7.3f      P_cal = %7.3f\n",G_std,P_std);
if(MM.uru_lost> 0)
fprintf(stderr,"SPHERE: M_no_lost = %7.3f\n",M_none);
fprintf(stderr,"SPHERE:       M_R = %7.3f\n",*M_R);
}

/*:71*/
#line 1198 "iad_calc.w"

/*75:*/
#line 1517 "iad_calc.w"

{
double r_cal,r_third,P_su,P_ss;
double r_first= 1;

if(MM.at_t==0){
r_cal= MM.rw_t;
r_third= MM.rw_t;
}
else if(MM.fraction_of_tu_in_mt> 0){
r_cal= MM.rstd_t;
r_third= MM.rstd_t;
}else{
r_cal= MM.rstd_t;
r_third= 0;
}

if(MM.baffle_t)r_first= MM.rw_t*(1-MM.at_t)+r_third*MM.at_t;

UT1_calc= UT1-Tu-MM.ut1_lost;
if(UT1_calc<0)UT1_calc= 0;

G= Gain(TRANSMISSION_SPHERE,MM,URU_calc,r_third);
G_std= Gain(TRANSMISSION_SPHERE,MM,0,r_cal);

P_su= r_third*Tu*MM.fraction_of_tu_in_mt;
P_ss= r_first*UT1_calc;

*M_T= (P_su+P_ss)*G/G_std;

if(Debug(DEBUG_SPHERE_GAIN)&&!CALCULATING_GRID){
double G_none,P_none,M_none;
G_none= Gain(TRANSMISSION_SPHERE,MM,URU,0.0);
P_none= (P_ss+P_su);
M_none= (P_su+P_ss)*G_none/G_std;
fprintf(stderr,"SPHERE: TRANSMISSION\n");
fprintf(stderr,"SPHERE:       baffle = %d\n",MM.baffle_t);
fprintf(stderr,"SPHERE:       T_u collected = %5.1f%%\n",MM.fraction_of_tu_in_mt*100);
fprintf(stderr,"SPHERE:       UR1 = %7.3f   UR1_calc = %7.3f\n",UR1,UR1_calc);
fprintf(stderr,"SPHERE:       URU = %7.3f   URU_calc = %7.3f\n",URU,URU_calc);
fprintf(stderr,"SPHERE:       UT1 = %7.3f   UT1_calc = %7.3f\n",UT1,UT1_calc);
fprintf(stderr,"SPHERE:       T_u = %7.3f\n",Tu);
fprintf(stderr,"SPHERE:         G = %7.3f          P = %7.3f\n",G,P_su+P_ss);
if(MM.uru_lost> 0)
fprintf(stderr,"SPHERE: G_no_lost = %7.3f  P_no_lost = %7.3f\n",G_none,P_none);
fprintf(stderr,"SPHERE:     G_cal = %7.3f      P_cal = %7.3f\n",G_std,1.0);
fprintf(stderr,"SPHERE:   r_third = %7.3f      r_cal = %7.3f\n",r_third,r_cal);
fprintf(stderr,"SPHERE:   r_first = %7.3f\n",r_first);
fprintf(stderr,"SPHERE:       Psu = %7.3f        Pss = %7.3f\n",P_su,P_ss);
if(MM.uru_lost> 0)
fprintf(stderr,"SPHERE: M_no_lost = %7.3f\n",M_none);
fprintf(stderr,"SPHERE:       M_T = %7.3f\n",*M_T);
fprintf(stderr,"\n");
}
}

/*:75*/
#line 1199 "iad_calc.w"

}break;

case 2:
/*79:*/
#line 1616 "iad_calc.w"

{
double R_0,T_0;
R_0= Two_Sphere_R(MM,0,0,0,0);
T_0= Two_Sphere_T(MM,0,0,0,0);

*M_R= MM.rstd_r*(Two_Sphere_R(MM,UR1_calc,URU_calc,UT1_calc,UTU_calc)-R_0)/
(Two_Sphere_R(MM,MM.rstd_r,MM.rstd_r,0,0)-R_0);
*M_T= (Two_Sphere_T(MM,UR1_calc,URU_calc,UT1_calc,UTU_calc)-T_0)/
(Two_Sphere_T(MM,0,0,1,1)-T_0);
}

/*:79*/
#line 1203 "iad_calc.w"

break;

default:
fprintf(stderr,"Bad number of spheres = %d\n",MM.num_spheres);
exit(EXIT_FAILURE);
}

/*80:*/
#line 1633 "iad_calc.w"


if(RR.search==FIND_A||RR.search==FIND_G||RR.search==FIND_B||
RR.search==FIND_Bs||RR.search==FIND_Ba){
/*81:*/
#line 1649 "iad_calc.w"


if(MM.m_t> 0){
if(RR.metric==RELATIVE)
*dev= fabs(MM.m_t-*M_T)/(MM.m_t+ABIT);
else
*dev= fabs(MM.m_t-*M_T);
}else{
if(RR.metric==RELATIVE)
*dev= fabs(MM.m_r-*M_R)/(MM.m_r+ABIT);
else
*dev= fabs(MM.m_r-*M_R);
}

/*:81*/
#line 1637 "iad_calc.w"

}else{
/*82:*/
#line 1668 "iad_calc.w"


if(RR.metric==RELATIVE){
if(MM.m_t> ABIT)
*dev= T_TRUST_FACTOR*fabs(MM.m_t-*M_T)/(UTU_calc+ABIT);
if(RR.default_a!=0){
*dev+= fabs(MM.m_r-*M_R)/(URU_calc+ABIT);
}
}else{
*dev= T_TRUST_FACTOR*fabs(MM.m_t-*M_T);
if(RR.default_a!=0)
*dev+= fabs(MM.m_r-*M_R);
}


/*:82*/
#line 1639 "iad_calc.w"

}

/*:80*/
#line 1211 "iad_calc.w"

/*83:*/
#line 1687 "iad_calc.w"

if((Debug(DEBUG_ITERATIONS)&&!CALCULATING_GRID)||
(Debug(DEBUG_GRID_CALC)&&CALCULATING_GRID)){
fprintf(stderr,"%10.5f %10.4f %10.5f |",RR.slab.a,RR.slab.b,RR.slab.g);
fprintf(stderr," %10.5f %10.5f |",MM.m_r,*M_R);
fprintf(stderr," %10.5f %10.5f |",MM.m_t,*M_T);
fprintf(stderr,"%10.3f\n",*dev);
}

/*:83*/
#line 1212 "iad_calc.w"

}

/*:61*/
#line 91 "iad_calc.w"

/*59:*/
#line 1121 "iad_calc.w"

/*58:*/
#line 1118 "iad_calc.w"

double Calculate_Grid_Distance(int i,int j)

/*:58*/
#line 1122 "iad_calc.w"

{
double ur1,ut1,uru,utu,Ru,Tu,b,dev,LR,LT;

if(Debug(DEBUG_GRID_CALC)&&i==0&&j==0){
fprintf(stderr,"+   i   j ");
fprintf(stderr,"      a         b          g     |");
fprintf(stderr,"     M_R        grid   |");
fprintf(stderr,"     M_T        grid   |  distance\n");
}

if(Debug(DEBUG_GRID_CALC))
fprintf(stderr,"g %3d %3d ",i,j);

b= The_Grid[GRID_SIZE*i+j][B_COLUMN];
ur1= The_Grid[GRID_SIZE*i+j][UR1_COLUMN];
ut1= The_Grid[GRID_SIZE*i+j][UT1_COLUMN];
uru= The_Grid[GRID_SIZE*i+j][URU_COLUMN];
utu= The_Grid[GRID_SIZE*i+j][UTU_COLUMN];
RR.slab.a= The_Grid[GRID_SIZE*i+j][A_COLUMN];
RR.slab.b= The_Grid[GRID_SIZE*i+j][B_COLUMN];
RR.slab.g= The_Grid[GRID_SIZE*i+j][G_COLUMN];
RR.slab.g= The_Grid[GRID_SIZE*i+j][G_COLUMN];

if(Debug(DEBUG_MC)){
MM.ur1_lost= The_Grid[GRID_SIZE*i+j][UR1_LOST_COLUMN];
MM.ut1_lost= The_Grid[GRID_SIZE*i+j][UT1_LOST_COLUMN];
MM.uru_lost= The_Grid[GRID_SIZE*i+j][URU_LOST_COLUMN];
MM.utu_lost= The_Grid[GRID_SIZE*i+j][UTU_LOST_COLUMN];
}

Sp_mu_RT_Flip(MM.flip_sample,
RR.slab.n_top_slide,RR.slab.n_slab,RR.slab.n_bottom_slide,
RR.slab.b_top_slide,b,RR.slab.b_bottom_slide,
RR.slab.cos_angle,&Ru,&Tu);

CALCULATING_GRID= 1;
Calculate_Distance_With_Corrections(ur1,ut1,Ru,Tu,uru,utu,&LR,&LT,&dev);
CALCULATING_GRID= 0;

return dev;
}

/*:59*/
#line 92 "iad_calc.w"

/*57:*/
#line 1096 "iad_calc.w"

/*56:*/
#line 1093 "iad_calc.w"

void Calculate_Distance(double*M_R,double*M_T,double*deviation)

/*:56*/
#line 1097 "iad_calc.w"

{
double Ru,Tu,ur1,ut1,uru,utu;

if(RR.slab.b<=1e-6)
RR.slab.b= 1e-6;

RT_Flip(MM.flip_sample,RR.method.quad_pts,&RR.slab,&ur1,&ut1,&uru,&utu);

Sp_mu_RT_Flip(MM.flip_sample,
RR.slab.n_top_slide,RR.slab.n_slab,RR.slab.n_bottom_slide,
RR.slab.b_top_slide,RR.slab.b,RR.slab.b_bottom_slide,
RR.slab.cos_angle,&Ru,&Tu);

if((!CALCULATING_GRID&&Debug(DEBUG_ITERATIONS))||
(CALCULATING_GRID&&Debug(DEBUG_GRID_CALC)))
fprintf(stderr,"        ");

Calculate_Distance_With_Corrections(ur1,ut1,Ru,Tu,uru,utu,M_R,M_T,deviation);
}

/*:57*/
#line 93 "iad_calc.w"

/*33:*/
#line 574 "iad_calc.w"

/*32:*/
#line 571 "iad_calc.w"

void abg_distance(double a,double b,double g,guess_type*guess)

/*:32*/
#line 575 "iad_calc.w"

{
double m_r,m_t,distance;
struct measure_type old_mm;
struct invert_type old_rr;

Get_Calc_State(&old_mm,&old_rr);

RR.slab.a= a;
RR.slab.b= b;
RR.slab.g= g;

Calculate_Distance(&m_r,&m_t,&distance);

Set_Calc_State(old_mm,old_rr);

guess->a= a;
guess->b= b;
guess->g= g;
guess->distance= distance;
}

/*:33*/
#line 94 "iad_calc.w"


/*85:*/
#line 1699 "iad_calc.w"

/*84:*/
#line 1696 "iad_calc.w"

double Find_AG_fn(double x[])

/*:84*/
#line 1700 "iad_calc.w"

{
double m_r,m_t,deviation;
RR.slab.a= acalc2a(x[1]);
RR.slab.g= gcalc2g(x[2]);
Calculate_Distance(&m_r,&m_t,&deviation);
return deviation;
}

/*:85*/
#line 96 "iad_calc.w"

/*87:*/
#line 1712 "iad_calc.w"

/*86:*/
#line 1709 "iad_calc.w"

double Find_AB_fn(double x[])

/*:86*/
#line 1713 "iad_calc.w"

{
double m_r,m_t,deviation;
RR.slab.a= acalc2a(x[1]);
RR.slab.b= bcalc2b(x[2]);
Calculate_Distance(&m_r,&m_t,&deviation);
return deviation;
}

/*:87*/
#line 97 "iad_calc.w"

/*89:*/
#line 1729 "iad_calc.w"

/*88:*/
#line 1722 "iad_calc.w"

double Find_Ba_fn(double x)

/*:88*/
#line 1730 "iad_calc.w"

{
double m_r,m_t,deviation,ba,bs;

bs= RR.slab.b;
ba= bcalc2b(x);
RR.slab.b= ba+bs;
RR.slab.a= bs/(ba+bs);

Calculate_Distance(&m_r,&m_t,&deviation);

RR.slab.b= bs;
return deviation;
}

/*:89*/
#line 98 "iad_calc.w"

/*91:*/
#line 1751 "iad_calc.w"

/*90:*/
#line 1748 "iad_calc.w"

double Find_Bs_fn(double x)

/*:90*/
#line 1752 "iad_calc.w"

{
double m_r,m_t,deviation,ba,bs;

ba= RR.slab.b;
bs= bcalc2b(x);
RR.slab.b= ba+bs;
RR.slab.a= bs/(ba+bs);

Calculate_Distance(&m_r,&m_t,&deviation);

RR.slab.b= ba;
return deviation;
}

/*:91*/
#line 99 "iad_calc.w"

/*93:*/
#line 1770 "iad_calc.w"

/*92:*/
#line 1767 "iad_calc.w"

double Find_A_fn(double x)

/*:92*/
#line 1771 "iad_calc.w"

{
double m_r,m_t,deviation;
RR.slab.a= acalc2a(x);
Calculate_Distance(&m_r,&m_t,&deviation);
return deviation;
}

/*:93*/
#line 100 "iad_calc.w"

/*95:*/
#line 1782 "iad_calc.w"

/*94:*/
#line 1779 "iad_calc.w"

double Find_B_fn(double x)

/*:94*/
#line 1783 "iad_calc.w"

{
double m_r,m_t,deviation;
RR.slab.b= bcalc2b(x);
Calculate_Distance(&m_r,&m_t,&deviation);
return deviation;
}

/*:95*/
#line 101 "iad_calc.w"

/*97:*/
#line 1794 "iad_calc.w"

/*96:*/
#line 1791 "iad_calc.w"

double Find_G_fn(double x)

/*:96*/
#line 1795 "iad_calc.w"

{
double m_r,m_t,deviation;
RR.slab.g= gcalc2g(x);
Calculate_Distance(&m_r,&m_t,&deviation);
return deviation;
}

/*:97*/
#line 102 "iad_calc.w"

/*99:*/
#line 1806 "iad_calc.w"

/*98:*/
#line 1803 "iad_calc.w"

double Find_BG_fn(double x[])

/*:98*/
#line 1807 "iad_calc.w"

{
double m_r,m_t,deviation;
RR.slab.b= bcalc2b(x[1]);
RR.slab.g= gcalc2g(x[2]);
RR.slab.a= RR.default_a;
Calculate_Distance(&m_r,&m_t,&deviation);
return deviation;
}

/*:99*/
#line 103 "iad_calc.w"

/*101:*/
#line 1826 "iad_calc.w"

/*100:*/
#line 1823 "iad_calc.w"

double Find_BaG_fn(double x[])

/*:100*/
#line 1827 "iad_calc.w"

{
double m_r,m_t,deviation;

RR.slab.b= bcalc2b(x[1])+RR.default_bs;
if(RR.slab.b<=0)
RR.slab.a= 0;
else
RR.slab.a= RR.default_bs/RR.slab.b;

RR.slab.g= gcalc2g(x[2]);

Calculate_Distance(&m_r,&m_t,&deviation);
return deviation;
}

/*:101*/
#line 104 "iad_calc.w"

/*103:*/
#line 1846 "iad_calc.w"

/*102:*/
#line 1843 "iad_calc.w"

double Find_BsG_fn(double x[])

/*:102*/
#line 1847 "iad_calc.w"

{
double m_r,m_t,deviation;

RR.slab.b= bcalc2b(x[1])+RR.default_ba;
if(RR.slab.b<=0)
RR.slab.a= 0;
else
RR.slab.a= 1.0-RR.default_ba/RR.slab.b;

RR.slab.g= gcalc2g(x[2]);
Calculate_Distance(&m_r,&m_t,&deviation);
return deviation;
}

/*:103*/
#line 105 "iad_calc.w"

/*105:*/
#line 1870 "iad_calc.w"

/*104:*/
#line 1867 "iad_calc.w"

double maxloss(double f)

/*:104*/
#line 1871 "iad_calc.w"

{
struct measure_type m_old;
struct invert_type r_old;
double m_r,m_t,deviation;

Get_Calc_State(&m_old,&r_old);

RR.slab.a= 1.0;
MM.ur1_lost*= f;
MM.ut1_lost*= f;

Calculate_Distance(&m_r,&m_t,&deviation);

Set_Calc_State(m_old,r_old);
deviation= ((MM.m_r+MM.m_t)-(m_r+m_t));

return deviation;
}

/*:105*/
#line 106 "iad_calc.w"

/*107:*/
#line 1900 "iad_calc.w"

/*106:*/
#line 1896 "iad_calc.w"

void Max_Light_Loss(struct measure_type m,struct invert_type r,
double*ur1_loss,double*ut1_loss)

/*:106*/
#line 1901 "iad_calc.w"

{
struct measure_type m_old;
struct invert_type r_old;

*ur1_loss= m.ur1_lost;
*ut1_loss= m.ut1_lost;

if(Debug(DEBUG_LOST_LIGHT))
fprintf(stderr,"\nlost before ur1=%7.5f, ut1=%7.5f\n",*ur1_loss,*ut1_loss);

Get_Calc_State(&m_old,&r_old);

Set_Calc_State(m,r);

if(maxloss(1.0)*maxloss(0.0)<0){
double frac;
frac= zbrent(maxloss,0.00,1.0,0.001);

*ur1_loss= m.ur1_lost*frac;
*ut1_loss= m.ut1_lost*frac;
}

Set_Calc_State(m_old,r_old);
if(Debug(DEBUG_LOST_LIGHT))
fprintf(stderr,"lost after  ur1=%7.5f, ut1=%7.5f\n",*ur1_loss,*ut1_loss);
}

/*:107*/
#line 107 "iad_calc.w"


/*:1*/
