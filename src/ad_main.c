/*2:*/
#line 18 "../src/ad_main.w"


#include <stdio.h> 
#include <string.h> 
#include <stdlib.h> 
#include <math.h> 
#include <unistd.h> 
#include <errno.h> 

#include "ad_globl.h"
#include "ad_prime.h"
#include "ad_cone.h"
#include "version.h"

/*10:*/
#line 242 "../src/ad_main.w"


static void print_version(void)
{
fprintf(stdout,"ad %s\n",Version);
fprintf(stdout,"Copyright 1993-2025 Scott Prahl, scott.prahl@oit.edu\n");
fprintf(stdout,"          (see Applied Optics, 32:559-568, 1993)\n\n");
fprintf(stdout,"This is free software; see the source for copying conditions.\n");
fprintf(stdout,"There is no warranty; not even for MERCHANTABILITY or FITNESS.\n");
fprintf(stdout,"FOR A PARTICULAR PURPOSE.\n");
}

/*:10*/
#line 32 "../src/ad_main.w"

/*11:*/
#line 254 "../src/ad_main.w"

static void print_usage(void)
{
fprintf(stdout,"ad %s\n\n",Version);
fprintf(stdout,"ad finds the reflection and transmission from optical properties\n\n");
fprintf(stdout,"Usage:  ad [options] input\n\n");
fprintf(stdout,"Options:\n");
fprintf(stdout,"  -a #             albedo (0-1)\n");
fprintf(stdout,"  -b #             optical thickness (>0)\n");
fprintf(stdout,"  -g #             scattering anisotropy (-1 to 1)\n");
fprintf(stdout,"  -h               display help\n");
fprintf(stdout,"  -i theta         oblique incidence at angle theta\n");
fprintf(stdout,"  -m               machine readable output\n");
fprintf(stdout,"  -n #             specify index of refraction of slab\n");
fprintf(stdout,"  -o filename      explicitly specify filename for output\n");
fprintf(stdout,"  -q #             quadrature points 4, 8, 16, 32\n");
fprintf(stdout,"  -s #             specify index of refraction of slide\n");
fprintf(stdout,"  -v               version information\n");
fprintf(stdout,"Examples:\n");
fprintf(stdout,"  ad data                        UR1, UT1, URU, UTU in data.rt\n");
fprintf(stdout,"  ad -m data                     data.rt in machine readable format\n");
fprintf(stdout,"  ad data -o out.txt             out.txt is the \n");
fprintf(stdout,"  ad -a 0.3                      a=0.3, b=inf, g=0.0, n=1.0\n");
fprintf(stdout,"  ad -a 0.3 -b 0.4               a=0.3, b=0.4, g=0.0, n=1.0\n");
fprintf(stdout,"  ad -a 0.3 -b 0.4 -g 0.5        a=0.3, b=0.4, g=0.5, n=1.0\n");
fprintf(stdout,"  ad -a 0.3 -b 0.4 -n 1.5        a=0.3, b=0.4, g=0.0, n=1.5\n\n");
fprintf(stdout,"inputfile has lines of the form:\n");
fprintf(stdout,"    a b g nslab ntopslide nbottomlslide btopslide bbottomslide q\n");
fprintf(stdout,"where:\n");
fprintf(stdout,"    1) a = albedo\n");
fprintf(stdout,"    2) b = optical thickness\n");
fprintf(stdout,"    3) g = anisotropy\n");
fprintf(stdout,"    4) nslab = index of refraction of slab\n");
fprintf(stdout,"    5) ntopslide = index of refraction of glass slide on top\n");
fprintf(stdout,"    6) nbottomslide = index of refraction of glass slide on bottom\n");
fprintf(stdout,"       negative is diffuse, i.e., -0.7 => 0.7 Lambertian reflection\n");
fprintf(stdout,"    7) btopslide = optical depth of top slide (for absorbing slides)\n");
fprintf(stdout,"    8) bbottomslide = optical depth of bottom slide (for absorbing slides)\n");
fprintf(stdout,"    9) q = number of quadrature points\n\n");
fprintf(stdout,"Report bugs to <scott.prahl@oit.edu>\n\n");
}

/*:11*/
#line 33 "../src/ad_main.w"

/*12:*/
#line 297 "../src/ad_main.w"

static char*strdup_together(char*s,char*t)
{
char*both;

if(s==NULL){
if(t==NULL)return NULL;
return strdup(t);
}

if(t==NULL)
return strdup(s);

both= malloc(strlen(s)+strlen(t)+1);
if(both==NULL)
fprintf(stderr,"Could not allocate memory for both strings.\n");

strcpy(both,s);
strcat(both,t);
return both;
}

/*:12*/
#line 34 "../src/ad_main.w"

/*13:*/
#line 320 "../src/ad_main.w"


static double my_strtod(const char*str)
{
char*endptr;
errno= 0;

double val= strtod(str,&endptr);

if(endptr==str){

fprintf(stderr,"Error: No conversion could be performed for `%s`.\n",str);
exit(EXIT_FAILURE);
}

if(*endptr!='\0'){

printf("Partial conversion: converted value = %f, remaining string = %s\n",val,endptr);
exit(EXIT_FAILURE);
}

if(errno==ERANGE){

printf("Error: The value is out of range of double.\n");
exit(EXIT_FAILURE);
}

return val;
}

/*:13*/
#line 35 "../src/ad_main.w"

/*15:*/
#line 359 "../src/ad_main.w"

static void validate_slab(struct AD_slab_type slab,int nstreams,int machine)
{
if(slab.a<0||slab.a> 1){
fprintf(stderr,"Bad Albedo a=%f\n",slab.a);
exit(EXIT_FAILURE);
}

if(slab.b<0){
fprintf(stderr,"Bad Optical Thickness b=%f\n",slab.b);
exit(EXIT_FAILURE);
}

if(slab.g<=-1||slab.g>=1){
fprintf(stderr,"Bad Anisotropy g=%f\n",slab.g);
exit(EXIT_FAILURE);
}

if(slab.n_slab<0||slab.n_slab> 10){
fprintf(stderr,"Bad Slab Index n=%f\n",slab.n_slab);
exit(EXIT_FAILURE);
}

if(slab.n_top_slide<1||slab.n_top_slide> 10){
fprintf(stderr,"Bad Top Slide Index n=%f\n",slab.n_top_slide);
exit(EXIT_FAILURE);
}

if(slab.n_bottom_slide> 0){
if(slab.n_bottom_slide<1||slab.n_bottom_slide> 10){
fprintf(stderr,"Bad Bottom Slide Index n=%f\n",slab.n_bottom_slide);
exit(EXIT_FAILURE);
}
}else{
if(slab.n_bottom_slide<-1){
fprintf(stderr,"Bad Lambertian Bottom Reflection n=%f\n",slab.n_bottom_slide);
exit(EXIT_FAILURE);
}
}

if(slab.b_top_slide<0||slab.b_top_slide> 10){
fprintf(stderr,"Bad Top Slide Optical Thickness b=%f\n",slab.b_top_slide);
exit(EXIT_FAILURE);
}

if(slab.b_bottom_slide<0||slab.b_bottom_slide> 10){
fprintf(stderr,"Bad Bottom Slide Optical Thickness b=%f\n",slab.b_bottom_slide);
exit(EXIT_FAILURE);
}

if(nstreams<4||nstreams%4!=0){
fprintf(stderr,"Bad Number of Quadrature Points npts=%d\n",nstreams);
fprintf(stderr,"Should be a multiple of four!\n");
exit(EXIT_FAILURE);
}
}/*:15*/
#line 36 "../src/ad_main.w"


int main(int argc,char**argv){

/*3:*/
#line 64 "../src/ad_main.w"

struct AD_slab_type slab;
int nstreams= 24;
double anisotropy= 0;
double albedo= 0.5;
double index_of_refraction= 1.0;
double index_of_slide1= 1.0;
double index_of_slide2= 1.0;
double optical_thickness= 100;
char*g_out_name= NULL;
double g_incident_cosine= 1;
int machine_readable_output= 0;
double R1,T1,URU,UTU;

/*:3*/
#line 40 "../src/ad_main.w"


if(argc==1){
print_usage();
exit(EXIT_FAILURE);
}

/*7:*/
#line 139 "../src/ad_main.w"

{
int c;
double x;

while((c= getopt(argc,argv,"hvma:b:g:i:n:o:q:s:t:"))!=-1){
switch(c){

case'a':
albedo= my_strtod(optarg);
break;

case'i':
x= my_strtod(optarg);
if(x<0||x> 90)
fprintf(stderr,"Incident angle must be between 0 and 90 degrees\n");
else
g_incident_cosine= cos(x*M_PI/180.0);
break;

case'o':
g_out_name= strdup(optarg);
break;

case'n':
index_of_refraction= my_strtod(optarg);
break;

case's':
index_of_slide1= my_strtod(optarg);
index_of_slide2= index_of_slide1;
break;

case't':
index_of_slide2= my_strtod(optarg);
break;

case'm':
machine_readable_output= 1;
break;

case'q':
nstreams= (int)my_strtod(optarg);
break;

case'b':
optical_thickness= my_strtod(optarg);
break;

case'g':
anisotropy= my_strtod(optarg);
break;

case'v':
print_version();
exit(EXIT_SUCCESS);

default:
fprintf(stderr,"unknown option '%c'\n",c);


case'h':
print_usage();
exit(EXIT_SUCCESS);
}
}

argc-= optind;
argv+= optind;
}


/*:7*/
#line 47 "../src/ad_main.w"


if(argc>=1){
/*8:*/
#line 213 "../src/ad_main.w"

if(argc> 1){
fprintf(stderr,"Only a single file can be processed at a time\n");
fprintf(stderr,"try 'apply ad file1 file2 ... fileN'\n");
exit(EXIT_FAILURE);
}

if(argc==1&&strcmp(argv[0],"-")!=0){

if(freopen(argv[0],"r",stdin)==NULL){
fprintf(stderr,"Could not open file '%s'\n",argv[0]);
exit(EXIT_FAILURE);
}

if(g_out_name==NULL)
g_out_name= strdup_together(argv[0],".rt");
}

/*:8*/
#line 50 "../src/ad_main.w"

/*9:*/
#line 233 "../src/ad_main.w"


if(g_out_name!=NULL){
if(freopen(g_out_name,"w",stdout)==NULL){
fprintf(stderr,"Could not open file <%s> for output",g_out_name);
exit(EXIT_FAILURE);
}
}

/*:9*/
#line 51 "../src/ad_main.w"

while(feof(stdin)==0){
slab.phase_function= HENYEY_GREENSTEIN;
/*5:*/
#line 97 "../src/ad_main.w"

{
int fileflag;
fileflag= scanf("%lf",&slab.a);
slab.cos_angle= g_incident_cosine;

if(fileflag!=EOF)
fileflag= scanf("%lf",&slab.b);
if(fileflag!=EOF)
fileflag= scanf("%lf",&slab.g);
if(fileflag!=EOF)
fileflag= scanf("%lf",&slab.n_slab);
if(fileflag!=EOF)
fileflag= scanf("%lf",&slab.n_top_slide);
if(fileflag!=EOF)
fileflag= scanf("%lf",&slab.n_bottom_slide);
if(fileflag!=EOF)
fileflag= scanf("%lf",&slab.b_top_slide);
if(fileflag!=EOF)
fileflag= scanf("%lf",&slab.b_bottom_slide);
if(fileflag!=EOF)
fileflag= scanf("%d",&nstreams);
}

/*:5*/
#line 54 "../src/ad_main.w"

/*6:*/
#line 121 "../src/ad_main.w"


validate_slab(slab,nstreams,machine_readable_output);
RT(nstreams,&slab,&R1,&T1,&URU,&UTU);

if(machine_readable_output)
printf("%9.5f \t%9.5f \t%9.5f \t%9.5f\n",R1,T1,URU,UTU);
else{
printf("UR1 = Total Reflection   for Normal  Illumination\n");
printf("UT1 = Total Transmission for Normal  Illumination\n");
printf("URU = Total Reflection   for Diffuse Illumination\n");
printf("UTU = Total Transmission for Diffuse Illumination\n\n");
printf("   UR1    \t   UT1    \t   URU    \t   UTU\n");
printf("%9.5f \t%9.5f \t%9.5f \t%9.5f\n",R1,T1,URU,UTU);
}

/*:6*/
#line 55 "../src/ad_main.w"

}
}else{
/*4:*/
#line 83 "../src/ad_main.w"


slab.phase_function= HENYEY_GREENSTEIN;
slab.a= albedo;
slab.b= optical_thickness;
slab.g= anisotropy;
slab.n_slab= index_of_refraction;
slab.n_top_slide= index_of_slide1;
slab.n_bottom_slide= index_of_slide2;
slab.b_top_slide= 0.0;
slab.b_bottom_slide= 0.0;
slab.cos_angle= g_incident_cosine;

/*:4*/
#line 58 "../src/ad_main.w"

/*6:*/
#line 121 "../src/ad_main.w"


validate_slab(slab,nstreams,machine_readable_output);
RT(nstreams,&slab,&R1,&T1,&URU,&UTU);

if(machine_readable_output)
printf("%9.5f \t%9.5f \t%9.5f \t%9.5f\n",R1,T1,URU,UTU);
else{
printf("UR1 = Total Reflection   for Normal  Illumination\n");
printf("UT1 = Total Transmission for Normal  Illumination\n");
printf("URU = Total Reflection   for Diffuse Illumination\n");
printf("UTU = Total Transmission for Diffuse Illumination\n\n");
printf("   UR1    \t   UT1    \t   URU    \t   UTU\n");
printf("%9.5f \t%9.5f \t%9.5f \t%9.5f\n",R1,T1,URU,UTU);
}

/*:6*/
#line 59 "../src/ad_main.w"

}
return EXIT_SUCCESS;
}

/*:2*/
