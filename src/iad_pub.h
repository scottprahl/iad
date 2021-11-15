/*2:*/
#line 43 "iad_pub.w"

/*4:*/
#line 71 "iad_pub.w"

void Inverse_RT(struct measure_type m,struct invert_type*r)

/*:4*/
#line 44 "iad_pub.w"
;
/*9:*/
#line 154 "iad_pub.w"

int measure_OK(struct measure_type m,struct invert_type r)

/*:9*/
#line 45 "iad_pub.w"
;
/*16:*/
#line 335 "iad_pub.w"

search_type determine_search(struct measure_type m,struct invert_type r)

/*:16*/
#line 46 "iad_pub.w"
;
/*20:*/
#line 506 "iad_pub.w"

void Initialize_Result(struct measure_type m,struct invert_type*r)

/*:20*/
#line 47 "iad_pub.w"
;
/*26:*/
#line 577 "iad_pub.w"

void ez_Inverse_RT(double n,double nslide,double UR1,double UT1,double Tc,
double*a,double*b,double*g,int*error)

/*:26*/
#line 48 "iad_pub.w"
;
/*28:*/
#line 618 "iad_pub.w"

void Initialize_Measure(struct measure_type*m)

/*:28*/
#line 49 "iad_pub.w"
;
/*37:*/
#line 818 "iad_pub.w"

void Calculate_MR_MT(struct measure_type m,
struct invert_type r,
int include_MC,
double*M_R,
double*M_T)

/*:37*/
#line 50 "iad_pub.w"
;
/*41:*/
#line 893 "iad_pub.w"

int MinMax_MR_MT(struct measure_type m,
struct invert_type r)

/*:41*/
#line 51 "iad_pub.w"
;
/*39:*/
#line 849 "iad_pub.w"

void Calculate_Minimum_MR(struct measure_type m,
struct invert_type r,double*mr,double*mt)

/*:39*/
#line 52 "iad_pub.w"
;
/*43:*/
#line 931 "iad_pub.w"

void Spheres_Inverse_RT2(double*sample,
double*illumination,
double*sphere_r,
double*sphere_t,
double*analysis,
double*measurement,
double*a,
double*b,
double*g)

/*:43*/
#line 53 "iad_pub.w"
;


/*:2*/
