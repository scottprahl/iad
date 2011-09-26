/*2:*/
#line 41 "iad_pub.w"

/*4:*/
#line 66 "iad_pub.w"

void Inverse_RT(struct measure_type m,struct invert_type*r)

/*:4*/
#line 42 "iad_pub.w"
;
/*9:*/
#line 148 "iad_pub.w"

int measure_OK(struct measure_type m,struct invert_type r)

/*:9*/
#line 43 "iad_pub.w"
;
/*16:*/
#line 312 "iad_pub.w"

search_type determine_search(struct measure_type m,struct invert_type r)

/*:16*/
#line 44 "iad_pub.w"
;
/*20:*/
#line 467 "iad_pub.w"

void Initialize_Result(struct measure_type m,struct invert_type*r)

/*:20*/
#line 45 "iad_pub.w"
;
/*26:*/
#line 538 "iad_pub.w"

void ez_Inverse_RT(double n,double nslide,double UR1,double UT1,double Tc,
double*a,double*b,double*g,int*error)

/*:26*/
#line 46 "iad_pub.w"
;
/*28:*/
#line 581 "iad_pub.w"

void Initialize_Measure(struct measure_type*m)

/*:28*/
#line 47 "iad_pub.w"
;
/*37:*/
#line 781 "iad_pub.w"

void Calculate_MR_MT(struct measure_type m,
struct invert_type r,
int include_MC,
double*M_R,
double*M_T)

/*:37*/
#line 48 "iad_pub.w"
;
/*41:*/
#line 852 "iad_pub.w"

int MinMax_MR_MT(struct measure_type m,
struct invert_type r)

/*:41*/
#line 49 "iad_pub.w"
;
/*39:*/
#line 812 "iad_pub.w"

void Calculate_Minimum_MR(struct measure_type m,
struct invert_type r,double*mr,double*mt)

/*:39*/
#line 50 "iad_pub.w"
;

/*:2*/
