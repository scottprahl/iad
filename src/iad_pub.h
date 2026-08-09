/*2:*/
#line 43 "iad_pub.w"

/*4:*/
#line 69 "iad_pub.w"

void Inverse_RT(struct measure_type m,struct invert_type*r)

/*:4*/
#line 44 "iad_pub.w"
;
/*12:*/
#line 229 "iad_pub.w"

int measure_OK(struct measure_type m,struct invert_type r,int flag_bad)

/*:12*/
#line 45 "iad_pub.w"
;
/*21:*/
#line 482 "iad_pub.w"

search_type determine_search(struct measure_type m,struct invert_type r)

/*:21*/
#line 46 "iad_pub.w"
;
/*25:*/
#line 688 "iad_pub.w"

void Initialize_Result(struct measure_type m,struct invert_type*r,int overwrite_defaults)

/*:25*/
#line 47 "iad_pub.w"
;
/*31:*/
#line 765 "iad_pub.w"

void ez_Inverse_RT(double n,double nslide,double UR1,double UT1,double Tu,
double*a,double*b,double*g,int*error)

/*:31*/
#line 48 "iad_pub.w"
;
/*33:*/
#line 806 "iad_pub.w"

void Initialize_Measure(struct measure_type*m)

/*:33*/
#line 49 "iad_pub.w"
;
/*42:*/
#line 1007 "iad_pub.w"

void Calculate_MR_MT(struct measure_type m,
struct invert_type r,
int include_MC,
int include_spheres,
double*M_R,
double*M_T)

/*:42*/
#line 50 "iad_pub.w"
;
/*46:*/
#line 1099 "iad_pub.w"

int MinMax_MR_MT(struct measure_type m,
struct invert_type r)

/*:46*/
#line 51 "iad_pub.w"
;
/*44:*/
#line 1054 "iad_pub.w"

void Calculate_Minimum_MR(struct measure_type m,
struct invert_type r,double*mr,double*mt)

/*:44*/
#line 52 "iad_pub.w"
;


/*:2*/
