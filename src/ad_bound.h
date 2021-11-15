/*2:*/
#line 28 "ad_bound.w"

#define TOP_BOUNDARY 0
#define BOTTOM_BOUNDARY 1 \


#line 29 "ad_bound.w"

/*4:*/
#line 51 "ad_bound.w"

void Init_Boundary(struct AD_slab_type slab,int n,
double*R01,double*R10,double*T01,double*T10,
char boundary)

/*:4*/
#line 30 "ad_bound.w"
;
/*7:*/
#line 83 "ad_bound.w"

void Boundary_RT(double n_i,double n_g,double n_t,int n,double b,
double*R,double*T)

/*:7*/
#line 31 "ad_bound.w"
;
/*15:*/
#line 204 "ad_bound.w"

void Add_Top(int n,double*R01,double*R10,double*T01,double*T10,
double**R12,double**R21,double**T12,double**T21,
double**R02,double**R20,double**T02,double**T20,
double**atemp,double**btemp)

/*:15*/
#line 32 "ad_bound.w"
;
/*17:*/
#line 229 "ad_bound.w"

void Add_Bottom(int n,double**R01,double**R10,double**T01,double**T10,
double*R12,double*R21,double*T12,double*T21,
double**R02,double**R20,double**T02,double**T20,
double**atemp,double**btemp)

/*:17*/
#line 33 "ad_bound.w"
;
/*19:*/
#line 296 "ad_bound.w"

void Add_Slides(int n,double*R01,double*R10,double*T01,double*T10,
double**R,double**T,
double**R_total,double**T_total,
double**atemp,double**btemp)

/*:19*/
#line 34 "ad_bound.w"
;
/*21:*/
#line 360 "ad_bound.w"

void Sp_RT(int n,struct AD_slab_type slab,double*ur1,double*ut1,double*uru,double*utu)

/*:21*/
#line 35 "ad_bound.w"
;

/*:2*/
