/*2:*/
#line 40 "iad_util.w"

/*3:*/
#line 77 "iad_util.w"

double What_Is_B(struct AD_slab_type slab,double Tc)

/*:3*/
#line 41 "iad_util.w"
;
/*9:*/
#line 218 "iad_util.w"

void Estimate_RT(struct measure_type m,struct invert_type r,double*rt,double*tt,
double*rd,double*rc,double*td,double*tc)

/*:9*/
#line 42 "iad_util.w"
;
/*15:*/
#line 313 "iad_util.w"

double a2acalc(double a)

/*:15*/
#line 43 "iad_util.w"
;
/*17:*/
#line 342 "iad_util.w"

double acalc2a(double acalc)

/*:17*/
#line 44 "iad_util.w"
;
/*19:*/
#line 367 "iad_util.w"

double g2gcalc(double g)

/*:19*/
#line 45 "iad_util.w"
;
/*21:*/
#line 385 "iad_util.w"

double gcalc2g(double gcalc)

/*:21*/
#line 46 "iad_util.w"
;
/*23:*/
#line 404 "iad_util.w"

double b2bcalc(double b)

/*:23*/
#line 47 "iad_util.w"
;
/*25:*/
#line 434 "iad_util.w"

double bcalc2b(double bcalc)

/*:25*/
#line 48 "iad_util.w"
;
/*27:*/
#line 449 "iad_util.w"

void twoprime(double a,double b,double g,double*ap,double*bp)

/*:27*/
#line 49 "iad_util.w"
;
/*29:*/
#line 469 "iad_util.w"

void twounprime(double ap,double bp,double g,double*a,double*b)

/*:29*/
#line 50 "iad_util.w"
;
/*31:*/
#line 488 "iad_util.w"

void abgg2ab(double a1,double b1,double g1,double g2,double*a2,double*b2)

/*:31*/
#line 51 "iad_util.w"
;
/*33:*/
#line 509 "iad_util.w"

void abgb2ag(double a1,double b1,double b2,double*a2,double*g2)

/*:33*/
#line 52 "iad_util.w"
;
/*40:*/
#line 610 "iad_util.w"

void quick_guess(struct measure_type m,struct invert_type r,double*a,double*b,double*g)

/*:40*/
#line 53 "iad_util.w"
;
/*53:*/
#line 750 "iad_util.w"

void Set_Debugging(unsigned long debug_level)

/*:53*/
#line 54 "iad_util.w"
;
/*55:*/
#line 761 "iad_util.w"

int Debug(unsigned long mask)

/*:55*/
#line 55 "iad_util.w"
;
/*57:*/
#line 775 "iad_util.w"

void Print_Invert_Type(struct invert_type r)

/*:57*/
#line 56 "iad_util.w"
;
/*59:*/
#line 795 "iad_util.w"

void Print_Measure_Type(struct measure_type m)

/*:59*/
#line 57 "iad_util.w"
;


/*:2*/
