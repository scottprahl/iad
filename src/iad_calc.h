/*2:*/
#line 98 "iad_calc.w"

/*5:*/
#line 228 "iad_calc.w"

double Gain(int sphere,struct measure_type m,double URU)

/*:5*/
#line 99 "iad_calc.w"
;
/*7:*/
#line 267 "iad_calc.w"

double Gain_11(struct measure_type m,double URU,double tdiffuse)

/*:7*/
#line 100 "iad_calc.w"
;
/*9:*/
#line 292 "iad_calc.w"

double Gain_22(struct measure_type m,double URU,double tdiffuse)

/*:9*/
#line 101 "iad_calc.w"
;
/*11:*/
#line 336 "iad_calc.w"

double Two_Sphere_R(struct measure_type m,
double UR1,double URU,double UT1,double UTU)

/*:11*/
#line 102 "iad_calc.w"
;
/*13:*/
#line 373 "iad_calc.w"

double Two_Sphere_T(struct measure_type m,
double UR1,double URU,double UT1,double UTU)

/*:13*/
#line 103 "iad_calc.w"
;
/*16:*/
#line 416 "iad_calc.w"

void Set_Calc_State(struct measure_type m,struct invert_type r)

/*:16*/
#line 104 "iad_calc.w"
;
/*18:*/
#line 433 "iad_calc.w"

void Get_Calc_State(struct measure_type*m,struct invert_type*r)

/*:18*/
#line 105 "iad_calc.w"
;
/*20:*/
#line 446 "iad_calc.w"

boolean_type Same_Calc_State(struct measure_type m,struct invert_type r)

/*:20*/
#line 106 "iad_calc.w"
;
/*26:*/
#line 514 "iad_calc.w"

boolean_type Valid_Grid(struct measure_type m,search_type s)

/*:26*/
#line 107 "iad_calc.w"
;
/*22:*/
#line 471 "iad_calc.w"

void Allocate_Grid(search_type s)

/*:22*/
#line 108 "iad_calc.w"
;
/*53:*/
#line 998 "iad_calc.w"

void Fill_Grid(struct measure_type m,struct invert_type r)

/*:53*/
#line 109 "iad_calc.w"
;
/*34:*/
#line 620 "iad_calc.w"

void Near_Grid_Points(double r,double t,search_type s,int*i_min,int*j_min)

/*:34*/
#line 110 "iad_calc.w"
;
/*24:*/
#line 485 "iad_calc.w"

void Grid_ABG(int i,int j,guess_type*guess)

/*:24*/
#line 111 "iad_calc.w"
;
/*69:*/
#line 1352 "iad_calc.w"

double Find_AG_fn(double x[])

/*:69*/
#line 112 "iad_calc.w"
;
/*71:*/
#line 1365 "iad_calc.w"

double Find_AB_fn(double x[])

/*:71*/
#line 113 "iad_calc.w"
;
/*73:*/
#line 1378 "iad_calc.w"

double Find_Ba_fn(double x)

/*:73*/
#line 114 "iad_calc.w"
;
/*75:*/
#line 1404 "iad_calc.w"

double Find_Bs_fn(double x)

/*:75*/
#line 115 "iad_calc.w"
;
/*77:*/
#line 1423 "iad_calc.w"

double Find_A_fn(double x)

/*:77*/
#line 116 "iad_calc.w"
;
/*79:*/
#line 1435 "iad_calc.w"

double Find_B_fn(double x)

/*:79*/
#line 117 "iad_calc.w"
;
/*81:*/
#line 1447 "iad_calc.w"

double Find_G_fn(double x)

/*:81*/
#line 118 "iad_calc.w"
;
/*83:*/
#line 1459 "iad_calc.w"

double Find_BG_fn(double x[])

/*:83*/
#line 119 "iad_calc.w"
;
/*87:*/
#line 1499 "iad_calc.w"

double Find_BsG_fn(double x[])

/*:87*/
#line 120 "iad_calc.w"
;
/*85:*/
#line 1479 "iad_calc.w"

double Find_BaG_fn(double x[])

/*:85*/
#line 121 "iad_calc.w"
;
/*47:*/
#line 890 "iad_calc.w"

void Fill_BG_Grid(struct measure_type m,struct invert_type r)

/*:47*/
#line 122 "iad_calc.w"
;
/*51:*/
#line 965 "iad_calc.w"

void Fill_BsG_Grid(struct measure_type m,struct invert_type r)

/*:51*/
#line 123 "iad_calc.w"
;
/*49:*/
#line 926 "iad_calc.w"

void Fill_BaG_Grid(struct measure_type m,struct invert_type r)

/*:49*/
#line 124 "iad_calc.w"
;
/*59:*/
#line 1116 "iad_calc.w"

void Calculate_Distance_With_Corrections(
double UR1,double UT1,
double Rc,double Tc,
double URU,double UTU,
double*M_R,double*M_T,double*dev)

/*:59*/
#line 125 "iad_calc.w"
;
/*55:*/
#line 1043 "iad_calc.w"

void Calculate_Distance(double*M_R,double*M_T,double*deviation)

/*:55*/
#line 126 "iad_calc.w"
;
/*57:*/
#line 1072 "iad_calc.w"

double Calculate_Grid_Distance(int i,int j)

/*:57*/
#line 127 "iad_calc.w"
;
/*32:*/
#line 587 "iad_calc.w"

void abg_distance(double a,double b,double g,guess_type*guess)

/*:32*/
#line 128 "iad_calc.w"
;
/*89:*/
#line 1523 "iad_calc.w"

double maxloss(double f)

/*:89*/
#line 129 "iad_calc.w"
;
/*91:*/
#line 1552 "iad_calc.w"

void Max_Light_Loss(struct measure_type m,struct invert_type r,
double*ur1_loss,double*ut1_loss)

/*:91*/
#line 130 "iad_calc.w"
;


/*:2*/
