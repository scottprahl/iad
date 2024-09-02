/*2:*/
#line 111 "iad_calc.w"

/*5:*/
#line 181 "iad_calc.w"

double Gain(int sphere,struct measure_type m,double uru_sample,double uru_third)

/*:5*/
#line 112 "iad_calc.w"
;
/*7:*/
#line 221 "iad_calc.w"

double Gain_11(struct measure_type m,double URU,double tdiffuse)

/*:7*/
#line 113 "iad_calc.w"
;
/*9:*/
#line 246 "iad_calc.w"

double Gain_22(struct measure_type m,double URU,double tdiffuse)

/*:9*/
#line 114 "iad_calc.w"
;
/*11:*/
#line 290 "iad_calc.w"

double Two_Sphere_R(struct measure_type m,
double UR1,double URU,double UT1,double UTU)

/*:11*/
#line 115 "iad_calc.w"
;
/*13:*/
#line 327 "iad_calc.w"

double Two_Sphere_T(struct measure_type m,
double UR1,double URU,double UT1,double UTU)

/*:13*/
#line 116 "iad_calc.w"
;
/*16:*/
#line 370 "iad_calc.w"

void Set_Calc_State(struct measure_type m,struct invert_type r)

/*:16*/
#line 117 "iad_calc.w"
;
/*18:*/
#line 383 "iad_calc.w"

void Get_Calc_State(struct measure_type*m,struct invert_type*r)

/*:18*/
#line 118 "iad_calc.w"
;
/*20:*/
#line 396 "iad_calc.w"

boolean_type Same_Calc_State(struct measure_type m,struct invert_type r)

/*:20*/
#line 119 "iad_calc.w"
;
/*26:*/
#line 473 "iad_calc.w"

boolean_type Valid_Grid(struct measure_type m,struct invert_type r)

/*:26*/
#line 120 "iad_calc.w"
;
/*22:*/
#line 421 "iad_calc.w"

void Allocate_Grid(search_type s)

/*:22*/
#line 121 "iad_calc.w"
;
/*54:*/
#line 1053 "iad_calc.w"

void Fill_Grid(struct measure_type m,struct invert_type r,int force_new)

/*:54*/
#line 122 "iad_calc.w"
;
/*34:*/
#line 604 "iad_calc.w"

void Near_Grid_Points(double r,double t,search_type s,int*i_min,int*j_min)

/*:34*/
#line 123 "iad_calc.w"
;
/*24:*/
#line 436 "iad_calc.w"

void Grid_ABG(int i,int j,guess_type*guess)

/*:24*/
#line 124 "iad_calc.w"
;
/*84:*/
#line 1696 "iad_calc.w"

double Find_AG_fn(double x[])

/*:84*/
#line 125 "iad_calc.w"
;
/*86:*/
#line 1709 "iad_calc.w"

double Find_AB_fn(double x[])

/*:86*/
#line 126 "iad_calc.w"
;
/*88:*/
#line 1722 "iad_calc.w"

double Find_Ba_fn(double x)

/*:88*/
#line 127 "iad_calc.w"
;
/*90:*/
#line 1748 "iad_calc.w"

double Find_Bs_fn(double x)

/*:90*/
#line 128 "iad_calc.w"
;
/*92:*/
#line 1767 "iad_calc.w"

double Find_A_fn(double x)

/*:92*/
#line 129 "iad_calc.w"
;
/*94:*/
#line 1779 "iad_calc.w"

double Find_B_fn(double x)

/*:94*/
#line 130 "iad_calc.w"
;
/*96:*/
#line 1791 "iad_calc.w"

double Find_G_fn(double x)

/*:96*/
#line 131 "iad_calc.w"
;
/*98:*/
#line 1803 "iad_calc.w"

double Find_BG_fn(double x[])

/*:98*/
#line 132 "iad_calc.w"
;
/*102:*/
#line 1843 "iad_calc.w"

double Find_BsG_fn(double x[])

/*:102*/
#line 133 "iad_calc.w"
;
/*100:*/
#line 1823 "iad_calc.w"

double Find_BaG_fn(double x[])

/*:100*/
#line 134 "iad_calc.w"
;
/*48:*/
#line 901 "iad_calc.w"

void Fill_BG_Grid(struct measure_type m,struct invert_type r)

/*:48*/
#line 135 "iad_calc.w"
;
/*52:*/
#line 1001 "iad_calc.w"

void Fill_BsG_Grid(struct measure_type m,struct invert_type r)

/*:52*/
#line 136 "iad_calc.w"
;
/*50:*/
#line 945 "iad_calc.w"

void Fill_BaG_Grid(struct measure_type m,struct invert_type r)

/*:50*/
#line 137 "iad_calc.w"
;
/*60:*/
#line 1177 "iad_calc.w"

void Calculate_Distance_With_Corrections(
double UR1,double UT1,
double Ru,double Tu,
double URU,double UTU,
double*M_R,double*M_T,double*dev)

/*:60*/
#line 138 "iad_calc.w"
;
/*56:*/
#line 1093 "iad_calc.w"

void Calculate_Distance(double*M_R,double*M_T,double*deviation)

/*:56*/
#line 139 "iad_calc.w"
;
/*58:*/
#line 1118 "iad_calc.w"

double Calculate_Grid_Distance(int i,int j)

/*:58*/
#line 140 "iad_calc.w"
;
/*32:*/
#line 571 "iad_calc.w"

void abg_distance(double a,double b,double g,guess_type*guess)

/*:32*/
#line 141 "iad_calc.w"
;
/*104:*/
#line 1867 "iad_calc.w"

double maxloss(double f)

/*:104*/
#line 142 "iad_calc.w"
;
/*106:*/
#line 1896 "iad_calc.w"

void Max_Light_Loss(struct measure_type m,struct invert_type r,
double*ur1_loss,double*ut1_loss)

/*:106*/
#line 143 "iad_calc.w"
;
/*36:*/
#line 653 "iad_calc.w"

void RT_Flip(int flip,int n,struct AD_slab_type*slab,double*UR1,double*UT1,
double*URU,double*UTU)

/*:36*/
#line 144 "iad_calc.w"
;


/*:2*/
