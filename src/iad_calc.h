

double Gain (int sphere, struct measure_type m, double URU);

double Gain_11 (struct measure_type m, double URU, double tdiffuse);

double Gain_22 (struct measure_type m, double URU, double tdiffuse);

double Two_Sphere_R (struct measure_type m,
		     double UR1, double URU, double UT1, double UTU);

double Two_Sphere_T (struct measure_type m,
		     double UR1, double URU, double UT1, double UTU);

void Set_Calc_State (struct measure_type m, struct invert_type r);

void Get_Calc_State (struct measure_type *m, struct invert_type *r);

boolean_type Same_Calc_State (struct measure_type m, struct invert_type r);

boolean_type Valid_Grid (struct measure_type m, search_type s);

void Allocate_Grid (search_type s);

void Fill_Grid (struct measure_type m, struct invert_type r);

void Near_Grid_Points (double r, double t, search_type s, int *i_min,
		       int *j_min);

void Grid_ABG (int i, int j, guess_type * guess);

double Find_AG_fn (double x[]);

double Find_AB_fn (double x[]);

double Find_Ba_fn (double x);

double Find_Bs_fn (double x);

double Find_A_fn (double x);

double Find_B_fn (double x);

double Find_G_fn (double x);

double Find_BG_fn (double x[]);

double Find_BsG_fn (double x[]);

double Find_BaG_fn (double x[]);

void Fill_BG_Grid (struct measure_type m, struct invert_type r);

void Fill_BsG_Grid (struct measure_type m, struct invert_type r);

void Fill_BaG_Grid (struct measure_type m, struct invert_type r);

void Calculate_Distance_With_Corrections (double UR1, double UT1,
					  double Rc, double Tc,
					  double URU, double UTU,
					  double *M_R, double *M_T,
					  double *dev);

void Calculate_Distance (double *M_R, double *M_T, double *deviation);

double Calculate_Grid_Distance (int i, int j);

void abg_distance (double a, double b, double g, guess_type * guess);

double maxloss (double f);

void Max_Light_Loss (struct measure_type m, struct invert_type r,
		     double *ur1_loss, double *ut1_loss);
