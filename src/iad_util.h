

double What_Is_B (struct AD_slab_type slab, double Tc);

void Estimate_RT (struct measure_type m, struct invert_type r, double *rt,
		  double *tt, double *rd, double *rc, double *td, double *tc);

double a2acalc (double a);

double acalc2a (double acalc);

double g2gcalc (double g);

double gcalc2g (double gcalc);

double b2bcalc (double b);

double bcalc2b (double bcalc);

void twoprime (double a, double b, double g, double *ap, double *bp);

void twounprime (double ap, double bp, double g, double *a, double *b);

void abgg2ab (double a1, double b1, double g1, double g2, double *a2,
	      double *b2);

void abgb2ag (double a1, double b1, double b2, double *a2, double *g2);

void quick_guess (struct measure_type m, struct invert_type r, double *a,
		  double *b, double *g);

void Set_Debugging (unsigned long debug_level);

int Debug (unsigned long mask);

void Print_Invert_Type (struct invert_type r);

void Print_Measure_Type (struct measure_type m);
