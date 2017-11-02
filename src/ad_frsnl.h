

double Cos_Critical_Angle (double ni, double nt);

double Cos_Snell (double n_i, double mu_i, double n_t);

void Absorbing_Glass_RT (double n_i, double n_g, double n_t, double mu_i,
			 double b, double *r, double *t);

void Sp_mu_RT (double n_top, double n_slab, double n_bottom,
	       double tau_top, double tau_slab, double tau_bottom, double mu,
	       double *r, double *t);

void Sp_mu_RT_Flip (int flip, double n_top, double n_slab, double n_bottom,
		    double tau_top, double tau_slab, double tau_bottom,
		    double mu, double *r, double *t);

double Diffuse_Glass_R (double nair, double nslide, double nslab);

double Glass (double n_i, double n_g, double n_t, double mu_i);
