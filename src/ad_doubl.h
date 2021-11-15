

void Add (int n,
	  double **R01, double **R10, double **T01, double **T10,
	  double **R12, double **R21, double **T12, double **T21,
	  double **R02, double **R20, double **T02, double **T20);

void Add_With_Sources (int n,
		       double **R01, double **R10, double **T01, double **T10,
		       double **J01, double **J10, double **R12, double **R21,
		       double **T12, double **T21, double **J12, double **J21,
		       double **R02, double **R20, double **T02, double **T20,
		       double **J02, double **J20);

void Add_Homogeneous (int n,
		      double **R01, double **T01,
		      double **R12, double **T12, double **R02, double **T02);

void Double_Once (int n, double **R, double **T);

void Double_Until (int n, double **r, double **t, double start, double end);

void Double_Until_Infinite (int n, double **r, double **t);

void Between (int n,
	      double **R01, double **R10, double **T01, double **T10,
	      double **R12, double **R21, double **T12, double **T21,
	      double **Lup, double **Ldown);
