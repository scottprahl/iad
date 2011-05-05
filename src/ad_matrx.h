

void Copy_Matrix (int n, double **A, double **B);

void One_Minus (int n, double **A);

void Transpose_Matrix (int n, double **a);

void Diagonal_To_Matrix (int n, double *Diag, double **Mat);

void Right_Diagonal_Multiply (int n, double **A, double *B, double **C);

void Left_Diagonal_Multiply (int n, double *A, double **B, double **C);

void Matrix_Multiply (int n, double **A, double **B, double **C);

void Matrix_Sum (int n, double **A, double **B, double **C);

void Solve (int n, double **A, double *B, int *ipvt);

void Decomp (int n, double **A, double *condition, int *ipvt);

void Matrix_Inverse (int n, double **A, double **Ainv);

void Left_Inverse_Multiply (int n, double **D, double **C, double **A);

void Right_Inverse_Multiply (int n, double **D, double **C, double **A);
