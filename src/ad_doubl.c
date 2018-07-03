
#include <math.h>
#include <float.h>
#include "nr_util.h"
#include "ad_matrx.h"
#include "ad_globl.h"
#include "ad_doubl.h"


static void
Star_Multiply (int n, double **A, double **B, double **C)
{
  Right_Diagonal_Multiply (n, A, twoaw, C);
  Matrix_Multiply (n, C, B, C);
}



static void
Star_One_Minus (int n, double **A)
{
  int i, j;

  for (i = 1; i <= n; i++)
    {
      for (j = 1; j <= n; j++)
	A[i][j] *= -1;
      A[i][i] += 1.0 / twoaw[i];
    }

}




static void
Basic_Add_Layers (int n,
		  double **R10, double **T01,
		  double **R12, double **R21, double **T12, double **T21,
		  double **R20, double **T02, double **a, double **b)
{
  Star_Multiply (n, R10, R12, a);
  Star_One_Minus (n, a);
  Left_Inverse_Multiply (n, a, T12, b);



  Matrix_Multiply (n, b, R10, a);

  Star_Multiply (n, a, T21, a);

  Matrix_Sum (n, R21, a, R20);
  Copy_Matrix (n, T01, a);
  Matrix_Multiply (n, b, a, T02);
}




static void
Basic_Add_Layers_With_Sources (int n,
			       double **R10, double **T01,
			       double **R12, double **R21, double **T12,
			       double **T21, double **R20, double **T02,
			       double **J01, double **J12, double **J21,
			       double **J02, double **a, double **b)
{
  Star_Multiply (n, R10, R12, a);
  Star_One_Minus (n, a);
  Left_Inverse_Multiply (n, a, T12, b);



  Matrix_Multiply (n, b, R10, a);

  Star_Multiply (n, a, T21, a);

  Matrix_Sum (n, R21, a, R20);
  Copy_Matrix (n, T01, a);
  Matrix_Multiply (n, b, a, T02);

  Star_Multiply (n, R10, J21, a);
  Matrix_Sum (n, J01, a, a);
  Matrix_Multiply (n, b, a, J02);


  Matrix_Sum (n, J02, J12, J02);
}




void
Add (int n,
     double **R01, double **R10, double **T01, double **T10,
     double **R12, double **R21, double **T12, double **T21,
     double **R02, double **R20, double **T02, double **T20)
{

  double **a, **b;

  a = dmatrix (1, n, 1, n);
  b = dmatrix (1, n, 1, n);



  Basic_Add_Layers (n, R10, T01, R12, R21, T12, T21, R20, T02, a, b);
  Basic_Add_Layers (n, R12, T21, R10, R01, T10, T01, R02, T20, a, b);


  free_dmatrix (a, 1, n, 1, n);
  free_dmatrix (b, 1, n, 1, n);
}




void
Add_With_Sources (int n,
		  double **R01, double **R10, double **T01, double **T10,
		  double **J01, double **J10, double **R12, double **R21,
		  double **T12, double **T21, double **J12, double **J21,
		  double **R02, double **R20, double **T02, double **T20,
		  double **J02, double **J20)
{

  double **a, **b;

  a = dmatrix (1, n, 1, n);
  b = dmatrix (1, n, 1, n);



  Basic_Add_Layers_With_Sources (n, R10, T01, R12, R21, T12, T21, R20, T02,
				 J01, J12, J21, J02, a, b);
  Basic_Add_Layers_With_Sources (n, R12, T21, R10, R01, T10, T01, R02, T20,
				 J21, J10, J01, J20, a, b);


  free_dmatrix (a, 1, n, 1, n);
  free_dmatrix (b, 1, n, 1, n);
}




void
Add_Homogeneous (int n,
		 double **R01, double **T01,
		 double **R12, double **T12, double **R02, double **T02)
{

  double **a, **b;

  a = dmatrix (1, n, 1, n);
  b = dmatrix (1, n, 1, n);



  Basic_Add_Layers (n, R01, T01, R12, R12, T12, T12, R02, T02, a, b);


  free_dmatrix (a, 1, n, 1, n);
  free_dmatrix (b, 1, n, 1, n);
}




void
Double_Once (int n, double **R, double **T)
{

  double **a, **b;

  a = dmatrix (1, n, 1, n);
  b = dmatrix (1, n, 1, n);


  Basic_Add_Layers (n, R, T, R, R, T, T, R, T, a, b);

  free_dmatrix (a, 1, n, 1, n);
  free_dmatrix (b, 1, n, 1, n);
}




void
Double_Until (int n, double **r, double **t, double start, double end)
{
  if (end == HUGE_VAL)
    {
      Double_Until_Infinite (n, r, t);
      return;
    }

  {

    double **a, **b;

    a = dmatrix (1, n, 1, n);
    b = dmatrix (1, n, 1, n);


    while (fabs (end - start) > 0.00001 && end > start)
      {
	Basic_Add_Layers (n, r, t, r, r, t, t, r, t, a, b);
	start *= 2;
      }

    free_dmatrix (a, 1, n, 1, n);
    free_dmatrix (b, 1, n, 1, n);
  }
}




void
Double_Until_Infinite (int n, double **r, double **t)
{
  double oldutu, UTU, UT1;


  double **a, **b;

  a = dmatrix (1, n, 1, n);
  b = dmatrix (1, n, 1, n);


  UTU = 0.0;
  do
    {
      oldutu = UTU;
      Basic_Add_Layers (n, r, t, r, r, t, t, r, t, a, b);
      URU_and_UR1 (n, 1.0, t, &UTU, &UT1);
    }
  while (fabs (UTU - oldutu) >= 0.000001);


  free_dmatrix (a, 1, n, 1, n);
  free_dmatrix (b, 1, n, 1, n);

}




void
Between (int n,
	 double **R01, double **R10, double **T01, double **T10,
	 double **R12, double **R21, double **T12, double **T21,
	 double **Lup, double **Ldown)
{

  double **a, **b;

  a = dmatrix (1, n, 1, n);
  b = dmatrix (1, n, 1, n);



  Star_Multiply (n, R10, R12, a);
  Star_One_Minus (n, a);
  Right_Inverse_Multiply (n, a, T01, Ldown);

  Star_Multiply (n, R12, R10, a);
  Star_One_Minus (n, a);
  Right_Inverse_Multiply (n, a, R12, b);
  Star_Multiply (n, b, T01, Lup);


  free_dmatrix (a, 1, n, 1, n);
  free_dmatrix (b, 1, n, 1, n);
}
