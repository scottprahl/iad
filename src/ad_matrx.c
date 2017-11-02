

#include<stddef.h>
#include<math.h>
#include"ad_globl.h"
#include"ad_matrx.h"
#include"nr_util.h"



void
Copy_Matrix (int n, double **A, double **B)
{
  double *a_ptr, *b_ptr, *a_last;

  a_last = &A[n][n];
  a_ptr = &A[1][1];
  b_ptr = &B[1][1];

  while (a_ptr <= a_last)
    *b_ptr++ = *a_ptr++;
}




void
One_Minus (int n, double **A)
{
  int i, j;

  for (i = 1; i <= n; i++)
    {
      for (j = 1; j <= n; j++)
	A[i][j] *= -1;
      A[i][i] += 1.0;
    }

}




void
Transpose_Matrix (int n, double **a)
{
  int i, j;
  double swap;

  for (i = 1; i <= n; i++)
    {
      for (j = i + 1; j <= n; j++)
	{
	  swap = a[i][j];
	  a[i][j] = a[j][i];
	  a[j][i] = swap;
	}
    }
}




void
Diagonal_To_Matrix (int n, double *Diag, double **Mat)
{
  int i, j;

  for (i = 1; i <= n; i++)
    {
      for (j = 1; j <= n; j++)
	Mat[i][j] = 0.0;
      Mat[i][i] = Diag[i];
    }
}




void
Right_Diagonal_Multiply (int n, double **A, double *B, double **C)
{
  int i, j;

  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      C[i][j] = A[i][j] * B[j];
}




void
Left_Diagonal_Multiply (int n, double *A, double **B, double **C)
{
  int i, j;

  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      C[i][j] = A[i] * B[i][j];
}




void
Matrix_Multiply (int n, double **A, double **B, double **C)
{

  double *a_ptr, *a_start;
  double *b_start, *b_last;
  double *c_start, *c_very_last, *c_ptr;
  double *D;
  double *d_start, *d_last;
  register double t, *d_ptr, *b_ptr;
  ptrdiff_t row;



  if (n <= 0)
    {
      AD_error ("Non-positive dimension passed to Matrix_Multiply");
    }
  else if (n == 1)
    {
      C[1][1] = A[1][1] * B[1][1];
      return;
    }



  D = dvector (1, n);



  a_start = &A[1][1];
  b_last = &B[n][1];
  row = &A[2][1] - a_start;
  c_very_last = &C[n][n];
  d_start = &D[1];
  d_last = &D[n];




  for (c_start = &C[1][1]; c_start <= c_very_last; c_start += row)
    {
      a_ptr = a_start;

      d_ptr = d_start;
      while (d_ptr <= d_last)
	*d_ptr++ = 0.0;



      for (b_start = &B[1][1]; b_start <= b_last; b_start += row)
	{
	  t = *a_ptr++;
	  b_ptr = b_start;
	  d_ptr = d_start;
	  while (d_ptr <= d_last)
	    *d_ptr++ += t * (*b_ptr++);
	}

      d_ptr = d_start;
      c_ptr = c_start;
      while (d_ptr <= d_last)
	*c_ptr++ = *d_ptr++;


      a_start += row;
    }



  free_dvector (D, 1, n);


}




void
Matrix_Sum (int n, double **A, double **B, double **C)
{
  int i, j;

  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      C[i][j] = A[i][j] + B[i][j];
}




void
Solve (int n, double **A, double *B, int *ipvt)
{
  int i, k, m;
  double t;


  for (k = 1; k < n; k++)
    {
      m = ipvt[k];
      t = B[m];
      B[m] = B[k];
      B[k] = t;
      for (i = k + 1; i <= n; i++)
	B[i] += A[i][k] * t;
    }



  for (k = n; k > 1; k--)
    {
      B[k] /= A[k][k];
      t = -B[k];
      for (i = 1; i < k; i++)
	B[i] += A[i][k] * t;
    }

  B[1] /= A[1][1];



}




void
Decomp (int n, double **A, double *condition, int *ipvt)
{
  double t, anorm;
  int i, j, k, m;


  ipvt[n] = 1;

  if (n == 1)
    {
      if (A[1][1] == 0)
	{
	  AD_error ("1 X 1 Matrix is Singular --- i.e. zero");
	  return;
	}
    }



  anorm = 0.0;
  for (j = 1; j <= n; j++)
    {
      t = 0.0;
      for (i = 1; i <= n; i++)
	t += fabs (A[i][j]);
      if (t > anorm)
	anorm = t;
    }



  for (k = 1; k < n; k++)
    {

      m = k;
      for (i = k + 1; i <= n; i++)
	if (fabs (A[i][k]) > fabs (A[m][k]))
	  m = i;

      ipvt[k] = m;
      if (m != k)
	ipvt[n] *= -1;
      t = A[m][k];
      A[m][k] = A[k][k];
      A[k][k] = t;


      if (t == 0)
	continue;



      for (i = k + 1; i <= n; i++)
	A[i][k] /= -t;



      for (j = k + 1; j <= n; j++)
	{
	  t = A[m][j];
	  A[m][j] = A[k][j];
	  A[k][j] = t;
	  if (t == 0)
	    continue;
	  for (i = k + 1; i <= n; i++)
	    A[i][j] += A[i][k] * t;
	}


    }




  *condition = 1.0;
  for (k = 1; k <= n; k++)
    {
      if (A[k][k] == 0.0)
	{
	  *condition = 1e32;
	  return;
	}
    }


}




void
Matrix_Inverse (int n, double **A, double **Ainv)
{
  int *ipvt;
  int i, j;
  double *work;
  double condition;

  ipvt = ivector (1, n);
  Decomp (n, A, &condition, ipvt);
  if (condition == (condition + 1) || condition == 1e32)
    {
      free_ivector (ipvt, 1, n);
      AD_error ("Singular Matrix ... failed in Inverse_Multiply\n");
    }
  work = dvector (1, n);
  for (i = 1; i <= n; i++)
    {
      for (j = 1; j <= n; j++)
	work[j] = 0.0;
      work[i] = 1.0;
      Solve (n, A, work, ipvt);
      for (j = 1; j <= n; j++)
	Ainv[j][i] = work[j];
    }

  free_dvector (work, 1, n);
  free_ivector (ipvt, 1, n);
}





void
Left_Inverse_Multiply (int n, double **D, double **C, double **A)
{
  int *ipvt;
  int i, j;
  double *work;
  double condition;

  Transpose_Matrix (n, D);
  ipvt = ivector (1, n);
  Decomp (n, D, &condition, ipvt);


  if (condition == (condition + 1) || condition == 1e32)
    {
      free_ivector (ipvt, 1, n);
      AD_error ("Singular Matrix ... failed in Left_Inverse_Multiply\n");
    }

  work = dvector (1, n);
  for (i = 1; i <= n; i++)
    {

      for (j = 1; j <= n; j++)
	work[j] = C[i][j];
      Solve (n, D, work, ipvt);
      for (j = 1; j <= n; j++)
	A[i][j] = work[j];
    }

  free_dvector (work, 1, n);
  free_ivector (ipvt, 1, n);
}




void
Right_Inverse_Multiply (int n, double **D, double **C, double **A)
{
  int *ipvt;
  int i, j;
  double *work;
  double condition;

  ipvt = ivector (1, n);
  Decomp (n, D, &condition, ipvt);


  if (condition == (condition + 1) || condition == 1e32)
    {
      free_ivector (ipvt, 1, n);
      AD_error ("Singular Matrix ... failed in Right_Inverse_Multiply\n");
    }

  work = dvector (1, n);
  for (i = 1; i <= n; i++)
    {

      for (j = 1; j <= n; j++)
	work[j] = C[j][i];
      Solve (n, D, work, ipvt);
      for (j = 1; j <= n; j++)
	A[j][i] = work[j];
    }

  free_dvector (work, 1, n);
  free_ivector (ipvt, 1, n);
}
