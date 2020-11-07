
#include <math.h>
#include <float.h>
#include <stdio.h>
#include "nr_util.h"
#include "ad_globl.h"
#include "ad_prime.h"
#include "ad_matrx.h"
#include "ad_prime.h"
#include "ad_layers.h"


static void
PrintTestResults (int test, int cas,
		  double aUR1, double aUT1, double aURU, double aUTU,
		  double bUR1, double bUT1, double bURU, double bUTU)
{
  printf ("\nTest:%d.%d\n", test, cas);
  printf ("            truth        layers\n");
  printf ("UR1     %10.5f    %10.5f\n", aUR1, bUR1);
  printf ("UT1     %10.5f    %10.5f\n", aUT1, bUT1);
  printf ("URU     %10.5f    %10.5f\n", aURU, bURU);
  printf ("UTU     %10.5f    %10.5f\n", aUTU, bUTU);
}



int
main (int argc, char **argv)
{
  double aUR1, aURU, aUT1, aUTU, bUR1, bURU, bUT1, bUTU;
  double dUR1, dUT1, dURU, dUTU, cUR1, cUT1, cURU, cUTU;

  double a[15], b[15], g[15];
  int i;
  struct AD_slab_type slab;
  int N = 32;



  a[0] = 0.0;
  b[0] = 0.1;
  g[0] = 0.875;
  RT_Layers (N, 1.0, 1.0, 1.0, 1, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (1, 1, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  ez_RT (N, 1.0, 1.0, 1.0, a[0], b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.0, 1.0, 1.0, 1, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (1, 1, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  a[0] = 0.0;
  b[0] = 0.1;
  g[0] = 0.875;
  ez_RT (N, 1.4, 1.0, 1.0, a[0], b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.4, 1.0, 1.0, 1, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (1, 2, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  a[0] = 0.0;
  b[0] = 0.1;
  g[0] = 0.875;
  ez_RT (N, 1.4, 1.5, 1.5, a[0], b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.4, 1.5, 1.5, 1, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (1, 3, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  a[0] = 0.0;
  b[0] = 0.1;
  g[0] = 0.875;
  ez_RT (N, 1.4, 1.5, 1.6, a[0], b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.4, 1.5, 1.6, 1, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (1, 4, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  a[0] = 0.5;
  b[0] = 0.1;
  g[0] = 0.875;
  ez_RT (N, 1.0, 1.0, 1.0, a[0], b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.0, 1.0, 1.0, 1, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (1, 5, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  a[0] = 0.5;
  b[0] = 0.1;
  g[0] = 0.875;
  ez_RT (N, 1.4, 1.0, 1.0, a[0], b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.4, 1.0, 1.0, 1, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (1, 6, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  a[0] = 0.5;
  b[0] = 0.1;
  g[0] = 0.875;
  ez_RT (N, 1.4, 1.5, 1.5, a[0], b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.4, 1.5, 1.5, 1, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (1, 7, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  a[0] = 0.5;
  b[0] = 0.1;
  g[0] = 0.875;
  ez_RT (N, 1.4, 1.5, 1.6, a[0], b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.4, 1.5, 1.6, 1, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (1, 8, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);




  for (i = 0; i < 2; i++)
    {
      a[i] = 0.5;
      b[i] = 0.05;
      g[i] = 0.875;
    }
  ez_RT (N, 1.0, 1.0, 1.0, a[0], 2 * b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.0, 1.0, 1.0, 2, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (2, 1, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  ez_RT (N, 1.4, 1.0, 1.0, a[0], 2 * b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.4, 1.0, 1.0, 2, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (2, 2, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  ez_RT (N, 1.4, 1.5, 1.5, a[0], 2 * b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.4, 1.5, 1.5, 2, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (2, 3, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  ez_RT (N, 1.4, 1.5, 1.6, a[0], 2 * b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.4, 1.5, 1.6, 2, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (2, 4, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  for (i = 0; i < 5; i++)
    {
      a[i] = 0.5;
      b[i] = 0.02;
      g[i] = 0.875;
    }
  ez_RT (N, 1.0, 1.0, 1.0, a[0], 5 * b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.0, 1.0, 1.0, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (2, 5, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  ez_RT (N, 1.4, 1.0, 1.0, a[0], 5 * b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.4, 1.0, 1.0, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (2, 6, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  ez_RT (N, 1.4, 1.5, 1.5, a[0], 5 * b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.4, 1.5, 1.5, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (2, 7, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  ez_RT (N, 1.4, 1.5, 1.6, a[0], 5 * b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.4, 1.5, 1.6, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (2, 8, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);




  for (i = 0; i < 5; i++)
    {
      a[i] = 0.5;
      b[i] = 0.02;
      g[i] = 0.875;
    }

  b[0] = 0.0;
  ez_RT (N, 1.0, 1.0, 1.0, a[0], 4 * b[1], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.0, 1.0, 1.0, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (3, 1, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  ez_RT (N, 1.4, 1.5, 1.6, a[0], 4 * b[1], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.4, 1.5, 1.6, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (3, 2, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  for (i = 0; i < 5; i++)
    {
      a[i] = 0.5;
      b[i] = 0.02;
      g[i] = 0.875;
    }

  b[4] = 0.0;
  ez_RT (N, 1.0, 1.0, 1.0, a[0], 4 * b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.0, 1.0, 1.0, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (3, 3, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  ez_RT (N, 1.4, 1.5, 1.6, a[0], 4 * b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.4, 1.5, 1.6, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (3, 4, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  for (i = 0; i < 5; i++)
    {
      a[i] = 0.5;
      b[i] = 0.02;
      g[i] = 0.875;
    }

  b[3] = 0.0;
  ez_RT (N, 1.0, 1.0, 1.0, a[0], 4 * b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.0, 1.0, 1.0, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (3, 5, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  ez_RT (N, 1.4, 1.5, 1.6, a[0], 4 * b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.4, 1.5, 1.6, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (3, 6, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  for (i = 0; i < 5; i++)
    {
      a[i] = 0.5;
      b[i] = 0.02;
      g[i] = 0.875;
    }

  b[2] = 0.0;
  b[4] = 0.0;
  ez_RT (N, 1.0, 1.0, 1.0, a[0], 3 * b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.0, 1.0, 1.0, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (3, 7, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  ez_RT (N, 1.4, 1.5, 1.6, a[0], 3 * b[0], g[0], &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.4, 1.5, 1.6, 5, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (3, 8, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);




  slab.n_slab = 1.0;
  slab.n_top_slide = 1.0;
  slab.n_bottom_slide = 1.0;
  slab.b_top_slide = 0.2;
  slab.b_bottom_slide = 0.1;
  slab.a = 0.9;
  slab.b = 2.0;
  slab.g = 0.0;
  slab.phase_function = HENYEY_GREENSTEIN;
  slab.cos_angle = 1.0;

  a[0] = 0.0;
  b[0] = 0.2;
  g[0] = 0.0;
  a[1] = 0.9;
  b[1] = 2.0;
  g[1] = 0.0;
  a[2] = 0.0;
  b[2] = 0.1;
  g[2] = 0.0;

  RT (N, &slab, &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.0, 1.0, 1.0, 3, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (4, 1, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  slab.n_slab = 1.4;
  slab.n_top_slide = 1.4;
  slab.n_bottom_slide = 1.4;
  RT (N, &slab, &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.4, 1.4, 1.4, 3, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (4, 2, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);




  slab.n_slab = 1.0;
  slab.n_top_slide = 1.0;
  slab.n_bottom_slide = 1.0;
  slab.b_top_slide = 0.02;
  slab.b_bottom_slide = 0.02;
  slab.a = 0.5;
  slab.b = 0.06;
  slab.g = 0.875;
  slab.phase_function = HENYEY_GREENSTEIN;
  slab.cos_angle = 1.0;


  for (i = 0; i < 7; i++)
    {
      a[i] = 0.5;
      b[i] = 0.02;
      g[i] = 0.875;
    }

  a[0] = 0.0;
  b[2] = 0.0;
  b[4] = 0.0;
  a[6] = 0.0;

  RT (N, &slab, &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.0, 1.0, 1.0, 7, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (5, 1, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);

  slab.n_slab = 1.4;
  slab.n_top_slide = 1.4;
  slab.n_bottom_slide = 1.4;
  RT (N, &slab, &aUR1, &aUT1, &aURU, &aUTU);
  RT_Layers (N, 1.4, 1.4, 1.4, 7, a, b, g, &bUR1, &bUT1, &bURU, &bUTU);
  PrintTestResults (5, 2, aUR1, aUT1, aURU, aUTU, bUR1, bUT1, bURU, bUTU);



  for (i = 0; i < 5; i++)
    {
      a[i] = 0.5;
      b[i] = 0.1;
      g[i] = 0.875;
    }

  a[0] = 0.1;
  a[1] = 0.4;
  a[2] = 0.5;
  a[3] = 0.3;
  b[4] = 3;
  a[4] = 0.99;

  RT_Layers_All (N, 1.4, 1.4, 1.4, 5, a, b, g,
		 &aUR1, &aUT1, &aURU, &aUTU, &bUR1, &bUT1, &bURU, &bUTU);

  for (i = 0; i < 5; i++)
    {
      a[i] = 0.5;
      b[i] = 0.1;
      g[i] = 0.875;
    }
  a[0] = 0.99;
  b[0] = 3;
  a[1] = 0.4;
  a[2] = 0.5;
  a[3] = 0.3;
  a[4] = 0.1;

  RT_Layers_All (N, 1.4, 1.4, 1.4, 5, a, b, g,
		 &cUR1, &cUT1, &cURU, &cUTU, &dUR1, &dUT1, &dURU, &dUTU);
  PrintTestResults (6, 1, aUR1, aUT1, aURU, aUTU, dUR1, dUT1, dURU, dUTU);
  PrintTestResults (6, 2, bUR1, bUT1, bURU, bUTU, cUR1, cUT1, cURU, cUTU);


  return 0;
}
