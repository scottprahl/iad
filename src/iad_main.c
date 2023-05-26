


#define _CRT_SECURE_NO_WARNINGS
#define _CRT_NONSTDC_NO_WARNINGS

#define NO_SLIDES                 0
#define ONE_SLIDE_ON_TOP          1
#define TWO_IDENTICAL_SLIDES      2
#define ONE_SLIDE_ON_BOTTOM       3
#define ONE_SLIDE_NEAR_SPHERE     4
#define ONE_SLIDE_NOT_NEAR_SPHERE 5

#define MR_IS_ONLY_RD        1
#define MT_IS_ONLY_TD        2
#define NO_UNSCATTERED_LIGHT 3

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <ctype.h>

#include "ad_globl.h"
#include "ad_prime.h"
#include "iad_type.h"
#include "iad_pub.h"
#include "iad_io.h"
#include "iad_calc.h"
#include "iad_util.h"
#include "mygetopt.h"
#include "version.h"
#include "mc_lost.h"
#include "ad_frsnl.h"





static void
print_version (void)
{
  fprintf (stderr, "iad %s\n", Version);
  fprintf (stderr, "Copyright 2020 Scott Prahl, scott.prahl@oit.edu\n");
  fprintf (stderr, "          (see Applied Optics, 32:559-568, 1993)\n\n");
  fprintf (stderr,
	   "This is free software; see the source for copying conditions.\n");
  fprintf (stderr,
	   "There is no warranty; not even for MERCHANTABILITY or FITNESS.\n");
  fprintf (stderr, "FOR A PARTICULAR PURPOSE.\n");
  exit (EXIT_SUCCESS);
}



static void
print_usage (void)
{
  fprintf (stderr, "iad %s\n\n", Version);
  fprintf (stderr, "iad finds optical properties from measurements\n\n");
  fprintf (stderr, "Usage:  iad [options] input\n\n");
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "  -1 '# # # # #'   reflection sphere parameters \n");
  fprintf (stderr, "  -2 '# # # # #'   transmission sphere parameters \n");
  fprintf (stderr,
	   "     'sphere d, sample d, entrance d, detector d, wall r'\n");
  fprintf (stderr, "  -a #             use this albedo \n");
  fprintf (stderr, "  -A #             use this absorption coefficient \n");
  fprintf (stderr, "  -b #             use this optical thickness \n");
  fprintf (stderr, "  -B #             beam diameter \n");
  fprintf (stderr, "  -c #             fraction of unscattered refl in MR\n");
  fprintf (stderr,
	   "  -C #             fraction of unscattered trans in MT\n");
  fprintf (stderr, "  -d #             thickness of sample \n");
  fprintf (stderr, "  -D #             thickness of slide \n");
  fprintf (stderr, "  -e #             error tolerance (default 0.0001) \n");
  fprintf (stderr, "  -E #             optical depth (=mua*D) for slides\n");
  fprintf (stderr,
	   "  -f #             allow a fraction 0.0-1.0 of light to hit sphere wall first\n");
  fprintf (stderr, "  -F #             use this scattering coefficient \n");
  fprintf (stderr,
	   "  -F 'P lambda0 mus0 gamma'   mus=mus0*(lambda/lambda0)^gamma\n");
  fprintf (stderr,
	   "  -F 'R lambda0 musp0 gamma'  musp=musp0*(lambda/lambda0)^gamma\n");
  fprintf (stderr, "  -g #             scattering anisotropy (default 0) \n");
  fprintf (stderr,
	   "  -G #             type of boundary '0', '2', 't', 'b', 'n', 'f' \n");
  fprintf (stderr,
	   "                   '0' or '2'                --- number of slides\n");
  fprintf (stderr, "                   't' (top) or 'b' (bottom) \
--- one slide that is hit by light first\n");
  fprintf (stderr, "                   'n' (near) or 'f' (far)   \
--- one slide position relative to sphere\n");
  fprintf (stderr, "  -h               display help\n");
  fprintf (stderr,
	   "  -i #             light is incident at this angle in degrees\n");
  fprintf (stderr, "  -M #             number of Monte Carlo iterations\n");
  fprintf (stderr,
	   "  -n #             specify index of refraction of slab\n");
  fprintf (stderr,
	   "  -N #             specify index of refraction of slides\n");
  fprintf (stderr,
	   "  -o filename      explicitly specify filename for output\n");
  fprintf (stderr,
	   "  -p #             # of Monte Carlo photons (default 100000)\n");
  fprintf (stderr,
	   "                   a negative number is max MC time in milliseconds\n");
  fprintf (stderr,
	   "  -q #             number of quadrature points (default=8)\n");
  fprintf (stderr, "  -r #             total reflection measurement\n");
  fprintf (stderr,
	   "  -R #             actual reflectance for 100%% measurement \n");
  fprintf (stderr, "  -S #             number of spheres used\n");
  fprintf (stderr, "  -t #             total transmission measurement\n");
  fprintf (stderr,
	   "  -T #             actual transmission for 100%% measurement \n");
  fprintf (stderr,
	   "  -u #             unscattered transmission measurement\n");
  fprintf (stderr, "  -v               version information\n");
  fprintf (stderr,
	   "  -V 0             verbosity low --- no output to stderr\n");
  fprintf (stderr, "  -V 1             verbosity moderate \n");
  fprintf (stderr, "  -V 2             verbosity high\n");
  fprintf (stderr, "  -x #             set debugging level\n");
  fprintf (stderr, "  -X               dual beam configuration\n");
  fprintf (stderr, "  -z               do forward calculation\n");
  fprintf (stderr, "Examples:\n");
  fprintf (stderr,
	   "  iad file.rxt              Results will be put in file.txt\n");
  fprintf (stderr, "  iad file                  Same as above\n");
  fprintf (stderr, "  iad -c 0.9 file.rxt       \
Assume M_R includes 90%% of unscattered reflectance\n");
  fprintf (stderr, "  iad -C 0.8 file.rxt       \
Assume M_T includes 80%% of unscattered transmittance\n");
  fprintf (stderr,
	   "  iad -e 0.0001 file.rxt    Better convergence to R & T values\n");
  fprintf (stderr,
	   "  iad -f 1.0 file.rxt       All light hits reflectance sphere wall first\n");
  fprintf (stderr, "  iad -o out file.rxt       Calculated values in out\n");
  fprintf (stderr,
	   "  iad -r 0.3                R_total=0.3, b=inf, find albedo\n");
  fprintf (stderr,
	   "  iad -r 0.3 -t 0.4         R_total=0.3, T_total=0.4, find a,b,g\n");
  fprintf (stderr,
	   "  iad -r 0.3 -t 0.4 -n 1.5  R_total=0.3, T_total=0.4, n=1.5, find a,b\n");
  fprintf (stderr,
	   "  iad -r 0.3 -t 0.4         R_total=0.3, T_total=0.4, find a,b\n");
  fprintf (stderr, "  iad -p 1000 file.rxt      Only 1000 photons\n");
  fprintf (stderr,
	   "  iad -p -100 file.rxt      Allow only 100ms per iteration\n");
  fprintf (stderr, "  iad -q 4 file.rxt         Four quadrature points\n");
  fprintf (stderr, "  iad -M 0 file.rxt         No MC    (iad)\n");
  fprintf (stderr,
	   "  iad -M 1 file.rxt         MC once  (iad -> MC -> iad)\n");
  fprintf (stderr,
	   "  iad -M 2 file.rxt         MC twice (iad -> MC -> iad -> MC -> iad)\n");
  fprintf (stderr, "  iad -M 0 -q 4 file.rxt    Fast and crude conversion\n");
  fprintf (stderr,
	   "  iad -G t file.rxt         One top slide with properties from file.rxt\n");
  fprintf (stderr,
	   "  iad -G b -N 1.5 -D 1 file Use 1 bottom slide with n=1.5 and thickness=1\n");
  fprintf (stderr,
	   "  iad -x   1 file.rxt       Show sphere and MC effects\n");
  fprintf (stderr, "  iad -x   2 file.rxt       DEBUG_GRID\n");
  fprintf (stderr, "  iad -x   4 file.rxt       DEBUG_ITERATIONS\n");
  fprintf (stderr, "  iad -x   8 file.rxt       DEBUG_LOST_LIGHT\n");
  fprintf (stderr, "  iad -x  16 file.rxt       DEBUG_SPHERE_EFFECTS\n");
  fprintf (stderr, "  iad -x  32 file.rxt       DEBUG_BEST_GUESS\n");
  fprintf (stderr, "  iad -x  64 file.rxt       DEBUG_EVERY_CALC\n");
  fprintf (stderr, "  iad -x 128 file.rxt       DEBUG_SEARCH\n");
  fprintf (stderr, "  iad -x 255 file.rxt       All debugging output\n");
  fprintf (stderr,
	   "  iad -X -i 8 file.rxt      Dual beam spectrometer with 8 degree incidence\n\n");
  fprintf (stderr,
	   "  iad -z -a 0.9 -b 1 -i 45  Forward calc assuming 45 degree incidence\n\n");
  fprintf (stderr, "  apply iad x.rxt y.rxt     Process multiple files\n\n");
  fprintf (stderr, "Report bugs to <scott.prahl@oit.edu>\n\n");
  exit (EXIT_SUCCESS);
}



static char *
strdup_together (char *s, char *t)
{
  char *both;

  if (s == NULL)
    {
      if (t == NULL)
	return NULL;
      return strdup (t);
    }

  if (t == NULL)
    return strdup (s);

  both = malloc (strlen (s) + strlen (t) + 1);
  if (both == NULL)
    fprintf (stderr, "Could not allocate memory for both strings.\n");

  strcpy (both, s);
  strcat (both, t);
  return both;
}




static double
seconds_elapsed (clock_t start_time)
{
  clock_t finish_time = clock ();
  return (double) (finish_time - start_time) / CLOCKS_PER_SEC;
}



static void
print_error_legend (void)
{
  fprintf (stderr,
	   "----------------- Sorry, but ... errors encountered ---------------\n");
  fprintf (stderr, "   *  ==> Success          ");
  fprintf (stderr, "  0-9 ==> Monte Carlo Iteration\n");
  fprintf (stderr, "   R  ==> M_R is too big   ");
  fprintf (stderr, "   r  ==> M_R is too small\n");
  fprintf (stderr, "   T  ==> M_T is too big   ");
  fprintf (stderr, "   t  ==> M_T is too small\n");
  fprintf (stderr, "   U  ==> M_U is too big   ");
  fprintf (stderr, "   u  ==> M_U is too small\n");
  fprintf (stderr, "   !  ==> M_R + M_T > 1    ");
  fprintf (stderr, "   +  ==> Did not converge\n\n");
}





static char
what_char (int err)
{
  if (err == IAD_NO_ERROR)
    return '*';
  if (err == IAD_TOO_MANY_ITERATIONS)
    return '+';
  if (err == IAD_MR_TOO_BIG)
    return 'R';
  if (err == IAD_MR_TOO_SMALL)
    return 'r';
  if (err == IAD_MT_TOO_BIG)
    return 'T';
  if (err == IAD_MT_TOO_SMALL)
    return 't';
  if (err == IAD_MU_TOO_BIG)
    return 'U';
  if (err == IAD_MU_TOO_SMALL)
    return 'u';
  if (err == IAD_TOO_MUCH_LIGHT)
    return '!';
  return '?';
}

static void
print_dot (clock_t start_time, int err, int count, int points,
	   int final, int verbosity, int *any_error)
{
  static int counter = 0;

  counter++;

  if (err != IAD_NO_ERROR)
    *any_error = err;

  if (verbosity == 0)
    return;

  if (final == 99)
    fprintf (stderr, "%c", what_char (err));
  else
    {
      counter--;
      fprintf (stderr, "%1d\b", final % 10);
    }

  if (final == 99)
    {
      if (counter % 50 == 0)
	{
	  double rate = (seconds_elapsed (start_time) / points);
	  fprintf (stderr, "  %3d done (%5.2f s/pt)\n", points, rate);
	}
      else if (counter % 10 == 0)
	fprintf (stderr, " ");
    }

  fflush (stderr);
}


static void
Calculate_Mua_Musp (struct measure_type m,
		    struct invert_type r, double *musp, double *mua)
{
  if (r.b == HUGE_VAL)
    {
      if (r.a <= 1e-5)
	{
	  *musp = 0.0;
	  *mua = 1.0;
	  return;
	}
      if (r.default_mus != UNINITIALIZED)
	{
	  *musp = r.default_mus * (1 - r.g);
	  *mua = r.default_mus / r.a - r.default_mus;
	  return;
	}
      if (r.default_mua != UNINITIALIZED)
	{
	  *musp = (r.default_mua / (1 - r.a) - r.default_mua) * (1 - r.g);
	  *mua = r.default_mua;
	  return;
	}

      *musp = 1.0 - r.g;
      *mua = (1.0 - r.a) / r.a;
      return;
    }

  *musp = r.a * r.b / m.slab_thickness * (1.0 - r.g);
  *mua = (1 - r.a) * r.b / m.slab_thickness;
}


static void
calculate_coefficients (struct measure_type m,
			struct invert_type r,
			double *LR, double *LT, double *musp, double *mua)
{
  double delta;
  *LR = 0;
  *LT = 0;
  Calculate_Distance (LR, LT, &delta);
  Calculate_Mua_Musp (m, r, musp, mua);
}




static int
parse_string_into_array (char *s, double *a, int n)
{
  char *t, *last, *r;
  int i = 0;
  t = s;
  last = s + strlen (s);

  while (t < last)
    {


      r = t;
      while (*r != ' ' && *r != '\0')
	r++;
      *r = '\0';


      if (sscanf (t, "%lf", &(a[i])) == 0)
	return 1;
      i++;


      if (i == n)
	return 0;


      t = r + 1;
    }
  return 1;
}



static void
print_results_header (FILE * fp)
{
  fprintf (fp,
	   "#     \tMeasured \t   M_R   \tMeasured \t   M_T   \tEstimated\tEstimated\tEstimated");
  if (Debug (DEBUG_LOST_LIGHT))
    fprintf (fp,
	     "\t  Lost   \t  Lost   \t  Lost   \t  Lost   \t   MC    \t   IAD   \t  Error  ");
  fprintf (fp, "\n");

  fprintf (fp,
	   "##wave\t   M_R   \t   fit   \t   M_T   \t   fit   \t  mu_a   \t  mu_s'  \t    g    ");
  if (Debug (DEBUG_LOST_LIGHT))
    fprintf (fp,
	     "\t   UR1   \t   URU   \t   UT1   \t   UTU   \t    #    \t    #    \t  State  ");
  fprintf (fp, "\n");

  fprintf (fp,
	   "# [nm]\t  [---]  \t  [---]  \t  [---]  \t  [---]  \t  1/mm   \t  1/mm   \t  [---]  ");
  if (Debug (DEBUG_LOST_LIGHT))
    fprintf (fp,
	     "\t  [---]  \t  [---]  \t  [---]  \t  [---]  \t  [---]  \t  [---]  \t  [---]  ");
  fprintf (fp, "\n");
}



void
print_optical_property_result (FILE * fp,
			       struct measure_type m,
			       struct invert_type r,
			       double LR,
			       double LT,
			       double mu_a,
			       double mu_sp, int mc_iter, int line)
{
  if (m.lambda != 0)
    fprintf (fp, "%6.1f\t", m.lambda);
  else
    fprintf (fp, "%6d\t", line);

  if (mu_a >= 200)
    mu_a = 199.9999;
  if (mu_sp >= 1000)
    mu_sp = 999.9999;

  fprintf (fp, "% 9.4f\t% 9.4f\t", m.m_r, LR);
  fprintf (fp, "% 9.4f\t% 9.4f\t", m.m_t, LT);
  fprintf (fp, "% 9.4f\t", mu_a);
  fprintf (fp, "% 9.4f\t", mu_sp);
  fprintf (fp, "% 9.4f\t", r.g);

  if (Debug (DEBUG_LOST_LIGHT))
    {
      fprintf (fp, "% 9.4f\t% 9.4f\t", m.ur1_lost, m.uru_lost);
      fprintf (fp, "% 9.4f\t% 9.4f\t", m.ut1_lost, m.utu_lost);
      fprintf (fp, " %2d  \t", mc_iter);
      fprintf (fp, " %4d\t", r.iterations);
    }

  fprintf (fp, "# %c \n", what_char (r.error));
  fflush (fp);
}



int
main (int argc, char **argv)
{

  struct measure_type m;
  struct invert_type r;
  char *g_out_name = NULL;
  char c;

  long n_photons = 100000;
  int MC_iterations = 19;
  int any_error = 0;
  int process_command_line = 0;
  int params = 0;

  int cl_quadrature_points = UNINITIALIZED;
  int cl_verbosity = 2;

  double cl_forward_calc = UNINITIALIZED;
  double cl_default_a = UNINITIALIZED;
  double cl_default_g = UNINITIALIZED;
  double cl_default_b = UNINITIALIZED;
  double cl_default_mua = UNINITIALIZED;
  double cl_default_mus = UNINITIALIZED;
  double cl_tolerance = UNINITIALIZED;
  double cl_slide_OD = UNINITIALIZED;

  double cl_cos_angle = UNINITIALIZED;
  double cl_beam_d = UNINITIALIZED;
  double cl_sample_d = UNINITIALIZED;
  double cl_sample_n = UNINITIALIZED;
  double cl_slide_d = UNINITIALIZED;
  double cl_slide_n = UNINITIALIZED;
  double cl_slides = UNINITIALIZED;
  double cl_default_fr = UNINITIALIZED;
  double cl_rstd_t = UNINITIALIZED;
  double cl_rstd_r = UNINITIALIZED;
  double cl_rc_fraction = UNINITIALIZED;
  double cl_tc_fraction = UNINITIALIZED;

  double cl_search = UNINITIALIZED;
  double cl_mus0 = UNINITIALIZED;
  double cl_musp0 = UNINITIALIZED;
  double cl_mus0_pwr = UNINITIALIZED;
  double cl_mus0_lambda = UNINITIALIZED;

  double cl_UR1 = UNINITIALIZED;
  double cl_UT1 = UNINITIALIZED;
  double cl_Tc = UNINITIALIZED;

  double cl_method = UNINITIALIZED;
  double cl_num_spheres = UNINITIALIZED;
  double cl_sphere_one[5] = { UNINITIALIZED, UNINITIALIZED, UNINITIALIZED,
    UNINITIALIZED, UNINITIALIZED
  };
  double cl_sphere_two[5] = { UNINITIALIZED, UNINITIALIZED, UNINITIALIZED,
    UNINITIALIZED, UNINITIALIZED
  };

  clock_t start_time = clock ();
  char command_line_options[] =
    "?1:2:a:A:b:B:c:C:d:D:e:E:f:F:g:G:hi:n:N:M:o:p:q:r:R:s:S:t:T:u:vV:x:Xz";




  while ((c = my_getopt (argc, argv, command_line_options)) != EOF)
    {
      int n;
      char cc;

      switch (c)
	{

	case '1':
	  parse_string_into_array (optarg, cl_sphere_one, 5);
	  break;

	case '2':
	  parse_string_into_array (optarg, cl_sphere_two, 5);
	  break;

	case 'a':
	  cl_default_a = strtod (optarg, NULL);
	  break;

	case 'A':
	  cl_default_mua = strtod (optarg, NULL);
	  break;

	case 'b':
	  cl_default_b = strtod (optarg, NULL);
	  break;

	case 'B':
	  cl_beam_d = strtod (optarg, NULL);
	  break;

	case 'c':
	  cl_rc_fraction = strtod (optarg, NULL);
	  if (cl_rc_fraction < 0.0 || cl_rc_fraction > 1.0)
	    {
	      fprintf (stderr,
		       "required: 0 <= fraction of unscattered refl. in M_R <= 1\n");
	      exit (EXIT_SUCCESS);
	    }
	  break;

	case 'C':
	  cl_tc_fraction = strtod (optarg, NULL);
	  if (cl_tc_fraction < 0.0 || cl_tc_fraction > 1.0)
	    {
	      fprintf (stderr,
		       "required: 0 <= fraction of unscattered trans. in M_T <= 1\n");
	      exit (EXIT_SUCCESS);
	    }
	  break;

	case 'd':
	  cl_sample_d = strtod (optarg, NULL);
	  break;

	case 'D':
	  cl_slide_d = strtod (optarg, NULL);
	  break;

	case 'e':
	  cl_tolerance = strtod (optarg, NULL);
	  break;

	case 'E':
	  cl_slide_OD = strtod (optarg, NULL);
	  break;

	case 'f':
	  cl_default_fr = strtod (optarg, NULL);
	  break;

	case 'F':

	  if (isdigit (optarg[0]))
	    {
	      cl_default_mus = strtod (optarg, NULL);
	      break;
	    }


	  n =
	    sscanf (optarg, "%c %lf %lf %lf", &cc, &cl_mus0_lambda, &cl_mus0,
		    &cl_mus0_pwr);

	  if (n != 4 || (cc != 'P' && cc != 'R'))
	    {
	      fprintf (stderr,
		       "Screwy argument for -F option.  Try something like\n");
	      fprintf (stderr, " -F 1.0              for mus =1.0\n");
	      fprintf (stderr,
		       " -F 'P 500 1.0 -1.3' for mus =1.0*(lambda/500)^(-1.3)\n");
	      fprintf (stderr,
		       " -F 'R 500 1.0 -1.3' for mus'=1.0*(lambda/500)^(-1.3)\n");
	      exit (EXIT_FAILURE);
	    }

	  if (cc == 'R' || cc == 'r')
	    {
	      cl_musp0 = cl_mus0;
	      cl_mus0 = UNINITIALIZED;
	    }
	  break;

	case 'g':
	  cl_default_g = strtod (optarg, NULL);
	  break;

	case 'G':
	  if (optarg[0] == '0')
	    cl_slides = NO_SLIDES;
	  else if (optarg[0] == '2')
	    cl_slides = TWO_IDENTICAL_SLIDES;
	  else if (optarg[0] == 't' || optarg[0] == 'T')
	    cl_slides = ONE_SLIDE_ON_TOP;
	  else if (optarg[0] == 'b' || optarg[0] == 'B')
	    cl_slides = ONE_SLIDE_ON_BOTTOM;
	  else if (optarg[0] == 'n' || optarg[0] == 'N')
	    cl_slides = ONE_SLIDE_NEAR_SPHERE;
	  else if (optarg[0] == 'f' || optarg[0] == 'F')
	    cl_slides = ONE_SLIDE_NOT_NEAR_SPHERE;
	  else
	    {
	      fprintf (stderr, "Argument for -G option must be \n");
	      fprintf (stderr,
		       "    't' --- light always hits top slide first\n");
	      fprintf (stderr,
		       "    'b' --- light always hits bottom slide first\n");
	      fprintf (stderr,
		       "    'n' --- slide always closest to sphere\n");
	      fprintf (stderr,
		       "    'f' --- slide always farthest from sphere\n");
	      exit (EXIT_FAILURE);
	    }
	  break;

	case 'i':
	  cl_cos_angle = strtod (optarg, NULL);
	  if (cl_cos_angle < 0 || cl_cos_angle > 90)
	    fprintf (stderr,
		     "Incident angle must be between 0 and 90 degrees\n");
	  else
	    cl_cos_angle = cos (cl_cos_angle * 3.1415926535 / 180.0);
	  break;

	case 'M':
	  MC_iterations = (int) strtod (optarg, NULL);
	  break;

	case 'n':
	  cl_sample_n = strtod (optarg, NULL);
	  break;

	case 'N':
	  cl_slide_n = strtod (optarg, NULL);
	  break;

	case 'o':
	  g_out_name = strdup (optarg);
	  break;

	case 'p':
	  n_photons = (int) strtod (optarg, NULL);
	  break;

	case 'q':
	  cl_quadrature_points = (int) strtod (optarg, NULL);
	  if (cl_quadrature_points % 4 != 0)
	    {
	      fprintf (stderr,
		       "Number of quadrature points must be a multiple of 4\n");
	      exit (EXIT_FAILURE);
	    }
	  if ((cl_cos_angle != UNINITIALIZED)
	      && (cl_quadrature_points % 12 != 0))
	    {
	      fprintf (stderr,
		       "Quadrature must be 12, 24, 36,... for oblique incidence\n");
	      exit (EXIT_FAILURE);
	    }
	  break;

	case 'r':
	  cl_UR1 = strtod (optarg, NULL);
	  process_command_line = 1;
	  break;

	case 'R':
	  cl_rstd_r = strtod (optarg, NULL);
	  break;

	case 's':
	  cl_search = (int) strtod (optarg, NULL);
	  break;

	case 'S':
	  cl_num_spheres = (int) strtod (optarg, NULL);
	  break;

	case 't':
	  cl_UT1 = strtod (optarg, NULL);
	  process_command_line = 1;
	  break;

	case 'T':
	  cl_rstd_t = strtod (optarg, NULL);
	  break;

	case 'u':
	  cl_Tc = strtod (optarg, NULL);
	  process_command_line = 1;
	  break;

	case 'v':
	  print_version ();
	  break;

	case 'V':
	  cl_verbosity = strtod (optarg, NULL);
	  break;

	case 'x':
	  Set_Debugging ((int) strtod (optarg, NULL));
	  break;

	case 'X':
	  cl_method = COMPARISON;
	  break;

	case 'z':
	  cl_forward_calc = 1;
	  process_command_line = 1;
	  break;

	default:
	  fprintf (stderr, "unknown option '%c'\n", c);


	case 'h':
	case '?':
	  print_usage ();
	  break;
	}
    }

  argc -= optind;
  argv += optind;



  Initialize_Measure (&m);


  if (cl_cos_angle != UNINITIALIZED)
    {
      m.slab_cos_angle = cl_cos_angle;
      if (cl_quadrature_points == UNINITIALIZED)
	cl_quadrature_points = 12;

      if (cl_quadrature_points != 12 * (cl_quadrature_points / 12))
	{
	  fprintf (stderr,
		   "If you use the -i option to specify an oblique incidence angle, then\n");
	  fprintf (stderr,
		   "the number of quadrature points must be a multiple of 12\n");
	  exit (EXIT_SUCCESS);
	}
    }

  if (cl_sample_n != UNINITIALIZED)
    m.slab_index = cl_sample_n;

  if (cl_slide_n != UNINITIALIZED)
    {
      m.slab_bottom_slide_index = cl_slide_n;
      m.slab_top_slide_index = cl_slide_n;
    }

  if (cl_slide_OD != UNINITIALIZED)
    {
      m.slab_bottom_slide_b = cl_slide_OD;
      m.slab_top_slide_b = cl_slide_OD;
    }

  if (cl_sample_d != UNINITIALIZED)
    m.slab_thickness = cl_sample_d;

  if (cl_beam_d != UNINITIALIZED)
    m.d_beam = cl_beam_d;

  if (cl_slide_d != UNINITIALIZED)
    {
      m.slab_bottom_slide_thickness = cl_slide_d;
      m.slab_top_slide_thickness = cl_slide_d;
    }

  if (cl_slides == NO_SLIDES)
    {
      m.slab_bottom_slide_index = 1.0;
      m.slab_bottom_slide_thickness = 0.0;
      m.slab_top_slide_index = 1.0;
      m.slab_top_slide_thickness = 0.0;
    }

  if (cl_slides == ONE_SLIDE_ON_TOP || cl_slides == ONE_SLIDE_NEAR_SPHERE)
    {
      m.slab_bottom_slide_index = 1.0;
      m.slab_bottom_slide_thickness = 0.0;
    }

  if (cl_slides == ONE_SLIDE_ON_BOTTOM ||
      cl_slides == ONE_SLIDE_NOT_NEAR_SPHERE)
    {
      m.slab_top_slide_index = 1.0;
      m.slab_top_slide_thickness = 0.0;
    }

  if (cl_slides == ONE_SLIDE_NEAR_SPHERE ||
      cl_slides == ONE_SLIDE_NOT_NEAR_SPHERE)
    m.flip_sample = 1;
  else
    m.flip_sample = 0;

  if (cl_method != UNINITIALIZED)
    m.method = (int) cl_method;

  if (cl_rstd_t != UNINITIALIZED)
    m.rstd_t = cl_rstd_t;

  if (cl_rstd_r != UNINITIALIZED)
    m.rstd_r = cl_rstd_r;

  if (cl_sphere_one[4] != UNINITIALIZED)
    {
      double d_sample_r, d_entrance_r, d_detector_r;

      m.d_sphere_r = cl_sphere_one[0];
      d_sample_r = cl_sphere_one[1];
      d_entrance_r = cl_sphere_one[2];
      d_detector_r = cl_sphere_one[3];
      m.rw_r = cl_sphere_one[4];

      m.as_r =
	(d_sample_r / m.d_sphere_r / 2) * (d_sample_r / m.d_sphere_r / 2);
      m.ae_r =
	(d_entrance_r / m.d_sphere_r / 2) * (d_entrance_r / m.d_sphere_r / 2);
      m.ad_r =
	(d_detector_r / m.d_sphere_r / 2) * (d_detector_r / m.d_sphere_r / 2);

      m.aw_r = 1.0 - m.as_r - m.ae_r - m.ad_r;

      m.d_sphere_t = m.d_sphere_r;
      m.as_t = m.as_r;
      m.ae_t = m.ae_r;
      m.ad_t = m.ad_r;
      m.aw_t = m.aw_r;
      m.rw_t = m.rw_r;

      if (cl_num_spheres == UNINITIALIZED)
	m.num_spheres = 1;
    }

  if (cl_sphere_two[4] != UNINITIALIZED)
    {
      double d_sample_t, d_entrance_t, d_detector_t;

      m.d_sphere_t = cl_sphere_two[0];
      d_sample_t = cl_sphere_two[1];
      d_entrance_t = cl_sphere_two[2];
      d_detector_t = cl_sphere_two[3];
      m.rw_t = cl_sphere_two[4];

      m.as_t =
	(d_sample_t / m.d_sphere_t / 2) * (d_sample_t / m.d_sphere_t / 2);
      m.ae_t =
	(d_entrance_t / m.d_sphere_t / 2) * (d_entrance_t / m.d_sphere_t / 2);
      m.ad_t =
	(d_detector_t / m.d_sphere_t / 2) * (d_detector_t / m.d_sphere_t / 2);
      m.aw_t = 1.0 - m.as_t - m.ae_t - m.ad_t;

      if (cl_num_spheres == UNINITIALIZED)
	m.num_spheres = 2;
    }

  if (cl_num_spheres != UNINITIALIZED)
    {
      m.num_spheres = (int) cl_num_spheres;
      if (m.num_spheres > 0 && m.method == UNKNOWN)
	m.method = SUBSTITUTION;
    }

  if (cl_rc_fraction != UNINITIALIZED)
    m.fraction_of_rc_in_mr = cl_rc_fraction;

  if (cl_tc_fraction != UNINITIALIZED)
    m.fraction_of_tc_in_mt = cl_tc_fraction;

  if (cl_UR1 != UNINITIALIZED)
    m.m_r = cl_UR1;

  if (cl_UT1 != UNINITIALIZED)
    m.m_t = cl_UT1;

  if (cl_Tc != UNINITIALIZED)
    m.m_u = cl_Tc;

  if (cl_default_fr != UNINITIALIZED)
    m.f_r = cl_default_fr;



  Initialize_Result (m, &r);


  if (cl_quadrature_points != UNINITIALIZED)
    r.method.quad_pts = cl_quadrature_points;
  else
    r.method.quad_pts = 8;

  if (cl_default_a != UNINITIALIZED)
    r.default_a = cl_default_a;

  if (cl_default_mua != UNINITIALIZED)
    {
      r.default_mua = cl_default_mua;
      if (cl_sample_d != UNINITIALIZED)
	r.default_ba = cl_default_mua * cl_sample_d;
      else
	r.default_ba = cl_default_mua * m.slab_thickness;
    }

  if (cl_default_b != UNINITIALIZED)
    r.default_b = cl_default_b;

  if (cl_default_g != UNINITIALIZED)
    r.default_g = cl_default_g;

  if (cl_tolerance != UNINITIALIZED)
    {
      r.tolerance = cl_tolerance;
      r.MC_tolerance = cl_tolerance;
    }

  if (cl_musp0 != UNINITIALIZED)
    cl_mus0 =
      (r.default_g !=
       UNINITIALIZED) ? cl_musp0 / (1 - r.default_g) : cl_musp0;

  if (cl_mus0 != UNINITIALIZED && m.lambda != 0)
    cl_default_mus = cl_mus0 * pow (m.lambda / cl_mus0_lambda, cl_mus0_pwr);

  if (cl_default_mus != UNINITIALIZED)
    {
      r.default_mus = cl_default_mus;
      if (cl_sample_d != UNINITIALIZED)
	r.default_bs = cl_default_mus * cl_sample_d;
      else
	r.default_bs = cl_default_mus * m.slab_thickness;
    }

  if (cl_search != UNINITIALIZED)
    r.search = cl_search;



  if (cl_forward_calc != UNINITIALIZED)
    {

      if (cl_default_a == UNINITIALIZED)
	{

	  if (cl_default_mus == UNINITIALIZED)
	    r.a = 0;
	  else if (cl_default_mua == UNINITIALIZED)
	    r.a = 1;
	  else
	    r.a = cl_default_mus / (cl_default_mua + cl_default_mus);

	}
      else
	r.a = cl_default_a;


      if (cl_default_b == UNINITIALIZED)
	{

	  if (cl_sample_d == UNINITIALIZED)
	    r.b = HUGE_VAL;

	  else if (r.a == 0)
	    {
	      if (cl_default_mua == UNINITIALIZED)
		r.b = HUGE_VAL;
	      else
		r.b = cl_default_mua * cl_sample_d;
	    }
	  else
	    {
	      if (cl_default_mus == UNINITIALIZED)
		r.b = HUGE_VAL;
	      else
		r.b = cl_default_mus / r.a * cl_sample_d;
	    }
	}
      else
	r.b = cl_default_b;


      if (cl_default_g == UNINITIALIZED)
	r.g = 0;
      else
	r.g = cl_default_g;


      r.slab.a = r.a;
      r.slab.b = r.b;
      r.slab.g = r.g;

      {
	double mu_sp, mu_a, m_r, m_t;
	Calculate_MR_MT (m, r, MC_iterations, &m_r, &m_t);
	Calculate_Mua_Musp (m, r, &mu_sp, &mu_a);
	if (cl_verbosity > 0)
	  {
	    Write_Header (m, r, -1);
	    print_results_header (stdout);
	  }
	print_optical_property_result (stdout, m, r, m_r, m_t, mu_a, mu_sp, 0,
				       0);
      }


      return EXIT_SUCCESS;
    }


  if (argc > 1)
    {
      fprintf (stderr, "Only a single file can be processed at a time\n");
      fprintf (stderr, "try 'apply iad file1 file2 ... fileN'\n");
      exit (EXIT_FAILURE);
    }

  if (argc == 1 && strcmp (argv[0], "-") != 0)
    {
      int n;
      char *base_name, *rt_name;
      base_name = strdup (argv[0]);
      n = (int) (strlen (base_name) - strlen (".rxt"));

      if (n > 0 && strstr (base_name + n, ".rxt") != NULL)
	base_name[n] = '\0';

      rt_name = strdup_together (base_name, ".rxt");

      if (freopen (argv[0], "r", stdin) == NULL &&
	  freopen (rt_name, "r", stdin) == NULL)
	{
	  fprintf (stderr, "Could not open either '%s' or '%s'\n",
		   argv[0], rt_name);
	  exit (EXIT_FAILURE);
	}

      if (g_out_name == NULL)
	g_out_name = strdup_together (base_name, ".txt");

      free (rt_name);
      free (base_name);
      process_command_line = 0;
    }

  if (g_out_name != NULL)
    {
      if (freopen (g_out_name, "w", stdout) == NULL)
	{
	  fprintf (stderr, "Could not open file '%s' for output\n",
		   g_out_name);
	  exit (EXIT_FAILURE);
	}
    }




  if (process_command_line)
    {


      m.num_measures = 3;
      if (m.m_t == 0)
	m.num_measures--;
      if (m.m_u == 0)
	m.num_measures--;
      params = m.num_measures;

      if (m.num_measures == 3)
	{

	  struct AD_slab_type s;
	  s.n_slab = m.slab_index;
	  s.n_top_slide = m.slab_top_slide_index;
	  s.n_bottom_slide = m.slab_bottom_slide_index;
	  s.b_top_slide = m.slab_top_slide_b;
	  s.b_bottom_slide = m.slab_bottom_slide_b;
	  s.cos_angle = m.slab_cos_angle;
	  cl_default_b = What_Is_B (s, m.m_u);
	}



      {

	static int rt_total = 0;
	static int mc_total = 0;
	int mc_iter = 0;

	double ur1 = 0;
	double ut1 = 0;
	double uru = 0;
	double utu = 0;
	double mu_a = 0;
	double mu_sp = 0;
	double LR = 0;
	double LT = 0;

	rt_total++;



	Initialize_Result (m, &r);


	if (cl_quadrature_points != UNINITIALIZED)
	  r.method.quad_pts = cl_quadrature_points;
	else
	  r.method.quad_pts = 8;

	if (cl_default_a != UNINITIALIZED)
	  r.default_a = cl_default_a;

	if (cl_default_mua != UNINITIALIZED)
	  {
	    r.default_mua = cl_default_mua;
	    if (cl_sample_d != UNINITIALIZED)
	      r.default_ba = cl_default_mua * cl_sample_d;
	    else
	      r.default_ba = cl_default_mua * m.slab_thickness;
	  }

	if (cl_default_b != UNINITIALIZED)
	  r.default_b = cl_default_b;

	if (cl_default_g != UNINITIALIZED)
	  r.default_g = cl_default_g;

	if (cl_tolerance != UNINITIALIZED)
	  {
	    r.tolerance = cl_tolerance;
	    r.MC_tolerance = cl_tolerance;
	  }

	if (cl_musp0 != UNINITIALIZED)
	  cl_mus0 =
	    (r.default_g !=
	     UNINITIALIZED) ? cl_musp0 / (1 - r.default_g) : cl_musp0;

	if (cl_mus0 != UNINITIALIZED && m.lambda != 0)
	  cl_default_mus =
	    cl_mus0 * pow (m.lambda / cl_mus0_lambda, cl_mus0_pwr);

	if (cl_default_mus != UNINITIALIZED)
	  {
	    r.default_mus = cl_default_mus;
	    if (cl_sample_d != UNINITIALIZED)
	      r.default_bs = cl_default_mus * cl_sample_d;
	    else
	      r.default_bs = cl_default_mus * m.slab_thickness;
	  }

	if (cl_search != UNINITIALIZED)
	  r.search = cl_search;



	if (cl_method == COMPARISON && m.d_sphere_r != 0 && m.as_r == 0)
	  {
	    fprintf (stderr,
		     "A dual-beam measurement is specified, but no port sizes.\n");
	    fprintf (stderr,
		     "You might forsake the -X option and use zero spheres (which gives\n");
	    fprintf (stderr,
		     "the same result except lost light is not taken into account).\n");
	    fprintf (stderr,
		     "Alternatively, bite the bullet and enter your sphere parameters,\n");
	    fprintf (stderr,
		     "with the knowledge that only the beam diameter and sample port\n");
	    fprintf (stderr, "diameter are worth obsessing over.\n");
	    exit (EXIT_SUCCESS);
	  }



	if (rt_total == 1 && cl_verbosity > 0)
	  {
	    Write_Header (m, r, params);
	    if (MC_iterations > 0)
	      {
		if (n_photons >= 0)
		  fprintf (stdout,
			   "#  Photons used to estimate lost light =   %ld\n",
			   n_photons);
		else
		  fprintf (stdout,
			   "#     Time used to estimate lost light =   %ld ms\n",
			   -n_photons);
	      }
	    else
	      fprintf (stdout,
		       "#  Photons used to estimate lost light =   0\n");

	    fprintf (stdout, "#\n");

	    print_results_header (stdout);
	  }



	Inverse_RT (m, &r);

	if (r.error == IAD_NO_ERROR)
	  {
	    calculate_coefficients (m, r, &LR, &LT, &mu_sp, &mu_a);



	    if (m.as_r != 0 && r.default_a != 0 && m.num_spheres > 0)
	      {
		double mu_sp_last = mu_sp;
		double mu_a_last = mu_a;

		if (Debug (DEBUG_LOST_LIGHT))
		  {
		    print_results_header (stderr);
		    print_optical_property_result (stderr, m, r, LR, LT, mu_a,
						   mu_sp, mc_iter, rt_total);
		  }

		while (mc_iter < MC_iterations)
		  {

		    MC_Lost (m, r, -1000, &ur1, &ut1, &uru, &utu,
			     &m.ur1_lost, &m.ut1_lost, &m.uru_lost,
			     &m.utu_lost);

		    mc_total++;
		    mc_iter++;

		    Inverse_RT (m, &r);
		    calculate_coefficients (m, r, &LR, &LT, &mu_sp, &mu_a);

		    if (fabs (mu_a_last - mu_a) / (mu_a + 0.0001) <
			r.MC_tolerance
			&& fabs (mu_sp_last - mu_sp) / (mu_sp + 0.0001) <
			r.MC_tolerance)
		      break;

		    mu_a_last = mu_a;
		    mu_sp_last = mu_sp;

		    if (Debug (DEBUG_LOST_LIGHT))
		      print_optical_property_result (stderr, m, r, LR, LT,
						     mu_a, mu_sp, mc_iter,
						     rt_total);
		    else
		      print_dot (start_time, r.error, mc_total, rt_total,
				 mc_iter, cl_verbosity, &any_error);

		    if (r.error != IAD_NO_ERROR)
		      break;
		  }
	      }


	  }
	print_optical_property_result (stdout, m, r, LR, LT, mu_a, mu_sp,
				       mc_iter, rt_total);


	if (Debug (DEBUG_LOST_LIGHT))
	  fprintf (stderr, "\n");
	else
	  print_dot (start_time, r.error, mc_total, rt_total, 99,
		     cl_verbosity, &any_error);
      }


      return EXIT_SUCCESS;
    }

  if (Read_Header (stdin, &m, &params) == 0)
    {
      start_time = clock ();
      while (Read_Data_Line (stdin, &m, params) == 0)
	{


	  if (cl_cos_angle != UNINITIALIZED)
	    {
	      m.slab_cos_angle = cl_cos_angle;
	      if (cl_quadrature_points == UNINITIALIZED)
		cl_quadrature_points = 12;

	      if (cl_quadrature_points != 12 * (cl_quadrature_points / 12))
		{
		  fprintf (stderr,
			   "If you use the -i option to specify an oblique incidence angle, then\n");
		  fprintf (stderr,
			   "the number of quadrature points must be a multiple of 12\n");
		  exit (EXIT_SUCCESS);
		}
	    }

	  if (cl_sample_n != UNINITIALIZED)
	    m.slab_index = cl_sample_n;

	  if (cl_slide_n != UNINITIALIZED)
	    {
	      m.slab_bottom_slide_index = cl_slide_n;
	      m.slab_top_slide_index = cl_slide_n;
	    }

	  if (cl_slide_OD != UNINITIALIZED)
	    {
	      m.slab_bottom_slide_b = cl_slide_OD;
	      m.slab_top_slide_b = cl_slide_OD;
	    }

	  if (cl_sample_d != UNINITIALIZED)
	    m.slab_thickness = cl_sample_d;

	  if (cl_beam_d != UNINITIALIZED)
	    m.d_beam = cl_beam_d;

	  if (cl_slide_d != UNINITIALIZED)
	    {
	      m.slab_bottom_slide_thickness = cl_slide_d;
	      m.slab_top_slide_thickness = cl_slide_d;
	    }

	  if (cl_slides == NO_SLIDES)
	    {
	      m.slab_bottom_slide_index = 1.0;
	      m.slab_bottom_slide_thickness = 0.0;
	      m.slab_top_slide_index = 1.0;
	      m.slab_top_slide_thickness = 0.0;
	    }

	  if (cl_slides == ONE_SLIDE_ON_TOP ||
	      cl_slides == ONE_SLIDE_NEAR_SPHERE)
	    {
	      m.slab_bottom_slide_index = 1.0;
	      m.slab_bottom_slide_thickness = 0.0;
	    }

	  if (cl_slides == ONE_SLIDE_ON_BOTTOM ||
	      cl_slides == ONE_SLIDE_NOT_NEAR_SPHERE)
	    {
	      m.slab_top_slide_index = 1.0;
	      m.slab_top_slide_thickness = 0.0;
	    }

	  if (cl_slides == ONE_SLIDE_NEAR_SPHERE ||
	      cl_slides == ONE_SLIDE_NOT_NEAR_SPHERE)
	    m.flip_sample = 1;
	  else
	    m.flip_sample = 0;

	  if (cl_method != UNINITIALIZED)
	    m.method = (int) cl_method;

	  if (cl_rstd_t != UNINITIALIZED)
	    m.rstd_t = cl_rstd_t;

	  if (cl_rstd_r != UNINITIALIZED)
	    m.rstd_r = cl_rstd_r;

	  if (cl_sphere_one[4] != UNINITIALIZED)
	    {
	      double d_sample_r, d_entrance_r, d_detector_r;

	      m.d_sphere_r = cl_sphere_one[0];
	      d_sample_r = cl_sphere_one[1];
	      d_entrance_r = cl_sphere_one[2];
	      d_detector_r = cl_sphere_one[3];
	      m.rw_r = cl_sphere_one[4];

	      m.as_r =
		(d_sample_r / m.d_sphere_r / 2) * (d_sample_r / m.d_sphere_r /
						   2);
	      m.ae_r =
		(d_entrance_r / m.d_sphere_r / 2) * (d_entrance_r /
						     m.d_sphere_r / 2);
	      m.ad_r =
		(d_detector_r / m.d_sphere_r / 2) * (d_detector_r /
						     m.d_sphere_r / 2);

	      m.aw_r = 1.0 - m.as_r - m.ae_r - m.ad_r;

	      m.d_sphere_t = m.d_sphere_r;
	      m.as_t = m.as_r;
	      m.ae_t = m.ae_r;
	      m.ad_t = m.ad_r;
	      m.aw_t = m.aw_r;
	      m.rw_t = m.rw_r;

	      if (cl_num_spheres == UNINITIALIZED)
		m.num_spheres = 1;
	    }

	  if (cl_sphere_two[4] != UNINITIALIZED)
	    {
	      double d_sample_t, d_entrance_t, d_detector_t;

	      m.d_sphere_t = cl_sphere_two[0];
	      d_sample_t = cl_sphere_two[1];
	      d_entrance_t = cl_sphere_two[2];
	      d_detector_t = cl_sphere_two[3];
	      m.rw_t = cl_sphere_two[4];

	      m.as_t =
		(d_sample_t / m.d_sphere_t / 2) * (d_sample_t / m.d_sphere_t /
						   2);
	      m.ae_t =
		(d_entrance_t / m.d_sphere_t / 2) * (d_entrance_t /
						     m.d_sphere_t / 2);
	      m.ad_t =
		(d_detector_t / m.d_sphere_t / 2) * (d_detector_t /
						     m.d_sphere_t / 2);
	      m.aw_t = 1.0 - m.as_t - m.ae_t - m.ad_t;

	      if (cl_num_spheres == UNINITIALIZED)
		m.num_spheres = 2;
	    }

	  if (cl_num_spheres != UNINITIALIZED)
	    {
	      m.num_spheres = (int) cl_num_spheres;
	      if (m.num_spheres > 0 && m.method == UNKNOWN)
		m.method = SUBSTITUTION;
	    }

	  if (cl_rc_fraction != UNINITIALIZED)
	    m.fraction_of_rc_in_mr = cl_rc_fraction;

	  if (cl_tc_fraction != UNINITIALIZED)
	    m.fraction_of_tc_in_mt = cl_tc_fraction;

	  if (cl_UR1 != UNINITIALIZED)
	    m.m_r = cl_UR1;

	  if (cl_UT1 != UNINITIALIZED)
	    m.m_t = cl_UT1;

	  if (cl_Tc != UNINITIALIZED)
	    m.m_u = cl_Tc;

	  if (cl_default_fr != UNINITIALIZED)
	    m.f_r = cl_default_fr;



	  {

	    static int rt_total = 0;
	    static int mc_total = 0;
	    int mc_iter = 0;

	    double ur1 = 0;
	    double ut1 = 0;
	    double uru = 0;
	    double utu = 0;
	    double mu_a = 0;
	    double mu_sp = 0;
	    double LR = 0;
	    double LT = 0;

	    rt_total++;



	    Initialize_Result (m, &r);


	    if (cl_quadrature_points != UNINITIALIZED)
	      r.method.quad_pts = cl_quadrature_points;
	    else
	      r.method.quad_pts = 8;

	    if (cl_default_a != UNINITIALIZED)
	      r.default_a = cl_default_a;

	    if (cl_default_mua != UNINITIALIZED)
	      {
		r.default_mua = cl_default_mua;
		if (cl_sample_d != UNINITIALIZED)
		  r.default_ba = cl_default_mua * cl_sample_d;
		else
		  r.default_ba = cl_default_mua * m.slab_thickness;
	      }

	    if (cl_default_b != UNINITIALIZED)
	      r.default_b = cl_default_b;

	    if (cl_default_g != UNINITIALIZED)
	      r.default_g = cl_default_g;

	    if (cl_tolerance != UNINITIALIZED)
	      {
		r.tolerance = cl_tolerance;
		r.MC_tolerance = cl_tolerance;
	      }

	    if (cl_musp0 != UNINITIALIZED)
	      cl_mus0 =
		(r.default_g !=
		 UNINITIALIZED) ? cl_musp0 / (1 - r.default_g) : cl_musp0;

	    if (cl_mus0 != UNINITIALIZED && m.lambda != 0)
	      cl_default_mus =
		cl_mus0 * pow (m.lambda / cl_mus0_lambda, cl_mus0_pwr);

	    if (cl_default_mus != UNINITIALIZED)
	      {
		r.default_mus = cl_default_mus;
		if (cl_sample_d != UNINITIALIZED)
		  r.default_bs = cl_default_mus * cl_sample_d;
		else
		  r.default_bs = cl_default_mus * m.slab_thickness;
	      }

	    if (cl_search != UNINITIALIZED)
	      r.search = cl_search;



	    if (cl_method == COMPARISON && m.d_sphere_r != 0 && m.as_r == 0)
	      {
		fprintf (stderr,
			 "A dual-beam measurement is specified, but no port sizes.\n");
		fprintf (stderr,
			 "You might forsake the -X option and use zero spheres (which gives\n");
		fprintf (stderr,
			 "the same result except lost light is not taken into account).\n");
		fprintf (stderr,
			 "Alternatively, bite the bullet and enter your sphere parameters,\n");
		fprintf (stderr,
			 "with the knowledge that only the beam diameter and sample port\n");
		fprintf (stderr, "diameter are worth obsessing over.\n");
		exit (EXIT_SUCCESS);
	      }



	    if (rt_total == 1 && cl_verbosity > 0)
	      {
		Write_Header (m, r, params);
		if (MC_iterations > 0)
		  {
		    if (n_photons >= 0)
		      fprintf (stdout,
			       "#  Photons used to estimate lost light =   %ld\n",
			       n_photons);
		    else
		      fprintf (stdout,
			       "#     Time used to estimate lost light =   %ld ms\n",
			       -n_photons);
		  }
		else
		  fprintf (stdout,
			   "#  Photons used to estimate lost light =   0\n");

		fprintf (stdout, "#\n");

		print_results_header (stdout);
	      }



	    Inverse_RT (m, &r);

	    if (r.error == IAD_NO_ERROR)
	      {
		calculate_coefficients (m, r, &LR, &LT, &mu_sp, &mu_a);



		if (m.as_r != 0 && r.default_a != 0 && m.num_spheres > 0)
		  {
		    double mu_sp_last = mu_sp;
		    double mu_a_last = mu_a;

		    if (Debug (DEBUG_LOST_LIGHT))
		      {
			print_results_header (stderr);
			print_optical_property_result (stderr, m, r, LR, LT,
						       mu_a, mu_sp, mc_iter,
						       rt_total);
		      }

		    while (mc_iter < MC_iterations)
		      {

			MC_Lost (m, r, -1000, &ur1, &ut1, &uru, &utu,
				 &m.ur1_lost, &m.ut1_lost, &m.uru_lost,
				 &m.utu_lost);

			mc_total++;
			mc_iter++;

			Inverse_RT (m, &r);
			calculate_coefficients (m, r, &LR, &LT, &mu_sp,
						&mu_a);

			if (fabs (mu_a_last - mu_a) / (mu_a + 0.0001) <
			    r.MC_tolerance
			    && fabs (mu_sp_last - mu_sp) / (mu_sp + 0.0001) <
			    r.MC_tolerance)
			  break;

			mu_a_last = mu_a;
			mu_sp_last = mu_sp;

			if (Debug (DEBUG_LOST_LIGHT))
			  print_optical_property_result (stderr, m, r, LR, LT,
							 mu_a, mu_sp, mc_iter,
							 rt_total);
			else
			  print_dot (start_time, r.error, mc_total, rt_total,
				     mc_iter, cl_verbosity, &any_error);

			if (r.error != IAD_NO_ERROR)
			  break;
		      }
		  }


	      }
	    print_optical_property_result (stdout, m, r, LR, LT, mu_a, mu_sp,
					   mc_iter, rt_total);


	    if (Debug (DEBUG_LOST_LIGHT))
	      fprintf (stderr, "\n");
	    else
	      print_dot (start_time, r.error, mc_total, rt_total, 99,
			 cl_verbosity, &any_error);
	  }


	}
    }
  if (cl_verbosity > 0)
    fprintf (stderr, "\n\n");
  if (any_error && cl_verbosity > 1)
    print_error_legend ();
  return EXIT_SUCCESS;
}
