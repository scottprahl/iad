
#define _CRT_SECURE_NO_WARNINGS
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "ad_globl.h"
#include "iad_type.h"
#include "iad_io.h"
#include "iad_pub.h"
#include "version.h"



int
skip_white (FILE * fp)
{
  int c = fgetc (fp);

  while (!feof (fp))
    {
      if (isspace (c))
	c = fgetc (fp);
      else if (c == '#')
	do
	  c = fgetc (fp);
	while (!feof (fp) && c != '\n' && c != '\r');
      else
	break;
    }

  if (feof (fp))
    return 1;

  ungetc (c, fp);
  return 0;
}




int
read_number (FILE * fp, double *x)
{
  if (skip_white (fp))
    return 1;

  if (fscanf (fp, "%lf", x))
    return 0;
  else
    return 1;
}




int
check_magic (FILE * fp)
{
  char magic[] = "IAD1";
  int i, c;

  for (i = 0; i < 4; i++)
    {
      c = fgetc (fp);
      if (feof (fp) || c != magic[i])
	{
	  fprintf (stderr,
		   "Sorry, but iad input files must begin with IAD1\n");
	  fprintf (stderr,
		   "       as the first four characters of the file.\n");
	  fprintf (stderr,
		   "       Perhaps you are using an old iad format?\n");
	  return 1;
	}
    }

  return 0;
}




int
Read_Header (FILE * fp, struct measure_type *m, int *params)
{
  double x;
  Initialize_Measure (m);
  if (check_magic (fp))
    return 1;
  if (read_number (fp, &m->slab_index))
    return 1;
  if (read_number (fp, &m->slab_top_slide_index))
    return 1;
  if (read_number (fp, &m->slab_thickness))
    return 1;
  if (read_number (fp, &m->slab_top_slide_thickness))
    return 1;
  if (read_number (fp, &m->d_beam))
    return 1;

  if (m->slab_top_slide_thickness == 0.0)
    m->slab_top_slide_index = 1.0;
  if (m->slab_top_slide_index == 1.0)
    m->slab_top_slide_thickness = 0.0;
  if (m->slab_top_slide_index == 0.0)
    {
      m->slab_top_slide_thickness = 0.0;
      m->slab_top_slide_index = 1.0;
    }

  m->slab_bottom_slide_index = m->slab_top_slide_index;
  m->slab_bottom_slide_thickness = m->slab_top_slide_thickness;

  if (read_number (fp, &m->rstd_r))
    return 1;

  if (read_number (fp, &x))
    return 1;
  m->num_spheres = (int) x;
  m->method = SUBSTITUTION;


  {
    double d_sample_r, d_entrance_r, d_detector_r;
    if (read_number (fp, &m->d_sphere_r))
      return 1;
    if (read_number (fp, &d_sample_r))
      return 1;
    if (read_number (fp, &d_entrance_r))
      return 1;
    if (read_number (fp, &d_detector_r))
      return 1;
    if (read_number (fp, &m->rw_r))
      return 1;

    m->as_r =
      (d_sample_r / m->d_sphere_r) * (d_sample_r / m->d_sphere_r) / 4.0;
    m->ae_r =
      (d_entrance_r / m->d_sphere_r) * (d_entrance_r / m->d_sphere_r) / 4.0;
    m->ad_r =
      (d_detector_r / m->d_sphere_r) * (d_detector_r / m->d_sphere_r) / 4.0;
    m->aw_r = 1.0 - m->as_r - m->ae_r - m->ad_r;
  }




  {
    double d_sample_t, d_entrance_t, d_detector_t;
    if (read_number (fp, &m->d_sphere_t))
      return 1;
    if (read_number (fp, &d_sample_t))
      return 1;
    if (read_number (fp, &d_entrance_t))
      return 1;
    if (read_number (fp, &d_detector_t))
      return 1;
    if (read_number (fp, &m->rw_t))
      return 1;

    m->as_t =
      (d_sample_t / m->d_sphere_t) * (d_sample_t / m->d_sphere_t) / 4.0;
    m->ae_t =
      (d_entrance_t / m->d_sphere_t) * (d_entrance_t / m->d_sphere_t) / 4.0;
    m->ad_t =
      (d_detector_t / m->d_sphere_t) * (d_detector_t / m->d_sphere_t) / 4.0;
    m->aw_t = 1.0 - m->as_t - m->ae_t - m->ad_t;
  }



  if (read_number (fp, &x))
    return 1;
  *params = (int) x;
  m->num_measures = (*params >= 3) ? 3 : *params;

  return 0;
}




void
Write_Header (struct measure_type m, struct invert_type r, int params)
{

  double xx;

  printf ("# Inverse Adding-Doubling %s \n", Version);
  printf ("# \n");
  printf ("#                        Beam diameter = %7.1f mm\n", m.d_beam);
  printf ("#                     Sample thickness = %7.3f mm\n",
	  m.slab_thickness);
  printf ("#                  Top slide thickness = %7.3f mm\n",
	  m.slab_top_slide_thickness);
  printf ("#               Bottom slide thickness = %7.3f mm\n",
	  m.slab_bottom_slide_thickness);
  printf ("#           Sample index of refraction = %7.4f\n", m.slab_index);
  printf ("#        Top slide index of refraction = %7.4f\n",
	  m.slab_top_slide_index);
  printf ("#     Bottom slide index of refraction = %7.4f\n",
	  m.slab_bottom_slide_index);



  printf ("# \n");




  printf ("#    Fraction unscattered refl. in M_R = %7.1f %%\n",
	  m.fraction_of_rc_in_mr * 100);
  printf ("#   Fraction unscattered trans. in M_T = %7.1f %%\n",
	  m.fraction_of_tc_in_mt * 100);
  printf ("# \n");



  printf ("# Reflection sphere\n");
  printf ("#                      sphere diameter = %7.1f mm\n",
	  m.d_sphere_r);
  printf ("#                 sample port diameter = %7.1f mm\n",
	  2 * m.d_sphere_r * sqrt (m.as_r));
  printf ("#               entrance port diameter = %7.1f mm\n",
	  2 * m.d_sphere_r * sqrt (m.ae_r));
  printf ("#               detector port diameter = %7.1f mm\n",
	  2 * m.d_sphere_r * sqrt (m.ad_r));
  printf ("#                     wall reflectance = %7.1f %%\n",
	  m.rw_r * 100);
  printf ("#                 standard reflectance = %7.1f %%\n",
	  m.rstd_r * 100);
  printf ("#                 detector reflectance = %7.1f %%\n",
	  m.rd_r * 100);
  printf ("#\n");



  printf ("# Transmission sphere\n");
  printf ("#                      sphere diameter = %7.1f mm\n",
	  m.d_sphere_t);
  printf ("#                 sample port diameter = %7.1f mm\n",
	  2 * m.d_sphere_r * sqrt (m.as_t));
  printf ("#               entrance port diameter = %7.1f mm\n",
	  2 * m.d_sphere_r * sqrt (m.ae_t));
  printf ("#               detector port diameter = %7.1f mm\n",
	  2 * m.d_sphere_r * sqrt (m.ad_t));
  printf ("#                     wall reflectance = %7.1f %%\n",
	  m.rw_t * 100);
  printf ("#               standard transmittance = %7.1f %%\n",
	  m.rstd_t * 100);
  printf ("#                 detector reflectance = %7.1f %%\n",
	  m.rd_t * 100);



  printf ("#\n");
  switch (params)
    {
    case -1:
      printf ("# No M_R or M_T -- forward calculation.\n");
      break;
    case 1:
      printf ("# Just M_R was measured");
      break;
    case 2:
      printf ("# M_R and M_T were measured");
      break;
    case 3:
      printf ("# M_R, M_T, and M_U were measured");
      break;
    case 4:
      printf ("# M_R, M_T, M_U, and r_w were measured");
      break;
    case 5:
      printf ("# M_R, M_T, M_U, r_w, and t_w were measured");
      break;
    case 6:
      printf ("# M_R, M_T, M_U, r_w, t_w, and r_std were measured");
      break;
    case 7:
      printf ("# M_R, M_T, M_U, r_w, t_w, r_std and t_std were measured");
      break;
    default:
      printf ("# Something went wrong ... measures should be 1 to 5!\n");
      break;
    }

  if (1 <= params && params <= 7)
    {
      if (m.flip_sample)
	printf (" (sample flipped) ");

      switch (m.method)
	{
	case UNKNOWN:
	  printf (" using an unknown method.\n");
	  break;
	case SUBSTITUTION:
	  printf (" using the substitution (single-beam) method.\n");
	  break;
	case COMPARISON:
	  printf (" using the comparison (dual-beam) method.\n");
	}
    }

  switch (m.num_spheres)
    {
    case 0:
      printf ("# No sphere corrections were used");
      break;

    case 1:
      printf ("# Single sphere corrections were used");
      break;

    case 2:
      printf ("# Double sphere corrections were used");
      break;
    }

  printf (" with light incident at %d degrees from the normal",
	  (int) (acos (m.slab_cos_angle) * 57.2958));
  printf (".\n");

  switch (r.search)
    {
    case FIND_AB:
      printf ("# The inverse routine varied the albedo and optical depth.\n");
      printf ("# \n");
      xx = (r.default_g != UNINITIALIZED) ? r.default_g : 0;
      printf ("# Default single scattering anisotropy = %7.3f \n", xx);
      break;
    case FIND_AG:
      printf ("# The inverse routine varied the albedo and anisotropy.\n");
      printf ("# \n");
      if (r.default_b != UNINITIALIZED)
	printf ("#                     Default (mu_t*d) = %7.3g\n",
		r.default_b);
      else
	printf ("# \n");
      break;
    case FIND_AUTO:
      printf ("# The inverse routine adapted to the input data.\n");
      printf ("# \n");
      printf ("# \n");
      break;
    case FIND_A:
      printf ("# The inverse routine varied only the albedo.\n");
      printf ("# \n");
      xx = (r.default_g != UNINITIALIZED) ? r.default_g : 0;
      printf ("# Default single scattering anisotropy is %7.3f ", xx);
      xx = (r.default_b != UNINITIALIZED) ? r.default_b : HUGE_VAL;
      printf (" and (mu_t*d) = %7.3g\n", xx);
      break;
    case FIND_B:
      printf ("# The inverse routine varied only the optical depth.\n");
      printf ("# \n");
      xx = (r.default_g != UNINITIALIZED) ? r.default_g : 0;
      printf ("# Default single scattering anisotropy is %7.3f ", xx);
      if (r.default_a != UNINITIALIZED)
	printf ("and default albedo = %7.3g\n", r.default_a);
      else
	printf ("\n");
      break;
    case FIND_Ba:
      printf ("# The inverse routine varied only the absorption.\n");
      printf ("# \n");
      xx = (r.default_bs != UNINITIALIZED) ? r.default_bs : 0;
      printf ("#                     Default (mu_s*d) = %7.3g\n", xx);
      break;
    case FIND_Bs:
      printf ("# The inverse routine varied only the scattering.\n");
      printf ("# \n");
      xx = (r.default_ba != UNINITIALIZED) ? r.default_ba : 0;
      printf ("#                     Default (mu_a*d) = %7.3g\n", xx);
      break;
    default:
      printf ("# \n");
      printf ("# \n");
      printf ("# \n");
      break;
    }

  printf ("#                 AD quadrature points = %3d\n",
	  r.method.quad_pts);
  printf ("#             AD tolerance for success = %9.5f\n", r.tolerance);
  printf ("#      MC tolerance for mu_a and mu_s' = %7.3f %%\n",
	  r.MC_tolerance);
}




int
Read_Data_Line (FILE * fp, struct measure_type *m, int params)
{
  if (read_number (fp, &m->m_r))
    return 1;
  if (m->m_r > 1)
    {
      m->lambda = m->m_r;
      if (read_number (fp, &m->m_r))
	return 1;
    }

  if (params == 1)
    return 0;

  if (read_number (fp, &m->m_t))
    return 1;
  if (params == 2)
    return 0;

  if (read_number (fp, &m->m_u))
    return 1;
  if (params == 3)
    return 0;

  if (read_number (fp, &m->rw_r))
    return 1;
  m->rw_t = m->rw_r;
  if (params == 4)
    return 0;

  if (read_number (fp, &m->rw_t))
    return 1;
  if (params == 5)
    return 0;

  if (read_number (fp, &m->rstd_r))
    return 1;
  if (params == 6)
    return 0;

  if (read_number (fp, &m->rstd_t))
    return 1;
  return 0;
}
