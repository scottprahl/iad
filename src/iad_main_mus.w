@** Main Program.

Here is a quick skeleton that I put together to show how the
inverse adding-doubling code works.  I have only cursorily tested
this.  If you find obvious bugs, they are probably real but should
not extend beyond this code snippet. 

All the output for this web file goes into \.{iad\_main.c}

@(iad_main_mus.c@>=

@<Include files for |main|@>@;

int main (int argc, char **argv){

@<Declare variables for |main|@>@;

if (Read_Header (&m, &r) == TRUE) {
		@<Process the header@>@;
		m.num_measures = 2;
		m.m_r=0.0;
		m.slab_thickness = 0.1;
		while (fp!=EOF)
		{
			fp = scanf("%lf %lf %lf %lf", &lambda, &r.mu_a, &m.default_g, &m.m_t);
			fp = readln(&line);
			@<Calculate and write optical properties@>@;
		}
}

@ @<Include files for |main|@>=
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "ad_globl.h"
#include "iad_type.h"
#include "iad_pub.h"
#include "iad_io.h"

@ @<Declare variables for |main|@>=
  struct measure_type m;
  struct invert_type r;
  int lines;
char *format1 = "    R    \t    T    \t    Tc    \t     a    \
\t     b    \t     g    \t    !/?\n";
  char *format2 = " %9.5f\t %9.5f\t %9.5f\t %9.5f\t %9.5f\t %9.5f\t %9c\n";
  char found = '?';
  int fp;
  double lambda;
  int line = 1;

  
@ @<Process the header@>=

      m.slab_thickness = 1;
      Initialize_Result (m, &r);
      Write_Header (m, r);

      lines = 1;
      printf (format1, m.m_r, m.m_t, m.m_u, r.a, r.b, r.g);

@ @<Calculate and write optical properties@>=
	{
	  r.search = FIND_mus;
	  Inverse_RT (m, &r);
	  if (r.found == TRUE)
	    found = '!';
	  else
	    found = '?';
	  printf (format2, m.m_r, m.m_t, m.m_u, r.a, r.b, r.g, found);
	  fflush (stdout);
	}
