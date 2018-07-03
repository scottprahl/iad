@*1 IAD Find.
March 1995.  Incorporated the |quick_guess| algorithm for low albedos.
@(iad_find.c@>=
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ad_globl.h"
#include "nr_util.h"
#include "nr_mnbrk.h"
#include "nr_brent.h"
#include "nr_amoeb.h"
#include "iad_type.h"
#include "iad_util.h"
#include "iad_calc.h"
#include "iad_find.h"
#define NUMBER_OF_GUESSES 10

guess_type guess[NUMBER_OF_GUESSES];

int compare_guesses (const void *p1, const void *p2)
{
	guess_type *g1 = (guess_type *) p1;
	guess_type *g2 = (guess_type *) p2;
	
	if (g1->distance < g2->distance)
		return -1;
	else if (g1->distance == g2->distance)
		return 0;
	else
		return 1;
}

	@<Definition for |U_Find_Ba|@>@;
	@<Definition for |U_Find_Bs|@>@;
	@<Definition for |U_Find_A|@>@;
	@<Definition for |U_Find_B|@>@;
	@<Definition for |U_Find_G|@>@;
	@<Definition for |U_Find_AG|@>@;
	@<Definition for |U_Find_AB|@>@;
	@<Definition for |U_Find_BG|@>@;
	@<Definition for |U_Find_BaG|@>@;
	@<Definition for |U_Find_BsG|@>@;

@ All the information that needs to be written to
the header file \.{iad\_find.h}. 
This eliminates the need to maintain a set of header files as 
well.

@(iad_find.h@>=
	@<Prototype for |U_Find_Ba|@>;
	@<Prototype for |U_Find_Bs|@>;
	@<Prototype for |U_Find_A|@>;
	@<Prototype for |U_Find_B|@>;
	@<Prototype for |U_Find_G|@>;
	@<Prototype for |U_Find_AG|@>;
	@<Prototype for |U_Find_AB|@>;
	@<Prototype for |U_Find_BG|@>;
	@<Prototype for |U_Find_BaG|@>;
	@<Prototype for |U_Find_BsG|@>;

@*2 Fixed Anisotropy.

This is the most common case.

@<Prototype for |U_Find_AB|@>=
	void U_Find_AB(struct measure_type m, struct invert_type *r)

@ @<Definition for |U_Find_AB|@>=
@<Prototype for |U_Find_AB|@> 
{
	@<Allocate local simplex variables@>@;

	if (Debug(DEBUG_SEARCH)) {
		fprintf(stderr,"In U_Find_AB");
		fprintf(stderr," (mu=%6.4f)",r->slab.cos_angle);
		if (r->default_g != UNINITIALIZED) 
			fprintf(stderr,"  default_g = %8.5f", r->default_g);
		fprintf(stderr,"\n");
	}
	
	r->slab.g = (r->default_g == UNINITIALIZED) ? 0 : r->default_g;
	Set_Calc_State(m, *r);
	
	@<Get the initial |a|, |b|, and |g|@>@;
	@<Initialize the nodes of the |a| and |b| simplex@>@;
	@<Evaluate the |a| and |b| simplex at the nodes@>@;
	amoeba(p, y, 2, r->tolerance, Find_AB_fn, &r->iterations);
	@<Choose the best node of the |a| and |b| simplex@>@;

	@<Free simplex data structures@>@;
	@<Put final values in result@>@;
}

@ To use the simplex algorithm, we need to vectors and a matrix.
@<Allocate local simplex variables@>=
	int i, i_best, j_best;
	double *x, *y, **p;
	
	x=dvector(1,2);
	y=dvector(1,3);
	p=dmatrix(1,3,1,2);

@ Just get the optimal optical properties to start the
search process. 

I had to add the line that tests to make sure the albedo is greater 
than 0.2 because the grid just does not work so well in this case.
The problem is that for low albedos there is really very little 
information about the anisotropy available.  This change was also
made in the analagous code for |a| and |b|.

@<Get the initial |a|, |b|, and |g|@>=
{ 	
/*	double a3,b3,g3;*/
	size_t  count = NUMBER_OF_GUESSES;
	/* distance to last result */
	abg_distance(r->slab.a, r->slab.b, r->slab.g, &(guess[0]));
	
	if (!Valid_Grid(m, r->search)) Fill_Grid(m,*r);	
	
	/* distance to nearest grid point */
	Near_Grid_Points(m.m_r,m.m_t,r->search, &i_best, &j_best);
	Grid_ABG(i_best  ,j_best  ,&(guess[1]));
	Grid_ABG(i_best+1,j_best  ,&(guess[2]));
	Grid_ABG(i_best-1,j_best  ,&(guess[3]));
	Grid_ABG(i_best  ,j_best+1,&(guess[4]));
	Grid_ABG(i_best  ,j_best-1,&(guess[5]));
	Grid_ABG(i_best+1,j_best+1,&(guess[6]));
	Grid_ABG(i_best-1,j_best-1,&(guess[7]));
	Grid_ABG(i_best+1,j_best-1,&(guess[8]));
	Grid_ABG(i_best-1,j_best+1,&(guess[9]));
		
	qsort((void *) guess, count, sizeof(guess_type), compare_guesses);
	
	if (Debug(DEBUG_BEST_GUESS)) {
		int k;
		fprintf(stderr, "after\n");
		for (k=0; k<=6; k++) {
			fprintf(stderr, "%3d  ", k);
			fprintf(stderr, "%10.5f ", guess[k].a);
			fprintf(stderr, "%10.5f ", guess[k].b);
			fprintf(stderr, "%10.5f ", guess[k].g);
			fprintf(stderr, "%10.5f\n", guess[k].distance);
		}
	}
}

@ @<Initialize the nodes of the |a| and |b| simplex@>=
	{
	int k,kk;
	
	p[1][1] = a2acalc(guess[0].a);
	p[1][2] = b2bcalc(guess[0].b);
	
	for (k=1;k<7;k++) {
		if (guess[0].a != guess[k].a) 
			break;
	}
	
	p[2][1] = a2acalc(guess[k].a);
	p[2][2] = b2bcalc(guess[k].b);
	
	for (kk=1;kk<7;kk++) {
		if (guess[0].b != guess[kk].b && guess[k].b != guess[kk].b) 
			break;
	}
	p[3][1] = a2acalc(guess[kk].a);	
	p[3][2] = b2bcalc(guess[kk].b);

	if (Debug(DEBUG_BEST_GUESS)) {
		fprintf(stderr, "guess 1");
		fprintf(stderr, "%10.5f ", guess[0].a);
		fprintf(stderr, "%10.5f ", guess[0].b);
		fprintf(stderr, "%10.5f ", guess[0].g);
		fprintf(stderr, "%10.5f\n", guess[0].distance);
		fprintf(stderr, "guess 2");
		fprintf(stderr, "%10.5f ", guess[k].a);
		fprintf(stderr, "%10.5f ", guess[k].b);
		fprintf(stderr, "%10.5f ", guess[k].g);
		fprintf(stderr, "%10.5f\n", guess[k].distance);
		fprintf(stderr, "guess 3");
		fprintf(stderr, "%10.5f ", guess[kk].a);
		fprintf(stderr, "%10.5f ", guess[kk].b);
		fprintf(stderr, "%10.5f ", guess[kk].g);
		fprintf(stderr, "%10.5f\n", guess[kk].distance);
	}
}

@ @<Evaluate the |a| and |b| simplex at the nodes@>=

	for (i = 1; i <= 3; i++) {
		x[1] = p[i][1];
		x[2] = p[i][2];
		y[i] = Find_AB_fn(x);
	}

@ @<Choose the best node of the |a| and |b| simplex@>=
    r->final_distance=10;
	for(i=1;i<=3;i++) {		
		if (y[i] < r->final_distance) {
			r->slab.a = acalc2a(p[i][1]);
			r->slab.b = bcalc2b(p[i][2]);
			r->final_distance = y[i];
		}
    }

@ @<Put final values in result@>=
	r->a = r->slab.a;
	r->b = r->slab.b;
	r->g = r->slab.g;
    r->found = (r->tolerance <= r->final_distance);

@ Since we allocated these puppies, we got to get rid of them.
@<Free simplex data structures@>=
	free_dvector(x,1,2);
	free_dvector(y,1,3);
	free_dmatrix(p,1,3,1,2);

@*2 Fixed Absorption and Anisotropy.
Typically, this routine is called when the absorption coefficient is
known, the anisotropy is known, and the physical thickness of the sample
is known.  This routine calculates the varies the scattering coefficient
until the measurements are matched.

This was written for Ted Moffitt to analyze some intralipid data.  We
wanted to know what the scattering coefficient of the Intralipid was and
made total transmission measurements through a sample with a fixed physical
thickness.  We did not make reflection measurements because the light source
diverged too much, and we could not make reflection measurements easily.

In retrospect, we could have made URU measurements by illuminating the wall
of the integrating sphere.  However, these diffuse type of measurements are
very difficult to make accurately.

This is tricky only because the value in |slab.b| is used to hold the
value of |ba| or $d \cdot \mu_a$ when the |Find_Bs_fn| is used.

@<Prototype for |U_Find_Bs|@>=
	void U_Find_Bs(struct measure_type m, struct invert_type *r)

@ @<Definition for |U_Find_Bs|@>=
@<Prototype for |U_Find_Bs|@> 
{
	double ax, bx, cx, fa, fb, fc, bs;
	
	if (Debug(DEBUG_SEARCH)) {
		fprintf(stderr,"In U_Find_Bs");
		fprintf(stderr," (mu=%6.4f)",r->slab.cos_angle);
		if (r->default_ba != UNINITIALIZED) 
			fprintf(stderr,"  default_ba = %8.5f", r->default_ba);
		if (r->default_g != UNINITIALIZED) 
			fprintf(stderr,"  default_g = %8.5f", r->default_g);
		fprintf(stderr,"\n");
	}

	r->slab.a = 0;
	r->slab.g = (r->default_g  == UNINITIALIZED) ? 0        : r->default_g;
	r->slab.b = (r->default_ba == UNINITIALIZED) ? HUGE_VAL : r->default_ba; 
	
	Set_Calc_State(m,*r);   /* store ba in RR.slab.b */
	
	ax = b2bcalc(0.1);      /* first try for bs */
	bx = b2bcalc(1.0);
	mnbrak(&ax,&bx,&cx,&fa,&fb,&fc, Find_Bs_fn);
	r->final_distance=brent(ax,bx,cx,Find_Bs_fn,r->tolerance,&bs);

	/* recover true values */
	r->slab.a = bcalc2b(bs)/(bcalc2b(bs)+r->slab.b);
	r->slab.b = bcalc2b(bs) + r->slab.b;   	
	Set_Calc_State(m,*r);
	
	@<Put final values in result@>@;
}

@*2 Fixed Absorption and Scattering.
Typically, this routine is called when the scattering coefficient is
known, the anisotropy is known, and the physical thickness of the sample
is known.  This routine calculates the varies the absorption coefficient
until the measurements are matched.

This is tricky only because the value in |slab.b| is used to hold the
value of |bs| or $d \cdot \mu_s$ when the |Find_Ba_fn| is used.

@<Prototype for |U_Find_Ba|@>=
	void U_Find_Ba(struct measure_type m, struct invert_type *r)

@ @<Definition for |U_Find_Ba|@>=
@<Prototype for |U_Find_Ba|@> 
{
	double ax, bx, cx, fa, fb, fc, ba;
	
	if (Debug(DEBUG_SEARCH)) {
		fprintf(stderr,"In U_Find_Bs");
		fprintf(stderr," (mu=%6.4f)",r->slab.cos_angle);
		if (r->default_bs != UNINITIALIZED) 
			fprintf(stderr,"  default_bs = %8.5f", r->default_bs);
		if (r->default_g != UNINITIALIZED) 
			fprintf(stderr,"  default_g = %8.5f", r->default_g);
		fprintf(stderr,"\n");
	}

	r->slab.a = 0;
	r->slab.g = (r->default_g  == UNINITIALIZED) ? 0        : r->default_g;
	r->slab.b = (r->default_bs == UNINITIALIZED) ? HUGE_VAL : r->default_bs; 
	
	Set_Calc_State(m,*r);   /* store bs in RR.slab.b */
	
	ax = b2bcalc(0.1);      /* first try for ba */
	bx = b2bcalc(1.0);
	mnbrak(&ax,&bx,&cx,&fa,&fb,&fc, Find_Ba_fn);
	r->final_distance=brent(ax,bx,cx,Find_Ba_fn,r->tolerance,&ba);
	
	/* recover true values */
	r->slab.a = (r->slab.b)/(bcalc2b(ba)+r->slab.b);
	r->slab.b = bcalc2b(ba) + r->slab.b;   /*actual value of b */	
	Set_Calc_State(m,*r);   
	
	@<Put final values in result@>@;
}

@*2 Fixed Optical Depth and Anisotropy.
Typically, this routine is called when the optical thickness is assumed
infinite.  However, it may also be called when the optical thickness is
assumed to be fixed at a particular value.  Typically the only reasonable
situation for this to occur is when the diffuse transmission is non-zero
but the collimated transmission is zero.  If this is the case then there is
no information in the collimated transmission measurement and there
is no sense even using it because the slab is not infinitely thick.

@<Prototype for |U_Find_A|@>=
	void U_Find_A(struct measure_type m, struct invert_type *r)

@ @<Definition for |U_Find_A|@>=
@<Prototype for |U_Find_A|@> {
	double Rt, Tt, Rd, Rc, Td, Tc;
	
	if (Debug(DEBUG_SEARCH)) {
		fprintf(stderr,"In U_Find_A");
		fprintf(stderr," (mu=%6.4f)",r->slab.cos_angle);
		if (r->default_b != UNINITIALIZED) 
			fprintf(stderr,"  default_b = %8.5f", r->default_b);
		if (r->default_g != UNINITIALIZED) 
			fprintf(stderr,"  default_g = %8.5f", r->default_g);
		fprintf(stderr,"\n");
	}
	
    Estimate_RT(m, *r, &Rt, &Tt, &Rd, &Rc, &Td, &Tc);

	r->slab.g = (r->default_g == UNINITIALIZED) ? 0        : r->default_g;
	r->slab.b = (r->default_b == UNINITIALIZED) ? HUGE_VAL : r->default_b;
	r->slab.a = 0.0;
	r->final_distance = 0.0;
	Set_Calc_State(m,*r);

	if (Rt > 0.99999) 
		r->final_distance = Find_A_fn(a2acalc(1.0));		
	else {
		double x, ax, bx, cx, fa, fb, fc;
	
		ax=a2acalc(0.3);
		bx=a2acalc(0.5);
	
		mnbrak(&ax,&bx,&cx,&fa,&fb,&fc, Find_A_fn);
		r->final_distance=brent(ax,bx,cx,Find_A_fn,r->tolerance,&x);
		r->slab.a = acalc2a(x);
	}
		
	@<Put final values in result@>@;
}

@*2 Fixed Optical Depth and Albedo.

@<Prototype for |U_Find_G|@>=
	void U_Find_G(struct measure_type m, struct invert_type *r)

@ @<Definition for |U_Find_G|@>=
@<Prototype for |U_Find_G|@> {
	double Rt, Tt, Rd, Rc, Td, Tc;
	
	if (Debug(DEBUG_SEARCH)) {
		fprintf(stderr,"In U_Find_G");
		fprintf(stderr," (mu=%6.4f)",r->slab.cos_angle);
		if (r->default_a != UNINITIALIZED) 
			fprintf(stderr,"  default_a = %8.5f", r->default_a);
		if (r->default_b != UNINITIALIZED) 
			fprintf(stderr,"  default_b = %8.5f", r->default_b);
		fprintf(stderr,"\n");
	}

	Estimate_RT(m, *r, &Rt, &Tt, &Rd, &Rc, &Td, &Tc);

	r->slab.a = (r->default_a == UNINITIALIZED) ? 0.5      : r->default_a;
	r->slab.b = (r->default_b == UNINITIALIZED) ? HUGE_VAL : r->default_b;
	r->slab.g = 0.0;
	r->final_distance = 0.0;
	Set_Calc_State(m,*r);

	if (Rd>0.0) {
		double x, ax, bx, cx, fa, fb, fc;
	
		ax=g2gcalc(-0.99);
		bx=g2gcalc(0.99);
	
		mnbrak(&ax,&bx,&cx,&fa,&fb,&fc, Find_G_fn);
		r->final_distance=brent(ax,bx,cx,Find_G_fn,r->tolerance,&x);
		r->slab.g = gcalc2g(x);
		Set_Calc_State(m,*r);
	}
		
	@<Put final values in result@>@;
}


@*2 Fixed Anisotropy and Albedo.
This routine can be called in three different situations: (1) the albedo is
zero, (2) the albedo is one, or (3) the albedo is fixed at a default value.
I calculate the individual reflections and transmissions to establish
which of these cases we happen to have.
 
@<Prototype for |U_Find_B|@>=
	void U_Find_B(struct measure_type m, struct invert_type *r)

@ @<Definition for |U_Find_B|@>=
@<Prototype for |U_Find_B|@> {
	double Rt, Tt, Rd, Rc, Td, Tc;
	
	if (Debug(DEBUG_SEARCH)) {
		fprintf(stderr,"In U_Find_B");
		fprintf(stderr," (mu=%6.4f)",r->slab.cos_angle);
		if (r->default_a != UNINITIALIZED) 
			fprintf(stderr,"  default_a = %8.5f", r->default_a);
		if (r->default_g != UNINITIALIZED) 
			fprintf(stderr,"  default_g = %8.5f", r->default_g);
		fprintf(stderr,"\n");
	}

	Estimate_RT(m, *r, &Rt, &Tt, &Rd, &Rc, &Td, &Tc);

	r->slab.g = (r->default_g == UNINITIALIZED) ? 0 : r->default_g;
	r->slab.a = (r->default_a == UNINITIALIZED) ? 0 : r->default_a;
	r->slab.b = 0.5;
	r->final_distance = 0.0;
	Set_Calc_State(m,*r);
	
	@<Iteratively solve for |b|@>@;
	
	@<Put final values in result@>@;

	if (Debug(DEBUG_SEARCH)) {
		fprintf(stderr,"In U_Find_B final (a,b,g) = ");
		fprintf(stderr,"(%8.5f,%8.5f,%8.5f)\n", r->a, r->b, r->g);
	}
}
	
@ This could be improved tremendously.  I just don't want to 
mess with it at the moment.
@<Iteratively solve for |b|@>=
{
	double x,ax,bx,cx,fa,fb,fc;
	
	ax=b2bcalc(0.1);
	bx=b2bcalc(10);
	
	mnbrak(&ax,&bx,&cx,&fa,&fb,&fc, Find_B_fn);	
	r->final_distance=brent(ax,bx,cx,Find_B_fn,r->tolerance,&x);
	r->slab.b = bcalc2b(x);
	Set_Calc_State(m,*r);
}


@*2 Fixed Optical Depth.

We can get here a couple of different ways.

First there can be three real measurements, i.e., $t_c$ is not zero,
in this case we want to fix $b$ based on the $t_c$ measurement.

Second, we can get here if a default value for $b$ has been set.

Otherwise, we really should not be here.  Just set $b=1$ and calculate
away.


@<Prototype for |U_Find_AG|@>=
	void U_Find_AG(struct measure_type m, struct invert_type *r)
	
@ @<Definition for |U_Find_AG|@>=
@<Prototype for |U_Find_AG|@> 
{
	@<Allocate local simplex variables@>@;

	if (Debug(DEBUG_SEARCH)) {
		fprintf(stderr,"In U_Find_AG");
		fprintf(stderr," (mu=%6.4f)",r->slab.cos_angle);
		if (r->default_b != UNINITIALIZED) 
			fprintf(stderr,"  default_b = %8.5f", r->default_b);
		fprintf(stderr,"\n");
	}
	
	if (m.num_measures==3) 
		r->slab.b = What_Is_B(r->slab, m.m_u);
	else if (r->default_b == UNINITIALIZED)
		r->slab.b = 1;
	else
		r->slab.b = r->default_b;

	Set_Calc_State(m, *r);
	@<Get the initial |a|, |b|, and |g|@>@;
	@<Initialize the nodes of the |a| and |g| simplex@>@;
	@<Evaluate the |a| and |g| simplex at the nodes@>@;
	amoeba(p, y, 2, r->tolerance, Find_AG_fn, &r->iterations);
	@<Choose the best node of the |a| and |g| simplex@>@;
	@<Free simplex data structures@>@;

	@<Put final values in result@>@;
}

@ @<Initialize the nodes of the |a| and |g| simplex@>=
	{
	int k,kk;
	
	p[1][1] = a2acalc(guess[0].a);
	p[1][2] = g2gcalc(guess[0].g);
	
	for (k=1;k<7;k++) {
		if (guess[0].a != guess[k].a) 
			break;
	}
	
	p[2][1] = a2acalc(guess[k].a);
	p[2][2] = g2gcalc(guess[k].g);
	
	for (kk=1;kk<7;kk++) {
		if (guess[0].g != guess[kk].g && guess[k].g != guess[kk].g) 
			break;
	}
	p[3][1] = a2acalc(guess[kk].a);	
	p[3][2] = g2gcalc(guess[kk].g);

	if (Debug(DEBUG_BEST_GUESS)) {
		fprintf(stderr, "guess 1");
		fprintf(stderr, "%10.5f ", guess[0].a);
		fprintf(stderr, "%10.5f ", guess[0].b);
		fprintf(stderr, "%10.5f ", guess[0].g);
		fprintf(stderr, "%10.5f\n", guess[0].distance);
		fprintf(stderr, "guess 2");
		fprintf(stderr, "%10.5f ", guess[k].a);
		fprintf(stderr, "%10.5f ", guess[k].b);
		fprintf(stderr, "%10.5f ", guess[k].g);
		fprintf(stderr, "%10.5f\n", guess[k].distance);
		fprintf(stderr, "guess 3");
		fprintf(stderr, "%10.5f ", guess[kk].a);
		fprintf(stderr, "%10.5f ", guess[kk].b);
		fprintf(stderr, "%10.5f ", guess[kk].g);
		fprintf(stderr, "%10.5f\n", guess[kk].distance);
	}
}

@ @<Evaluate the |a| and |g| simplex at the nodes@>=

  for (i = 1; i <= 3; i++) {
      x[1] = p[i][1];
      x[2] = p[i][2];
      y[i] = Find_AG_fn(x);
  }


@ Here we find the node of the simplex that gave the best
result and save that one.  At the same time we save the whole
simplex for later use if needed.

@<Choose the best node of the |a| and |g| simplex@>=
    r->final_distance=10;
	for(i=1;i<=3;i++) {
		if (y[i] < r->final_distance) {
			r->slab.a = acalc2a(p[i][1]);
			r->slab.g = gcalc2g(p[i][2]);
			r->final_distance = y[i];
		}
    }


@*2 Fixed Albedo.
Here the optical depth and the anisotropy are varied (for a fixed
albedo).

@<Prototype for |U_Find_BG|@>=
	void U_Find_BG(struct measure_type m, struct invert_type *r)
	
@ @<Definition for |U_Find_BG|@>=
@<Prototype for |U_Find_BG|@> 
{
	@<Allocate local simplex variables@>@;

	if (Debug(DEBUG_SEARCH)) {
		fprintf(stderr,"In U_Find_BG");
		fprintf(stderr," (mu=%6.4f)",r->slab.cos_angle);
		if (r->default_a != UNINITIALIZED) 
			fprintf(stderr,"  default_a = %8.5f", r->default_a);
		fprintf(stderr,"\n");
	}

	r->slab.a = (r->default_a == UNINITIALIZED) ? 0 : r->default_a;
	Set_Calc_State(m, *r);

	@<Get the initial |a|, |b|, and |g|@>@;
	@<Initialize the nodes of the |b| and |g| simplex@>@;
	@<Evaluate the |bg| simplex at the nodes@>@;
	amoeba(p, y, 2, r->tolerance, Find_BG_fn, &r->iterations);
	@<Choose the best node of the |b| and |g| simplex@>@;

	@<Free simplex data structures@>@;
	@<Put final values in result@>@;
}

@ A very simple start for variation of |b| and |g|.  This should
work fine for the cases in which the absorption or scattering
are fixed. 

@ @<Initialize the nodes of the |b| and |g| simplex@>=
	{
	int k,kk;
	
	p[1][1] = b2bcalc(guess[0].b);
	p[1][2] = g2gcalc(guess[0].g);
	
	for (k=1;k<7;k++) {
		if (guess[0].b != guess[k].b) 
			break;
	}
	
	p[2][1] = b2bcalc(guess[k].b);
	p[2][2] = g2gcalc(guess[k].g);
	
	for (kk=1;kk<7;kk++) {
		if (guess[0].g != guess[kk].g && guess[k].g != guess[kk].g) 
			break;
	}
	p[3][1] = b2bcalc(guess[kk].b);	
	p[3][2] = g2gcalc(guess[kk].g);

	if (Debug(DEBUG_BEST_GUESS)) {
		fprintf(stderr, "guess 1");
		fprintf(stderr, "%10.5f ", guess[0].a);
		fprintf(stderr, "%10.5f ", guess[0].b);
		fprintf(stderr, "%10.5f ", guess[0].g);
		fprintf(stderr, "%10.5f\n", guess[0].distance);
		fprintf(stderr, "guess 2");
		fprintf(stderr, "%10.5f ", guess[k].a);
		fprintf(stderr, "%10.5f ", guess[k].b);
		fprintf(stderr, "%10.5f ", guess[k].g);
		fprintf(stderr, "%10.5f\n", guess[k].distance);
		fprintf(stderr, "guess 3");
		fprintf(stderr, "%10.5f ", guess[kk].a);
		fprintf(stderr, "%10.5f ", guess[kk].b);
		fprintf(stderr, "%10.5f ", guess[kk].g);
		fprintf(stderr, "%10.5f\n", guess[kk].distance);
	}
}

@ @<Evaluate the |bg| simplex at the nodes@>=

  for (i = 1; i <= 3; i++) {
      x[1] = p[i][1];
      x[2] = p[i][2];
      y[i] = Find_BG_fn(x);
  }

@ Here we find the node of the simplex that gave the best
result and save that one.  At the same time we save the whole
simplex for later use if needed.

@<Choose the best node of the |b| and |g| simplex@>=
    r->final_distance=10;
	for(i=1;i<=3;i++) {
		if (y[i] < r->final_distance) {
			r->slab.b = bcalc2b(p[i][1]);
			r->slab.g = gcalc2g(p[i][2]);
			r->final_distance = y[i];
		}
    }

@*2 Fixed Scattering.
Here I assume that a constant $b_s$,
$$
b_s = \mu_s d
$$
where $d$ is the physical thickness of the sample and $\mu_s$ is
of course the absorption coefficient.  This is just like 
|U_Find_BG| except that $b_a=\mu_a d$ is varied instead of
$b$.  

@<Prototype for |U_Find_BaG|@>=
	void U_Find_BaG(struct measure_type m, struct invert_type *r)
	
@ @<Definition for |U_Find_BaG|@>=
@<Prototype for |U_Find_BaG|@> 
{	
	@<Allocate local simplex variables@>@;
	Set_Calc_State(m, *r);
	@<Get the initial |a|, |b|, and |g|@>@;
	@<Initialize the nodes of the |ba| and |g| simplex@>@;
	@<Evaluate the |BaG| simplex at the nodes@>@;
	amoeba(p, y, 2, r->tolerance, Find_BaG_fn, &r->iterations);
	@<Choose the best node of the |ba| and |g| simplex@>@;

	@<Free simplex data structures@>@;
	@<Put final values in result@>@;
}

@ @<Initialize the nodes of the |ba| and |g| simplex@>=

	if (guess[0].b > r->default_bs) {
		p[1][1] = b2bcalc(guess[0].b - r->default_bs);
		p[2][1] = b2bcalc(2 * (guess[0].b - r->default_bs));
		p[3][1] = p[1][1];	
	} else {
		p[1][1] = b2bcalc(0.0001);
		p[2][1] = b2bcalc(0.001);
		p[3][1] = p[1][1];	
	}
	
	p[1][2] = g2gcalc(guess[0].g);
	p[2][2] = p[1][2];
	p[3][2] = g2gcalc(0.9*guess[0].g+0.05);

@ @<Evaluate the |BaG| simplex at the nodes@>=

  for (i = 1; i <= 3; i++) {
      x[1] = p[i][1];
      x[2] = p[i][2];
      y[i] = Find_BaG_fn(x);
  }

@ Here we find the node of the simplex that gave the best
result and save that one.  At the same time we save the whole
simplex for later use if needed.

@<Choose the best node of the |ba| and |g| simplex@>=
    r->final_distance=10;
	for(i=1;i<=3;i++) {
		if (y[i] < r->final_distance) {
			r->slab.b = bcalc2b(p[i][1]) + r->default_bs;
			r->slab.a = r->default_bs/r->slab.b;
			r->slab.g = gcalc2g(p[i][2]);
			r->final_distance = y[i];
		}
    }

@*2 Fixed Absorption.
Here I assume that a constant $b_a$,
$$
b_a = \mu_a d
$$
where $d$ is the physical thickness of the sample and $\mu_a$ is
of course the absorption coefficient.  This is just like 
|U_Find_BG| except that $b_s=\mu_s d$ is varied instead of
$b$.

@<Prototype for |U_Find_BsG|@>=
	void U_Find_BsG(struct measure_type m, struct invert_type *r)
	
@ @<Definition for |U_Find_BsG|@>=
@<Prototype for |U_Find_BsG|@> 
{	
	@<Allocate local simplex variables@>@;

	if (Debug(DEBUG_SEARCH)) {
		fprintf(stderr,"In U_Find_BsG");
		fprintf(stderr," (mu=%6.4f)",r->slab.cos_angle);
		if (r->default_ba != UNINITIALIZED) 
			fprintf(stderr,"  default_ba = %8.5f", r->default_ba);
		fprintf(stderr,"\n");
	}

	Set_Calc_State(m, *r);
	@<Get the initial |a|, |b|, and |g|@>@;
	@<Initialize the nodes of the |bs| and |g| simplex@>@;
	@<Evaluate the |BsG| simplex at the nodes@>@;
	amoeba(p, y, 2, r->tolerance, Find_BsG_fn, &r->iterations);
	@<Choose the best node of the |bs| and |g| simplex@>@;

	@<Free simplex data structures@>@;
	@<Put final values in result@>@;
}

@ @<Initialize the nodes of the |bs| and |g| simplex@>=

	p[1][1] = b2bcalc(guess[0].b - r->default_ba);
	p[1][2] = g2gcalc(guess[0].g);

	p[2][1] = b2bcalc(2 * guess[0].b - 2 * r->default_ba);
	p[2][2] = p[1][2];

	p[3][1] = p[1][1];	
	p[3][2] = g2gcalc(0.9*guess[0].g+0.05);


@ @<Evaluate the |BsG| simplex at the nodes@>=

  for (i = 1; i <= 3; i++) {
      x[1] = p[i][1];
      x[2] = p[i][2];
      y[i] = Find_BsG_fn(x);
  }

@ @<Choose the best node of the |bs| and |g| simplex@>=
    r->final_distance=10;
	for(i=1;i<=3;i++) {
		if (y[i] < r->final_distance) {
			r->slab.b = bcalc2b(p[i][1]) + r->default_ba;
			r->slab.a = 1-r->default_ba/r->slab.b;
			r->slab.g = gcalc2g(p[i][2]);
			r->final_distance = y[i];
		}
    }
