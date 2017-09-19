/*
 * In article <MATT.93Jan10152553@physics2.berkeley.edu>,
 * matt@physics2.berkeley.edu (Matt Austern) writes:
 *
 * > Does anybody out there have a recommendation for a minimization > routine?
 * I'll be using it to minimize a function of two or three > variables; the
 * function is reasonably well behaved, but I am not able > to calculate the
 * derivatives very easily. > Code in C or C++ would be preferable.
 *
 * For small numbers of variables (N<12), I find that the Hooke and Jeeves
 * algorithm [1963] works quite well.  Its good points are
 *
 *  1* Doesn't require the user to supply derivatives
 *
 *  2* Will minimize discontinuous functions
 *
 *  3* (A consequence of *2*) Will minimize functions whose derivatives are
 * discontinuous
 *
 *  4* Very short subroutine (2700 bytes), easily transliterated from one
 * programming language to another
 *
 *  5* (A consequence of *2*) One of the very quickest ways of getting a problem
 * up and running.
 *
 *
 * The algorithm has VERY MUCH fallen into academic disfavour.  Principally this
 * is because it is a _heuristic_, and because no proofs are known of its
 * rates of convergence.  The darlings of the I-gotta-see-a-proof community
 * are called "Huang's class of Variable Metric" algorithms, and their
 * parentage is quasi-Newton.  Generally these require that the function have
 * continuous derivatives.  Often they require that the subroutine user
 * _supply_ the derivatives, leading to disadvantage 5*.  However, for
 * certain assumptions (including the cheerfully unrealistic one "for
 * starting points sufficiently close to the optimum"), Huang's class has
 * provably quadratic convergence.
 *
 * So, if microseconds are life-or-death, if you must have the very very very
 * fastest possible algorithm, prepare to use a quasi-Newton method.  In this
 * case I recommend BFGS or LM, the latter available from netlib (send index
 * from minpack).
 *
 * On the other hand, if you can accept a runtime of 8 seconds instead of 5, and
 * if you want avoid derivatives, and if you are using small numbers of
 * variables, then it's worth giving Hooke and Jeeves a try.
 *
 * Pragmatically speaking, quoting Fred Dryer, "It works for me."  I use it for
 * modeling semiconductor devices, especially short-channel MOSFETs.
 *
 * Technical note: because H&J accepts functions with discontinuous derivatives,
 * it mates beautifully with penalty functions to give quick and easy
 * _constrained_ minimization.  And you don't have to sweat the details of
 * slope-matching the penalties or barriers.
 *
 */
 
 /* hooke() has been rewritten so that it can be a drop in replacement for amoeba() */
 
#include <stdio.h>
#include <math.h>
#define VARS        (2)         /* max # of variables */
#define RBEGIN      (0.6)       /* stepsize geometric shrink   */
#define NMAX 5000

/* given a point, look for a better one nearby, one coord at a time */
double best_nearby(double (*funk) (double[]), int *evals,
                   double delta[VARS], double point[VARS], double prevbest, int nvars)
{
    double z[3];
    double minf, ftmp;
    int i;
    *evals = 0;
    
    minf = prevbest;
    for (i = 0; i < nvars; i++)
    	z[i+1] = point[i];
    	
    for (i = 0; i < nvars; i++) {
    	z[i+1] = point[i] + delta[i];
		ftmp = funk(z);
		(*evals)++;
		if (ftmp < minf)
			minf = ftmp;
		else {
			delta[i] = 0.0 - delta[i];
			z[i+1] = point[i] + delta[i];
			ftmp = funk(z);
			(*evals)++;
			if (ftmp < minf)
				minf = ftmp;
			else
				z[i+1] = point[i];
		}
    }
    for (i = 0; i < nvars; i++)
    	point[i] = z[i+1];
    return (minf);
}

void hooke(double **p, double y[], int nvars, double epsilon, double (*funk)(double []), int *nfunk)
{
    double newf, fbefore, steplength, tmp;
    double xbefore[VARS], newx[VARS], delta[VARS];
    double rho = RBEGIN;
    int i, keep, evals;

    for (i = 0; i < nvars; i++) {
    	xbefore[i] = p[1][i+1];
    	delta[i]   = p[2+i][i+1] - p[1][i+1];
    }
        
    fbefore = y[1];
    newf = fbefore;
    
	*nfunk = 0;
    steplength = rho;
    
    while ((*nfunk < NMAX) && (steplength > epsilon)) {

/*
 		fprintf(stderr, "\nAfter %5d funk() evals, funk(x) =  %.7f at\n", *nfunk, fbefore);
        for (j = 0; j < nvars; j++)
            fprintf(stderr, "   x[%2d] = %.7f\n", j, xbefore[j]);
*/

        /* find best new point, one coord at a time */
        for (i = 0; i < nvars; i++) 
            newx[i] = xbefore[i];
        
        newf = best_nearby(funk, &evals, delta, newx, fbefore, nvars);
        *nfunk = *nfunk + evals;
        
        /* if we made some improvements, pursue that direction */
        keep = 1;
        while ((newf < fbefore) && (keep == 1)) {

            for (i = 0; i < nvars; i++) {
    
                /* firstly, arrange the sign of delta[] */
                if (newx[i] <= xbefore[i])
                    delta[i] = 0.0 - fabs(delta[i]);
                else
                    delta[i] = fabs(delta[i]);
                    
                /* now, move further in this direction */
                tmp = xbefore[i];
                xbefore[i] = newx[i];
                newx[i] = newx[i] + newx[i] - tmp;
            }
            fbefore = newf;
            newf = best_nearby(funk, &evals, delta, newx, fbefore, nvars);
        	*nfunk = *nfunk + evals;
            
            /* if the further (optimistic) move was bad.... */
            if (newf >= fbefore)
                break;
            
            /* make sure that the differences between the new */
            /* and the old points are due to actual */
            /* displacements; beware of roundoff errors that */
            /* might cause newf < fbefore */
            keep = 0;
            for (i = 0; i < nvars; i++) {
                keep = 1;
                if (fabs(newx[i] - xbefore[i]) >
                    (0.5 * fabs(delta[i])))
                    break;
                else
                    keep = 0;
            }
        }
        if ((steplength >= epsilon) && (newf >= fbefore)) {
            steplength = steplength * rho;
            for (i = 0; i < nvars; i++) 
                delta[i] *= rho;
        }
    }
    for (i = 0; i < nvars; i++) 
        p[1][i+1] = xbefore[i];
    
    y[1] = fbefore;
/*	fprintf(stderr, "<%5d>\n", *nfunk); */
}
