/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1997--2000  Robert Gentleman, Ross Ihaka and the
 *                            R Development Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <search.h>

#include <Rinternals.h>
#include <R_ext/Random.h>
#include <Rmath.h> /* for rxxx functions, MATH_CHECK  */

/* Equal probability sampling; without-replacement case */
/* From R/src/main/random.c */

static void SampleNoReplace(int k, int n, int *y, int *x)
{
    int i, j;
    for (i = 0; i < n; i++)
	x[i] = i;
    for (i = 0; i < k; i++) {
	j = n * unif_rand();
	y[i] = x[j] + 1;
	x[j] = x[--n];
    }
}

/* Comparison function for arrays of integers
 * to be used with lsearch */

int compare_integers (const void *a, const void *b)
{
    const int *ia = (const int *) a;
    const int *ib = (const int *) b;
    return (*ia > *ib) - (*ia < *ib);
}

/*
 *  Sampford Sampling.
 *  Unequal probability sampling without replacement.
 *  Probability of sample is proportional to the product
 *  of the probabilities of its elements.
 *
 *  Modelled after algorithm described in page 263 of:
 *    William G. Cochran, "Sampling techniques",
 *    3rd Edition, Wiley 1977.
 */

static void ProbSampleSampford(int n, double *p, int *perm,
			       int nans, int *ans)
{
    double random, mass, totalmass;
    int i, j, nm1=n-1, result;
    unsigned int nmemb;
    double *p2 = (double *) malloc (n*sizeof(double));
    
    /* Record element identities */
    for (i = 0; i < n; i++) 
	perm[i] = i + 1;
    
    /* Sort probabilities into descending order */
    /* Order element identities in parallel */
    revsort(p, perm, n);
    
    /* Loop until non-repeated units */
    totalmass = 0;
    for (j = 0; j < n; j++) {
	totalmass += (p2[j] = p[j]/(1-nans*p[j]));
    }  
    do {
	random = unif_rand();
	mass = 0;
	for (j = 0; j < nm1; j++) {
	    mass += p[j];
	    if (random <= mass)
		break;
	}
	ans[0] = perm[j];
	for (i=1; i<nans; i++) {
	    random = totalmass*unif_rand();
	    mass = 0;
	    for (j = 0; j < nm1; j++) {
		mass += p2[j];
		if (random <= mass)
		    break;
	    }
	    nmemb = i;
	    result = (int) lsearch (&perm[j], ans, &nmemb, 
				    sizeof(int), compare_integers);
	    if (nmemb == i) break;
	}
    } while (i < nans);
}

/* sampfle - equal probability sampling with/without replacement. */
/* Implements sample(n, k, r) - choose k elements from 1 to n */
/* with/without replacement according to r. */
/* adapted from "do_sample" */
static void FixupProbSampford(double *p, int n, int k)
{
    double sum;
    int i, npos;
    npos = 0;
    sum = 0.;
    for (i = 0; i < n; i++) {
	if (!R_FINITE(p[i]))
	    error("NA in probability vector");
	if (p[i] < 0)
	    error("non-positive probability");
	if (p[i] > 0)
	    npos++;
	sum += p[i];
    }
    if (npos == 0 || (k > npos))
	error("insufficient positive probabilities");
    for (i = 0; i < n; i++)
	p[i] = p[i] / sum;
}

SEXP sampfle(SEXP args)
{
    SEXP x, y, prob;
    int i, k, n;
    args = CDR(args);
    n = asInteger(CAR(args)); args = CDR(args);
    k = asInteger(CAR(args)); args = CDR(args);
    prob = CAR(args);
    if (n == NA_INTEGER || n < 1)
	error("invalid first argument");
    if (k == NA_INTEGER || k < 0)
	error("invalid second argument");
    if (k > n)
	error("can't take a sample larger than the population\n with Sampford sampling");
    GetRNGstate();
    PROTECT(y = allocVector(INTSXP, k));
    if (!isNull(prob)) {
      prob = coerceVector(prob, REALSXP);
      if (NAMED(prob)) prob = duplicate(prob);
      PROTECT(prob);
      if (length(prob) != n)
	  error("incorrect number of probabilities");
      FixupProbSampford(REAL(prob), n, k);
      PROTECT(x = allocVector(INTSXP, n));
      for (i = 0; i < n; i++) {
	  if (k*(REAL(prob)[i])>=1) 
	      error("probabilities do not fulfil hypotheses");
      }
      ProbSampleSampford(n, REAL(prob), INTEGER(x), k, INTEGER(y));
      UNPROTECT(2);
    }
    else {
	x = allocVector(INTSXP, n);
	SampleNoReplace(k, n, INTEGER(y), INTEGER(x));
    }
    PutRNGstate();
    UNPROTECT(1);
    return y;
}
