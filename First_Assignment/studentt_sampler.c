#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int main(void) {
	const gsl_rng_type * T;
	gsl_rng * r;
	
	double nu = 3.0; /* The order for student-t dist */
	int i, n = 10;
	
	r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r,1); /* Sets the seed for generator r */
	
	for (i=0;i<n;i++) {
		unsigned int k = gsl_ran_poisson(r,nu);
		printf(" %u\n",k);
	}
	
	printf("\n");
	gsl_rng_free(r);
	return 0;
}