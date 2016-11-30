#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

/* A data type to specify which parameters will be stored */
struct param_type {
	double nu;
	double mu;
	double sigma;
};

double lklhood(double x, void *params) {
	/*
	 * Student-t disitribution.
	 */ 
	 
	/* caste params to the struct type desired */
	struct param_type *params_ptr = params;
	/* access the parameters. The arrow -> is used bc params_ptr is
	 * a pointer and needs to be dereference them */
	double nu  = params_ptr->nu;
	double mu = params_ptr->mu;
	double sigma = params_ptr->sigma;
	
	double arg = (x-mu)/sigma;
	
	return gsl_ran_tdist_pdf(arg,nu)/sigma;
}



int main(int argc, char *argv[]) {
	/* Initialization */
	const gsl_rng_type * T;
	gsl_rng * r;
	
	
	/* set iteration variables and the order of the student-t distribution
	 * from the command line arguments */
	int i, itr = atoi(argv[1]);
	
	/* parameters of student t distributions */
	double nu = atof(argv[2]); 
	double mu = atof(argv[3]);
	double sigma = atof(argv[4]);
	
	/* store the parameters in param_type struct */
	struct param_type params = {nu,mu,sigma};
	
	/* open text file for writing  and make sure it works*/
	
	FILE *f = fopen(argv[5], "w"); 
	if (f == NULL) {
    	printf("Error opening file!\n");
    	exit(1);
	}
	
	/* allocate memory for generator and set its seed */
	r = gsl_rng_alloc(gsl_rng_mt19937); 
	gsl_rng_set(r,1); 
	
	/* Start initial value */
	double x_cur = 1.0; 
	double proposal_sigma = atof(argv[6]);
	double alpha;
	double x_prop;
	int accept; /* keep track of acceptance rate */
	double u; /* to make decision of accept proposal or not */
	double accept_prob;
	
	/* Start the MCMC */
	for (i=0;i<itr;i++) {
		/* propose a new x */
		x_prop = gsl_ran_gaussian(r,proposal_sigma) + x_cur;
		
		/* Calculate acceptance probability */
		accept_prob = lklhood(x_prop, &params)/lklhood(x_cur, &params);
		alpha = GSL_MIN(1.0,accept_prob);
		
		/* Accept or not, decide */
		u = gsl_ran_flat(r,0.0,1.0);
		if (u < alpha) {
			x_cur = x_prop;
			accept = 1;
		}/* print to data file */
		else {
			accept = 0;
		}
		fprintf(f," %.5f %i\n",x_cur,accept);
	} 
	
	/* Clean up time */
	fclose(f);
	gsl_rng_free(r);
	
	return 0;
}