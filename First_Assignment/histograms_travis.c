/* 
 * gcc histograms_travis.c -lm -o histograms
 *   -makes marginalized pdfs from MCMC chains
 *	 -ignores first two columns as iteration and log-likelihood
 *
 * Handed down to me from Meg (Margaret) Millhouse
 * Last modified: 12/3/16
 */

/***************************  REQUIRED LIBRARIES  ***************************/


// Below is an example 
// plot 'datafile' u 1:2 w steps (1 and 2 first parameter, 3 and 4 for second etc...)



#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define NR_END 1
#define FREE_ARG char*

/*************  PROTOTYPE DECLARATIONS FOR INTERNAL FUNCTIONS  **************/

/* ============================  MAIN PROGRAM  ============================ */

int main(int argc, char* argv[])
{
	/* -------------  DECLARATION OF VARIABLES  ------------- */
	int i, j, k;
	int count, bin;
	int Nparam; // number of parameters
	int Nhist; 
	
	double norm;
	double max;
	double bw;
	double x;
	
	char filename[100];
	
	FILE* infile;
	FILE* outfile;
	
	// Check that there are the correct number of arguments passed 
	if(argc<3) { printf("usage: ./histograms chainfile Nparam [bin width]\n"); return 0;}
	
	Nparam	= atoi(argv[2]); // second argument is the number of parameters
    //	Nparam = 13;                              //this is just set up for the ppE runs
    
    // if bin_width is specified (4th argument) set it 
	if(argc==4) bw = (double)atof(argv[3]); 
	
	double Parameters[Nparam+2]; 
	double Psigma[Nparam+2];
	double Pwidth[Nparam+2]; 
	double Pmean[Nparam+2]; // store means of each parameter
	double Pmode[Nparam+2]; // store modes of each parameter
	double Pmax[Nparam+2]; // store the max of each parameter
	double Pmin[Nparam+2]; // store the min of each parameter
	double Pmle[Nparam+2];
	
	/* initialize arrays */
	for(i=0; i<Nparam; i++)
	{
		Pmean[i] = 0.0;
		Pmin[i] = 1.0e60;
		Pmax[i] = -1.0e60;
	}
	
	
	/*  open chain file and find max and min for each parameter  */
	sprintf(filename,"%s",argv[1]); 
	infile = fopen(filename,"r"); 
	
	
	count = 0;
	while(!feof(infile)) // while not at end of infile
	{
		for(i=0; i<Nparam; i++)
		{
			fscanf(infile,"%lf",  &Parameters[i]);
		}
        
		count++;
        
		for(i=0; i<Nparam; i++)
		{
			Pmean[i] += Parameters[i];
			if(Parameters[i] < Pmin[i])  Pmin[i] = Parameters[i];
			if(Parameters[i] > Pmax[i])  Pmax[i] = Parameters[i];
		}
		
	}
	rewind(infile);
	
	if(argc==3)
	{
		Nhist = 50;
		for(i=0; i<Nparam; i++)
		{
			Pwidth[i] = (Pmax[i]-Pmin[i])/(double)(Nhist - 1);
		}
	}
	if(argc==4)
	{
		for(i=0; i<Nparam; i++)
		{
			Pwidth[i] = bw;
			Nhist = (int)((Pmax[i]-Pmin[i])/(double)(Pwidth[i]));
		}
	}
	
	double histograms[Nparam][Nhist];
	for(i=0; i<Nparam; i++) for(j=0; j<Nhist; j++) histograms[i][j] = 0.0;
	
	for(i=0; i<Nparam; i++)
	{
		Pmean[i] /= (double)(count);
		Psigma[i] = 0.0;
	}
	
	/* calculate variance */
	while(!feof(infile))
	{
        //	fscanf(infile,"%d", &j);
		for(i=0; i<Nparam; i++)
		{
			fscanf(infile,"%lf",  &Parameters[i]);
		}
		
		for(i=0; i<Nparam; i++)
		{
			Psigma[i] += (Parameters[i]-Pmean[i])*(Parameters[i]-Pmean[i]);
		}
	}
	
	rewind(infile);
    
	/*  sort chain into bins for histograms  */
	norm = 1.0/((double)(count));
	
	while(!feof(infile))
	{
        //	fscanf(infile,"%d", &j);
		for(i=0; i<Nparam; i++)
		{
			fscanf(infile,"%lf",  &Parameters[i]);
		}
		
		
		for(i=0; i<Nparam; i++)
		{
			bin = (int)( (Parameters[i] - Pmin[i])/Pwidth[i] );
			if((bin >= 0 && bin < Nhist)) histograms[i][bin] += norm/Pwidth[i];
		}
	}
	fclose(infile);
	
    
	for(i=0; i<Nparam; i++)
	{
		max = -1.0e60;
		for(j=0; j<Nhist; j++)
		{
			if(histograms[i][j] > max)
			{
				max = histograms[i][j];
				bin = j;
			}
		}
		Pmode[i] = Pmin[i] + (double)bin*Pwidth[i];
	}
	
    /*  calculate std. deviation */
	printf("\n");
	printf("Parameter  |  mean  |  mode | sigma\n");
	for(i=0; i<Nparam; i++)
	{
		Psigma[i] /= (double)(count-1);
		Psigma[i] = sqrt(Psigma[i]);
		printf("%i  | %g | %g | %g\n", i+1, Pmean[i], Pmode[i], Psigma[i]);
	}
	printf("\n");
    
//    printf("sigma_Mc = %g\n",Psigma[3]);
//    printf("sigma_Mt = %g\n",Psigma[4]);
//    printf("sigma_tc = %g\n",Psigma[5]);
//    printf("sigma_chi1 = %g\n",Psigma[7]);
//    printf("sigma_chi2 = %g\n",Psigma[8]);
    
	
	/*  write histogram file  */
	sprintf(filename,"histograms_%s",argv[1]);
	outfile = fopen(filename,"w");
	fprintf(outfile,"# ");
	for(i=0; i<Nparam; i++)fprintf(outfile,"%.12g ",Pmode[i]);
	fprintf(outfile,"\n");
	for(i=0; i<Nhist; i++)
	{
		for(j=0;j<Nparam; j++)
		{
			x = Pmin[j] + (double)i*Pwidth[j];
			fprintf(outfile,"%.12e %e ", x, histograms[j][i]);
            //		if(j==12) printf("%.12g %.12g \n", x, histograms[j][i]);
		}
		fprintf(outfile,"\n");
	}
	fclose(outfile);
	
	printf("created file %s\n\n",filename);
	
	return 0;
}
