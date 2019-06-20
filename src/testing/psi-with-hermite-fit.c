#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include "../edited_gsl_src/gsl_multimin.h"
#include "../energy_src/nrutil.h"
#include "../headerfile.h"


int main(int argc, char **argv)
{
  void initialize_params_mpt(struct params *p,char **args);

  void initialize_param_vectors(struct params *p);

  void initialize_R_R_c_eta_delta(struct params *p);

  void initialize_file(FILE **output,char *path,char *fname,struct params p);

  double E_calc(struct params *p);

  void save_psivsr_with_hermite(FILE *psivsr,struct params *p);

  
  struct params p; 
  initialize_params(&p,argv);
  initialize_param_vectors(&p);

  initialize_R_R_c_eta_delta(&p);

  p.x_size = 4;

  FILE *psivsr;
  initialize_file(&psivsr,argv[1],"psivsr",p);

  double E;
  
  fprintf(psivsr,"# R = %e, R_c = %e, eta = %e, delta = %e\n",p.R,p.R_c,p.eta,p.delta);
  
  E = E_calc(&p);

  save_psivsr_with_hermite(psivsr,&p);
    

  free_vector(p.r,1,MAX_M);
  free_matrix(p.y,1,NE,1,MAX_M);
  free_vector(p.rf_fib,1,MAX_M);
  free_vector(p.z,1,MAX_M);
  free_matrix(p.s,1,NSI,1,NSJ);
  free_f3tensor(p.c,1,NCI,1,NCJ,1,MAX_M+1);

  fclose(psivsr);

  return 0;

}

void initialize_params_with_mpt(struct params *p,char **args)
{

  void print_only_some_params(struct params *p);
  
  sscanf(args[2],"%lf",&p->K33);
  sscanf(args[3],"%lf",&p->k24);
  sscanf(args[4],"%lf",&p->Lambda);
  sscanf(args[5],"%lf",&p->omega);
  sscanf(args[6],"%lf",&p->gamma_s);
  sscanf(args[7],"%lf",&p->Rguess);
  sscanf(args[8],"%lf",&p->R_cguess);
  sscanf(args[9],"%lf",&p->etaguess);
  sscanf(args[10],"%lf",&p->deltaguess);
  sscanf(args[11],"%d",&p->mpt);

  print_only_some_params(p);
  
  return;

}


void print_only_some_params(struct params *p)
{

  printf("\n\nparameter values for minimization:\n");

  printf("K33 = %e\n",p->K33);
  printf("k24 = %e\n",p->k24);
  printf("Lambda = %e\n",p->Lambda);
  printf("omega = %e\n",p->omega);
  printf("gamma_s = %e\n",p->gamma_s);


  printf("initial guesses for R, R_c, eta, and delta:\n");

  printf("initial R = %e\n",p->Rguess);
  printf("initial R_c = %e\n",p->R_cguess);
  printf("initial for eta = %e\n",p->etaguess);
  printf("initial for delta = %e\n",p->deltaguess);

  printf("mpt=%d\n",p->mpt);

  return;
}

void save_psivsr_with_hermite(FILE *psivsr,struct params *p)
{
  int i;

  for (i = 1; i <= p->mpt; i++) {
    fprintf(psivsr,"%13.6e\t%13.6e\t%13.6e\t%13.6e\n",p->r[i],p->y[1][i],p->y[2][i],
	    p->rf_fib[i]);
  }
  printf("psi(R) = %1.2e\n",p->y[1][p->mpt]);
  return;
}
