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
  void initialize_params(struct params *p,char **args);

  void initialize_param_vectors(struct params *p);

  void initialize_R_R_c_eta_delta(struct params *p);

  void initialize_file(FILE **output,char *path,char *fname,struct params p);

  double E_calc(struct params *p);

  void save_psivsr(FILE *psivsr,struct params *p);

  void save_psivsr_with_hermite(FILE *psivsr,struct params *p);
  


  struct params p; 
  initialize_params(&p,argv);
  initialize_param_vectors(&p);

  initialize_R_R_c_eta_delta(&p);

  p.x_size = 4;
  p.mpt = 1024+1;

  FILE *psivsr;
  initialize_file(&psivsr,argv[1],"psivsr",p);

  FILE *hermite_psivsr;
  initialize_file(&hermite_psivsr,argv[1],"hermite-psivsr",p);

  double E;
  
  fprintf(psivsr,"# R = %e, R_c = %e, eta = %e, delta = %e\n",p.R,p.R_c,p.eta,p.delta);

  fprintf(hermite_psivsr,"# R = %e, R_c = %e, eta = %e, delta = %e\n",p.R,p.R_c,p.eta,p.delta);

  E = E_calc(&p);

  save_psivsr(psivsr,&p);

  save_psivsr_with_hermite(hermite_psivsr,&p);

  free_vector(p.r,1,MAX_M);
  free_matrix(p.y,1,NE,1,MAX_M);
  free_matrix(p.s,1,NSI,1,NSJ);
  free_f3tensor(p.c,1,NCI,1,NCJ,1,MAX_M+1);

  fclose(psivsr);
  fclose(hermite_psivsr);

  return 0;

}


void save_psivsr_with_hermite(FILE *psivsr,struct params *p)
{

  void psi_interp(double x,double *r, double **y,double *psi,double *psiprime);

  double rf_fibril(double x,struct params *p);

  double R = p->R;

  double M = 4*(p->mpt-1)+1;

  double dr = R/(M-1);

  double x,psi,psiprime;

  for (int i = 1; i <= M; i++) {

    x = dr*(i-1);

    psi_interp(x,p->r,p->y,&psi,&psiprime);

    fprintf(psivsr,"%13.6e\t%13.6e\t%13.6e\t%13.6e\n",x,psi,psiprime,
	    rf_fibril(x,p));

  }

  return;
}
