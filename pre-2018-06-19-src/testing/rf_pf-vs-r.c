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

  double cubic_spline(double t,double *x,double *y, double *z);
  
  struct params p; 
  initialize_params(&p,argv);
  initialize_param_vectors(&p);

  initialize_R_R_c_eta_delta(&p);

  p.x_size = 4;

  FILE *rf_pf;
  FILE *splines;
  initialize_file(&rf_pf,argv[1],"rf_pf",p);
  initialize_file(&splines,argv[1],"splines",p);

  //  FILE *psivsr;
  //  initialize_file(&psivsr,argv[1],"psivsr",p);

  double E;

  fprintf(rf_pf,"# R = %e, R_c = %e, eta = %e, delta = %e\n",p.R,p.R_c,p.eta,p.delta);
    
  E = E_calc(&p);

  double t;

  double t0 = p.r[p.i_c];
  double tf = p.r[p.mpt];
  int N = 100000;
  double dt = (tf-t0)/(N-1);

  for (int i = p.i_c; i <= p.mpt; i++) {  
    
    fprintf(rf_pf,"%15.8e\t%15.8e\n",p.r[i],p.rf_fib[i]);

  }

  for (int i = 0; i < N; i++) {

    t = t0+i*dt;
    
    fprintf(splines,"%15.8e\t%15.8e\n",t,
	    cubic_spline(t,p.r+1,p.rf_fib+1,p.z+1));

  }

  free_vector(p.r,1,MAX_M);
  free_matrix(p.y,1,NE,1,MAX_M);
  free_vector(p.rf_fib,1,MAX_M);
  free_vector(p.z,1,MAX_M);
  free_matrix(p.s,1,NSI,1,NSJ);
  free_f3tensor(p.c,1,NCI,1,NCJ,1,MAX_M+1);

  fclose(splines);
  fclose(rf_pf);

  return 0;

}
