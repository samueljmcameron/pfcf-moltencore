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

  void save_psivsr(FILE *psivsr,struct params *p);

  void set_NAN(double *E,struct params *p);

  void save_observables(FILE *observables,double E,struct params *p);

  void reset_guess_vals(struct params *p);
  
  int full3var_driver(double *E,struct params *p,FILE *energy);

  double E_calc(struct params *p);

  double f_dumb(const gsl_vector *placeholder,void *ps);

  int deriv_xi(double (*f)(const gsl_vector *,void *),const gsl_vector *x,
	       int i,void *ps,double h,double *result,double *abserr);
  
  struct params p; 
  initialize_params(&p,argv);
  initialize_param_vectors(&p);

  initialize_R_R_c_eta_delta(&p);

  p.x_size = 4;

  FILE *observables;
  initialize_file(&observables,argv[1],"observables",p);

  //  FILE *psivsr;
  //  initialize_file(&psivsr,argv[1],"psivsr",p);

  double E;

  double dEdR_c;
  double abserr;
  double h;

  // allocate a placeholder vector, but no need to initialize it.
  gsl_vector *placeholder = gsl_vector_alloc(p.x_size);

  for (p.R_c = 0.00; p.R_c < 0.5; p.R_c += 0.05) {
    printf("%lf\n",p.R_c);
    E = E_calc(&p);
    deriv_xi(f_dumb,placeholder,3,&p,h,&dEdR_c,&abserr);
    fprintf(observables,"%15.8e\t%15.8e\t%15.8e\t%15.8e\n",p.R_c,E,dEdR_c,abserr);
  }

  //  save_psivsr(psivsr,&p);
    

  free_vector(p.r,1,MAX_M);
  free_matrix(p.y,1,NE,1,MAX_M);
  free_vector(p.r_cp,1,MAX_M);
  free_matrix(p.y_cp,1,NE,1,MAX_M);
  free_vector(p.rf_fib,1,MAX_M);
  free_vector(p.z,1,MAX_M);
  free_matrix(p.s,1,NSI,1,NSJ);
  free_f3tensor(p.c,1,NCI,1,NCJ,1,MAX_M+1);

  gsl_vector_free(placeholder);

  //  fclose(psivsr);
  fclose(observables);

  return 0;

}


double f_dumb(const gsl_vector *placeholder,void *ps)
// This function is just a dumb wrapper for the energy function. It returns
// the same value as the energy function REGARDLESS of the elements in the
// placeholder vector. It only requires that placeholder be of size p.x_size.
// 
// The sole purpose of this function is to be used as input to the numerical
// differentiation routine I have written (based off of the gsl one).
{
  
  double E_calc(struct params *p);


  double E;
  struct params *p = ps;
  
  
  E = E_calc(p);

  return E;
}


