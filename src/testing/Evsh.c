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
  void initialize_params_with_mpt(struct params *p,char **args);

  void initialize_param_vectors(struct params *p);

  void initialize_R_R_c_eta_delta(struct params *p);

  void initialize_file(FILE **output,char *path,char *fname,struct params p);

  void update_params(const gsl_vector *x,struct params *p);

  double E_calc(struct params *p);

  double f_dumb(const gsl_vector *placeholder,void *ps);

  int deriv_xi(double (*f)(const gsl_vector *,void *),const gsl_vector *x,
	       int i,void *ps,double h,double *result,double *abserr);

  void c_deriv(double (*f)(const gsl_vector *,void *),gsl_vector *x,int i,
	       void *ps,double h,double *result,double *abserr_round,
	       double *abserr_trunc);
  
  struct params p; 
  initialize_params_with_mpt(&p,argv);
  initialize_param_vectors(&p);

  initialize_R_R_c_eta_delta(&p);

  p.x_size = 4;

  char fpre[400];

  snprintf(fpre,sizeof(fpre),"Evsh-%d",p.mpt);
  
  FILE *observables;
  initialize_file(&observables,argv[1],fpre,p);

  double E;

  double dEdR_c;
  double abserr_round,abserr_trunc;

  
  gsl_vector *x = gsl_vector_alloc(p.x_size);

  gsl_vector_set(x,0,p.R);
  gsl_vector_set(x,1,p.R_c);
  gsl_vector_set(x,2,p.eta);
  gsl_vector_set(x,3,p.delta);

  fprintf(observables,"# R = %e, eta = %e, delta = %e\n",p.R,p.eta,p.delta);
  
  for (double h = 1e-12; h < 1e-2 ; h *= 2) {
    c_deriv(f_dumb,x,1,&p,h,&dEdR_c,&abserr_round,&abserr_trunc);
    update_params(x,&p);
    fprintf(observables,"%15.8e\t%15.8e\t%15.8e\t%15.8e\n",
	    h,dEdR_c,abserr_round,abserr_trunc);
  }
    

  free_vector(p.r,1,p.M0+1);
  free_matrix(p.y,1,NE,1,p.M0+1);
  free_matrix(p.s,1,NSI,1,NSJ);
  free_f3tensor(p.c,1,NCI,1,NCJ,1,p.M0+2);

  gsl_vector_free(x);

  //  fclose(psivsr);
  fclose(observables);

  return 0;

}


double f_dumb(const gsl_vector *x,void *ps)
{
  
  double E_calc(struct params *p);


  double E;
  struct params *p = ps;

  p->R = gsl_vector_get(x,0);
  p->R_c = gsl_vector_get(x,1);
  p->eta = gsl_vector_get(x,2);
  p->delta = gsl_vector_get(x,3);
  
  
  E = E_calc(p);

  return E;
}


void update_params(const gsl_vector *x,struct params *p)
{

  p->R = gsl_vector_get(x,0);
  p->R_c = gsl_vector_get(x,1);
  p->eta = gsl_vector_get(x,2);
  p->delta = gsl_vector_get(x,3);
  

  return;
}
