#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "headerfile.h"
#include "energy_src/nrutil.h"


void reset_guess_vals(struct params *p)
/*==============================================================================

  Purpose: This function is called if the scaling of R,eta,delta is poor (i.e.
  if scaledx values in the minimization procedure are consistently not between 
  -1 and 1 for at least one of the R,eta,delta components (usually it is R).
  The function assumes that the current values of R,eta,delta are okay (as 
  they've been optimized somewhat by the optimization procedure) and only the
  scaling itself needs to be reset.
  ============================================================================*/
{

  double current_R = p->R;
  double current_R_c = p->R_c;
  double current_eta = p->eta;
  double current_delta = p->delta;
  
  p->Rguess = current_R;
  p->Rupper = 1.5*p->Rguess;
  p->Rlower = 0.75*p->Rguess;

  p->R_cguess = current_R_c;
  p->R_cupper = p->R_cguess+0.01;
  p->R_clower = 0.75*p->R_cguess;

  if (current_eta == current_eta) { // if eta is not NAN
    p->etaguess = current_eta;
    p->etaupper = p->etaguess+0.1;
    p->etalower = p->etaguess-0.02;
  }
  if (current_delta == current_delta && fabs(current_delta)>DELTA_CLOSE_TO_ZERO) {
    // if delta is not NAN and not close to zero
    p->deltaguess = current_delta;
    p->deltaupper = 0.818;
    if (p->deltaguess < 0.81) {
      p->deltalower = 0.95*p->deltaguess;
    } else {
      p->deltalower = 0.81;
    }
  }
  
  return;
}


void save_observables(FILE *observables,double E,struct params *p)
{
  fprintf(observables,"%15.8e\t%15.8e\t%15.8e\t%15.8e\t%15.8e\t%15.8e\n",
	  E,p->R,p->R_c,p->eta,p->delta,p->y[1][p->mpt]);
  return;
}



void save_psivsr(FILE *psivsr,struct params *p)
{

  double rf_fibril(double x,struct params *p);

  int i;

  for (i = 1; i <= p->mpt; i++) {
    fprintf(psivsr,"%13.6e\t%13.6e\t%13.6e\t%13.6e\n",p->r[i],p->y[1][i],p->y[2][i],
	    rf_fibril(p->r[i],p));
  }
  printf("psi(R) = %1.2e\n",p->y[1][p->mpt]);
  return;
}



void set_NAN(double *E,struct params *p)
{

  p->R = sqrt(-1);
  p->R_c = sqrt(-1);
  p->eta = sqrt(-1);
  p->delta = sqrt(-1);
  *E = sqrt(-1);
  
  return;
}

void initialize_R_R_c_eta_delta(struct params *p)
{
  p->R = p->Rguess;
  p->R_c = p->R_cguess;
  p->eta = p->etaguess;
  p->delta = p->deltaguess;
  return;
}

void propagate_r(double *r, double R,int M0)
{

  double h = R/(M0-1);    // compute stepsize in r[1..mpt] 
  // reset r with evenly spaced points from 0 to R.
  for (int k=1;k <=M0; k++) r[k] = (k-1)*h; 
  return;
}

void initialize_param_vectors(struct params *p)
{

  void propagate_r(double *r, double R,int M0);
  void constantGuess(double *r,double **y,double val,int mpt);


  p->r = vector(1,p->M0+1);
  p->y = matrix(1,NE,1,p->M0+1);
  p->s = matrix(1,NSI,1,NSJ);
  p->c = f3tensor(1,NCI,1,NCJ,1,p->M0+2);

  propagate_r(p->r,p->R,p->M0);
  constantGuess(p->r,p->y,0.3,p->M0); //linear initial guess for psi(r)

  return;
}


void initialize_param_vectors_linearguess(struct params *p)
{

  void propagate_r(double *r, double R,int M0);
  void linearGuess(double *r, double **y, double initialSlope,int mpt);

  p->r = vector(1,MAX_M);
  p->y = matrix(1,NE,1,MAX_M);
  p->s = matrix(1,NSI,1,NSJ);
  p->c = f3tensor(1,NCI,1,NCJ,1,MAX_M+1);

  propagate_r(p->r,p->R,p->mpt);

  double slopeguess = M_PI/(4.0*p->Rguess);

  linearGuess(p->r,p->y,slopeguess,p->mpt); //linear initial guess for psi(r)

  return;
}

void initialize_params(struct params *p,char **args)
{

  void print_params(struct params *p);

  sscanf(args[2],"%lf",&p->K33);
  sscanf(args[3],"%lf",&p->k24);
  sscanf(args[4],"%lf",&p->Lambda);
  sscanf(args[5],"%lf",&p->omega);
  sscanf(args[6],"%lf",&p->gamma_s);
  sscanf(args[7],"%lf",&p->Rguess);
  sscanf(args[8],"%lf",&p->R_cguess);
  sscanf(args[9],"%lf",&p->etaguess);
  sscanf(args[10],"%lf",&p->deltaguess);
  sscanf(args[11],"%lf",&p->Rupper);
  sscanf(args[12],"%lf",&p->Rlower);
  sscanf(args[13],"%lf",&p->R_cupper);
  sscanf(args[14],"%lf",&p->R_clower);
  sscanf(args[15],"%lf",&p->etaupper);
  sscanf(args[16],"%lf",&p->etalower);
  sscanf(args[17],"%lf",&p->deltaupper);
  sscanf(args[18],"%lf",&p->deltalower);

  p->M0 = MAX_M-1;
  p->mpt = p->M0;

  print_params(p);

  return;

}

void initialize_params_withstrain(struct params *p,char **args,double *strain)
{

  void print_params(struct params *p);

  sscanf(args[2],"%lf",&p->K33);
  sscanf(args[3],"%lf",&p->k24);
  sscanf(args[4],"%lf",&p->Lambda);
  sscanf(args[5],"%lf",&p->omega);
  sscanf(args[6],"%lf",&p->gamma_s);
  sscanf(args[7],"%lf",&p->Rguess);
  sscanf(args[8],"%lf",&p->R_cguess);
  sscanf(args[9],"%lf",&p->etaguess);
  sscanf(args[10],"%lf",&p->deltaguess);
  sscanf(args[11],"%lf",&p->Rupper);
  sscanf(args[12],"%lf",&p->Rlower);
  sscanf(args[13],"%lf",&p->R_cupper);
  sscanf(args[14],"%lf",&p->R_clower);
  sscanf(args[15],"%lf",&p->etaupper);
  sscanf(args[16],"%lf",&p->etalower);
  sscanf(args[17],"%lf",&p->deltaupper);
  sscanf(args[18],"%lf",&p->deltalower);
  sscanf(args[19],"%lf",strain);

  p->M0 = MAX_M-1;
  p->mpt = p->M0;

  print_params(p);


  printf("strain = %e\n",*strain);

  return;

}

void print_params(struct params *p)
{

  printf("\n\nparameter values for minimization:\n");

  printf("K33 = %e\n",p->K33);
  printf("k24 = %e\n",p->k24);
  printf("Lambda = %e\n",p->Lambda);
  printf("omega = %e\n",p->omega);
  printf("gamma_s = %e\n",p->gamma_s);


  printf("initial guesses for R, R_c, eta, and delta:\n");

  printf("guess for R = %e\n",p->Rguess);
  printf("guess for R_c = %e\n",p->R_cguess);
  printf("guess for eta = %e\n",p->etaguess);
  printf("guess for delta = %e\n",p->deltaguess);


  printf("estimated ranges of R, eta, and delta:\n");

  printf("upper bound estimate of R = %e\n",p->Rupper);
  printf("lower bound estimate of R = %e\n",p->Rlower);
  printf("upper bound estimate of R_c = %e\n",p->R_cupper);
  printf("lower bound estimate of R_c = %e\n",p->R_clower);
  printf("upper bound estimate of eta = %e\n",p->etaupper);
  printf("lower bound estimate of eta = %e\n",p->etalower);
  printf("upper bound estimate of delta = %e\n",p->deltaupper);
  printf("lower bound estimate of delta = %e\n",p->deltalower);

  return;
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

  p->M0 = p->mpt;

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


void initialize_file(FILE **output,char *path,char *fname,struct params p)
{

  int num_chars = 400;
  char suffix[num_chars];
  char f[num_chars];

  snprintf(suffix,sizeof(suffix),"%1.4e_%1.4e_%1.4e_%1.4e_"
	   "%1.4e.txt",p.K33,p.k24,p.Lambda,p.omega,p.gamma_s);
  snprintf(f,sizeof(f),"%s_%s_%s",path,fname,suffix);
  *output = fopen(f,"w");

  return;
}

void psivsrstrain_filename(FILE **output,char *path,char *fname,struct params p,
			   double strain)
{

  int num_chars = 400;
  char suffix[num_chars];
  char f[num_chars];

  snprintf(suffix,sizeof(suffix),"%1.4e_%1.4e_%1.4e_%1.4e_"
	   "%1.4e_%1.4e.txt",p.K33,p.k24,p.Lambda,p.omega,p.gamma_s,
	   strain);
  snprintf(f,sizeof(f),"%s_%s_%s",path,fname,suffix);
  *output = fopen(f,"w");

  return;
}
