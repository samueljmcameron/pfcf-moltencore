/* Functions to calculate the cost function, which is the energy per unit
   volume E(x), for a given value of x. */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <gsl/gsl_vector.h>
#include "nrutil.h"
#include "../headerfile.h"

#define LRG_NBR 1e10
#define CPTR_ERR 1e-15



double E_calc(struct params *p)
/*==============================================================================

  Purpose: This function calculates E(x), by first solving the ODE for psi(r)
  using relaxation methods (with solvde_wrapper function), and then integrating
  with successful_E_count. This is all done with *mpt grid points in r and y.
  If successful_E_count fails to calculate E(x), then the whole calculation 
  (for psi(r) and E(x)) is retried  with 2*((*mpt)-1)+1 grid points. This is 
  continued until either E(x) is successfully calculated, or the maximum
  number of grid points max_mpt is reached, in which case the program exits to 
  the system.

  ------------------------------------------------------------------------------

  Parameters:

  p -- This struct has all of the constant parameter info (e.g. K33, k24), as
  well as the array info (r,y,rf_fib,etc).

  ------------------------------------------------------------------------------

  Returns: Returns the value of E(x) if the calculation is successful. If the 
  calculation is unsuccessful, the form of psi(r), dpsi/dr, and r*f_fibril(r)
  are saved in a file starting with "QROMB" and then an exit status exit(1) is
  invoked, exiting to the system.

  ============================================================================*/

{

  void solvde_wrapper(double scalv[],struct params *p,bool ignore_first_y);

  void propagate_r(double *r, double R,int M0);

  void add_R_c(struct params *p);
  
  bool successful_E_count(double *E,struct params *p);

  double h;
  double E;  
  double scalv[2+1];



  scalv[1] = .1;    // guess for magnitude of the psi values
  scalv[2] = 4.0;   // guess for magnitude of the psi' values



  if (p->R <= 0 || p->eta <= 0 || p->eta >= 8.0
      || fabs(p->delta) >= 1.0 || p->R_c >= p->R || p->R_c < -1e-8) {
    printf("p->R = %e, p->R_c = %e, p->eta = %e, p->delta = %e\n",p->R,p->R_c,p->eta,p->delta);
    printf("something is too big or less than zero, so returning failed calculation.\n");
    return FAILED_E;
  }

  propagate_r(p->r,p->R,p->M0);

  add_R_c(p);

  solvde_wrapper(scalv,p,false);

  if(successful_E_count(&E,p)) return E;

  return FAILED_E;
}


bool successful_E_count(double *E,struct params *p)
/*==============================================================================

  Purpose: Given the form of psi(r) (in the array y[1..2][1..mpt]), compute the
  integrand rf_fib[1..mpt], numerically integrate this integrand, and then
  compute E(x).

  ------------------------------------------------------------------------------

  Parameters:

  *E -- pointer to variable which will store the value of E(x), if the
  calculation is successful.

  p -- This struct has all of the constant parameter info (e.g. K33, k24), as
  well as the array info (r,y,rf_fib,etc).

  ------------------------------------------------------------------------------

  Returns: Returns true if the calculation of E(x) was successful, and stores
  E(x) in *E. Returns false if E(x) integral was not successful.

  ============================================================================*/  
  
{


  double rf_fibril(double x,struct params *p);

  double qromb(double (*func)(double,struct params *),double t0,double t1,
	       struct params *p,bool *failure);

  double E_bulk_surface(struct params *p,double psiR,
			double integration_2233b1);

  
  double integration_2233b1;

  bool failure;

  if (p->R_c > CPTR_ERR) {

    integration_2233b1 = qromb(rf_fibril,0,p->R_c-CPTR_ERR,p,&failure);
  
    if (failure) {
      *E = integration_2233b1;
      printf("failed to integrate frank free energy with 0<r<R_c at"
	     " (R,R_c,eta,delta) = (%e,%e,%e,%e)\n",
	     p->R,p->R_c,p->eta,p->delta);
      return false;
    }
  } else { integration_2233b1 = 0.0;}

  integration_2233b1 += qromb(rf_fibril,p->R_c,p->R,p,&failure);

  if (failure) {
    *E = integration_2233b1;
    printf("failed to integrate frank free energy with R_c<r<R at (R,R_c,eta,delta) = (%e,%e,%e,%e)\n",
	   p->R,p->R_c,p->eta,p->delta);
    return false;
  }

  *E = E_bulk_surface(p,p->y[1][p->mpt],integration_2233b1);


  return true;
}


double E_bulk_surface(struct params *p,double psiR,
		      double integration_2233b1)
/*==============================================================================

  Purpose: This function computes E(x) given that the integral terms have
  already been calculated and stored in integration_2233b1.

  ------------------------------------------------------------------------------

  Parameters:

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  psiR -- This is the surface twist of the fibrils, psi(R). In this code, this
  term would be y[1][mpt].

  integration_2233b1 -- This is the value of the integral that is computed
  elsewhere.

  ------------------------------------------------------------------------------

  Returns: Returns E(x) if the integral term integration_2233b1 is the true
  integral value.

  ============================================================================*/

{

  double E;

  // first calculate bulk energy per unit length
  E = 2.0/(p->R*p->R)*integration_2233b1; 

  // add density fluctuations term
  E += (p->delta*p->delta*p->omega*0.5
	*(0.75*p->delta*p->delta-1)*(1-p->R_c*p->R_c/(p->R*p->R)));

  // add surface term tension terms
  E += +0.5+1.0/p->R*(-(1+p->k24)*(sin(psiR)*sin(psiR))/p->R+2.0*p->gamma_s);  
  

  return E;
}



double rf_fibril(double x,struct params *p)
{

  void psi_interp(double x,double *r,double **y,double *psi,double *psiprime);

  if (x < CPTR_ERR) return 0;

  double ans;

  double psi,psi_p;

  psi_interp(x,p->r,p->y,&psi,&psi_p);

  double sin_psi = sin(psi);
  double sin_2psi = sin(2*psi);
  double cos_psi = cos(psi);

  

  ans = (-(psi_p+0.5*sin_2psi/x)+0.5*(psi_p+0.5*sin_2psi/x)
	 *(psi_p+0.5*sin_2psi/x)+0.5*p->K33*sin_psi*sin_psi*sin_psi
	 *sin_psi/(x*x));

  if (x >= p->R_c) {

    ans += (p->Lambda*p->delta*p->delta/4.0
	    *(4*M_PI*M_PI-p->eta*p->eta*cos_psi*cos_psi)
	    *(4*M_PI*M_PI-p->eta*p->eta*cos_psi*cos_psi));

  }

  return ans*x;
  
}


void add_R_c(struct params *p)
{

  int get_i_c(double *r,double R_c);

  int i_c = get_i_c(p->r,p->R_c);

  p->mpt = p->M0;

  if (fabs(p->r[i_c]-p->R_c) > CPTR_ERR) {

    p->mpt += 1;
    
    for (int i = p->mpt; i > i_c; i--) p->r[i] = p->r[i-1];

    p->r[i_c] = p->R_c;

    p->y[1][p->mpt] = p->y[1][p->mpt-1];
    p->y[2][p->mpt] = p->y[2][p->mpt-1];

  } else { printf(" close to r[%d] when R_c = %e\n",i_c,p->R_c);}

  return;
  
}


int get_i_c(double *r,double R_c)
{

  int i = 1;

  while (r[i] < R_c) i += 1;

  return i;
  
}

