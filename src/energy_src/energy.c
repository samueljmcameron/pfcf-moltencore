/* Functions to calculate the cost function, which is the energy per unit
   volume E(x), for a given value of x. */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include "nrutil.h"
#include "../headerfile.h"

#define LRG_NBR 1e10



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

  void solvde_wrapper(double scalv[],struct params *p,double h,
		      bool ignore_first_y);
  
  void propagate_r(double *r, double h,int mpt);
  
  int find_i0(double R0, double *r);
  
  bool successful_E_count(double *E,struct params *p);

  double h;
  double E;  
  double scalv[2+1];



  scalv[1] = .1;    // guess for magnitude of the psi values
  scalv[2] = 4.0;   // guess for magnitude of the psi' values



  if (p->R <= 0 || p->eta <= 0 || p->eta >= 8.0
      || fabs(p->delta) >= 1.0 || p->R0 >= p->R || p->R0 < -1e-8) {
    printf("p->R = %e, p->R0 = %e, p->eta = %e, p->delta = %e\n",p->R,p->R0,p->eta,p->delta);
    printf("something is too big or less than zero, so returning failed calculation.\n");
    return FAILED_E;
  }
  

  h = p->R/(p->mpt-1);    // compute stepsize in r[1..mpt] 

  propagate_r(p->r,h,p->mpt);

  p->i0 = find_i0(p->R0,p->r);
  
  solvde_wrapper(scalv,p,h,false);

  if(successful_E_count(&E,p)) return E;

  return FAILED_E;
}


void propagate_r(double *r, double h,int mpt)
{
  // only change r since psi, psi' are stored from last loop
  int k;
  for (k=1;k <=mpt; k++) r[k] = (k-1)*h; 
  return;
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

  void compute_rf2233(struct params *p);
  void compute_rfb(struct params *p);
  double E_bulk_surface(struct params *p,double psiR,
			double integration_2233b1);

  double cubic_spline(double t,double *x,double *y, double *z);
  void find_zs(double *x, double *y, double *z,int m,int N);
  double qromb(double (*func)(double,double *,double *,double *),double t0,double t1,
	       double *x, double *y,double *z,bool *failure);

  
  double integration_2233;
  double integration_b;

  bool failure;
  
  compute_rf2233(p);

  find_zs(p->r+1,p->rf_fib+1,p->z+1,0,p->mpt);
  
  integration_2233 = qromb(cubic_spline,0,p->R,p->r+1,p->rf_fib+1,p->z+1,&failure);
  
  if (failure) {
    *E = integration_2233;
    printf("failed to integrate frank free energy at (R,eta,delta) = (%e,%e,%e)\n",
	   p->R,p->eta,p->delta);
    return false;
  }

  compute_rfb(p);
  
  find_zs(p->r+1,p->rf_fib+1,p->z+1,p->i0-1,p->mpt);
  
  integration_b = qromb(cubic_spline,p->R0,p->R,p->r+1,p->rf_fib+1,p->z+1,&failure);

  if (failure) {
    *E = integration_b;
    printf("failed to integrate pf free energy at (R,eta,delta) = (%e,%e,%e)\n",
	   p->R,p->eta,p->delta);
    return false;
  }
  
  *E = E_bulk_surface(p,p->y[1][p->mpt],integration_2233+integration_b);


  return true;
}

int find_i0(double R0, double *r)
{

  int i0 = 1;

  while (r[i0]<=R0) i0 += 1;

  i0 -= 1;
  
  return (i0 < 1) ? 1 : i0;

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
	*(0.75*p->delta*p->delta-1)*(1-p->R0*p->R0/(p->R*p->R)));

  // add surface term tension terms
  E += +0.5+1.0/p->R*(-(1+p->k24)*(sin(psiR)*sin(psiR))/p->R+2.0*p->gamma_s);  
  

  return E;
}

void compute_rfb(struct params *p)
/*==============================================================================

  Purpose: This function computes r*f_fibril(r) for all ri in r[1..mpt], and
  stores it into the array rf_fib[1..mpt].

  ------------------------------------------------------------------------------

  Parameters:

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  ============================================================================*/
{
  double fb_r(struct params *p,double ri,double cos_yi);

  double cosy;
  
  for (int i = p->i0; i <= p->mpt; i++) {  // compute f_fibril*r

    cosy = cos(p->y[1][i]);

    p->rf_fib[i] = p->r[i]*fb_r(p,p->r[i],cosy);

  }

  return;
}



void compute_rf2233(struct params *p)
/*==============================================================================

  Purpose: This function computes r*f_fibril(r) for all ri in r[1..mpt], and
  stores it into the array rf_fib[1..mpt].

  ------------------------------------------------------------------------------

  Parameters:

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  ============================================================================*/
{

  double f2233_r(struct params *p,double ri,double sin_yi,
		   double sin_2yi,double cos_yi,double yi_p);



  double siny, sin2y,cosy;
  
  p->rf_fib[1] = 0;
  
  for (int i = 2; i <= p->mpt; i++) {  // compute f_fibril*r

    siny = sin(p->y[1][i]);

    sin2y = sin(2*p->y[1][i]);

    cosy = cos(p->y[1][i]);

    p->rf_fib[i] = p->r[i]*f2233_r(p,p->r[i],siny,sin2y,cosy,
				   p->y[2][i]);

  }

  return;
}


double fb_r(struct params *p,double ri,double cos_yi)
/*==============================================================================

  Purpose: This function computes the integrand of the energy functional (which
  is the energy density r*f_fibril(r) in the model) at the point ri.

  ------------------------------------------------------------------------------

  Parameters:

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  ri -- This is the current grid point where the function is being calculated
  at.

  sin_yi, sin_2yi, cos_yi -- These are just shortcut variables for calculating
  sin(psi(ri)), sin(2*psi(ri)), and cos(psi(ri)) (used to save time instead of
  evaluating trig functions multiple times).

  yi_p -- This is the value dpsi/dr evaluated at ri.

  ------------------------------------------------------------------------------

  Returns: Returns r*f_fibril(r) at the point ri.

  ============================================================================*/
{

  double ans;

  ans = (p->Lambda*p->delta*p->delta/4.0
	 *(4*M_PI*M_PI-p->eta*p->eta*cos_yi*cos_yi)
	 *(4*M_PI*M_PI-p->eta*p->eta*cos_yi*cos_yi));

  return ans;

}

double f2233_r(struct params *p,double ri,double sin_yi,
	       double sin_2yi,double cos_yi,double yi_p)
/*==============================================================================

  Purpose: This function computes the integrand of the energy functional (which
  is the energy density r*f_fibril(r) in the model) at the point ri.

  ------------------------------------------------------------------------------

  Parameters:

  p -- This struct has all of the constant parameter info (e.g. K33, k24).

  ri -- This is the current grid point where the function is being calculated
  at.

  sin_yi, sin_2yi, cos_yi -- These are just shortcut variables for calculating
  sin(psi(ri)), sin(2*psi(ri)), and cos(psi(ri)) (used to save time instead of
  evaluating trig functions multiple times).

  yi_p -- This is the value dpsi/dr evaluated at ri.

  ------------------------------------------------------------------------------

  Returns: Returns r*f_fibril(r) at the point ri.

  ============================================================================*/
{

  double ans;

  ans = (-(yi_p+0.5*sin_2yi/ri)+0.5*(yi_p+0.5*sin_2yi/ri)
	 *(yi_p+0.5*sin_2yi/ri)+0.5*p->K33*sin_yi*sin_yi*sin_yi
	 *sin_yi/(ri*ri));
  
  return ans;

}
