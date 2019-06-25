#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../headerfile.h"

#define CPTR_ERR 1e-15


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

  if (p->R_c < x) {
    ans += (p->Lambda*p->delta*p->delta/4.0
	    *(4*M_PI*M_PI-p->eta*p->eta*cos_psi*cos_psi)
	    *(4*M_PI*M_PI-p->eta*p->eta*cos_psi*cos_psi));

  }

  return ans*x;
  
}
