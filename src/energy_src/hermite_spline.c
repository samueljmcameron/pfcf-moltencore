#include <stdio.h>
#include <stdlib.h>
#include "../headerfile.h"

#define CMPT_ERR 1e-15

#define NO_ADD_I -1

void psi_interp(double x,double *r, double **y,double *psi,double *psiprime)
{

  double hermite_poly(double x, double x_1,double x_2,double y_1,double y_2,
		      double yp_1,double yp_2);

  double hermite_deriv(double x, double x_1,double x_2,double y_1,double y_2,
		       double yp_1,double yp_2);

  int i = 1;

  while (r[i]<x) i += 1;

  *psi =  hermite_poly(x,r[i],r[i+1],y[1][i],y[1][i+1],
		       y[2][i],y[2][i+1]);

  *psiprime =  hermite_deriv(x,r[i],r[i+1],y[1][i],y[1][i+1],
			     y[2][i],y[2][i+1]);

  return;
  
}

double hermite_poly(double x, double x_1,double x_2,double y_1,double y_2,
		    double yp_1,double yp_2)
{

  double t = (x-x_1)/(x_2-x_1);

  if (fabs(t)<CMPT_ERR) printf("t is zero\n");

  double h00 = 2*t*t*t-3*t*t+1;

  double h10 = t*t*t-2*t*t+t;

  double h01 = -2*t*t*t+3*t*t;

  double h11 = t*t*t-t*t;

  return h00*y_1 + h10*(x_2-x_1)*yp_1 + h01*y_2 + h11*(x_2-x_1)*yp_2;

}


double hermite_deriv(double x, double x_1,double x_2,double y_1,double y_2,
		     double yp_1,double yp_2)
{

  double t = (x-x_1)/(x_2-x_1);

  double h00p = 6*t*t-6*t;

  double h10p = 3*t*t-4*t+1;

  double h01p = -6*t*t+6*t;

  double h11p = 3*t*t-2*t;

  return h00p*y_1/(x_2-x_1) + h10p*yp_1+h01p*y_2/(x_2-x_1) + h11p*yp_2;

}
