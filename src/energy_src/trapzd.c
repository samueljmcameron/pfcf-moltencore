#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "../headerfile.h"

#define MACHINE_ERR 1e-14

double trapzd(double (*func)(double,struct params *),double t0,
	      double t1, int n,struct params *p, bool *failure)
/*This routine computes the nth stage of refinement of an extended trapezoidal rule. func is input as a pointer to the function to be integrated between limits a and b, also input. When called with n=1, the routine returns the crudest estimate of integral_a^b(f(x)dx). Subsequent calls with n=2,3,... (in that sequential order) will improve the accuracy by adding 2n-2 additional interior points. */
{
  double t,tnm,sum,del;
  static double s;
  int it,j,index,spacing;

  if (n == 1) {
    return (s=0.5*(t1-t0)*((*func)(t0,p)+(*func)(t1,p)));
  } else {
    
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    
    tnm=it; // equals 1 if n = 2
    
    del=(t1-t0)/tnm; //This is the spacing of the points to be added.

    if (del < 2*MACHINE_ERR) {
      printf("Stepsize of integration smaller than discrete data spacing.\n");
      *failure = true;
      return s;
    }
    
    t = t0 + 0.5*del;
    
    for (sum=0.0,j=1;j<=it;j++,t += del) sum += (*func)(t,p);
    
    s=0.5*(s+(t1-t0)*sum/tnm); //This replaces s by its refined value.

    *failure = false;
    return s;
  }
}
