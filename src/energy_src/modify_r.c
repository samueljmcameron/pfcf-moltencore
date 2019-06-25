#include <stdio.h>
#include <stdlib.h>
#include "../headerfile.h"

#define CPTR_ERR 1e-15
#define NO_ADD_I -1


void add_R_c(double *r,double R_c, int *N)
{

  int get_i_c(double *r,double R_c);

  int i_c = get_i_c(r,R_c);

  if (fabs(r[i]-R_c) > CPTR_ERR) {

    *N += 1;
    
    for (int i = *N; i > i_c; i--) r[i] = r[i-1];

    r[i_c] = R_c;

  }

  return;
  
}


int get_i_c(double *r,double R_c)
{

  int i = 1;

  while (r[i] < R_c) i += 1;

  return i;
  
}

