#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include <stdbool.h>

#define STARTX   0.0
#define ENDX     M_PI



int main()
{

  void set_x_y(double *x, double *y, int m, int N);
  double cubic_spline(double t,double *x,double *y, double *z);
  void find_zs(double *x, double *y, double *z,int m,int N);
  double qromb(double (*func)(double,double *,double *,double *),double t0,double t1,
	       double *x, double *y, double *z,double tol,bool *failure);


  void savefile(char name[],double *x, double *y,int N);
  
  int N = 100;

  int m = N/4;

  int tlength = 2*N;

  double *x = malloc(sizeof(double)*(N+1));

  double *y = malloc(sizeof(double)*(N+1));

  double *z = malloc(sizeof(double)*(N+1));

  double *t = malloc(sizeof(double)*tlength);

  double *f = malloc(sizeof(double)*tlength);
  
  set_x_y(x+1,y+1,m,N);

  find_zs(x+1,y+1,z,m,N);

  bool failure = true;

  printf("%d\n",failure);

  double integral = qromb(cubic_spline,x[m+1],x[N],x+1,y+1,z,1e-15,&failure);

  printf("%lf\n",integral);
  printf("%d\n",failure);

  printf("%lf\n",sin(x[N])-sin(x[m+1]));
  
  savefile("discrete.txt",x+1,y+1,N);

  double dt = (x[N]-x[m+1])/(tlength-1);

  for (int i = 0; i < tlength; i++) {

    t[i] = x[m+1]+dt*i;

    f[i] = cos(t[i]);

  }

  savefile("spline.txt",t,f,tlength);

  
  free(x);
  free(y);
  free(z);
  free(t);
  free(f);
  return 0;

}


void savefile(char name[],double *x, double *y,int N)
{

  FILE *out = fopen(name,"w");

  for (int i = 0; i < N; i++) {

    fprintf(out,"%lf\t%lf\n",x[i],y[i]);
  }

  fclose(out);
  
  return;

}

void set_x_y(double *x, double *y, int m, int N)
{

  double x0 = STARTX;
  double xf = ENDX;

  double dx = (xf-x0)/(N-1);

  
  for (int i = 0; i < m; i++) {

    x[i] = x0+i*dx;
    
    y[i] = 0.0;

  }

  for (int i = m; i < N; i++) {

    x[i] = x0+i*dx;
    
    y[i] = cos(x[i]);

  }

  return;

}

double cubic_spline(double t,double *x,double *y, double *z)
{
  double s_of_t(double t, double *x,double *y, double *z, int i);

  int i = 0;

  while (x[i] < t) i++;

  if (fabs(t-x[i])<1e-6) return y[i];
  
  return s_of_t(t,x,y,z,i);
  
}

double s_of_t(double t, double *x, double *y,double *z, int i)
{

  double h_i = x[i]-x[i-1];
  
  double ans = z[i]*(t-x[i-1])*(t-x[i-1])*(t-x[i-1]);

  ans += z[i-1]*(x[i]-t)*(x[i]-t)*(x[i]-t);

  ans /= (6.0*h_i);

  ans += (y[i]/h_i-z[i]*h_i/6.0)*(t-x[i-1]);

  ans += (y[i-1]/h_i-h_i*z[i-1]/6.0)*(x[i]-t);

  return ans;

}

void find_zs(double *x, double *y, double *z,int m,int N)
{

  void set_matrix_eqn(double *a,double *D,double *b, double *x, double *y,
		      int N);

  void tridiag(double *a, double *D, double *b,double *ans, int N);

  double *a = malloc(sizeof(double)*(N-m-2));

  double *D = malloc(sizeof(double)*(N-m-1));

  double *b = malloc(sizeof(double)*(N-m-1));

  set_matrix_eqn(a,D,b,x+m,y+m,N-m-1);
  
  z[m] = z[N-1] = 0.0;
  
  tridiag(a,D,b,z+1+m,N-m-2);


  free(a);
  free(D);
  free(b);

  return;
  
}


void set_matrix_eqn(double *a,double *D,double *b, double *x, double *y,
		    int M)
{

  double h_i = x[1]-x[0];
  double b_i = 1/h_i*(y[1]-y[0]);

  double h_im1 = h_i;
  double b_im1 = b_i;

  double v_i,u_i; // initialized below

  for (int i = 1; i < M; i++) {

    h_i = x[i+1]-x[i];

    b_i = 1.0/h_i*(y[i+1]-y[i]);

    v_i = 2*(h_im1+h_i);

    u_i = 6*(b_i-b_im1);


    if (i<M-1) a[i-1] = h_i;

    D[i-1] = v_i;

    b[i-1] = u_i;

    h_im1 = h_i;

    b_im1 = b_i;

  }

  return;
  
}

void tridiag(double *a, double *D, double *b,double *ans, int N)
{

  a[0] = a[0]/D[0];
  b[0] = b[0]/D[0];

  for (int i = 1; i < N-1; i++) {

    D[i] = D[i]-D[i-1]*a[i-1]*a[i-1];
    a[i] = a[i]/D[i];
    b[i] = (b[i]-b[i-1]*D[i-1]*a[i-1])/D[i];
  }
  
  D[N-1] = D[N-1]-D[N-2]*a[N-2]*a[N-2];
  b[N-1] = (b[N-1]-b[N-2]*D[N-2]*a[N-2])/D[N-1];

  ans[N-1] = b[N-1];

  
  for (int i = N-1; i > 0; i--) {

    ans[i-1] = b[i-1]-a[i-1]*ans[i];

  }

  return;

}
