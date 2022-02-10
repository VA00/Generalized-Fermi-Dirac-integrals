/*
A. Odrzywolek, andrzej.odrzywolek@uj.edu.pl
*/
#include "../fermidirac.h"

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <stdio.h>
#include <float.h>
#include <quadmath.h>

//without leading g !
double g_derivative(double x, double theta, int n)
{
  double g2 = 1.0+0.5*x*theta;
  double g  = sqrt(g2);
  double fac = 0.25*x/g2;

  if(n==0) return 1.0;
  //if(n==0) return g; 

  return -(2.0*n-3.0)*fac*g_derivative(x, theta, n-1);

}

long double g_derivative_long(long double x, long double theta, int n)
{
  long double g2 = 1.0L+0.5L*x*theta;
  long double g  = sqrtl(g2);
  long double fac = 0.25L*x/g2;

  if(n==0) return 1.0L;
  //if(n==0) return g; 

  return -(2.0L*n-3.0L)*fac*g_derivative_long(x, theta, n-1);

}

__float128 g_derivative_quad(__float128 x, __float128 theta, int n)
{
  __float128 g2 = 1.0q+0.5q*x*theta;
  __float128 g  = sqrtq(g2);
  __float128 fac = 0.25q*x/g2;

  if(n==0) return 1.0q;
  //if(n==0) return g; 

  return -(2.0q*n-3.0q)*fac*g_derivative_quad(x, theta, n-1);

}

//without leading g !
void g_derivative_vector(double x, double theta, double dg[DERIVATIVE_MATRIX_SIZE])
{
  double g2 = 1.0+0.5*x*theta;
  double fac = 0.25*x/g2;
  int i;
  
  for(i=0;i<DERIVATIVE_MATRIX_SIZE;i++)
   {
    dg[i] = ((i%2) ? 1.0 : -1.0)*fac2(2*i-3)*power_squaring(fac,i);
   }

  return;

}

/* Linear recurrence required to compute general mixed partial derivatives D[\[Eta]^k Sqrt[1 + \[Eta] \[Theta]/2], {\[Theta], n}, {\[Eta], m}]*/
double r(int i, double k, double z, int n, int m)
{
  double c1,c2;

  if(i==0) return 0.0;
  if(i==1) return pow(1.0+z,0.5-n)*binom(k+n,m);
     
  

  //c1 = -(((-2. + i - m)*(-5. + 2.*i + 2.*n)*z)/(2.*(-1. + i)*(-1. + i + k - m + n)*(1. + z))) ;
  
  //c2 = ((2. + 2.*m - 2.*n + 12.*z + 7.*m*z - 6.*n*z - 2.*m*n*z - 2.*k*(1. + z) + i*i*(2. + 4.*z) + 
    //i*(-4. - 2.*m + 2.*n - 13.*z - 4.*m*z + 4.*n*z + 2.*k*(1. + z))))/(2.*(-1. + i)*(-1. + i + k - m + n)*(1. + z));


   c1 = (2 + 2*m - 2*n + 12*z + 7*m*z - 6*n*z - 2*m*n*z - 2*k*(1 + z) + i*i*(2 + 4*z) + i*(-4 - 2*m + 2*n - 13*z - 4*m*z + 4*n*z + 2*k*(1 + z)))/
   (2.*(-1 + i)*(-1 + i + k - m + n)*(1 + z)); 

   if(i==2) return c1*r(1,k, z, n, m); 

   c2 = -0.5*((-2 + i - m)*(-5 + 2*i + 2*n)*z)/((-1 + i)*(-1 + i + k - m + n)*(1 + z));
  

   printf("DBG1 c1 c2 :\ti=%d\tn=%d\tm=%d\t%lf %lf\n", i, n, m,  c1,c2);

  return c2*r(i-2, k, z, n, m) + c1*r(i-1, k, z, n, m);
              
   
  
}

/* OBSOLETE: new explicit formula for derivatives found
/* I give up for now. (-1 + i + k - m + n) cause division by zero if k is an integer. It cancel,
e.g  for i=3, n=0, m=2 we get k(k-1)(k-2) from factorial_power divided by k,
k(k-1)(k-2)/k = (k-1)(k-2) */
double r2(int i, double k, double z, int n, int m)
{
  double c1,c2;

  if(i==0) return 0.0;
  if(i==1) return pow(1.0+z,0.5-n)*factorial_power(k+n,m);
     
  

  //c1 = -(((-2. + i - m)*(-5. + 2.*i + 2.*n)*z)/(2.*(-1. + i)*(-1. + i + k - m + n)*(1. + z))) ;
  
  //c2 = ((2. + 2.*m - 2.*n + 12.*z + 7.*m*z - 6.*n*z - 2.*m*n*z - 2.*k*(1. + z) + i*i*(2. + 4.*z) + 
    //i*(-4. - 2.*m + 2.*n - 13.*z - 4.*m*z + 4.*n*z + 2.*k*(1. + z))))/(2.*(-1. + i)*(-1. + i + k - m + n)*(1. + z));


   c1 = (2 + 2*m - 2*n + 12*z + 7*m*z - 6*n*z - 2*m*n*z - 2*k*(1 + z) + i*i*(2 + 4*z) + i*(-4 - 2*m + 2*n - 13*z - 4*m*z + 4*n*z + 2*k*(1 + z)))/
   (2.*(-1 + i)*(-1 + i + k - m + n)*(1 + z)); 

   if(i==2) return c1*r(1,k, z, n, m); 

   c2 = -0.5*((-2 + i - m)*(-5 + 2*i + 2*n)*z)/((-1 + i)*(-1 + i + k - m + n)*(1 + z));
  

   printf("DBG2 c1 c2 :\ti=%d\tn=%d\tm=%d\t%lf %lf\n", i, n, m,  c1,c2);

  return c2*r(i-2, k, z, n, m) + c1*r(i-1, k, z, n, m);
              
   
  
}


/* Functions below are integrated with so-called DoubleExponential or Tanh-Sinh quadrature.
 * 
 * Some references:
 * 
 * Mori, Masatake (2005), "Discovery of the double exponential transformation and its developments", 
 * Publications of the Research Institute for Mathematical Sciences 41 (4): 897–935, 
 * doi:10.2977/prims/1145474600, 
 * ISSN 0034-5318
 * http://www.kurims.kyoto-u.ac.jp/~okamoto/paper/Publ_RIMS_DE/41-4-38.pdf, eq. (4.17)
 * 
 * See also: http://en.wikipedia.org/wiki/Tanh-sinh_quadrature and references therein.
 * 
 */



/*

		SECTION FOR RELATIVISTIC Fermi-Dirac integrals (F-function)
        INCLUDING DERIVATIVES UP TO THIRD ORDER
 
        UNFORTUNATELY, vector length is 10, neither AVX or AVX512 aligned. 
        Are all of the 3-rd order derivatives required?

*/



void integrandF_derivatives(const double t, const double k, const double eta, const double theta,        
       double * integrand, double *integrand_deta, double *integrand_deta2,
       double * integrand_dtheta, double *integrand_dtheta2, double *integrand_deta_dtheta,
       double * integrand_dtheta3, double *integrand_dtheta2_deta, double *integrand_dtheta_deta2, double *integrand_deta3)

{
  
  double x,dx,exp_t,s,g,g2,z,f;
  //double s1,s2,s3;

  exp_t  = exp(-t); //this might be faster, THX Karol U.
  x      = exp(  t - exp_t ); /* Masatake Mori, eq. (4.17) */
  dx     = 1.0+exp_t; /* in this case x is adsorbed in integrand, and x^k -> x^(k+1) */
  g2  = 1.0+ 0.5*theta*x;
  g = sqrt(g2);
  z = 0.25*x/g2;
  s = sigmoid(eta-x); // if using machine precison sigmoid is equal 1.0 implementation should handle this  
    /* possible optimization for eta derivatives
    s1 = s*(1.0-s);    s2 = s1*(1.0-2.0*s);    s3 = s1*(1.0-6.0*s1);
    */
  
  if(x-eta<-log(DBL_EPSILON)) // if using machine precison we are unable to add 1.0 to exp(), then approximation is optimal
    {
	
	f = exp( (k+1.0)*(t - exp_t) );
    f = f*g*s*dx;
	
    *integrand                  =   f;
	*integrand_deta             =   f*(1.0-s); 
	*integrand_deta2            =   f*(1.0-s)*(1.0-2.0*s); 
	*integrand_dtheta           =   f*z;
	*integrand_dtheta2          =  -f*z*z;
	*integrand_deta_dtheta      =   f*z*(1.0-s);
	*integrand_dtheta3          =   f*3.0*z*z*z;
	*integrand_dtheta2_deta     =  -f*z*z*(1.0-s);
	*integrand_dtheta_deta2     =   f*z*(1.0+s*(2.0*s-3.0));
	*integrand_deta3            =   f*(1.0 + s*(-7.0 + (12.0 - 6.0*s)*s));
	
	}
  else
    {
    //sigma = exp(eta-x) sigmoid adsorbed into exp, to avoid 0*infinity mess 
	f = exp((k+1.0)*(t - exp_t) + eta - x );
	f = f*g*dx;

    *integrand                  =   f;
	*integrand_deta             =   f*(1.0-s); 
	*integrand_deta2            =   f*(1.0-s)*(1.0-2.0*s); 
	*integrand_dtheta           =   f*z;
	*integrand_dtheta2          =  -f*z*z;
	*integrand_deta_dtheta      =   f*z*(1.0-s);
	*integrand_dtheta3          =   f*3.0*z*z*z;
	*integrand_dtheta2_deta     =  -f*z*z*(1.0-s);
	*integrand_dtheta_deta2     =   f*z*(1.0+s*(2.0*s-3.0));
	*integrand_deta3            =   f*(1.0 + s*(-7.0 + (12.0 - 6.0*s)*s));
    }


  
  
}

void integrandF_derivatives_v2(const double t, const double k, const double eta, const double theta, double integrand[10])
{
  
  double x,dx,exp_t,s,g,g2,z,f;
  //double s1,s2,s3;

  exp_t  = exp(-t); //this might be faster, THX Karol U.
  x      = exp(  t - exp_t ); /* Masatake Mori, eq. (4.17) */
  dx     = 1.0+exp_t; /* in this case x is adsorbed in integrand, and x^k -> x^(k+1) */
  g2  = 1.0+ 0.5*theta*x;
  g = sqrt(g2);
  z = 0.25*x/g2;
  s = sigmoid(eta-x); // if using machine precison sigmoid is equal 1.0 implementation should handle this  
    /* possible optimization for eta derivatives
    s1 = s*(1.0-s);    s2 = s1*(1.0-2.0*s);    s3 = s1*(1.0-6.0*s1);
    */
  
  if(x-eta<-log(DBL_EPSILON)) // if using machine precison we are able to add 1.0 to exp() in sigmoid
    {
	
	f = exp( (k+1.0)*(t - exp_t) );
    f = f*g*s*dx;
	
    integrand[0]  =   f;
	integrand[1]  =   f*(1.0-s); 
	integrand[2]  =   f*(1.0-s)*(1.0-2.0*s); 
	integrand[3]  =   f*z;
	integrand[4]  =  -f*z*z;
	integrand[5]  =   f*z*(1.0-s);
	integrand[6]  =   f*3.0*z*z*z;
	integrand[7]  =  -f*z*z*(1.0-s);
	integrand[8]  =   f*z*(1.0+s*(2.0*s-3.0));
	integrand[9]  =   f*(1.0 + s*(-7.0 + (12.0 - 6.0*s)*s));
	
	}
  else // if using machine precison we are UNABLE to add 1.0 to exp() in sigmoid
    {
    //sigma = exp(eta-x) sigmoid adsorbed into exp, to avoid 0*infinity mess 
	f = exp((k+1.0)*(t - exp_t) + eta - x );
	f = f*g*dx;

    integrand[0]  =   f;
	integrand[1]  =   f*(1.0-s); 
	integrand[2]  =   f*(1.0-s)*(1.0-2.0*s); 
	integrand[3]  =   f*z;
	integrand[4]  =  -f*z*z;
	integrand[5]  =   f*z*(1.0-s);
	integrand[6]  =   f*3.0*z*z*z;
	integrand[7]  =  -f*z*z*(1.0-s);
	integrand[8]  =   f*z*(1.0+s*(2.0*s-3.0));
	integrand[9]  =   f*(1.0 + s*(-7.0 + (12.0 - 6.0*s)*s));
    }


  
  
}


void integrandF_derivatives_v3(const double t, const double k, const double eta, const double theta, double integrand[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE])
{
    //double result[10];
    double ds[DERIVATIVE_MATRIX_SIZE],dg[DERIVATIVE_MATRIX_SIZE];
    double x,dx,exp_t,s,g,g2,z,f;
    int m,n;

  exp_t  = exp(-t); //this might be faster, THX Karol U.
  x      = exp(  t - exp_t ); /* Masatake Mori, eq. (4.17) */
  dx     = 1.0+exp_t; /* in this case x is adsorbed in integrand, and x^k -> x^(k+1) */
  g2  = 1.0+ 0.5*theta*x;
  g = sqrt(g2);
  z = 0.25*x/g2;
  s = sigmoid(eta-x);
  
  if(x-eta<-log(DBL_EPSILON)) // if using machine precison we are able to add 1.0 to exp() in sigmoid
    {
	
	f = exp( (k+1.0)*(t - exp_t) );
    f = f*g*s*dx;
	
	}
  else // if using machine precison we are UNABLE to add 1.0 to exp() in sigmoid
    {
    //sigma = exp(eta-x) sigmoid adsorbed into exp, to avoid 0*infinity mess 
	f = exp((k+1.0)*(t - exp_t) + eta - x );
	f = f*g*dx;
    }

    //calling eta and theta derivatives
    sigmoid_derivative_polynomial_vector(s, ds);
    g_derivative_vector(x, theta, dg);
    
    for(m=0;m<DERIVATIVE_MATRIX_SIZE;m++)
      for(n=0;n<DERIVATIVE_MATRIX_SIZE;n++)
       {
         if(n+m>DERIVATIVE_MAX_ORDER) continue;
         integrand[m][n] = f*ds[m]*dg[n];
       }
     
/*
    integrandF_derivatives_v2(t, k, eta, theta, result);

    integrand[0][0] = result[0];
    integrand[1][0] = result[1];
    integrand[2][0] = result[2];
    integrand[3][0] = result[9];
    integrand[0][1] = result[3];
    integrand[1][1] = result[5];
    integrand[2][1] = result[8];
    integrand[3][1] = -1.0;
    integrand[0][2] = result[4];
    integrand[1][2] = result[7];
    integrand[2][2] = -1.0;
    integrand[3][2] = -1.0;
    integrand[0][3] = result[6];
    integrand[1][3] = -1.0;
    integrand[2][3] = -1.0;
    integrand[3][3] = -1.0;
*/
  return;
}


double integrandF_derivatives_m_n(const double t, const double k, const double eta, const double theta, const int m, const int n)
{
    
    double ds,dg;
    double x,dx,exp_t,s,g,g2,z,f, integrand;
    

  exp_t  = exp(-t); //this might be faster, THX Karol U.
  x      = exp(  t - exp_t ); /* Masatake Mori, eq. (4.17) */
  dx     = 1.0+exp_t; /* in this case x is adsorbed in integrand, and x^k -> x^(k+1) */
  g2  = 1.0+ 0.5*theta*x;
  g = sqrt(g2);
  z = 0.25*x/g2;
  s = sigmoid(eta-x);
  
  if(x-eta<-log(DBL_EPSILON)) // if using machine precison we are able to add 1.0 to exp() in sigmoid
    {
	
	f = exp( (k+1.0)*(t - exp_t) );
    f = f*g*s*dx;
	
	}
  else // if using machine precison we are UNABLE to add 1.0 to exp() in sigmoid
    {
    //sigma = exp(eta-x) sigmoid adsorbed into exp, to avoid 0*infinity mess 
	f = exp((k+1.0)*(t - exp_t) + eta - x );
	f = f*g*dx;
    }


    ds=sigmoid_derivative_polynomial(s, m);
    dg=g_derivative(x, theta, n);
    
    integrand = f*ds*dg;

  return integrand;
}


long double integrandF_derivatives_m_n_long(const long double t, const long double k, const long double eta, const long double theta, const int m, const int n)
{
    
    long double ds,dg;
    long double x,dx,exp_t,s,g,g2,z,f, integrand;
    

  exp_t  = expl(-t); //this might be faster, THX Karol U.
  x      = expl(  t - exp_t ); /* Masatake Mori, eq. (4.17) */
  dx     = 1.0L+exp_t; /* in this case x is adsorbed in integrand, and x^k -> x^(k+1) */
  g2  = 1.0L + 0.5L*theta*x;
  g = sqrtl(g2);
  z = 0.25L*x/g2;
  s = sigmoid_long(eta-x);
  
  if(x-eta<-logl(LDBL_EPSILON)) // if using machine precison we are able to add 1.0 to exp() in sigmoid
    {
	
	f = expl( (k+1.0L)*(t - exp_t) );
    f = f*g*s*dx;
	
	}
  else // if using machine precison we are UNABLE to add 1.0 to exp() in sigmoid
    {
    //sigma = exp(eta-x) sigmoid adsorbed into exp, to avoid 0*infinity mess 
	f = expl((k+1.0L)*(t - exp_t) + eta - x );
	f = f*g*dx;
    }


    ds=sigmoid_derivative_polynomial_long(s, m);
    dg=g_derivative_long(x, theta, n);
    
    integrand = f*ds*dg;

  return integrand;
}

__float128 integrandF_derivatives_m_n_quad(const __float128 t, const __float128 k, const __float128 eta, const __float128 theta, const int m, const int n)
{
    
    __float128 ds,dg;
    __float128 x,dx,exp_t,s,g,g2,z,f, integrand;
    

  exp_t  = expq(-t); //this might be faster, THX Karol U.
  x      = expq(  t - exp_t ); /* Masatake Mori, eq. (4.17) */
  dx     = 1.0q+exp_t; /* in this case x is adsorbed in integrand, and x^k -> x^(k+1) */
  g2  = 1.0q + 0.5q*theta*x;
  g = sqrtq(g2);
  z = 0.25q*x/g2;
  s = sigmoid_quad(eta-x);
  
  if(x-eta<-logq(FLT128_EPSILON)) // if using machine precison we are able to add 1.0 to exp() in sigmoid
    {
	
	f = expq( (k+1.0q)*(t - exp_t) );
    f = f*g*s*dx;
	
	}
  else // if using machine precison we are UNABLE to add 1.0 to exp() in sigmoid
    {
    //sigma = exp(eta-x) sigmoid adsorbed into exp, to avoid 0*infinity mess 
	f = expq((k+1.0q)*(t - exp_t) + eta - x );
	f = f*g*dx;
    }


    ds=sigmoid_derivative_polynomial_quad(s, m);
    dg=g_derivative_quad(x, theta, n);
    
    integrand = f*ds*dg;

  return integrand;
}

/* vector 10-component version */

void Ffermi_estimate_derivatives(double h, double last_result[10], double k, double eta, double theta, double new_result[10])
{
  
  int step,i,j;
  double sum_Left_old[10] ={ [ 0 ... 9 ] = 0.0 }, sum_Right_old[10]={ [ 0 ... 9 ] = 0.0 };
  double sum_Left_new[10] ={ [ 0 ... 9 ] = 0.0 }, sum_Right_new[10]={ [ 0 ... 9 ] = 0.0 };
  double old_result[10], integrand[10];
  
  
  if(last_result[0]<0.0) /* Negative value means first iteration*/
  {
    step=1;
    integrandF_derivatives_v2(0.0, k, eta, theta, integrand);
    for(j=0;j<10;j++) old_result[j] = 2.0*h*integrand[j];
  }
  else
  {
    step=2;
    for(j=0;j<10;j++) old_result[j] = last_result[j];//Is this necessary? old_result===last_result?
  }
  
  /* integral for 0 < t < Infinity  */
  
  //sum_Right_old = 0.0;
  //sum_Right_new = 0.0;
  
  
  i=1;

  do
  {
	for(j=0;j<10;j++) sum_Right_old[j] = sum_Right_new[j];
    integrandF_derivatives_v2(h*i, k, eta, theta,integrand);
    for(j=0;j<10;j++) sum_Right_new[j] = sum_Right_old[j] + integrand[j];
    
	i = i + step;
  }
  while  ( sum_Right_old[0]<sum_Right_new[0] ); //floating point fixed-point method on first vector component!

  /* integral for -Infinity < t <0  */
  
  //sum_Left_old = 0.0;
  //sum_Left_new = 0.0;
  
  
  i=-1;
  do
  {
	for(j=0;j<10;j++) sum_Left_old[j] = sum_Left_new[j];
    integrandF_derivatives_v2(h*i, k, eta, theta,integrand);
    for(j=0;j<10;j++) sum_Left_new[j] = sum_Left_old[j] + integrand[j];
    i = i - step;
  }
  while  (sum_Left_old[0]<sum_Left_new[0]);
  
  
    for(j=0;j<10;j++) new_result[j] = h*(sum_Left_new[j]  + sum_Right_new[j]) + 0.5*old_result[j];


}



void Ffermi_dblexp_derivatives(const double k, const double eta, const double theta,
  const double precision, const int recursion_limit, double result[10])
{
  
  double old[10]={ [ 0 ... 9 ] = -1.0 }; //Setting old to -1.0 cause Ffermi_estimate_derivatives to restart at the first call
  double new[10]={ [ 0 ... 9 ] =  0.0 };
  double h=0.5;
  int j;
  
  //if(k<=-1.0) return nan("NaN"); /* not converging for k <= -1 */

  Ffermi_estimate_derivatives(h, old, k, eta, theta, new);

  for(j=0;j<10;j++) old[j] = 0.0;/*Is this necessary? Two lines below we reset old to new which is zero anyway. 
         Except precision goal is achieved at first run (in theory possible, if one modify code and set e.g h=0.125 or less what MIGHT may have sense in future
   optimization ) */ 

  
  while( fabs(old[0]-new[0])>precision*fabs(new[0]) && h>pow(2.0,-recursion_limit))
  {
    for(j=0;j<10;j++) old[j]=new[j];
    h=0.5*h;
    Ffermi_estimate_derivatives(h, old, k, eta, theta, new);
  }

  for(j=0;j<10;j++) result[j]=new[j];
    
}

/* Arbitrary order derivative version */

void Ffermi_estimate_derivatives_matrix(double h, double last_result[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE], double k, double eta, double theta, double new_result[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE])
{
  
  int step,i,m,n;
  /* I hope I understand correctly https://gcc.gnu.org/onlinedocs/gcc-4.1.2/gcc/Designated-Inits.html */
  double sum_Left_old[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE] = { [0][0]  = 0.0 }, sum_Right_old[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE] = { [0][0]  = 0.0 };
  double sum_Left_new[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE] = { [0][0]  = 0.0 }, sum_Right_new[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE] = { [0][0]  = 0.0 };
  double old_result[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE], integrand[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE];
  
  
  if(last_result[0][0]<0.0) /* Negative value means first iteration*/
  {
    step=1;
    integrandF_derivatives_v3(0.0, k, eta, theta, integrand);
    for(m=0;m<DERIVATIVE_MATRIX_SIZE;m++) 
     for(n=0;n<DERIVATIVE_MATRIX_SIZE;n++) 
      old_result[m][n] = 2.0*h*integrand[m][n];
  }
  else
  {
    step=2;
    for(m=0;m<DERIVATIVE_MATRIX_SIZE;m++) 
     for(n=0;n<DERIVATIVE_MATRIX_SIZE;n++) 
      old_result[m][n] = last_result[m][n];//Is this necessary? old_result===last_result?
  }
  
  /* integral for 0 < t < Infinity  */
  
  //sum_Right_old = 0.0;
  //sum_Right_new = 0.0;
  
  
  i=1;

  do
  {
    for(m=0;m<DERIVATIVE_MATRIX_SIZE;m++) 
     for(n=0;n<DERIVATIVE_MATRIX_SIZE;n++) 
      sum_Right_old[m][n] = sum_Right_new[m][n];

    integrandF_derivatives_v3(h*i, k, eta, theta,integrand);

    for(m=0;m<DERIVATIVE_MATRIX_SIZE;m++) 
     for(n=0;n<DERIVATIVE_MATRIX_SIZE;n++) 
      sum_Right_new[m][n] = sum_Right_old[m][n] + integrand[m][n];
    
	i = i + step;
  }
  while  ( sum_Right_old[0][0]<sum_Right_new[0][0] ); //floating point fixed-point method on first matrix component!

  /* integral for -Infinity < t <0  */
  
  //sum_Left_old = 0.0;
  //sum_Left_new = 0.0;
  
  
  i=-1;
  do
  {
    for(m=0;m<DERIVATIVE_MATRIX_SIZE;m++) 
     for(n=0;n<DERIVATIVE_MATRIX_SIZE;n++) 
      sum_Left_old[m][n] = sum_Left_new[m][n];

    integrandF_derivatives_v3(h*i, k, eta, theta,integrand);

    for(m=0;m<DERIVATIVE_MATRIX_SIZE;m++) 
     for(n=0;n<DERIVATIVE_MATRIX_SIZE;n++) 
 
      sum_Left_new[m][n] = sum_Left_old[m][n] + integrand[m][n];
    i = i - step;
  }
  while  (sum_Left_old[0][0]<sum_Left_new[0][0]);
  
  
    for(m=0;m<DERIVATIVE_MATRIX_SIZE;m++) 
     for(n=0;n<DERIVATIVE_MATRIX_SIZE;n++) 
       new_result[m][n] = h*(sum_Left_new[m][n]  + sum_Right_new[m][n]) + 0.5*old_result[m][n];


}

/* Arbitrary order derivative version, one derivative per call */

double Ffermi_estimate_derivatives_m_n(double h, double last_result, double k, double eta, double theta, int m, int n)
{
  

  
  double sum_Left_old = 0.0 , sum_Right_old = 0.0 ;
  double sum_Left_new = 0.0 , sum_Right_new = 0.0 ;
  double old_result, new_result, integrand;

  int step,i, i_peak=0;
  double peak_position=0.0;
  if(eta>1024.0) peak_position = log(eta);//+1.0/eta-1.0/eta/eta; Precise formula is Log[eta] + W(1/eta), where W is LambertW==ProductLog
  
  
  if(last_result==0.0) /* ZERO value means first iteration, this should work also for derivatives, which can be negative ! Maybe nan/inf will be more appropriate*/
  {
    step=1;
    integrand = integrandF_derivatives_m_n(peak_position, k, eta, theta, m, n);
    old_result = 2.0*h*integrand;
  }
  else
  {
    step=2;
    old_result = last_result;//Is this necessary? old_result===last_result?
  }
  
  /* integral for 0 < t < Infinity  */
  
  sum_Right_old = 0.0;
  sum_Right_new = 0.0;
  
  //if(eta>64.0) i_peak = (int) ceil(log(eta)/h); 
  
  i=1;
  //i=i_peak+1;

  do
  {
    sum_Right_old = sum_Right_new;

    integrand = integrandF_derivatives_m_n(peak_position+h*i, k, eta, theta, m, n);

      sum_Right_new = sum_Right_old + integrand;
    
	i = i + step;
  }
  while  ( (sum_Right_old!=sum_Right_new) /*|| (h*i<=peak_position)*/ ); 

  /* integral for -Infinity < t <0  */
  
  sum_Left_old = 0.0;
  sum_Left_new = 0.0;
  
  
  i=-1;
  //i=i_peak-1;

  do
  {
    sum_Left_old = sum_Left_new;

    integrand = integrandF_derivatives_m_n(peak_position+h*i, k, eta, theta, m, n);

    sum_Left_new = sum_Left_old + integrand;
    i = i - step;
  }
  while  (sum_Left_old!=sum_Left_new); 
  
  
       new_result = h*(sum_Left_new  + sum_Right_new) + 0.5*old_result;

  return new_result;


}


long double Ffermi_estimate_derivatives_m_n_long(long double h, long double last_result, long double k, long double eta, long double theta, int m, int n)
{
  

  
  long double sum_Left_old = 0.0L , sum_Right_old = 0.0L ;
  long double sum_Left_new = 0.0L , sum_Right_new = 0.0L ;
  long double old_result, new_result, integrand;

  int step,i, i_peak=0;
  long double peak_position=0.0L;
  if(eta>1024.0L) peak_position = logl(eta);//+1.0/eta-1.0/eta/eta; Precise formula is Log[eta] + W(1/eta), where W is LambertW==ProductLog
  
  
  if(last_result==0.0L) /* ZERO value means first iteration, this should work also for derivatives, which can be negative ! Maybe nan/inf will be more appropriate*/
  {
    step=1;
    integrand = integrandF_derivatives_m_n_long(peak_position, k, eta, theta, m, n);
    old_result = 2.0L*h*integrand;
  }
  else
  {
    step=2;
    old_result = last_result;//Is this necessary? old_result===last_result?
  }
  
  /* integral for 0 < t < Infinity  */
  
  sum_Right_old = 0.0L;
  sum_Right_new = 0.0L;
  
  //if(eta>64.0) i_peak = (int) ceil(log(eta)/h); 
  
  i=1;
  //i=i_peak+1;

  do
  {
    sum_Right_old = sum_Right_new;

    integrand = integrandF_derivatives_m_n_long(peak_position+h*i, k, eta, theta, m, n);

      sum_Right_new = sum_Right_old + integrand;
    
	i = i + step;
  }
  while  ( (sum_Right_old!=sum_Right_new) /*|| (h*i<=peak_position)*/ ); 

  /* integral for -Infinity < t <0  */
  
  sum_Left_old = 0.0L;
  sum_Left_new = 0.0L;
  
  
  i=-1;
  //i=i_peak-1;

  do
  {
    sum_Left_old = sum_Left_new;

    integrand = integrandF_derivatives_m_n_long(peak_position+h*i, k, eta, theta, m, n);

    sum_Left_new = sum_Left_old + integrand;
    i = i - step;
  }
  while  (sum_Left_old!=sum_Left_new); 
  
  
       new_result = h*(sum_Left_new  + sum_Right_new) + 0.5L*old_result;

  return new_result;


}

__float128 Ffermi_estimate_derivatives_m_n_quad(__float128 h, __float128 last_result, __float128 k, __float128 eta, __float128 theta, int m, int n)
{
  

  
  __float128 sum_Left_old = 0.0q , sum_Right_old = 0.0q ;
  __float128 sum_Left_new = 0.0q , sum_Right_new = 0.0q ;
  __float128 old_result, new_result, integrand;

  int step,i, i_peak=0;
  __float128 peak_position=0.0q;
  if(eta>32.0q) peak_position = logq(eta);//+1.0/eta-1.0/eta/eta; Precise formula is Log[eta] + W(1/eta), where W is LambertW==ProductLog
  
  
  if(last_result==0.0q) /* ZERO value means first iteration, this should work also for derivatives, which can be negative ! Maybe nan/inf will be more appropriate*/
  {
    step=1;
    integrand = integrandF_derivatives_m_n_quad(peak_position, k, eta, theta, m, n);
    old_result = 2.0q*h*integrand;
  }
  else
  {
    step=2;
    old_result = last_result;//Is this necessary? old_result===last_result?
  }
  
  /* integral for 0 < t < Infinity  */
  
  sum_Right_old = 0.0q;
  sum_Right_new = 0.0q;
  
  //if(eta>64.0) i_peak = (int) ceil(log(eta)/h); 
  
  i=1;
  //i=i_peak+1;

  do
  {
    sum_Right_old = sum_Right_new;

    integrand = integrandF_derivatives_m_n_quad(peak_position+h*i, k, eta, theta, m, n);

      sum_Right_new = sum_Right_old + integrand;
    
	i = i + step;
  }
  while  ( (sum_Right_old!=sum_Right_new) /*|| (h*i<=peak_position)*/ ); 

  /* integral for -Infinity < t <0  */
  
  sum_Left_old = 0.0q;
  sum_Left_new = 0.0q;
  
  
  i=-1;
  //i=i_peak-1;

  do
  {
    sum_Left_old = sum_Left_new;

    integrand = integrandF_derivatives_m_n_quad(peak_position+h*i, k, eta, theta, m, n);

    sum_Left_new = sum_Left_old + integrand;
    i = i - step;
  }
  while  (sum_Left_old!=sum_Left_new); 
  
  
       new_result = h*(sum_Left_new  + sum_Right_new) + 0.5q*old_result;

  return new_result;


}



void Ffermi_dblexp_derivatives_matrix(const double k, const double eta, const double theta,
  const double precision, const int recursion_limit, double result[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE])
{
 
  /* I hope I understand correctly https://gcc.gnu.org/onlinedocs/gcc-4.1.2/gcc/Designated-Inits.html */ 
  double old[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE]={ [0][0] = -1.0 }; //Setting old to -1.0 cause Ffermi_estimate_derivatives to restart at the first call
  double new[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE]={ [0][0] =  0.0 };
  double h=0.5; //initial dbl. exp. step
  int m,n;
  
  //if(k<=-1.0) return nan("NaN"); /* not converging for k <= -1 */

  Ffermi_estimate_derivatives_matrix(h, old, k, eta, theta, new);

    for(m=0;m<DERIVATIVE_MATRIX_SIZE;m++) 
     for(n=0;n<DERIVATIVE_MATRIX_SIZE;n++)  
      old[m][n] = 0.0;/*Is this necessary? Two lines below we reset old to new which is zero anyway. 
         Except precision goal is achieved at first run (in theory possible, if one modify code and set e.g h=0.125 or less what MIGHT may have sense in future
   optimization ) */ 

  
  while( fabs(old[0][0]-new[0][0])>precision*fabs(new[0][0]) && h>pow(2.0,-recursion_limit))
  {
    for(m=0;m<DERIVATIVE_MATRIX_SIZE;m++) 
     for(n=0;n<DERIVATIVE_MATRIX_SIZE;n++)  
      old[m][n]=new[m][n];
    
    h=0.5*h;
    Ffermi_estimate_derivatives_matrix(h, old, k, eta, theta, new);
  }

    for(m=0;m<DERIVATIVE_MATRIX_SIZE;m++) 
     for(n=0;n<DERIVATIVE_MATRIX_SIZE;n++)  
      result[m][n]=new[m][n];
    
}


double Ffermi_dblexp_derivatives_m_n(const double k, const double eta, const double theta, const int m, const int n,
  const double precision, const int recursion_limit)
{
 
  double old =  0.0 ; //Setting old to 0.0 cause Ffermi_estimate_derivatives to restart at the first call
  double new =  0.0 ;
  double h=0.5; //initial dbl. exp. step
  //if(eta>4.0) h = log(eta); // this force quadrature to hit peak of the transformed integrand, peak position for sigma(eta-x) is log(eta)+W(1/eta) = log(eta)+1/eta-1/eta^2+....
  
  
  
  //if(k<=-1.0) return nan("NaN"); /* not converging for k <= -1 */

  new = Ffermi_estimate_derivatives_m_n(h, old, k, eta, theta, m, n);

      old = 0.0;/*Is this necessary? Two lines below we reset old to new which is zero anyway. 
         Except precision goal is achieved at first run (in theory possible, if one modify code and set e.g h=0.125 or less what MIGHT may have sense in future
   optimization ) */ 

  
  while( fabs(old-new)>precision*fabs(new) && h>pow(2.0,-recursion_limit))
  {
    old=new;
    h=0.5*h;
    new = Ffermi_estimate_derivatives_m_n(h, old, k, eta, theta, m, n);
  }

    return new;
    
}

long double Ffermi_dblexp_derivatives_m_n_long(const long double k, const long double eta, const long double theta, const int m, const int n,
  const long double precision, const int recursion_limit)
{
 
  long double old =  0.0L ; //Setting old to 0.0 cause Ffermi_estimate_derivatives to restart at the first call
  long double new =  0.0L ;
  long double h=0.5L; //initial dbl. exp. step
  //if(eta>4.0) h = log(eta); // this force quadrature to hit peak of the transformed integrand, peak position for sigma(eta-x) is log(eta)+W(1/eta) = log(eta)+1/eta-1/eta^2+....
  
  
  
  //if(k<=-1.0) return nan("NaN"); /* not converging for k <= -1 */

  new = Ffermi_estimate_derivatives_m_n_long(h, old, k, eta, theta, m, n);

      old = 0.0L; 

  
  while( fabsl(old-new)>precision*fabsl(new) && h>powl(2.0L,-recursion_limit))
  {
    old=new;
    h=0.5L*h;
    new = Ffermi_estimate_derivatives_m_n_long(h, old, k, eta, theta, m, n);
  }

    return new;
    
}


__float128 Ffermi_dblexp_derivatives_m_n_quad(const __float128 k, const __float128 eta, const __float128 theta, const int m, const int n,
  const __float128 precision, const int recursion_limit)
{
 
  __float128 old =  0.0q ; //Setting old to 0.0 cause Ffermi_estimate_derivatives to restart at the first call
  __float128 new =  0.0q ;
  __float128 h=0.5q; //initial dbl. exp. step
  //if(eta>4.0) h = log(eta); // this force quadrature to hit peak of the transformed integrand, peak position for sigma(eta-x) is log(eta)+W(1/eta) = log(eta)+1/eta-1/eta^2+....
  
  
  
  //if(k<=-1.0) return nan("NaN"); /* not converging for k <= -1 */

  new = Ffermi_estimate_derivatives_m_n_quad(h, old, k, eta, theta, m, n);

      old = 0.0q; 

  
  while( fabsq(old-new)>precision*fabsq(new) && h>powq(2.0q,-recursion_limit))
  {
    old=new;
    h=0.5q*h;
    new = Ffermi_estimate_derivatives_m_n_quad(h, old, k, eta, theta, m, n);
  }

    return new;
    
}


/* 

i-th Sommerfeld term for partial derivative D^n/Dtheta^n D^m/Deta^m is:

2 DirichletEta[2 i] Derivative[i][f][\[Eta]]], {\[Theta], n}, {\[Eta], m}]

NOTE: formula for eta derivatives, apart from 2 DirichletEta[2i] term, is simply shifted formula for i-th
term. Therefore, once we have computed 1-st expansion for third derivative, we already have
3-rd order expansion for first derivative !

*/


/* Formula below computes D[eta^k Sqrt[1+eta*theta/2],{eta,m},{theta,n}] */
void sommerfeld_derivatives(const double k, const double eta, const double theta, double D[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE])
{
    #include "factorial.h"
    double sign;
    double eta_k;
    double z1 = 1.0 + 0.5*eta*theta;
    double sum=0.0;
    int i,m,n;

//FIXME: below is a very unoptimized code
    for(m=0;m<=3;m++)
     for(n=0;n<=3;n++)
      {
       sign = (n%2) ? 1.0 : -1.0; // (-1)^(n+1) 
       eta_k = pow(eta,k - m + n);
       sum=0.0;
       for(i=0;i<=m;i++)
         sum = sum + binom(m,i)*fac2(2*n + 2*m - 3 - 2*i)*pow(2.0, i - m - 2*n)*pochhammer(k + 1.5 - m, i)*pow(z1,0.5 + i - m - n);
       
       D[m][n] = sign*eta_k*sum;
      }


}

/* Formula below computes D[eta^k Sqrt[1+eta*theta/2],{eta,m},{theta,n}] */
double sommerfeld_derivatives_m_n(const double k, const double eta, const double theta, const int m, const int n)
{
    #include "factorial.h"
    double sign;
    double eta_k;
    double z1 = 1.0 + 0.5*eta*theta;
    double sum=0.0;
    int i;

       sign = ((n%2)==0) ? -1.0 : 1.0; // (-1)^(n+1);
       eta_k = pow(eta,k - m + n);
       for(i=0;i<=m;i++)
         sum = sum + binom(m,i)*fac2(2*n + 2*m - 3 - 2*i)*pow(2.0, i - m - 2*n)*pochhammer(k + 1.5 - m, i)*pow(z1,0.5 + i - m - n);
       
       return sign*eta_k*sum;


}

__float128 sommerfeld_derivatives_m_n_quad(const __float128 k, const __float128 eta, const __float128 theta, const int m, const int n)
{
    #include "factorial_quad.h"
    __float128 sign;
    __float128 eta_k;
    __float128 z1 = 1.0q + 0.5q*eta*theta;
    __float128 sum=0.0q;
    int i;

       sign = ((n%2)==0) ? -1.0q : 1.0q; // (-1)^(n+1);
       eta_k = powq(eta,k - m + n);
       for(i=0;i<=m;i++)
         sum = sum + binom_quad(m,i)*fac2_quad(2*n + 2*m - 3 - 2*i)*powq(2.0, i - m - 2*n)*pochhammer_quad(k + 1.5 - m, i)*powq(z1,0.5 + i - m - n);
       
       return sign*eta_k*sum;


}

/* TODO: error control not implemented ! */

void Ffermi_sommerfeld_derivatives(const double k, const double eta, const double theta, const double precision, const int SERIES_TERMS_MAX, double result[10])
{
	double z = -0.5*eta*theta,eta_k=pow(eta,k);
    double z1=1.0-z;
    double sqrt_1z = sqrt(z1);
    double S[DERIVATIVE_MATRIX_SIZE], derivatives[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE];
    int n,m; // order of partial derivatives with respect to theta and eta, respectively
    #include "factorial.h"

	int i,j;
    double derivative;
    /* S[z_] := Hypergeometric2F1[-1/2, 1 + k, 2 + k, z] */
    sommerfeld_leading_term_derivatives(k,z,S);

    /* Leading term */
    /* NOTE/FIXME: eta derivatives DO NOT require 2F1 function, and can be computed using std. math */

	result[0] =  eta_k*eta/(1.0+k)*S[0];
    //result[1] =  eta_k*S[0]-eta_k*eta*theta*S[1]/(2.0+2.0*k); 
    result[1] =  eta_k*sqrt_1z;
    //result[2] =  eta_k/eta*k*S[0]-eta_k*theta*S[1]+eta_k*eta*theta*theta*S[2]/(4.0+4.0*k);
    result[2] =  result[1]/eta*(0.5+k-0.5/z1);
    result[3] = -eta_k*eta*eta*S[1]/(2.0+2.0*k);
    result[4] =  eta_k*eta*eta*eta*S[2]/(4.0+4.0*k);
    result[5] =  eta_k*eta*(eta*theta*S[2]-(4.0+2.0*k)*S[1])/(4.0+4.0*k);
    //result[5] =  eta_k*eta*0.25/sqrt_1z;
    result[6] = -eta_k*eta*eta*eta*eta*S[3]/(8.0+8.0*k);
    result[7] =  eta_k*eta*eta*(2.0*(k+3.0)*S[2] - eta*theta*S[3])/8.0/(1.0+k);
    //result[7] =  eta_k*eta*eta/16.0/pow(1.0-z,1.5);
    result[8] = -eta_k*(4.0*(2.0+3.0*k+k*k)*S[1]+eta*theta*((-8.0-4.0*k)*S[2]+eta*theta*S[3]))/(8.0+8.0*k);
    result[9] =  k*(k-1.0)*eta_k/(eta*eta)*S[0]
                    -1.5*k*eta_k/eta*theta*S[1]
                   +0.75*eta_k*theta*theta*S[2]
  -eta_k*eta*theta*theta*theta/(8.0+8.0*k)*S[3];

	if(SERIES_TERMS_MAX<1) return;

    //Compute partial derivatives at i-th Sommerfeld expansion order
    i=1;
    for(n=0;n<=3;n++)
      for(m=0;m<=3;m++)
        {
          if(m+n>DERIVATIVE_MAX_ORDER) continue; //we do not need higher order derivatives for now
          
          derivatives[n][m] = sommerfeld_derivatives_m_n(k, eta, theta, m+2*i-1, n);


        }

    

    result[0] = result[0] + 2.0*etaTBL_odd[i]*derivatives[0][0];//+2.0*etaTBL[2*(i+1)]*derivatives[0][2];
    result[1] = result[1] + 2.0*etaTBL_odd[i]*derivatives[0][1];//+2.0*etaTBL[2*(i+1)]*derivatives[0][3];
    result[2] = result[2] + 2.0*etaTBL_odd[i]*derivatives[0][2];
    result[3] = result[3] + 2.0*etaTBL_odd[i]*derivatives[1][0];
    result[4] = result[4] + 2.0*etaTBL_odd[i]*derivatives[2][0];
    result[5] = result[5] + 2.0*etaTBL_odd[i]*derivatives[1][1];//+2.0*etaTBL[2*(i+1)]*derivatives[1][3];
    result[6] = result[6] + 2.0*etaTBL_odd[i]*derivatives[3][0];
    result[7] = result[7] + 2.0*etaTBL_odd[i]*derivatives[2][1];
    result[8] = result[8] + 2.0*etaTBL_odd[i]*derivatives[1][2];
    result[9] = result[9] + 2.0*etaTBL_odd[i]*derivatives[0][3];



	if(SERIES_TERMS_MAX<=1) return;

    for(i=2;i<=SERIES_TERMS_MAX;i++)
     {
      for(n=0;n<=3;n++)
        for(m=0;m<=3;m++)
          {
            if(m+n>DERIVATIVE_MAX_ORDER) continue; //we do not need higher order derivatives for now
            if(m+n<=1){ derivatives[n][m] = derivatives[n][m+2]; continue;}  //re-use already computed eta derivatives
            derivatives[n][m] = sommerfeld_derivatives_m_n(k, eta, theta, m+2*i-1, n);
          }
      result[0] = result[0] + 2.0*etaTBL_odd[i]*derivatives[0][0];
      //result[0] = result[0] /* Already computed */            +2.0*etaTBL[2*(i+1)]*derivatives[0][2];
      result[1] = result[1] + 2.0*etaTBL_odd[i]*derivatives[0][1];
      //result[1] = result[1] /* Already computed */            +2.0*etaTBL[2*(i+1)]*derivatives[0][3];
      result[2] = result[2] + 2.0*etaTBL_odd[i]*derivatives[0][2];
      result[3] = result[3] + 2.0*etaTBL_odd[i]*derivatives[1][0];
      result[4] = result[4] + 2.0*etaTBL_odd[i]*derivatives[2][0];
      result[5] = result[5] + 2.0*etaTBL_odd[i]*derivatives[1][1];
      //result[5] = result[5] /* Already computed */          +2.0*etaTBL[2*(i+1)]*derivatives[1][3];
      result[6] = result[6] + 2.0*etaTBL_odd[i]*derivatives[3][0];
      result[7] = result[7] + 2.0*etaTBL_odd[i]*derivatives[2][1];
      result[8] = result[8] + 2.0*etaTBL_odd[i]*derivatives[1][2];
      result[9] = result[9] + 2.0*etaTBL_odd[i]*derivatives[0][3]; 
     }

}

/* TODO: error control not implemented ! */

void Ffermi_sommerfeld_derivatives_matrix(const double k, const double eta, const double theta, const double precision, const int SERIES_TERMS_MAX, double result[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE])
{
	double z = -0.5*eta*theta,eta_k=pow(eta,k);
    double z1=1.0-z;
    double sqrt_1z = sqrt(z1);
    double S[DERIVATIVE_MATRIX_SIZE], derivatives[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE];
    int n,m; // order of partial derivatives with respect to theta and eta, respectively
    #include "factorial.h"

	int i,j;
    double derivative;


    /* S[z_] := Hypergeometric2F1[-1/2, 1 + k, 2 + k, z] */
    //sommerfeld_leading_term_derivatives(k,z,S);

    sommerfeld_leading_term_derivatives_matrix(k, eta, theta, result);

/*
	result[0][0] =  eta_k*eta/(1.0+k)*S[0];
    result[1][0] =  eta_k*sqrt_1z;
    result[2][0] =  result[1][0]/eta*(0.5+k-0.5/z1);
    result[0][1] = -eta_k*eta*eta*S[1]/(2.0+2.0*k);
    result[0][2] =  eta_k*eta*eta*eta*S[2]/(4.0+4.0*k);
    result[1][1] =  eta_k*eta*(eta*theta*S[2]-(4.0+2.0*k)*S[1])/(4.0+4.0*k);
    result[0][3] = -eta_k*eta*eta*eta*eta*S[3]/(8.0+8.0*k);
    result[1][2] =  eta_k*eta*eta*(2.0*(k+3.0)*S[2] - eta*theta*S[3])/8.0/(1.0+k);
    result[2][1] = -eta_k*(4.0*(2.0+3.0*k+k*k)*S[1]+eta*theta*((-8.0-4.0*k)*S[2]+eta*theta*S[3]))/(8.0+8.0*k);
    result[3][0] =  k*(k-1.0)*eta_k/(eta*eta)*S[0]
                    -1.5*k*eta_k/eta*theta*S[1]
                   +0.75*eta_k*theta*theta*S[2]
  -eta_k*eta*theta*theta*theta/(8.0+8.0*k)*S[3];
*/
	if(SERIES_TERMS_MAX<1) return;

    //Compute partial derivatives at i-th Sommerfeld expansion order
    //FIXME: some/most(?) of them are already computed above! 
    i=1;
    for(n=0;n<=3;n++)
      for(m=0;m<=3;m++)
        {
          if(m+n>DERIVATIVE_MAX_ORDER) continue; //we do not need higher order derivatives for now
          
          derivatives[m][n] = sommerfeld_derivatives_m_n(k, eta, theta, m+2*i-1, n);


        }

    
    for(m=0;m<DERIVATIVE_MATRIX_SIZE;m++)
     for(n=0;n<DERIVATIVE_MATRIX_SIZE;n++)
      result[m][n] = result[m][n] + 2.0*etaTBL_odd[i]*derivatives[m][n];



	if(SERIES_TERMS_MAX<=1) return;

    for(i=2;i<=SERIES_TERMS_MAX;i++)
     {
      for(n=0;n<=3;n++)
        for(m=0;m<=3;m++)
          {
            if(m+n>DERIVATIVE_MAX_ORDER) continue; //we do not need higher order derivatives for now
            if(m+n<=1){ derivatives[n][m] = derivatives[n][m+2]; continue;}  //re-use already computed eta derivatives
            derivatives[m][n] = sommerfeld_derivatives_m_n(k, eta, theta, m+2*i-1, n);
            result[m][n] = result[m][n] + 2.0*etaTBL_odd[i]*derivatives[m][n];
  
          }
     }

}



double Ffermi_sommerfeld_derivatives_m_n(const double k, const double eta, const double theta, const int m, const int n, const double precision, const int SERIES_TERMS_MAX)
{
	double z = -0.5*eta*theta,eta_k=pow(eta,k);
    double z1=1.0-z;
    double sqrt_1z = sqrt(z1);
    double result = 0.0,sum_new,sum_old;
    
    #include "factorial.h"

	int i,j;
    double derivative;


    /* S[z_] := Hypergeometric2F1[-1/2, 1 + k, 2 + k, z] */
    //sommerfeld_leading_term_derivatives(k,z,S);

    
    result = sommerfeld_leading_term_derivatives_m_n(k, eta, theta, m, n);
	if(SERIES_TERMS_MAX<1) return result; 

    //Compute partial derivatives at i-th Sommerfeld expansion order
    //FIXME: some/most(?) of them are already computed above! 

/*
    i=1;
    if(m+n>DERIVATIVE_MAX_ORDER) return -1.0; //FIXME: return nan
          
    derivative = sommerfeld_derivatives_m_n(k, eta, theta, m+2*i-1, n);

    result = result + 2.0*etaTBL_odd[i]*derivative;



	if(SERIES_TERMS_MAX<=1) return result;


    for(i=2;i<=SERIES_TERMS_MAX;i++)
     {
       derivative = sommerfeld_derivatives_m_n(k, eta, theta, m+2*i-1, n);
       result = result + 2.0*etaTBL_odd[i]*derivative;
     }
    return result; 
*/     
    i=1;
    sum_new = result; 
	do
	{
		sum_old = sum_new;
        derivative = sommerfeld_derivatives_m_n(k, eta, theta, m+2*i-1, n);
		sum_new = sum_old + 2.0*etaTBL_odd[i]*derivative;
		i++;
	}
    while ( ( (precision>0.0) ? fabs(sum_old-sum_new)>=precision*sum_new : sum_old!=sum_new )  && i<SERIES_TERMS_MAX );



    return sum_new;
     

}


__float128 Ffermi_sommerfeld_derivatives_m_n_quad(const __float128 k, const __float128 eta, const __float128 theta, const int m, const int n, const __float128 precision, const int SERIES_TERMS_MAX)
{
	__float128 z = -0.5q*eta*theta,eta_k=powq(eta,k);
    __float128 z1=1.0q-z;
    __float128 sqrt_1z = sqrtq(z1);
    __float128 result = 0.0q,sum_new,sum_old;
    
    #include "factorial_quad.h"

	int i,j;
    __float128 derivative;


    /* S[z_] := Hypergeometric2F1[-1/2, 1 + k, 2 + k, z] */
    //sommerfeld_leading_term_derivatives(k,z,S);

    
    result = sommerfeld_leading_term_derivatives_m_n_quad(k, eta, theta, m, n);
	if(SERIES_TERMS_MAX<1) return result; 
    
    i=1;
    sum_new = result; 
	do
	{
		sum_old = sum_new;
        derivative = sommerfeld_derivatives_m_n_quad(k, eta, theta, m+2*i-1, n);
		sum_new = sum_old + 2.0q*etaTBL_odd_quad[i]*derivative;
		i++;
	}
    while ( ( (precision>0.0q) ? fabs(sum_old-sum_new)>=precision*sum_new : sum_old!=sum_new )  && i<SERIES_TERMS_MAX );



    return sum_new;
     

}


void Ffermi_derivatives(const double k, const double eta, const double theta, double result[10])
{
   
   
   if( eta>2048.0) //original 56000.0
    {
	  Ffermi_sommerfeld_derivatives(k, eta, theta, PRECISION_GOAL, 2, result);
    }
  else
    {
      Ffermi_dblexp_derivatives(k,eta,theta,PRECISION_GOAL, MAX_REFINE, result);
    }
    
    return;
}

void Ffermi_derivatives_matrix(const double k, const double eta, const double theta, double FD[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE])
{
   

   if( eta>512.0) 
    {
	  Ffermi_sommerfeld_derivatives_matrix(k, eta, theta, PRECISION_GOAL, 5, FD);
    }
  else
    {
      Ffermi_dblexp_derivatives_matrix(k,eta,theta,PRECISION_GOAL, MAX_REFINE, FD);
    }
    
    return;
}


double Ffermi_derivatives_m_n(const double k, const double eta, const double theta, const int m, const int n)
{
    /* FIXME: improve initializer to work with 3-rd order derivs, but remain general*/
    //double eta[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE] = { [0][0] = 8192.0 };
    double eta_s[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE] = {{2048.0, 4096.0, 4096.0, 4096.0},
                         {8192.0,  256.0,   64.0,   64.0},
                         {  64.0,   64.0,   64.0,   64.0},
                         {  64.0,   64.0,   64.0,   64.0}
                         };
  
   
    
   if( eta>eta_s[m][n]) 
    {
	  return Ffermi_sommerfeld_derivatives_m_n(k, eta, theta, m, n, PRECISION_GOAL, 2); //Sommerfeld order might be m,n dependent or adaptive
    }
  else
    {
      return Ffermi_dblexp_derivatives_m_n(k,eta,theta, m, n, PRECISION_GOAL, MAX_REFINE);
    }
    
    
}


__float128 Ffermi_derivatives_m_n_quad(const __float128 k, const __float128 eta, const __float128 theta, const int m, const int n)
{
    /* FIXME: improve initializer to work with 3-rd order derivs, but remain general*/
    //double eta[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE] = { [0][0] = 8192.0 };
    __float128 eta_s[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE] = {{2048.0q, 4096.0q, 4096.0q, 4096.0q},
                              {8192.0q,  256.0q,   64.0q,   64.0q},
                              {  64.0q,   64.0q,   64.0q,   64.0q},
                              {  64.0q,   64.0q,   64.0q,   64.0q}
                              };
  
   
    
   if( eta>eta_s[m][n]) 
    {
	  return Ffermi_sommerfeld_derivatives_m_n_quad(k, eta, theta, m, n, PRECISION_GOAL_QUAD, 6); //Sommerfeld order might be m,n dependent or adaptive
    }
  else
    {
      return Ffermi_dblexp_derivatives_m_n_quad(k,eta,theta, m, n, PRECISION_GOAL_QUAD, MAX_REFINE);
    }
    
    
}
