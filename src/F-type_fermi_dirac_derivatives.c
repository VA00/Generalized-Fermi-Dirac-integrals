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
#define DEBUG 0
/* MAX_REFINE limit recursion depth for FFermi.
Note, that convergence might be slow, and
using very large large MAX_REFINE>16
together with high PRECISION settings
close to DBL_EPSILON results in extremely
slow computations.
*/
#define MAX_REFINE 16 // more than 16 result is significant slow-down
//#define PRECISION sqrt(DBL_EPSILON) // convergence is exponential, so in theory this is enough
#define PRECISION 8*DBL_EPSILON   // down to 2*DBL_EPSILON seem harmless, 1*DBL_EPSILON cause problems
//#define PRECISION pow(DBL_EPSILON,0.6666)
#define KAHAN 0 // Enable https://en.wikipedia.org/wiki/Kahan_summation_algorithm ; usually this has no significant effect, but results might be not identical, and computation slow
#define TGAMMA_MAX 170.62437695630272081244437878577 // FindInstance[LogGamma[k + 1] == Log[2^1024] tgamma overflow


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



void Ffermi_value_derivatives(const double k, const double eta, const double theta,
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


/* 

i-th Sommerfeld term for partial derivative D^n/Dtheta^n D^m/Deta^m is:

2 DirichletEta[2 i] Derivative[i][f][\[Eta]]], {\[Theta], n}, {\[Eta], m}]

NOTE: formula for eta derivatives, apart from 2 DirichletEta[2i] term, is simply shifted formula for i-th
term. Therefore, once we have computed 1-st expansion for third derivative, we already have
3-rd order expansion for first derivative !

*/


/* Formula below computes D[eta^k Sqrt[1+eta*theta/2],{eta,m},{theta,n}] */
void sommerfeld_derivatives(const double k, const double eta, const double theta, double D[4][4])
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
       sign = (n%2) ? 1.0 : -1.0; // (-1)^(n+1);
       eta_k = pow(eta,k - m + n);
       sum=0.0;
       for(i=0;i<=m;i++)
         sum = sum + binom(m,i)*fac2(2*n + 2*m - 3 - 2*i)*pow(2.0, i - m - 2*n)*pochhammer(k + 1.5 - m, i)*pow(z1,0.5 + i - m - n);
       
       D[m][n] = sign*eta_k*sum;
      }


}

/* Formula below computes D[eta^k Sqrt[1+eta*theta/2],{eta,m},{theta,n}] */
double sommerfeld_derivatives_m_n(const double k, const double eta, const double theta, int m, int n)
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

/* TODO: error control not implemented ! */

void Ffermi_sommerfeld_derivatives(const double k, const double eta, const double theta, const double precision, const int SERIES_TERMS_MAX, double result[10])
{
	double z = -0.5*eta*theta,eta_k=pow(eta,k);
    double z1=1.0-z;
    double sqrt_1z = sqrt(z1);
    double S[4], derivatives[4][4];
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
          if(m+n>3) continue; //we do not need higher order derivatives for now
          
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
            if(m+n>3) continue; //we do not need higher order derivatives for now
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