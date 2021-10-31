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
  
  double old[10]={ [ 0 ... 9 ] = -1.0 }; //Setting old to -1.0 cause Ffermi_estimate_derivatives to restart
  double new[10]={ [ 0 ... 9 ] =  0.0 };
  double h=0.5;
  int j;
  
  //if(k<=-1.0) return nan("NaN"); /* not converging for k <= -1 */

  Ffermi_estimate_derivatives(h, old, k, eta, theta, new);

  for(j=0;j<10;j++) old[j] = 0.0;//Is this necessary? 
  
  while( fabs(old[0]-new[0])>precision*fabs(new[0]) && h>pow(2.0,-recursion_limit))
  {
    for(j=0;j<10;j++) old[j]=new[j];
    h=0.5*h;
    Ffermi_estimate_derivatives(h, old, k, eta, theta, new);
  }

  for(j=0;j<10;j++) result[j]=new[j];
    
}