/*
A. Odrzywolek, AOdrzywolek 
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


*/



double integrandF(const double t, const double k, const double eta, const double theta)
{
  
  double x,dx,integrand,result,factor;
  
  /* if(t>-6.5)
  * 
  * this is min t=-9.3, for which exp(t-exp(-t)) 
    is still smaller than LDBL_MIN defined in <float.h>. For DBL_MIN it is t>-6.5, but 
    using proper (unsafe?) coding  modern CPU's do calculations internally in long double format anyway. 
    NOTE: obsolete, see comment below where optimal coding for integrand is described.
    
    */

  
  x = exp(  t - exp(-t) ); /* Masatake Mori, eq. (4.17) */
  //if( (eta>k) && (k>0) ) x = eta*exp(t-exp(-t)); else x = exp(t-exp(-t));
  
  
  

    
  //dx = x*(1 + exp(-t) ); /* dx/dt */
  dx = 1.0+exp(-t); /* in this case x is adsorbed in integrand, and x^k -> x^(k+1) */
  
  
  if(x-eta<-log(DBL_EPSILON)) // if using machine precison we are unable to add 1.0 to exp(), then approximation is optimal
    {
	factor = 1.0/(1.0+exp(x-eta) );
	integrand = exp( (k+1.0)*(t - exp(-t)) );
	//integrand = pow(x,k+1.0);
	integrand = integrand*sqrt(1.0+0.5*theta*x)*factor;
	}
  else
    {
    //factor = exp(eta-x) adsorbed into exp, to avoid 0*infinity mess 
	integrand = exp((k+1.0)*(t - exp(-t)) + eta - x );
	integrand = integrand*sqrt(1.0+ 0.5*theta*x);
    }
  
  /* NOTE:
   * 
   * if we use:
   * 
   * integrand = pow(x,k+1.0)*sqrt(1.0+ 0.5*theta*x)*factor;
   * 
   * then:
   * 
   * a) precision is lost, beacuse  x is double, while exp((k+1.0)*(t - exp(-t)) ) 
   * is internally handled as long double (96 bit)
   * b) if k<0 we lost advantage of postponed underflow ( k+1 << 1 in such a case )
   * 
  */
  

  
#if DEBUG  
  printf("DEBUG300: factor = %.20Lf, x=%.20Lf, dx=%.20Lf, integrand=%.20Lf, return = %.20Lf \t test= %.20Lf \n",factor, x,dx,integrand, (integrand*dx),test);
#endif  
  
  result = integrand*dx;
  
  
  return result;
  
  
}

long double integrandF_long(const long double t, const long double k, const long double eta, const long double theta)
{
  
  long double x,dx,integrand,result,factor;
  
  
  
  
  //const double lambda = M_E*100.5;//scaling factor

 /* if(t>-6.5)
  * 
  * this is min t=-9.3, for which exp(t-exp(-t)) 
    is still smaller than LDBL_MIN defined in <float.h>. For DBL_MIN it is t>-6.5, but 
    using proper (unsafe?) coding  modern CPU's do calculations internally in long double format anyway. 
    NOTE: obsolete, see comment below where optimal coding for integrand is described.
    
    */

  
  x = expl(  t - expl(-t) ); /* Masatake Mori, eq. (4.17) */
  
  
  //dx = x*(1 + exp(-t) ); /* dx/dt */
  dx = 1.0L+exp(-t); /* in this case x is adsorbed in integrand, and x^k -> x^(k+1) */
  
  
  if(x-eta<-logl(LDBL_EPSILON)) // if using machine precison we are unable to add 1.0 to exp(), then approximation is optimal
    {
	factor = 1.0L/(1.0L+expl(x-eta) );
	//integrand = expl( (kL+1.0L)*(tL - expl(-tL)) );
	integrand = powl(x,k+1.0L);
	integrand = integrand*sqrtl(1.0L+0.5L*theta*x)*factor;
	}
  else
    {
    //factor = exp(eta-x) adsorbed into exp, to avoid 0*infinity mess 
	integrand = expl((k+1.0L)*(t - expl(-t)) + eta - x );
	integrand = integrand*sqrtl(1.0L+ 0.5L*theta*x);
    }
  
  result = integrand*dx;
  
  
  return result;
  
  
}



double Ffermi_estimate(double h, double last_result, double k, double eta, double theta)
{
  
  int step,i;
  double sum_Left_old, sum_Right_old;
  double sum_Left_new, sum_Right_new;
  double old_result, new_result;
  #if KAHAN
  double c=0.0,t,y; // https://en.wikipedia.org/wiki/Kahan_summation_algorithm
  #endif 
  
  
  if(last_result<0.0) /* Negative value means first iteration*/
  {
    step=1;
    old_result = 2.0*h*integrandF(0.0, k, eta, theta);
  }
  else
  {
    step=2;
    old_result = last_result;
  }

#if DEBUG 
    
  printf("DEBUG2: old=%e,\tlast=%e\n",old_result,last_result);  

#endif

  
  /* integral for 0 < t < Infinity  */
  
  sum_Right_old = 0.0;
  sum_Right_new = 0.0;
  
  
  i=1;
  /* possible vectorization, but loop step must be known at compile time!
  #pragma omp simd
  #pragma ivdep
  for(i=1;i<=16;i+=2)
  {
    sum_Right_new += integrandF(h*i, k, eta, theta);
  }
  */
  do
  {
	sum_Right_old = sum_Right_new;
    #if KAHAN	
	y = integrandF(h*i, k, eta, theta) - c;
	t = sum_Right_new + y;
	c = (t-sum_Right_new) - y;
	sum_Right_new = t;
    #else 
    sum_Right_new = sum_Right_old + integrandF(h*i, k, eta, theta);
    //sum_Right_new = sum_Right_old + integrandF(h*i, k, eta, theta);
    #endif	
	i = i + step;
  }
  while  ( sum_Right_old<sum_Right_new ); //floating point fixed-point method

  /* integral for -Infinity < t <0  */
  
  sum_Left_old = 0.0;
  sum_Left_new = 0.0;
  #if KAHAN
  c = 0.0;
  #endif
  
  
  i=-1;
  do
  {
	sum_Left_old = sum_Left_new;
    #if KAHAN	
	y = integrandF(h*i, k, eta, theta) - c;
	t = sum_Left_new + y;
	c = (t-sum_Left_new) - y;
	sum_Left_new = t;
    #else 
    sum_Left_new = sum_Left_old + integrandF(h*i, k, eta, theta);
    #endif	
    i = i - step;
  }
  while  (sum_Left_old<sum_Left_new);
  
  
  new_result = h*(sum_Left_new  + sum_Right_new) + 0.5*old_result;

  return new_result;
}


long double Ffermi_estimate_long(long double h, long double last_result, long double k, long double eta, long double theta)
{
  
  int step,i;
  long double sum_Left_old, sum_Right_old;
  long double sum_Left_new, sum_Right_new;
  long double old_result, new_result;
  
  
  if(last_result<0.0L) /* Negative value means first iteration*/
  {
    step=1;
    old_result = 2.0L*h*integrandF_long(0.0L, k, eta, theta);
  }
  else
  {
    step=2;
    old_result = last_result;
  }
  
  /* integral for 0 < t < Infinity  */
  
  sum_Right_old = 0.0;
  sum_Right_new = 0.0;
  
  
  i=1;

  do
  {
	sum_Right_old = sum_Right_new;
    sum_Right_new = sum_Right_old + integrandF_long(h*i, k, eta, theta);
	i = i + step;
  }
  while  ( sum_Right_old<sum_Right_new ); //floating point fixed-point method

  /* integral for -Infinity < t <0  */
  
  sum_Left_old = 0.0;
  sum_Left_new = 0.0;
  
  
  i=-1;
  do
  {
	sum_Left_old = sum_Left_new;
    sum_Left_new = sum_Left_old + integrandF_long(h*i, k, eta, theta);
    i = i - step;
  }
  while  (sum_Left_old<sum_Left_new);
  
  
  new_result = h*(sum_Left_new  + sum_Right_new) + 0.5L*old_result;

  return new_result;
}

double Ffermi_value(const double k, const double eta, const double theta,
  const double precision, const int recursion_limit)
{
  
  double old=0.0, new=0.0, h=0.5;
  
  if(k<=-1.0) return nan("NaN"); /* not converging for k <= -1 */

#if DEBUG 
    
  printf("DEBUG0: h=%lf,\tval=%e\n",h,new);  

#endif

  
  old = 0.0;
  new = Ffermi_estimate(h, -1.0, k, eta, theta);

#if DEBUG 
    
  printf("DEBUG1: h=%lf,\tval=%e\n",h,new);  

#endif


  
  while( fabs(old-new)>precision*fabs(new) && h>pow(2.0,-recursion_limit))
  {
    old=new;
    h=0.5*h;
    new = Ffermi_estimate(h, old, k, eta, theta);

#if DEBUG 
    
  printf("DEBUG4: h=%lf,\tval=%e\n",h,new);  

#endif
    
  }
  
  
  return new;  
    
}


long double Ffermi_value_long(const long double k, const long double eta, const long double theta,  const long double precision, const int recursion_limit)
{
  
  long double old=0.0L, new=0.0L, h=0.5L;
  
  if(k<=-1.0L) return nan("NaN"); /* not converging for k <= -1 */

  
  old = 0.0L;
  new = Ffermi_estimate_long(h, -1.0L, k, eta, theta);

  while( fabsl(old-new)>precision*fabsl(new) && h>powl(2.0L,-recursion_limit))
  {
    old=new;
    h=0.5L*h;
    new = Ffermi_estimate_long(h, old, k, eta, theta);
  }
  
  
  return new;  
    
}


/* TODO: error control not implemented ! */
double Ffermi_sommerfeld(const double k, const double eta, const double theta, const double precision, const int SERIES_TERMS_MAX)
{
	double leading_term, derivative,asymptotic_terms=0.0;
	int i,j;
	double etaTBL[12] = {0.50000000000000000000000000000000, \
                         0.69314718055994530941723212145818, \
                         0.82246703342411321823620758332301, \
                         0.90154267736969571404980362113359, \
                         0.94703282949724591757650323447352, \
                         0.97211977044690930593565514355347, \
                         0.98555109129743510409843924448495, \
                         0.99259381992283028267042571313339, \
                         0.99623300185264789922728926008280, \
                         0.99809429754160533076778303185260, \
                         0.99903950759827156563922184569934, \
                         0.99951714349806075414409417482869};
	
	//leading_term = pow(eta,1.0+k)/(1.0+k)*hyp2f1(-0.5,1.0+k,2.0+k,-0.5*eta*theta);
	leading_term = pow(eta,1.0+k)/(1.0+k)*sommerfeld_leading_term(k,-0.5*eta*theta);
	
	if(SERIES_TERMS_MAX==0) return leading_term;
	if(SERIES_TERMS_MAX==1) return leading_term 
	+ M_PI*M_PI/6.0*(pow(eta,k)*theta/4.0/sqrt(1.0+theta*eta/2.0)+k*pow(eta,k-1.0)*sqrt(1.0+theta*eta/2.0));

	
	
	for(i=1;i<=SERIES_TERMS_MAX;i++)
	{
	  derivative = 0.0;
	  for(j=0;j<=2*i-1;j++)
		derivative = derivative + binom(2*i-1,j)*tgamma(1.5)*tgamma(1.0+k)/tgamma(1.5-j)/tgamma(2.0+k-2.0*i+j)
	                             *pow(0.5*theta,j)*pow(1.0+0.5*theta*eta,0.5-j)*pow(eta,1.0-2.0*i+j+k);
	
	  if(i>5)
	   asymptotic_terms = asymptotic_terms + derivative*dirichlet_eta(2.0*i,DBL_EPSILON,64);
      else
	   asymptotic_terms = asymptotic_terms + derivative*etaTBL[2*i];
	  
	}  
	
	
		return leading_term + 2.0*asymptotic_terms;
	
}

double Ffermi_series_neg(const double k, const double eta, const double theta, const double precision, const int SERIES_TERMS_MAX)
{

  double sum_old=0.0, sum_new=0.0,x;
  int i=0;
  
  x=2.0/theta;
  
  do
   {
	i++;
	sum_old = sum_new;
  
    sum_new += ( i % 2 == 0 ) ? exp(i*eta)*U(k,i*x) : -exp(i*eta)*U(k,i*x);
  
  }
   while( ( (precision>0) ? fabs(sum_old-sum_new) >= precision*sum_new: sum_old!=sum_new ) && i<SERIES_TERMS_MAX );
     
  
   return -sum_new*tgamma(1.0+k)*pow(x,1.0+k);

}

double Ffermi_series_sqrt_a(const double k, const double eta, const double theta, const double precision, const int SERIES_TERMS_MAX)
{
#include "factorial.h"
  int i;
  double sum_old=0.0, sum_new=0.0;
  
  //for(i=0;i<SERIES_TERMS_MAX;i++) sum = sum + Ffermi_complete(k+i,eta)*pow(0.5*theta,i)*binom12[i];
  i=0;
  do
   {
	sum_old = sum_new;
    sum_new += Ffermi_complete(k+i,eta)*pow(0.5*theta,i)*binom12[i];
    i++; 
   }
  while( ( (precision>0) ? fabs(sum_old-sum_new) >= precision*sum_new: sum_old!=sum_new ) && i<SERIES_TERMS_MAX );
  
  //printf("\nDBG:\t%e\t%d\n",theta,i);
   
  return sum_new;

}


double Ffermi_series_sqrt_b(const double k, const double eta, const double theta, const double precision, const int SERIES_TERMS_MAX)
{

#include "factorial.h"
  int i;
  double sum=0.0;
  
  for(i=0;( (i<SERIES_TERMS_MAX) && (k+0.5-i>-1.0) );i++) sum = sum + Ffermi_complete(k-i+0.5,eta)*pow(0.5*theta,0.5-i)*binom12[i];

  //printf("\nDBG:\t%e\t%d\n",theta,i);

  return sum;
  
}


double Ffermi(const double k, const double eta, const double theta)
{
  #if 0
  if( fmax(1.0+k-log(DBL_EPSILON),eta+1.0+k-log(DBL_EPSILON))*theta<sqrt(DBL_EPSILON) ) 
  {
	  /* special case for tiny theta relative to 1 and eta */
	  printf("SPECIAL\t");
	  return Ffermi_series_sqrt(k, eta, theta);
  }
  #endif
  
  if( eta>56000.0) 
	  return Ffermi_sommerfeld(k, eta, theta, DBL_EPSILON, 32);
  else if( (eta<0.0) && (k>25.0) && (theta>=1.0) )
	 return Ffermi_series_neg(k, eta, theta, DBL_EPSILON, 32);
  else
	 return Ffermi_value(k,eta,theta,PRECISION, MAX_REFINE);
}

long double Ffermi_long(const long double k, const long double eta, const long double theta)
{
  return Ffermi_value_long(k,eta,theta,PRECISION, MAX_REFINE);
}


