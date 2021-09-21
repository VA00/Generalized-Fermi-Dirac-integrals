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
 

		SECTION FOR RELATIVISTIC Fermi-Dirac integrals (G-function)


*/




double integrandG(const double t, const double n, const double alpha, const double beta)
{
  
  double x, dx, integrand, result, factor;
  
  //if(t<-6.5) return( 0.0 );
   
  x = exp(t - exp(-t) ); /* Masatake Mori, eq. (4.17) */
  
  
  dx = x*(1.0 + exp(-t) ); /* dx/dt */
  
  if( (x+alpha-beta) < -log(DBL_EPSILON) )
  {
    factor = 1.0/(1.0+exp(x+alpha-beta) );
	integrand = pow(1.0+x/alpha,2.0*n+1.0)*sqrt(x*x+2.0*x*alpha)*factor;
  }
  else
  {
    //factor = exp(beta-alpha-x);
	integrand = exp((2.0*n+1.0)*log(1.0+x/alpha)+beta-alpha-x)*sqrt(x*x+2.0*x*alpha);
  }
  
  
  

  result = integrand*dx/alpha/alpha;
  
#if DEBUG
   printf("t=%lf x=%e dx=%e fac=%e int=%e res=%e\n",t, x, dx, factor,integrand,result);
#endif  
  
  return result;
  
  
}


/* Implementation of the sample algorithm from Mathematica 8 manual */

/*
IRuleEstimate[F_, h_, oldSum_: None] :=
  Block[{step, i, temp, s1, s2},
   If[oldSum === None, step = 1, step = 2];
   i = 1;
   s1 = FixedPoint[(temp = F[i*h]; i += step; N[temp + #1]) &, 
     0];
   i = -1;
   s2 = FixedPoint[(temp = F[i*h]; i -= step; N[temp + #1]) &, 
     0];
   If[oldSum === None, h (s1 + F[0] + s2), h (s1 + s2) + oldSum/2]
   ];
*/

double Gfermi_estimate(double h, double last_result, double n, double alpha, double beta)
{
  
  int step,i;
  double sum_Left_old, sum_Right_old;
  double sum_Left_new, sum_Right_new;
  double old_result, new_result;
  
  
  if(last_result<0.0) /* Negative value means first iteration*/
  {
    step=1;
    old_result = 2.0*h*integrandG(0.0, n ,alpha, beta);
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
    sum_Right_new = sum_Right_old + integrandG(h*i,n,alpha,beta);
    i = i + step;
  }
  while  (sum_Right_old!=sum_Right_new);
  

  /* integral for -Infinity < t <0  */
  
  sum_Left_old = 0.0;
  sum_Left_new = 0.0;
  
  
  i=-1;
  do
  {
    sum_Left_old = sum_Left_new;
    sum_Left_new = sum_Left_old + integrandG(h*i,n,alpha,beta);
    i = i - step;
  }
  while  (sum_Left_old!=sum_Left_new);
  
  /* fail - safe integration in reversed order
     required if integrandG is 0.0 for t surrounding t=0
  */
  if( (sum_Left_new<=DBL_MIN) && (sum_Right_new<=DBL_MIN) ) 
    {  i=(int) 4.0/h-1.0;
	  do
	  {
		sum_Right_new += integrandG(h*i,n,alpha,beta);
		sum_Left_new += integrandG(-h*i,n,alpha,beta);
		i = i - step;
		#if DEBUG
		printf("%d\n",i);
		#endif
	  }
	  while( i>=1 ) ;
 

	};
   
   
   
  new_result = h*(sum_Left_new  + sum_Right_new) + 0.5*old_result;
  
  return new_result;
}


double Gfermi(const double n, const double alpha, const double beta)
{
  
  double old=0.0, new=0.0, h=0.5;
  

#if DEBUG
    printf("Entering GFermi\n");
#endif  
  
  old = 0.0;
  new = Gfermi_estimate(h, -1.0, n, alpha, beta);
#if DEBUG
    printf("BEFORE WHILE: old = %e, new=%e\n",old,new);
#endif     

  while( fabs(old-new)>PRECISION*fabs(new) && h>pow(2.0,-MAX_REFINE))
  {
    old=new;
    h=0.5*h;
    new = Gfermi_estimate(h, old, n, alpha, beta);
#if DEBUG
    printf("old = %e, new=%e\n",old,new);
#endif    
  }
  
#if DEBUG
    printf("new=%e alpha=%lf n=%lf\n",new, alpha, 2.0*n+3.0);
#endif     
  return new;  
    
}


double Gp(const double n, const double alpha, const double beta)
{
  return Gfermi(n,alpha,-beta);
}


double Gm(const double n, const double alpha, const double beta)
{
  return Gfermi(n,alpha,beta);
}

