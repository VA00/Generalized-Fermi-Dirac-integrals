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
#define DENORMALS 1


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
 

		SECTION FOR Bessel K


*/


double integrandK(const double t, const double nu, const double x)
{
  return 0.5*exp(-x*cosh(t)-nu*t);
}




double BesselK_estimate(double h, double last_result, double nu, double x)
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
    old_result = 2.0*h*integrandK(0.0, nu, x);
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
    sum_Right_new = sum_Right_old + h*integrandK(h*i, nu, x);
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
    sum_Left_new = sum_Left_old + h*integrandK(h*i, nu, x); //This way numeric overflow is postponed
    i = i - step;
  }
  while  (sum_Left_old<sum_Left_new);
  
  
  //new_result = h*(sum_Left_new  + sum_Right_new) + 0.5*old_result;
  new_result = (sum_Left_new  + sum_Right_new) + 0.5*old_result; //This way numeric overflow is postponed

  return new_result;
}


double BesselK_dbl_exp(const double nu, const double x,  const double precision, const int recursion_limit)
{
  
  double old=0.0, new=0.0, h=0.5;
  
	  
  old = 0.0;
  new = BesselK_estimate(h, -1.0, nu, x);

 
  while( ( (precision>0.0) ? fabs(old-new)>precision*fabs(new) : old!=new ) && h>pow(2.0,-recursion_limit))
  {
    old=new;
    h=0.5*h;
    new = BesselK_estimate(h, old, nu, x);
  }
  
  
  return new;  
    
}



double BesselK(const double nu, const double x)
{
	#ifdef DENORMALS
	if(x>741.36161901346672018165612143975) return 0.0;
	#else
	if(x>705.34286787390965314924968043042) return 0.0;
    #endif
	
	return BesselK_dbl_exp(nu,x,0.0,3);
	
}
