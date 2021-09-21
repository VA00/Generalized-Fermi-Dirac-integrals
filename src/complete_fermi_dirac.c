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
 

		SECTION FOR COMPLETE Fermi-Dirac integrals 


*/


double integrandF_complete(const double t, const double k, const double eta)
{
  
  double x,dx,integrand,result,factor,scale;
  
  //scale = log(M_E + exp(eta));
  //if(eta>k+1.0) scale = eta; else scale = 1.0;
  scale  = 1.0;
  //if( (eta>k+1.0) && (k>0) ) scale = M_E*eta; else scale=1.0;
    
  x = scale*exp(  t - exp(-t) ); /* Masatake Mori, eq. (4.17) */
  
  //dx = x*(1 + exp(-t) ); /* dx/dt */
  dx = 1.0+exp(-t); /* in this case x is adsorbed in integrand, and x^k -> x^(k+1) */
 
  if(x-eta<-log(DBL_EPSILON)) // if using machine precison we are unable to add 1.0 to exp(), then approximation is optimal
    {
	factor = 1.0/(1.0+exp(x-eta) );
	integrand = exp( (k+1.0)*(t - exp(-t) + log(scale)) );
	//integrand = pow(x,k+1.0);
	integrand = integrand*factor;
	}
  else
    {
    //factor = exp(eta-x) adsorbed into exp, to avoid 0*infinity mess 
	integrand = exp((k+1.0)*(t - exp(-t) + log(scale)) + eta - x );
	}


  result = integrand*dx;
  
  
  return result;
  
  
}

double integrandF_complete_internal_long(const double t, const double k, const double eta)
{
  
  long double tL = (long double) t, kL = (long double) k, etaL = (long double) eta;
  long double x,dx,integrand,result,factor;
    
  x = expl(  tL - expl(-tL) ); /* Masatake Mori, eq. (4.17) */
  
  //dx = x*(1 + exp(-t) ); /* dx/dt */
  dx = 1.0L+expl(-tL); /* in this case x is adsorbed in integrand, and x^k -> x^(k+1) */
 
  if(x-etaL<-logl(LDBL_EPSILON)) // if using machine precison we are unable to add 1.0 to exp(), then approximation is optimal
    {
	factor = 1.0L/(1.0L+expl(x-etaL) );
	integrand = expl( (kL+1.0L)*(tL - expl(-tL) ) );
	integrand = integrand*factor;
	}
  else
    {
    //factor = exp(eta-x) adsorbed into exp, to avoid 0*infinity mess 
	integrand = expl((kL+1.0L)*(tL - exp(-tL) ) + etaL - x );
	}

  result = (double) integrand*dx;
  
  
  return result;
  
  
}


long double integrandF_complete_long(const long double tL, const long double kL, const long double etaL)
{
  
  long double x,dx,integrand,result,factor;
    
  x = expl(  tL - expl(-tL) ); /* Masatake Mori, eq. (4.17) */
  
  //dx = x*(1 + exp(-t) ); /* dx/dt */
  dx = 1.0L+expl(-tL); /* in this case x is adsorbed in integrand, and x^k -> x^(k+1) */
 
  if(x-etaL<-logl(LDBL_EPSILON)) // if using machine precison we are unable to add 1.0 to exp(), then approximation is optimal
    {
	factor = 1.0L/(1.0L+expl(x-etaL) );
	integrand = expl( (kL+1.0L)*(tL - expl(-tL) ) );
	integrand = integrand*factor;
	}
  else
    {
    //factor = exp(eta-x) adsorbed into exp, to avoid 0*infinity mess 
	integrand = expl((kL+1.0L)*(tL - expl(-tL) ) + etaL - x );
	}

  result = integrand*dx;
  
  
  return result;
  
  
}

double Ffermi_complete_estimate(double h, double last_result, double k, double eta)
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
    old_result = 2.0*h*integrandF_complete(0.0, k, eta);
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
  do
  {
	sum_Right_old = sum_Right_new;
    #if KAHAN	
	y = h*integrandF_complete(h*i, k, eta) - c;
	t = sum_Right_new + y;
	c = (t-sum_Right_new) - y;
	sum_Right_new = t;
    #else 
    //sum_Right_new = sum_Right_old + integrandF_complete(h*i, k, eta);
    sum_Right_new = sum_Right_old + h*integrandF_complete(h*i, k, eta); //This way numeric overflow is postponed
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
	y = h*integrandF_complete(h*i, k, eta) - c;
	t = sum_Left_new + y;
	c = (t-sum_Left_new) - y;
	sum_Left_new = t;
    #else 
    //sum_Left_new = sum_Left_old + integrandF_complete(h*i, k, eta);
    sum_Left_new = sum_Left_old + h*integrandF_complete(h*i, k, eta); //This way numeric overflow is postponed
    #endif	
    i = i - step;
  }
  while  (sum_Left_old<sum_Left_new);
  
  
  //new_result = h*(sum_Left_new  + sum_Right_new) + 0.5*old_result;
  new_result = (sum_Left_new  + sum_Right_new) + 0.5*old_result; //This way numeric overflow is postponed

  return new_result;
}


long double Ffermi_complete_estimate_long(long double hL, long double last_result_long, long double kL, long double etaL)
{
  
  int step,i;
  long double sum_Left_old, sum_Right_old;
  long double sum_Left_new, sum_Right_new;
  long double old_result, new_result;
  
  
  if(last_result_long<0.0L) /* Negative value means first iteration*/
  {
    step=1;
    old_result = 2.0L*hL*integrandF_complete_long(0.0L, kL, etaL);
  }
  else
  {
    step=2;
    old_result = last_result_long;
  }


  
  /* integral for 0 < t < Infinity  */
  
  sum_Right_old = 0.0L;
  sum_Right_new = 0.0L;
  
  
  i=1;
  do
  {
	sum_Right_old = sum_Right_new;
    sum_Right_new = sum_Right_old + integrandF_complete_long(hL*i, kL, etaL);
	i = i + step;
  }
  while  ( sum_Right_old!=sum_Right_new ); //floating point fixed-point method

  /* integral for -Infinity < t <0  */
  
  sum_Left_old = 0.0L;
  sum_Left_new = 0.0L;
  
  
  i=-1;
  
  //while(hL*i > ( -logl(-logl(DBL_MIN)/(1.0L+kL))) ) i--;
  
  //i=i+1; //odd shift - skip some near-zero points
  //printf("\n\n\n%d\t%Lf\t%Lf\n\n\n",i,hL, i*hL);
  
  do
  {
	sum_Left_old = sum_Left_new;
    sum_Left_new = sum_Left_old + integrandF_complete_long(hL*i, kL, etaL);
	//printf("%.32Lf\t%.32Lf\t%.32Lf\t%d\n",sum_Left_old,sum_Left_new,sum_Left_old-sum_Left_new,i);
    i = i - step;
    //i = i + step;
  }
  while  (sum_Left_old!=sum_Left_new);
  //while  (i<0);
  
  
  new_result = hL*(sum_Left_new  + sum_Right_new) + 0.5L*old_result;

  return new_result;
}

double Ffermi_complete_series_polylog(const double k, const double eta, const double precision, const int SERIES_TERMS_MAX)
{
	double sum_old=0.0,sum_new=0.0;
	int i=0;
	
   do
   {
	i++;
	sum_old = sum_new;
    //sum_new += ( i % 2 == 0 ) ? exp(i*eta)/pow(i,k+1.0) : -exp(i*eta)/pow(i,k+1.0);
    sum_new += ( i % 2 == 0 ) ? exp(i*eta)*pow(i,-k-1.0) : -exp(i*eta)*pow(i,-k-1.0);
   }
   while( ( (precision>0) ? fabs(sum_old-sum_new) >= precision*sum_new: sum_old!=sum_new ) && i<SERIES_TERMS_MAX );
   
   
  
   return -sum_old*tgamma( 1.0+k);
	
}


double Ffermi_complete_series_asympt(const double k, const double eta, const double precision, const int SERIES_TERMS_MAX)
{
	
   const double coeff[32] = {1.0000000000000000000000000000000, \
1.6449340668482264364724151666460, 1.8940656589944918351530064689470, \
1.9711021825948702081968784889699, 1.9924660037052957984545785201656, \
1.9980790151965431312784436913987, 1.9995153702877163817063593574255, \
1.9998783406919594363419083845108, 1.9999695284298122128833655499228, \
1.9999923757392202269593784528232, 1.9999980932231630442301016851155, \
1.9999995232264616450957944098862, 1.9999998807977847892567228062804, \
1.9999999701984639931375753236117, 1.9999999925495068002174550553403, \
1.9999999981373645629079572545591, 1.9999999995343397919029816456306, \
1.9999999998835847980906318477607, 1.9999999999708961828677695279102, \
1.9999999999927240438663375110039, 1.9999999999981810107609577561821, \
1.9999999999995452526673917954700, 1.9999999999998863131643093067280, \
1.9999999999999715782907952543909, 1.9999999999999928945726674721731, \
1.9999999999999982236431633856591, 1.9999999999999995559107904594830, \
1.9999999999999998889776975718783, 1.9999999999999999722444243881926, \
1.9999999999999999930611060965174, 1.9999999999999999982652765240704, \
1.9999999999999999995663191310110}; // Remaining coefficients assumed to be 2.0
   double sum,c,prod=1.0;
   int i=1,j;
   
   sum = pow(eta,k+1.0)/(k+1.0) + M_PI*M_PI/6.0*k*pow(eta,k-1.0);// + (k*k-3.0*k+2.0)*7.0/360.0*M_PI*M_PI*M_PI*M_PI*k*pow(eta,k-3.0);
   
   if( k<0 )
     
    while( i < 6 )
     {
	   i++;
	   prod=1.0;
	   for(j=0;j<=2*i-2;j++) prod=prod*(k-j);
	   sum = sum + coeff[i]*pow(eta,k+1.0-2.0*i)*prod;
	   //printf("DEBUG:i=%d\tk+1=%e\t%lf,%lf,%lf,%lf,%e\n",i,k+1.0,coeff[i],pow(eta,k+1.0-2.0*i),tgamma(k+1.0),tgamma(k+2.0-2.0*i), prod);
	 }
 
   else
   
     while( (i<SERIES_TERMS_MAX) && ( (k+2.0-2.0*i)>0.0) )
     {
	   
	   if(k<=30.0) c = coeff[i]; else c = 2.0;
	   
	   i++;
	   prod=1.0;
	   for(j=0;j<=2*i-2;j++) prod=prod*(k-j);
	   sum = sum + coeff[i]*pow(eta,k+1.0-2.0*i)*prod;//tgamma(k+1.0)/tgamma(k+2.0-2.0*i);
	   
	 }
   

   return sum;
	
}

double Ffermi_complete_series_zeta(const double k, const double eta, const double precision, const int SERIES_TERMS_MAX)
{


	double sum_old=0.0,sum_new=0.0,term;
	int i=0;
	
	do
   {
	
	sum_old = sum_new;
	term = (-1.0+pow(2.0,i-k))*zeta3(1.0+k-i,2.0*DBL_EPSILON,16)*pow(eta,i)/tgamma(i+1.0);
    sum_new += term;
	
	i++;
   }
   while( ( (precision>0.0) ? fabs(sum_old-sum_new) >= precision*sum_new : sum_old!=sum_new ) && i<SERIES_TERMS_MAX );
	
   return -sum_new*tgamma(1.+ k);
	
}

double Ffermi_complete_dbl_exp(const double k, const double eta,  const double precision, const int recursion_limit)
{
  
  double old=0.0, new=0.0, h=0.5;
  
	  
  old = 0.0;
  new = Ffermi_complete_estimate(h, -1.0, k, eta);

 
  while( ( (precision>0.0) ? fabs(old-new)>precision*fabs(new) : old!=new ) && h>pow(2.0,-recursion_limit))
  {
    old=new;
    h=0.5*h;
    new = Ffermi_complete_estimate(h, old, k, eta);
  }
  
  
  return new;  
    
}

double Ffermi_complete(const double k, const double eta)
{
  
  double old=0.0, new=0.0, h=0.5;
  
  if(k<=-1.0) return nan("NaN"); /* not converging for k <= -1 */
  //if( k<=(-1.0+1.0*DBL_EPSILON) ) return 1.0/(1.0+exp(eta))/(1.0+k);
  if(k==0.0) 
  { 
	if(eta>-log(DBL_EPSILON) )
        return eta;		
	else if (eta<log(DBL_EPSILON) )
		return exp(eta);
	else
		return log1p(exp(eta));
  }
  
  if( (fabs(k) < FLT_EPSILON) && (fabs(eta) < FLT_EPSILON) )
	  return M_LN2 + 0.5*eta - 0.5*k*M_LN2*M_LN2; // this is enough for small k,eta
  
  if( (eta==0.0) && (k>32) ) return tgamma(1.0+k);
 
  
  if( (eta>=143.229) && ( log(eta)*(k+1.0) - log(k+1.0)  > log(DBL_MAX) ) ) return INFINITY;
  if( (eta < 143.229) && (eta>log(DBL_MAX)-lgamma(k+1.0)) ) return INFINITY;
  if( (eta<log(nextafter(0.0,DBL_MIN))-lgamma(k+1.0)-log(1.0-exp(eta)/pow(2.0,1.0+k))) && (eta<-700.0) ) // allow denormals below DBL_MIN  
  //if( eta<log(DBL_MIN)-lgamma(k+1.0)-log(1.0-exp(eta)/pow(2.0,1.0+k)) ) // flush denormals to 0.0
  {
      //printf("DBG: CONDITION ZERO\n");	  
	  return 0.0;
  }	  
  
  if( (k>=128.0) && (eta <= 8.0) && (eta>0.0) ) return Ffermi_complete_series_polylog(k, eta, 2.0*DBL_EPSILON, 4);
  if( (k>=20.0) && (eta<= 4.5) && (eta>-0.5)  ) return Ffermi_complete_series_zeta(k, eta, 0.0, 32); //MAX series terms
  if( (k>=20.0) && (eta<=-0.5) && (eta>=-8.0)  ) return Ffermi_complete_series_polylog(k, eta, DBL_EPSILON, 24); //MAX series terms

  if(eta<-8.0) 
  {
	  if( (k<TGAMMA_MAX) && (eta>log(DBL_MIN) ) ) 
	  {
		  //printf("DBG: SERIES POLYLOG\tk=%g\teta=%g\n",k,eta); 
		  //return exp(eta)*(1.0-exp(eta)/pow(2.0,k+1.0))*tgamma(k+1.0); /* PolyLog series is optimal for small eta */
		  return Ffermi_complete_series_polylog(k, eta, DBL_EPSILON, 12);//MAX 8 series terms 
	  }
	  else
		  return (double) expl( ( (long double) eta)  +  lgammal(  ((long double) k)  +1.0L) );
		  
  }
   
  if( (eta>=512.0) ) return Ffermi_complete_series_asympt(k, eta, 0.0, 32.0);
	  
   
  
  
  return Ffermi_complete_dbl_exp(k, eta, 128.0*DBL_EPSILON, 16);  
    
}

long double Ffermi_complete_value_long(const long double kL, const long double etaL,
  const long double precision, const int recursion_limit)
{
  
  long double old=0.0L, new=0.0L, hL=0.5L;
  
  if(kL<=-1.0L) return nanl("NaN"); /* not converging for k <= -1 */
  if( (etaL>=143.229L) && ( logl(etaL)*(kL+1.0L) - logl(kL+1.0L)  > logl(LDBL_MAX) ) ) return INFINITY;
  if( (etaL < 143.229L) && (etaL>logl(LDBL_MAX)-lgammal(kL+1.0L)) ) return INFINITY;
  if( etaL<logl(LDBL_MIN)-lgammal(kL+1.0L) ) return 0.0L;
  

  if(etaL<=-32.0L) 
  {
	  if( (kL<TGAMMA_MAX) && (etaL>logl(LDBL_MIN) ) ) 
		  return expl(etaL)*(1.0L-expl(etaL)/powl(2.0L,kL+1.0L))*tgammal(kL+1.0L); /* PolyLog series is optimal for small eta */
	  else
		  return expl( etaL  +  lgammal( kL+1.0L) );
		  
  }
  
  if( (etaL>=4096.0L) && (kL<8.0L) ) 
	  return powl(etaL,kL+1.0L)/(kL+1.0L)
				+M_PI*M_PI/6.0L*kL*powl(etaL,kL-1.0L)
				+(kL*kL-3.0L*kL+2.0L)*7.0L/360.0L*M_PI*M_PI*M_PI*M_PI*kL*powl(etaL,kL-3.0L); /* PolyLog asympt series is optimal for very large eta */
	  
  old = 0.0L;
  new = Ffermi_complete_estimate_long(hL, -1.0L, kL, etaL);
 
  while( fabsl(old-new)>precision*fabsl(new) && hL>powl(2.0L,-recursion_limit))
  {
    old=new;
    hL=0.5L*hL;
    new = Ffermi_complete_estimate_long(hL, old, kL, etaL);
  }
  
  
  return new;  
    
}

