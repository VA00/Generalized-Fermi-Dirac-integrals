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
#include "factorial.h" //pre-calculated factorials

// Pochhammer[a,n]
double pochhammer(double a, int n)
{
   double prod=1.0;
   int j;
   
   for(j=0;j<n;j++) prod=prod*(a+j);
   
   return prod;
	
}

// FactorialPower[a,n] == Pochhammer[a - n + 1, n]
double factorial_power(double a, int n)
{
   double prod=1.0;
   int j;
   
   for(j=0;j<n;j++) prod=prod*(a-j);
   
   return prod;
	
}
		
/*
Exponentiation by squaring
*/


double power_squaring(double x, int n)
{
  double p=1.0;
  
  while(n!=0)
  {
    if(n&1) p = p*x;
    x = x*x;
    n=n>>1;
  }
  
  return p;
}

/* Std. hypergeometric series, converging rapidly for |z|<0.5, and slowly for |z|<0.8 */

long double hyp2f1_series_long(long double a, long double b, long double c, long double z, long double precision, int SERIES_TERMS_MAX)
{
	long double sum_old=0.0L, sum_new=0.0L;
	int i=0;

	//for(i=0;i<=num;i++)	series = series + tgamma(a+i)*tgamma(b+i)/tgamma(c+i)/tgamma(1.0+i)*pow(z,i);

    if(z==0.0L) return 1.0L;
    
	do
	{
		sum_old=sum_new;
		sum_new = sum_old + tgammal(a+i)*tgammal(b+i)/tgammal(c+i)/tgammal(1.0+i)*powl(z,i);
		i++;
	}
    while ( ( (precision>0.0L) ? fabsl(sum_old-sum_new)>=precision*sum_new : sum_old!=sum_new )  && i<SERIES_TERMS_MAX );

	return sum_new*tgammal(c)/tgammal(a)/tgammal(b);
}



double hyp2f1_series(double a, double b, double c, double z, double precision, int SERIES_TERMS_MAX)
{
	double sum_old=0.0, sum_new=0.0;
	int i=0;

    if(z==0.0) return 1.0;
	if(b==c)  return pow(1.0-z,-a); 
	if(a==c)  return pow(1.0-z,-b); 
	
    
	do
	{
		sum_old = sum_new;
		sum_new = sum_old + pochhammer(a,i)*pochhammer(b,i)/pochhammer(c,i)/tgamma(1.0+i)*pow(z,i);
		i++;
	}
    while ( ( (precision>0.0) ? fabs(sum_old-sum_new)>=precision*sum_new : sum_old!=sum_new )  && i<SERIES_TERMS_MAX );

	return sum_new;
}




double hyp2f1_reflection0(double a, double b, double c, double z)
{
	return hyp2f1_series(a,b,c,z,0.0,64);
}

double hyp2f1_reflection1(double a, double b, double c, double z)
{
	
   // printf("%e %e %e %e %e %e \n", tgamma(b-a),tgamma(c),tgamma(c-a),tgamma(b),hyp2f1_series(a,c-b,a-b+1.0,1.0/(1.0-z),0.0,64),pow(1.0-z,a));
  //  printf("%e %e %e %e %e %e \n", tgamma(a-b),tgamma(c),tgamma(c-b),tgamma(a),hyp2f1_series(b,c-a,b-a+1.0,1.0/(1.0-z),0.0,64),pow(1.0-z,b));

	
	return  tgamma(b-a)*tgamma(c)/tgamma(c-a)/tgamma(b)*hyp2f1_series(a,c-b,a-b+1.0,1.0/(1.0-z),0.0,64)/pow(1.0-z,a)
	       +tgamma(a-b)*tgamma(c)/tgamma(c-b)/tgamma(a)*hyp2f1_series(b,c-a,b-a+1.0,1.0/(1.0-z),0.0,64)/pow(1.0-z,b)
			;
}

double hyp2f1_reflection2(double a, double b, double c, double z)
{
	return hyp2f1_series(a,c-b,c,z/(z-1.0),0.0,64)/pow(1.0-z,a);
}

double hyp2f1_reflection4(double a, double b, double c, double z)
{
	
	return  tgamma(c-b-a)*tgamma(c)/tgamma(c-a)/tgamma(c-b)*hyp2f1_series(a,b,a+b-c+1.0,1.0-z,0.0,64)
	       +
		    tgamma(a+b-c)*tgamma(c)/tgamma(b)/tgamma(a)*hyp2f1_series(c-a,c-b,c-b-a+1.0,1.0-z,0.0,64)*pow(1.0-z,c-a-b);
}

double hyp2f1(double a, double b, double c, double z)
{
	if(z<-1.62) 
		return hyp2f1_reflection1(a,b,c,z);
	else if (z<-0.5) 
		return hyp2f1_reflection2(a,b,c,z);
	else if (z<0.58)
		return hyp2f1_reflection0(a,b,c,z);
	else
		return hyp2f1_reflection4(a,b,c,z);
}

double hyp2f1_series_fd(double k, double z, double precision, int SERIES_TERMS_MAX)
{
	double sum_old=0.0, sum_new=0.0;
	int i=0;

    if(z==0.0) return 1.0;
	
    sum_new = 1.0-(1.0+k)/(2.0+k)*0.5*z;
	
	i=2;
	
	do
	{
		sum_old = sum_new;
		sum_new = sum_old + (1.0 + k)/(1.0 + k + i)*pow(z,i)*hypfac[i];
						 
		i++;
	}
    while ( ( (precision>0.0) ? fabs(sum_old-sum_new)>=precision*sum_new : sum_old!=sum_new )  && i<SERIES_TERMS_MAX );

	return sum_new;
}



/*
double hyp2f1_series_fd(double k, double z, double precision, int SERIES_TERMS_MAX)
{
	double sum_old=0.0, sum_new=0.0;
	int i=0;

    if(z==0.0) return 1.0;
	
    sum_new = 1.0-(1.0+k)/(2.0+k)*0.5*z;
	
	i=2;
	
	do
	{
		sum_old = sum_new;
		sum_new = sum_old - 0.5*(1.0 + k)/(1.0 + k + i)*factorial[2*i-3]/
                         factorial[i-2]/factorial[i]*pow(z,i)*pow(2.0,3.0-2.0*i);
						 
		i++;
	}
    while ( ( (precision>0.0) ? fabs(sum_old-sum_new)>=precision*sum_new : sum_old!=sum_new )  && i<SERIES_TERMS_MAX );

	return sum_new;
}
*/


/*
double hyp2f1_series_fd(double k, double z, double precision, int SERIES_TERMS_MAX)
{
	double sum_old=0.0, sum_new=0.0;
	int i=0;

    if(z==0.0) return 1.0;
	
 
	do
	{
		sum_old = sum_new;
		sum_new = sum_old - 0.5*(1.0 + k)/(1.0 + k + i)*tgamma(i-0.5)/
                         tgamma(i+1.0)*pow(z,i)/sqrt(M_PI);
						 
		i++;
	}
    while ( ( (precision>0.0) ? fabs(sum_old-sum_new)>=precision*sum_new : sum_old!=sum_new )  && i<SERIES_TERMS_MAX );

	return sum_new;
}
*/


double hyp2f1_reflection1_fd(double k, double z)
{
	double gamma_ratio;
	
	
	if(k<168.0) 
		gamma_ratio = tgamma(2.0 + k)*tgamma(-1.5 - k);
	else
		gamma_ratio = M_PI/sin((-1.5 - k)*M_PI)/sqrt(k); // this looks sufficient approximation
	

	return -pow(-z, -1.0 - k)/2.0/sqrt(M_PI)*gamma_ratio 
	+  2.0*(1.0 + k)
	*hyp2f1_series(-0.5, 1.0, -0.5 - k, 1.0/(1.0 - z),0.0,64)/(3.0 + 2.0*k)*sqrt(1.0 - z);
		
}


double hyp2f1_reflection1_fd_long(double k, double z)
{
	long double gamma_ratio;
	long double kL = (long double) k;
	long double zL = (long double) z;
	
	if(kL<168.0L) 
		gamma_ratio = tgammal(2.0L + kL)*tgammal(-1.5L - kL	);
	else
		gamma_ratio = M_PI/sinl((-1.5L - kL)*M_PI)/sqrtl(kL); // this looks sufficient approximation
	

	return -powl(-zL, -1.0L - kL)/2.0/sqrtl(M_PI)*gamma_ratio 
	+  2.0L*(1.0L + kL)/(3.0L + 2.0L*kL)*sqrtl(1.0L - zL)*( (long double) hyp2f1_series(-0.5, 1.0, -0.5 - k, 1.0/(1.0 - z),0.0,64));
	
	
}

double hyp2f1_reflection2_fd(double k, double z)
{
	
	return hyp2f1_series(-0.5,1.0,2.0+k,z/(z-1.0),0.0,64)*sqrt(1.0-z);

}

double hyp2f1_reflection4_fd(double k, double z)
{
	return  tgamma(0.5)*tgamma(2.0+k)/tgamma(1.5+k)/tgamma(1.0)*hyp2f1_series(0.5,1.0+k,0.5,1.0-z,0.0,64)
	       +
		    tgamma(-0.5)*tgamma(2.0+k)/tgamma(1.0+k)/tgamma(0.5)*hyp2f1_series(1.5+k,1.0,1.5,1.0-z,0.0,64)*pow(1.0-z,0.5);
}

/* unstable recursion ! use up to k<32 */
double recursion_half_frac_k(double k, double x)
{
	double recursion_half_frac_k_tmp;
	
	if(k>-0.5)
	  recursion_half_frac_k_tmp = 2.0*(1.0+k)/(3.0+2.0*k)/x*( recursion_half_frac_k(k-1.0,x) - pow(1.0-x,1.5) );
    else
	  if(x<0.0) 
			recursion_half_frac_k_tmp =  sqrt(1.0 - x)/2.0 + asinh(sqrt(-x))/(2.0*sqrt(-x));
	  else if (x>0.0 && x<=1.0)
			recursion_half_frac_k_tmp = sqrt(1.0 - x)/2.0 + asin(sqrt(x))/(2.0*sqrt(x));
	  else  recursion_half_frac_k_tmp = 1.0;
	  
	return  recursion_half_frac_k_tmp;  
}


double recursion_int_k(double k, double x)
{
	double recursion_int_k_tmp;
	
	if(k>-1.0)
	  recursion_int_k_tmp = 2.0*(1.0+k)/(3.0+2.0*k)/x*( recursion_int_k(k-1.0,x) - pow(1.0-x,1.5) );
    else
	  recursion_int_k_tmp = 1.0;
	  
	return  recursion_int_k_tmp;  
}


double sommerfeld_leading_term(double k, double x)
{
	
	if( (k-floor(k)==0.5) && (k<=64.0) && (x<-0.5) ) 
		return recursion_half_frac_k(k, x);
	else if ( (k-floor(k)==0.0) && (k<=32.0) && (x<-0.5) )
		return recursion_int_k(k, x);
	else
	{
	
	if(x<-1.62) //-1.62
		return hyp2f1_reflection1_fd(k,x);
	else if (x<-0.5) 
		return hyp2f1_reflection2_fd(k,x);
		//return hyp2f1_reflection2(-0.5,1.0+k,2.0+k,x);
	else if (x<0.58)
		//return hyp2f1_series(-0.5,1.0+k,2.0+k,x,0.0,64);
		return hyp2f1_series_fd(k,x,0.0,64);
	else
		return hyp2f1_reflection4(-0.5,1.0+k,2.0+k,x);
	}
}


void sommerfeld_leading_term_derivatives(double k, double z, double result[4])
{
	double s = sqrt(1.0-z);

	result[0] = sommerfeld_leading_term(k,z);
	//The first d/dz derivative of the 2F1(-0.5, 1+k,2+,z) hypergeometric function
    result[1] = (1.0+k)/z*(s - result[0]);
    //The second d2/dz2 derivative
	result[2] = (s + k*s + result[1]*(4.0 + k*(2.0 - 2.0*z) - 4.0*z))/(z*(-2.0 + 2.0*z));
    ////The third d3/dz3 derivative
	result[3] = (-s - k*s + result[2]*(-12.0 + (24.0 - 12.0*z)*z + k*(-4.0 + (8.0 - 4.0*z)*z)))/(z*(4.0 + z*(-8.0 + 4.0*z)));
	

}

/*
double sommerfeld_leading_term(double k, double z)
{
	if(z<-1.62) 
		return hyp2f1_reflection1(-0.5,1.0+k,2.0+k,z);
	else if (z<-0.5) 
		return hyp2f1_reflection2(-0.5,1.0+k,2.0+k,z);
	else if (z<0.58)
		return hyp2f1_reflection0(-0.5,1.0+k,2.0+k,z);
	else
		return hyp2f1_reflection4(-0.5,1.0+k,2.0+k,z);
}
*/


double sommerfeld_leading_term_int(double k, double x)
{
	
	double t[64] = {0.00034747913211393027154718782718184,0.0018299416140223603265377496618004,0.0044933142616278396303088082783484,0.0083318730576870215343503489215844,0.013336586105044518129073246323864,0.019495600173973140540692939051421,0.026794312570798591968759254326368,0.035215413934030212089254922720387,0.044738931460748597121809665995835,0.055342277002442947073297980863574,0.067000300922953590119608307464921,0.079685351873709818624154227652063,0.093367342438601220129038330956848,0.10801382052832929619488973739312,0.12359004636973405169406811255715,0.14005907491419458657552989108403,0.15738184347288337871822081448431,0.1755172643726713300711193840033,0.19442232241380337487557351449073,0.21405217689868298285806094167041,0.23436026799005272717099304822777,0.25529842714647352126073684648904,0.27681699137326795600752614262054,0.29886492101800419815211661436992,0.3213899208311659420247786924769,0.34433856400489452192124365071992,0.36765641889561629181301791374499,0.39128817812999645792517562550559,0.41517778978800359098134318512587,0.43926859035193972276481176825388,0.46350343910610048027522852902983,0.48782485366828778374552207857314,0.51217514633171221625447792142686,0.53649656089389951972477147097017,0.56073140964806027723518823174612,0.58482221021199640901865681487413,0.60871182187000354207482437449441,0.63234358110438370818698208625501,0.65566143599510547807875634928008,0.6786100791688340579752213075231,0.70113507898199580184788338563008,0.72318300862673204399247385737946,0.74470157285352647873926315351096,0.76563973200994727282900695177223,0.78594782310131701714193905832959,0.80557767758619662512442648550927,0.8244827356273286699288806159967,0.84261815652711662128177918551569,0.85994092508580541342447010891597,0.87640995363026594830593188744285,0.89198617947167070380511026260688,0.90663265756139877987096166904315,0.92031464812629018137584577234794,0.93299969907704640988039169253508,0.94465772299755705292670201913643,0.95526106853925140287819033400416,0.96478458606596978791074507727961,0.97320568742920140803124074567363,0.98050439982602685945930706094858,0.98666341389495548187092675367614,0.99166812694231297846564965107842,0.99550668573837216036969119172165,0.9981700583859776396734622503382,0.99965252086788606972845281217282};
	
	double w[64] = {0.000891640360848216473648039572486,0.0020735166302812338176437678642757,0.0032522289844891814280586801999906,0.0044233799131819738615154573298653,0.005584069730065564409295246509604,0.006731523948359321299030383342978,0.00786301523801235966098299764877,0.00897585788784867154252265100056,0.010067411576765104686170158364272,0.011135086904191627079649165192078,0.012176351284355436669088775204534,0.013188734857527329335845896312613,0.014169836307129741613755652600119,0.015117328536201239433987029909774,0.016028964177425776792733752173949,0.016902580918570804695782741055363,0.017736106628441191905346573357623,0.018527564270120023020207550904792,0.019275076589307814564481248473404,0.019976870566360170693328463064168,0.020631281621311764305078148736819,0.021236757561826794503669883954409,0.021791862264661726688413930486869,0.022295279081878281530067355015472,0.022745813963709072239885498485635,0.023142398290657208647976624616131,0.023484091408105008662663142877291,0.023770082857415154331141103472112,0.023999694298229153864063089935673,0.024172381117401478584884763579009,0.024287733720751713467399533391989,0.024345478504569860191682695367375,0.024345478504569860191682695367375,0.024287733720751713467399533391989,0.024172381117401478584884763579009,0.023999694298229153864063089935673,0.023770082857415154331141103472112,0.023484091408105008662663142877291,0.023142398290657208647976624616131,0.022745813963709072239885498485635,0.022295279081878281530067355015472,0.021791862264661726688413930486869,0.021236757561826794503669883954409,0.020631281621311764305078148736819,0.019976870566360170693328463064168,0.019275076589307814564481248473404,0.018527564270120023020207550904792,0.017736106628441191905346573357623,0.016902580918570804695782741055363,0.016028964177425776792733752173949,0.015117328536201239433987029909774,0.014169836307129741613755652600119,0.013188734857527329335845896312613,0.012176351284355436669088775204534,0.011135086904191627079649165192078,0.010067411576765104686170158364272,0.00897585788784867154252265100056,0.00786301523801235966098299764877,0.006731523948359321299030383342978,0.005584069730065564409295246509604,0.0044233799131819738615154573298653,0.0032522289844891814280586801999906,0.0020735166302812338176437678642757,0.000891640360848216473648039572486};
	
	double sum=0.0;
	int i;
	
	for(i=0;i<64;i++)
		sum = sum + w[i]*sqrt(1.0 - x*t[i]*t[i])*pow(t[i],2.0*k+1.0); 
	
	return sum*(2.0+2.0*k);
	
}

double hyp1f1_series(double a, double b, double z, double precision, int SERIES_TERMS_MAX)
{
	#include "factorial.h"
	double sum_old=0.0, sum_new=0.0, prod1=1.0, prod2=1.0;
	int i=0,j;

    if(z==0.0) return 1.0;
    
	do
	{
		sum_old=sum_new;
		prod1=1.0;
		for(j=0;j<i;j++) prod1=prod1*(a+j);
		prod2=1.0;
		for(j=0;j<i;j++) prod2=prod2*(b+j);
		sum_new = sum_old + prod1/prod2/factorial[i]*pow(z,i);
		i++;
	}
    while ( ( (precision>0.0) ? fabs(sum_old-sum_new)>=precision*sum_new : sum_old!=sum_new )  && i<SERIES_TERMS_MAX );

	return sum_new;
}

double hyp1f1(double a, double b, double z)
{
	return hyp1f1_series(a,b,z,0.0,64);
}

double hypU_asympt(double a, double b, double z, const double precision, const int SERIES_TERMS_MAX)
{
	#include "factorial.h"
	double sum_old,sum_new, prod1=1.0, prod2=1.0;
	double term_old, term_new=1.0;
	int i,j;
	

	i=0;
	sum_new=0.0;
	do
	{
		term_old=term_new;
		sum_old=sum_new;
		prod1=1.0;
		for(j=0;j<i;j++) prod1=prod1*(a+j);
		prod2=1.0;
		for(j=0;j<i;j++) prod2=prod2*(a-b+1.0+j);
	    //term_new = ( i % 2 == 0 ) ?  pow(z,-i)/factorial[i]*prod1*prod2 : -pow(z,-i)/factorial[i]*prod1*prod2;
	    //term_new = ( i % 2 == 0 ) ?  pow(z,-i-a)/factorial[i]*pochhammer(a,i)*pochhammer(a-b+1.0,i) : -pow(z,-i-a)/factorial[i]*pochhammer(a,i)*pochhammer(a-b+1.0,i);
	    term_new = ( i % 2 == 0 ) ?  1.0/power_squaring(z,i)/factorial[i]*pochhammer(a,i)*pochhammer(a-b+1.0,i) : -1.0/power_squaring(z,i)/factorial[i]*pochhammer(a,i)*pochhammer(a-b+1.0,i);
		
				//if(round(z)==z) printf("DBG: %d\t%.18e\n", i, term_new);
				
		
		sum_new += term_new;
        if(fabs(term_new) > 1.0*fabs(term_old))  break; //primitive divergence detection
		i++;
	}
	while ( ( (precision>0.0) ? fabs(sum_old-sum_new)>=precision*sum_new : sum_old!=sum_new )  && i<SERIES_TERMS_MAX );
 
/*
	sum_new=0.0;
	do
	{
	    term_new = ( i % 2 == 0 ) ?  1.0/power_squaring(z,i)/factorial[i]*pochhammer(a,i)*pochhammer(a-b+1.0,i) : -1.0/power_squaring(z,i)/factorial[i]*pochhammer(a,i)*pochhammer(a-b+1.0,i);
		 printf("DBG: %d\t%.18e\n", i, term_new);
		sum_new += term_new;
		i--;
	}
	while (i>=0);
*/
 
	return sum_new*pow(z,-a);
	//return sum_old;
}

double hypU_series(const double a, const double b, const double z, const double precision, const int SERIES_TERMS_MAX)
{
    
   return tgamma(1.0 - b)/tgamma(a - b + 1.0)*hyp1f1_series(a,                 b, z, precision, SERIES_TERMS_MAX) 
              + 
          tgamma(b - 1.0)/tgamma(a)          *hyp1f1_series(a - b + 1.0, 2.0 - b, z, precision, SERIES_TERMS_MAX)*pow(z,1.0 - b);
			
}



double hypU(double a, double b, double z, const double precision, const int SERIES_TERMS_MAX)
{
	if(fabs(z)<16.0)
		return hypU_series(a,b,z,precision, SERIES_TERMS_MAX);
	else
		return hypU_asympt(a,b,z,precision, SERIES_TERMS_MAX);
   
}

/* HypergeometricU[1 + k, 5/2 + k, x] */
double U_half_frac(double k, double x)
{
	
	if(k==-0.5) return exp(0.5*x)/2.0/sqrt(M_PI)*( BesselK(0.0,0.5*x) + BesselK(1.0,0.5*x) );
	if(k==0.5) return exp(0.5*x)/x/sqrt(M_PI)*BesselK(1.0,0.5*x);
		
	return ( U_half_frac(k-2.0, x) + (0.5+k-x)*U_half_frac(k-1.0, x) )/k/x;
}



/* This is particular case of the hypergeometric U(1+k, 5/2+k,2*j/theta) required 
for Fermi-Dirac integrals */
double U(double k, double x)
{
	double eps = 4.0*DBL_EPSILON;
	
	if( (round(2.0*k)==(2.0*k)) && (round(k)!=k) && (x<16.0) ) // Detect half-fractional k (!)
	  return U_half_frac(k,x); //half-fractional case
	  //return 0.0; //half-fractional case
    else 
	  return hypU(1.0+k, 2.5+k, x, eps, 64);
	
}


