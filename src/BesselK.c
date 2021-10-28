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

/*
 

		SECTION FOR Bessel I0 and I1

https://www.advanpix.com/2015/11/11/rational-approximations-for-the-modified-bessel-function-of-the-first-kind-i0-computations-double-precision/

*/


double BesselI0(const double z)
{ 
  const double A = 7.75; //Splitting point
  /* I0(x) */
  /* https://www.advanpix.com/wp-content/uploads/2015/11/I0_coefficients.txt */
/* [0,7.75): */
  const double P1[17] = {
  1.0000000000000000000000801e+00,
  2.4999999999999999999629693e-01,
  2.7777777777777777805664954e-02,
  1.7361111111111110294015271e-03,
  6.9444444444444568581891535e-05,
  1.9290123456788994104574754e-06,
  3.9367598891475388547279760e-08,
  6.1511873265092916275099070e-10,
  7.5940584360755226536109511e-12,
  7.5940582595094190098755663e-14,
  6.2760839879536225394314453e-16,
  4.3583591008893599099577755e-18,
  2.5791926805873898803749321e-20,
  1.3141332422663039834197910e-22,
  5.9203280572170548134753422e-25,
  2.0732014503197852176921968e-27,
  1.1497640034400735733456400e-29};
  
  /* [7.75,Inf): */
  const double P2[23] = {
   3.9894228040143265335649948e-01,
   4.9867785050353992900698488e-02,
   2.8050628884163787533196746e-02,
   2.9219501690198775910219311e-02,
   4.4718622769244715693031735e-02,
   9.4085204199017869159183831e-02,
  -1.0699095472110916094973951e-01,
   2.2725199603010833194037016e+01,
  -1.0026890180180668595066918e+03,
   3.1275740782277570164423916e+04,
  -5.9355022509673600842060002e+05,
   2.6092888649549172879282592e+06,
   2.3518420447411254516178388e+08,
  -8.9270060370015930749184222e+09,
   1.8592340458074104721496236e+11,
  -2.6632742974569782078420204e+12,
   2.7752144774934763122129261e+13,
  -2.1323049786724612220362154e+14,
   1.1989242681178569338129044e+15,
  -4.8049082153027457378879746e+15,
   1.3012646806421079076251950e+16,
  -2.1363029690365351606041265e+16,
   1.6069467093441596329340754e+16};
   
   double x, y, P22, P16;
   int i;


  x = fabs(z); //I0 is odd-functrion
  
  if( x>=A )
  {
    y=1.0/x;

    P22 = y*P2[22];
    for(i=21;i>0;i--) P22 = y*(P2[i]+P22);
    P22 = P22+P2[0];

    return exp(x)/sqrt(x)*P22;
  }
  else
  {
    y= 0.25*x*x;

    P16 = y*P1[16];
    for(i=15;i>0;i--) P16 = y*(P1[i]+P16);
    P16 = P16+P1[0];

    return y*P16+1.0;
  }

}


double BesselI1(const double z)
{
  const double A = 7.75; //Splitting point
/* I1(x) https://www.advanpix.com/wp-content/uploads/2015/11/I1/I1_coefficients.txt */

//[0,7.75):
const double P1[14] = {
8.3333333333333333311567967e-02,
6.9444444444444450632369337e-03,
3.4722222222221933634809047e-04,
1.1574074074079326676719210e-05,
2.7557319223490964726712181e-07,
4.9209498642000488498034902e-09,
6.8346524852208360284643288e-11,
7.5940608652484019265663823e-13,
6.9036483611746228414130722e-15,
5.2305536429160706017626195e-17,
3.3486060280590464196327437e-19,
1.8645262719811753663834319e-21,
7.9611250107842314599760659e-24,
5.3251032089995165438568695e-26};

//[7.75,Inf):
const double P2[23] = {
 3.9894228040143270388374079e-01,
-1.4960335515072058522575487e-01,
-4.6751048269476797374239762e-02,
-4.0907267094886972971863462e-02,
-5.7501487840859800117669379e-02,
-1.1428156617865937773864845e-01,
 6.7988447242260666801129937e-02,
-2.2694203870019250176636896e+01,
 9.7548286270114208672947525e+02,
-2.9286459257939415083570152e+04,
 4.9934855620495985742805154e+05,
 5.7682364160056137069002930e+05,
-3.1576840778898356890175020e+08,
 1.0484906321376589515223174e+10,
-2.0918193917759394367113655e+11,
 2.9320804098307168426392082e+12,
-3.0147278411132255281401004e+13,
 2.2950466603697814797615042e+14,
-1.2816007548999035598180100e+15,
 5.1086996139908353110844064e+15,
-1.3774917783425787550429723e+16,
 2.2531580094188348024267027e+16,
-1.6895178303473738478791245e+16};

  double x,y,P22,P13;
  int i;

  x = fabs(z);

  if( x>=A )
  {
    y=1.0/x;

    P22 = y*P2[22];
    for(i=21;i>0;i--) P22 = y*(P2[i]+P22);
    P22 = P22+P2[0];

    return exp(x)/sqrt(x)*P22*x/z;
  }
  else
  {
    y= 0.25*x*x;

    P13 = y*P1[12];
    for(i=11;i>0;i--) P13 = y*(P1[i]+P13);
    P13 = P13+P1[0];

    return 0.5*z*(1.0+y+y*y*P13);
  }


}
/*
 



		SECTION FOR Bessel K0 and K1

https://www.advanpix.com/2015/11/25/rational-approximations-for-the-modified-bessel-function-of-the-second-kind-k0-for-computations-with-double-precision/


*/


double BesselK0(const double x)
{
/* https://www.advanpix.com/wp-content/uploads/2015/11/K0/K0_coefficients.txt */
/* K0(x) */

/* [0,1): */
  const double P1[8] = {
  1.1593151565841244842077226e-01,
  2.7898287891460317300886539e-01,
  2.5248929932161220559969776e-02,
  8.4603509072136578707676406e-04,
  1.4914719243067801775856150e-05,
  1.6271068931224552553548933e-07,
  1.2082660336282566759313543e-09,
  6.6117104672254184399933971e-12};
  
  const double P2[7] = {
  1.0000000000000000044974165e+00,
  2.4999999999999822316775454e-01,
  2.7777777777892149148858521e-02,
  1.7361111083544590676709592e-03,
  6.9444476047072424198677755e-05,
  1.9288265756466775034067979e-06,
  3.9908220583262192851839992e-08};
  
  /* [1,Inf): */
  const double P3[22] = {
   1.0694678222191263215918328e-01,
   9.0753360415683846760792445e-01,
   1.7215172959695072045669045e+00,
  -1.7172089076875257095489749e-01,
   7.3154750356991229825958019e-02,
  -5.4975286232097852780866385e-02,
   5.7217703802970844746230694e-02,
  -7.2884177844363453190380429e-02,
   1.0443967655783544973080767e-01,
  -1.5741597553317349976818516e-01,
   2.3582486699296814538802637e-01,
  -3.3484166783257765115562496e-01,
   4.3328524890855568555069622e-01,
  -4.9470375304462431447923425e-01,
   4.8474122247422388055091847e-01,
  -3.9725799556374477699937953e-01,
   2.6507653322930767914034592e-01,
  -1.3951265948137254924254912e-01,
   5.5500667358490463548729700e-02,
  -1.5636955694760495736676521e-02,
   2.7741514506299244078981715e-03,
  -2.3261089001545715929104236e-04};
  
  const double Q3[3] = {
  8.5331186362410449871043129e-02,
  7.3477344946182065340442326e-01,
  1.4594189037511445958046540e+00};
  
  double y, yy, K0;
  double P21,Q2,P7;
  int i;



  if(x>=1.0) 
   {
    y  = 1.0/x;
    yy = 1.0;
    
    
    /* Polynomial evaluation */
    /*
    P21=0.0;
    for(i=0;i<22;i++){ 
      P21 = P21 + P3[i]*yy;
      yy=yy*y;
    }
    */
    
    /* Horner form */
    
    P21 = y*P3[21];
    for(i=20;i>0;i--) P21 = y*(P3[i]+P21);
    P21 = P21+P3[0];
    
    
    
    //Q2 = Q3[0] + Q3[1]*y + Q3[2]*y*y;
    Q2 = Q3[0] + y*(Q3[1] + Q3[2]*y); //HornerForm
  
    return exp(-x)*sqrt(y)*P21/Q2; 
   }
  else
   {
    y=x*x;
    
    P7 = y*P1[7];
    for(i=6;i>0;i--) P7 = y*(P1[i]+P7);
    P7 = P7+ P1[0];  
   
    return P7-log(x)*BesselI0(x);
  
   }

}


double BesselK1(const double x)
{
/* https://www.advanpix.com/wp-content/uploads/2015/11/K1/K1_coefficients.txt */
/* K1(x) */

/* [0,1): - K1 */
const double P1[9] = {
-3.0796575782920622440538935e-01,
-8.5370719728650778045782736e-02,
-4.6421827664715603298154971e-03,
-1.1253607036630425931072996e-04,
-1.5592887702110907110292728e-06,
-1.4030163679125934402498239e-08,
-8.8718998640336832196558868e-11,
-4.1614323580221539328960335e-13,
-1.5261293392975541707230366e-15};

/* [0,1): - I1 */
const double P2[6] = {
8.3333333333333325191635191e-02,
6.9444444444467956461838830e-03,
3.4722222211230452695165215e-04,
1.1574075952009842696580084e-05,
2.7555870002088181016676934e-07,
4.9724386164128529514040614e-09};

/* [1,Inf): */
const double P3[23] = {
 1.0234817795732426171122752e-01,
 9.4576473594736724815742878e-01,
 2.1876721356881381470401990e+00,
 6.0143447861316538915034873e-01,
-1.3961391456741388991743381e-01,
 8.8229427272346799004782764e-02,
-8.5494054051512748665954180e-02,
 1.0617946033429943924055318e-01,
-1.5284482951051872048173726e-01,
 2.3707700686462639842005570e-01,
-3.7345723872158017497895685e-01,
 5.6874783855986054797640277e-01,
-8.0418742944483208700659463e-01,
 1.0215105768084562101457969e+00,
-1.1342221242815914077805587e+00,
 1.0746932686976675016706662e+00,
-8.4904532475797772009120500e-01,
 5.4542251056566299656460363e-01,
-2.7630896752209862007904214e-01,
 1.0585982409547307546052147e-01,
-2.8751691985417886721803220e-02,
 4.9233441525877381700355793e-03,
-3.9900679319457222207987456e-04};

const double Q3[3] = {
8.1662031018453173425764707e-02,
7.2398781933228355889996920e-01,
1.4835841581744134589980018e+00};
  
  double y, yy;
  double P22,Q2,P8;
  int i;



  if(x>=1.0) 
   {
    y  = 1.0/x;
    yy = 1.0;
    
    
    /* Polynomial evaluation */
    /*
    P21=0.0;
    for(i=0;i<22;i++){ 
      P21 = P21 + P3[i]*yy;
      yy=yy*y;
    }
    */
    
    /* Horner form */
    
    P22 = y*P3[22];
    for(i=21;i>0;i--) P22 = y*(P3[i]+P22);
    P22 = P22+P3[0];
    
    
    
    //Q2 = Q3[0] + Q3[1]*y + Q3[2]*y*y;
    Q2 = Q3[0] + y*(Q3[1] + Q3[2]*y); //HornerForm
  
    return exp(-x)*sqrt(y)*P22/Q2; 
   }
  else
   {
    y=x*x;
    
    P8 = y*P1[8];
    for(i=7;i>0;i--) P8 = y*(P1[i]+P8);
    P8 = P8+ P1[0];  
   
    return x*P8+log(x)*BesselI1(x)+1.0/x;
  
   }

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
	
    if(nu==0.0) return BesselK0(x);
	return BesselK_dbl_exp(nu,x,0.0,3);
	
}
