/*
A. Odrzywolek, andrzej.odrzywolek@uj.edu.pl, 01-06-2020
*/
#include "../fermidirac.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <stdio.h>
#include <float.h>

/*


		SECTION FOR FAST VERSIONS


*/




void fixedFfermi_derivatives(const double k, const double eta, const double theta,
       const double h, double hMIN, double hMAX, 
       double * F, double *dF_deta, double *d2F_deta2,
       double * dF_dtheta, double *d2F_dtheta2, double *d2F_dtheta_deta)
{
  int ii,num;
  double t, x, dx, f, integral=0.0, integral_deta=0.0, integral_deta2=0.0, integral_dtheta=0.0, integral_dtheta2=0.0, integral_deta_dtheta=0.0;
  double xst,dxst,factor, denom, denomi; //AUX vars similar to used by FXT code
  
  num = (int) -hMIN/h;
  hMIN = -num*h;
  
  num = (int) hMAX/h;
  hMAX = num*h;
#if DEBUG  
  printf("fixedFfermi: %d %lf %lf\n", num, hMIN,hMAX);
#endif  
  
  num = (int) (hMAX-hMIN)/h;

  for(ii=0;ii<=num;ii++)
  {
    t = hMIN + ii*h;
    //t = hMAX - ii*h;

    x    = exp(t-exp(-t));
    dx   = 1.0 +exp(-t); // This is NOT D[x,t], D[x,t] = x*dx !
	xst  = 1.0+ 0.5*theta*x;
	dxst = sqrt(xst);
	factor = exp(x-eta);
	denom  = 1.0+factor;
    denomi = 1.0/denom;	
	f      = 1.0/(1.0+factor)*exp(   (k+1.0) * (t - exp(-t))   )*dxst * dx;
    integral             +=   f;
	integral_deta        +=   f*factor*denomi; 
	integral_deta2       +=   f*factor*denomi*denomi*(factor-1.0); 
	integral_dtheta      +=   f*0.25*x/xst;
	integral_dtheta2     +=  -f*0.25*x/xst*x/xst*0.25;
	integral_deta_dtheta +=   f*0.25*x/xst*factor*denomi;
  }
  
  *F                 = h*integral;
  *dF_deta           = h*integral_deta;
  *d2F_deta2         = h*integral_deta2;
  *dF_dtheta         = h*integral_dtheta;
  *d2F_dtheta2       = h*integral_dtheta2;
  *d2F_dtheta_deta   = h*integral_deta_dtheta;

}

void fixedFfermi_derivatives_v2(const double k, const double eta, const double theta,
       const double h, double hMIN, double hMAX, 
       double * F, double *dF_deta, double *d2F_deta2,
       double * dF_dtheta, double *d2F_dtheta2, double *d2F_dtheta_deta,
       double * d3F_dtheta3, double *d3F_dtheta2_deta, double *d3F_dtheta_deta3, double *d3F_deta3,
	   int D_MAX)
{
  int ii,num;
  double t, x, dx, f, integral=0.0, integral_deta=0.0, integral_deta2=0.0, integral_dtheta=0.0, integral_dtheta2=0.0, integral_deta_dtheta=0.0;
  double xst,dxst,factor, denom, denomi; //AUX vars similar to used by FXT code
  
  num = (int) -hMIN/h;
  hMIN = -num*h;
  
  num = (int) hMAX/h;
  hMAX = num*h;
#if DEBUG  
  printf("fixedFfermi: %d %lf %lf\n", num, hMIN,hMAX);
#endif  
  
  num = (int) (hMAX-hMIN)/h;

  for(ii=0;ii<=num;ii++)
  {
    t = hMIN + ii*h;
    //t = hMAX - ii*h;

    x    = exp(t-exp(-t));
    dx   = 1.0 +exp(-t); // This is NOT D[x,t], D[x,t] = x*dx !
	xst  = 1.0+ 0.5*theta*x;
	dxst = sqrt(xst);
	factor = exp(x-eta);
	denom  = 1.0+factor;
    denomi = 1.0/denom;	
	f      = 1.0/(1.0+factor)*exp(   (k+1.0) * (t - exp(-t))   )*dxst * dx;
    integral             +=   f;
	if(D_MAX>0) integral_deta        +=   f*factor*denomi; 
	if(D_MAX>1) integral_deta2       +=   f*factor*denomi*denomi*(factor-1.0); 
	if(D_MAX>0) integral_dtheta      +=   f*0.25*x/xst;
	if(D_MAX>1) integral_dtheta2     +=  -f*0.25*x/xst*x/xst*0.25;
	if(D_MAX>1) integral_deta_dtheta +=   f*0.25*x/xst*factor*denomi;
  }
  
  *F                 = h*integral;
  *dF_deta           = h*integral_deta;
  *d2F_deta2         = h*integral_deta2;
  *dF_dtheta         = h*integral_dtheta;
  *d2F_dtheta2       = h*integral_dtheta2;
  *d2F_dtheta_deta   = h*integral_deta_dtheta;

}


double fixedFfermi(const double k, const double eta, const double theta,
		   const double h, double hMIN, double hMAX)
{
  int ii,num;
  double t, x, dx, integral=0.0;
  double c=0.0,tmp,y; // https://en.wikipedia.org/wiki/Kahan_summation_algorithm

  
  num = (int) -hMIN/h;
  hMIN = -num*h;
  
  num = (int) hMAX/h;
  hMAX = num*h;
#if DEBUG  
  printf("fixedFfermi: %d %lf %lf\n", num, hMIN,hMAX);
#endif  
  
  num = (int) (hMAX-hMIN)/h;
//  #pragma omp simd
//  #pragma ivdep
  for(ii=0;ii<=num;ii++)
  {
    t = hMIN + ii*h;
    //t = hMAX -	 ii*h;

    x = exp(t-exp(-t));
    dx = 1.0 +exp(-t);
    integral+= 1.0/(1.0+exp(x)*exp(-eta) )*exp(   (k+1.0) * (t - exp(-t))   )*sqrt(1.0+ 0.5*theta*x) * dx;

  //integral+= 1.0/(1.0+exp(x-eta) )*pow(x,k+1.0)*sqrt(1.0+ 0.5*theta*x) * dx;


/*
    x = exp(t-exp(-t));
    dx = 1.0 +exp(-t);
	y = 1.0/(1.0+exp(x)*exp(-eta) )
	            *exp(   (k+1.0) * (t - exp(-t))   )
                *sqrt(1.0+ 0.5*theta*x) * dx 
				- c;
	tmp = integral + y;
	c = (tmp-integral) - y;
	integral=tmp; */
  }

  return h*integral;  

}

long double fixedFfermi_long(const long double k, const long double eta, const long double theta,
		   const long double h, const long double hMIN, const long double hMAX)
{
  int ii,num;
  long double t, x, dx, integral=0.0L;
  long double h_MIN,h_MAX;

  
  num = (int) -hMIN/h;
  h_MIN = -num*h;
  
  num = (int) hMAX/h;
  h_MAX = num*h;
  

  num = (int) (h_MAX-h_MIN)/h;
  #pragma omp simd
  #pragma ivdep
  for(ii=0;ii<=num;ii++)
  {
    t = h_MIN + ii*h;
    //t = h_MAX -	 ii*h;
    x = expl(t-expl(-t));
    dx = 1.0L +expl(-t);
    integral+= 1.0L/(1.0L+expl(x)*expl(-eta) )
	            *expl(   (k+1.0L) * (t - expl(-t))   )
                *sqrtl(1.0L+ 0.5L*theta*x) * dx;
              
 
  }
  return h*integral;  

}

float fixedFfermif(const float k, const float eta, const float theta,
		   const float h, float hMIN, float hMAX)
{
  int ii,num;
  float t, x, dx, fd, integral=0.0f;

  
  num = (int) -hMIN/h;
  hMIN = -num*h;
  
  num = (int) hMAX/h;
  hMAX = num*h;
  num = (int) (hMAX-hMIN)/h;

    #pragma omp simd
    #pragma ivdep
    for(ii=0;ii<=num;ii++)
    {
      t = hMIN + ii*h;
      if( (eta>k) && (k>0) ) x = eta*expf(t-expf(-t)); else x = expf(t-expf(-t));
	  dx = 1.0f +expf(-t);
	  fd = 1.0f/(1.0f+expf(x-eta) );
      integral+= fd*powf(x,k+1)*sqrtf(1.0f+ 0.5f*theta*x)*dx;
	}
  
  return h*integral;  

}

double gaussFfermi(const double k, const double eta, const double theta)
{
      const double xg_leg[4] =  {0.18343464249564980493947614236018398066675781291297378231718847369920\
44742215421141160682237111233537, \
0.52553240991632898581773904918924634904196424312039285775085709927245\
48207685612725239614001936319821, \
0.79666647741362673959155393647583043683717173161596483207017029503921\
73056764730921471519272957259390, \
0.96028985649753623168356086856947299042823523430145203827163977737242\
48977434192844394389592633122683
      };

      const double wg_leg[4] =   {0.36268378337836198296515044927719561219414603989433054052482306756668\
67347239066773243660420848285096, \
0.31370664587788728733796220198660131326032899900273493769026394507495\
62719421734969616980762339285560, \
0.22238103445337447054435599442624088443013087005124956472590928929361\
6814570449040853653142377197928, \
0.10122853629037625915253135430996219011539409105168495705900369806474\
0178763470784860282739304045007};
      const int num=4; 
      
      const double d   = 3.3609e0 ;
      const double sg  = 9.1186e-2 ;
      const double a1  = 6.7774e0 ;
      const double b1  = 1.1418e0 ;
      const double c1  = 2.9826e0 ;
      const double a2  = 3.7601e0 ;
      const double b2  = 9.3719e-2 ;
      const double c2  = 2.1063e-2 ;
      const double d2  = 3.1084e1 ;
      const double e2  = 1.0056e0 ;
      const double a3  = 7.5669e0 ;
      const double b3  = 1.1695e0 ;
      const double c3  = 7.5416e-1 ;
      const double d3  = 6.6558e0 ;
      const double e3  =-1.2819e-1 ;
      
      double x1,x2,x3,s1,s2,s3,s12,xi,xi2,eta1;
      
      double fd, res1,res2,res3,res4;
      
//   definition of xi:
      eta1=sg*(eta-d);
      if (eta1<50.0) 
        xi=log(1.0+exp(eta1))/sg;
      else
        xi=eta-d;
      
      xi2=xi*xi;

//   definition of the x_i:
      x1=(a1  +b1*xi+c1*   xi2)/(1.0+c1*xi);
      x2=(a2  +b2*xi+c2*d2*xi2)/(1.0+e2*xi+c2*   xi2);
      x3=(a3  +b3*xi+c3*d3*xi2)/(1.0+e3*xi+c3*   xi2);

//   breakpoints:
      s1=x1-x2;
      s2=x1;
      s3=x1+x3;
      s12=sqrt(s1);
      
      //first integral of four: int_0^Sqrt[s1] f(x^2) 2x dx
      double center, hlfrun, result, x, z, f, a, b;
      int i;
      
      a=0.0;
      b=s12;
      
      center   = 0.5* (a+b);
      hlfrun   = 0.5* (b-a);
      result   = 0.0;
      
#pragma omp simd
#pragma ivdep
      for(i=0;i<num;i++)
      {
        x = center + hlfrun*xg_leg[i];
	z = x*x;
        f = pow(z,k)*sqrt(1.0+0.5*z*theta)/(1.0+exp(z-eta))*2.0*x;
        result   += f*wg_leg[i];
      }
#pragma omp simd            
#pragma ivdep
      for(i=0;i<num;i++)
      {
        x = center - hlfrun*xg_leg[i];
	z = x*x;
        f = pow(z,k)*sqrt(1.0+0.5*z*theta)/(1.0+exp(z-eta))*2.0*x;
        result   += f*wg_leg[i];
      }
      
      res1   = result * hlfrun;
      
      //integral 2 of 4: int_s1^s3 f(x) dx
      a=s1;
      b=s2;
      
      center   = 0.5* (a+b);
      hlfrun   = 0.5* (b-a);
      result   = 0.0;
      
#pragma omp simd
#pragma ivdep
      for(i=0;i<num;i++)
      {
        x = center + hlfrun*xg_leg[i];
        f = pow(x,k)*sqrt(1.0+0.5*x*theta)/(1.0+exp(x-eta));
        result   += f*wg_leg[i];
      }
#pragma omp simd            
#pragma ivdep
      for(i=0;i<num;i++)
      {
        x = center - hlfrun*xg_leg[i];
        f = pow(x,k)*sqrt(1.0+0.5*x*theta)/(1.0+exp(x-eta));
        result   += f*wg_leg[i];
      }
      
      res2   = result * hlfrun;

            //integral 3 of 4: int_s2^s3 f(x) dx
      a=s2;
      b=s3;
      
      center   = 0.5* (a+b);
      hlfrun   = 0.5* (b-a);
      result   = 0.0;
      
#pragma omp simd
#pragma ivdep
      for(i=0;i<num;i++)
      {
        x = center + hlfrun*xg_leg[i];
        f = pow(x,k)*sqrt(1.0+0.5*x*theta)/(1.0+exp(x-eta));
        result   += f*wg_leg[i];
      }
#pragma omp simd            
#pragma ivdep
      for(i=0;i<num;i++)
      {
        x = center - hlfrun*xg_leg[i];
        f = pow(x,k)*sqrt(1.0+0.5*x*theta)/(1.0+exp(x-eta));
        result   += f*wg_leg[i];
      }
      
      res3   = result * hlfrun;

      const double xg_lag[10] = { 1.37793470540492430830772505652711188e-1,
       7.29454549503170498160373121676078781e-1 ,
       1.80834290174031604823292007575060883e0 ,
       3.40143369785489951448253222140839067e0 ,
       5.55249614006380363241755848686876285e0 ,
       8.33015274676449670023876719727452218e0 ,
       1.18437858379000655649185389191416139e1 ,
       1.62792578313781020995326539358336223e1 ,
       2.19965858119807619512770901955944939e1 ,
       2.99206970122738915599087933407991951e1 };

       const double wg_lag[10] = {3.54009738606996308762226891442067608e-1,
       8.31902301043580738109829658127849577e-1 ,
       1.33028856174932817875279219439399369e0 ,
       1.86306390311113098976398873548246693e0 ,
       2.45025555808301016607269373165752256e0 ,
       3.12276415513518249615081826331455472e0 ,
       3.93415269556152109865581245924823077e0 ,
       4.99241487219302310201148565243315445e0 ,
       6.57220248513080297518766871037611234e0 ,
       9.78469584037463069477008663871859813e0 };
       

      


            //integral 4 of 4: int_s3^\inf f(x) dx
      a=s3;
      b=1.0;
      result   = 0.0;

#pragma omp simd
      for(i=0;i<10;i++)
      {
       x = a+b*xg_lag[i];
       f = pow(x,k)*sqrt(1.0+0.5*x*theta)/(1.0+exp(x-eta));
       result   +=  f*wg_lag[i];
      } 
      
      res4   = result*b;  

#if DEBUG
      printf("gaussFfermi: %.18e %.18e %.18e %.18e\n",res1,res2,res3,res4);
#endif      
      
  return res1+res2+res3+res4;
      
}


float gaussFfermif(const float k, const float eta, const float theta)
{
      const float xg_leg[5] =  {1.48874338981631210884826001129719984e-1f,
        4.33395394129247190799265943165784162e-1f,
         6.79409568299024406234327365114873575e-1f,
         8.65063366688984510732096688423493048e-1f,
         9.73906528517171720077964012084452053e-1f
      };

      const float wg_leg[5] =   {2.95524224714752870173892994651338329e-1f,
      2.69266719309996355091226921569469352e-1f,
      2.19086362515982043995534934228163192e-1f,
      1.49451349150580593145776339657697332e-1f,
      6.66713443086881375935688098933317928e-2f
      };
      
      const float d   = 3.3609e0 ;
      const float sg  = 9.1186e-2 ;
      const float a1  = 6.7774e0 ;
      const float b1  = 1.1418e0 ;
      const float c1  = 2.9826e0 ;
      const float a2  = 3.7601e0 ;
      const float b2  = 9.3719e-2 ;
      const float c2  = 2.1063e-2 ;
      const float d2  = 3.1084e1 ;
      const float e2  = 1.0056e0 ;
      const float a3  = 7.5669e0 ;
      const float b3  = 1.1695e0 ;
      const float c3  = 7.5416e-1 ;
      const float d3  = 6.6558e0 ;
      const float e3  =-1.2819e-1 ;
      
      float x1,x2,x3,s1,s2,s3,s12,xi,xi2,eta1;
      
      float fd, res1,res2,res3,res4;
      
//   definition of xi:
      eta1=sg*(eta-d);
      if (eta1<50.0f) 
        xi=logf(1.0f+expf(eta1))/sg;
      else
        xi=eta-d;
      
      xi2=xi*xi;

//   definition of the x_i:
      x1=(a1  +b1*xi+c1*   xi2)/(1.0f+c1*xi);
      x2=(a2  +b2*xi+c2*d2*xi2)/(1.0f+e2*xi+c2*   xi2);
      x3=(a3  +b3*xi+c3*d3*xi2)/(1.0f+e3*xi+c3*   xi2);

//   breakpoints:
      s1=x1-x2;
      s2=x1;
      s3=x1+x3;
      s12=sqrtf(s1);
      
      //first integral of four: int_0^Sqrt[s1] f(x^2) 2x dx
      float center, hlfrun, result, x, z, f, a, b;
      int i;
      
      a=0.0f;
      b=s12;
      
      center   = 0.5f* (a+b);
      hlfrun   = 0.5f* (b-a);
      result   = 0.0f;
      
#pragma omp simd
#pragma ivdep
      for(i=0;i<5;i++)
      {
        x = center + hlfrun*xg_leg[i];
	z = x*x;
        f = powf(z,k)*sqrtf(1.0f+0.5f*z*theta)/(1.0f+expf(z-eta))*2.0f*x;
        result   += f*wg_leg[i];
      }
#pragma omp simd            
#pragma ivdep
      for(i=0;i<5;i++)
      {
        x = center - hlfrun*xg_leg[i];
	z = x*x;
        f = powf(z,k)*sqrtf(1.0f+0.5f*z*theta)/(1.0f+expf(z-eta))*2.0f*x;
        result   += f*wg_leg[i];
      }
      
      res1   = result * hlfrun;
      
      //integral 2 of 4: int_s1^s3 f(x) dx
      a=s1;
      b=s2;
      
      center   = 0.5f* (a+b);
      hlfrun   = 0.5f* (b-a);
      result   = 0.0f;
      
#pragma omp simd
#pragma ivdep
      for(i=0;i<5;i++)
      {
        x = center + hlfrun*xg_leg[i];
        f = powf(x,k)*sqrtf(1.0f+0.5f*x*theta)/(1.0f+expf(x-eta));
        result   += f*wg_leg[i];
      }
#pragma omp simd            
#pragma ivdep
      for(i=0;i<5;i++)
      {
        x = center - hlfrun*xg_leg[i];
        f = powf(x,k)*sqrtf(1.0f+0.5f*x*theta)/(1.0f+expf(x-eta));
        result   += f*wg_leg[i];
      }
      
      res2   = result * hlfrun;

            //integral 3 of 4: int_s2^s3 f(x) dx
      a=s2;
      b=s3;
      
      center   = 0.5f* (a+b);
      hlfrun   = 0.5f* (b-a);
      result   = 0.0f;
      
#pragma omp simd
#pragma ivdep
      for(i=0;i<5;i++)
      {
        x = center + hlfrun*xg_leg[i];
        f = powf(x,k)*sqrtf(1.0f+0.5f*x*theta)/(1.0f+expf(x-eta));
        result   += f*wg_leg[i];
      }
#pragma omp simd            
#pragma ivdep
      for(i=0;i<5;i++)
      {
        x = center - hlfrun*xg_leg[i];
        f = powf(x,k)*sqrtf(1.0f+0.5f*x*theta)/(1.0f+expf(x-eta));
        result   += f*wg_leg[i];
      }
      
      res3   = result * hlfrun;

      const float xg_lag[10] = { 1.37793470540492430830772505652711188e-1f,
       7.29454549503170498160373121676078781e-1f ,
       1.80834290174031604823292007575060883e0f ,
       3.40143369785489951448253222140839067e0f ,
       5.55249614006380363241755848686876285e0f ,
       8.33015274676449670023876719727452218e0f ,
       1.18437858379000655649185389191416139e1f ,
       1.62792578313781020995326539358336223e1f ,
       2.19965858119807619512770901955944939e1f ,
       2.99206970122738915599087933407991951e1f };

       const float wg_lag[10] = {3.54009738606996308762226891442067608e-1f,
       8.31902301043580738109829658127849577e-1f ,
       1.33028856174932817875279219439399369e0f ,
       1.86306390311113098976398873548246693e0f ,
       2.45025555808301016607269373165752256e0f ,
       3.12276415513518249615081826331455472e0f ,
       3.93415269556152109865581245924823077e0f ,
       4.99241487219302310201148565243315445e0f ,
       6.57220248513080297518766871037611234e0f ,
       9.78469584037463069477008663871859813e0f };
       

      



      result   = 0.0f;
      a=s3;
      b=1.0f;

#pragma omp simd
      for(i=0;i<10;i++)
      {
       x = a+b*xg_lag[i];
       f = powf(x,k)*sqrtf(1.0f+0.5f*x*theta)/(1.0f+expf(x-eta));
       result   +=  f*wg_lag[i];
      } 
      
      res4   = result*b;  
      
  return res1+res2+res3+res4;
      
}

double quickFfermi0(const double k, const double eta, const double theta)
{
   
  int ii;
  double t, x, dx, integral=0.0;
  #include "tbl_double.h"  
  t=-4.0;
  #pragma omp simd
  #pragma ivdep
  for(ii=0;ii<num;ii++)
  {
    x = exp(t-exp(-t));
    dx = 1.0 +exp(-t);
    integral+= 1.0/(1.0+exp(x)*exp(-eta) )
	            *exp(   (k+1.0) * (t - exp(-t))   )
                *sqrt(1.0+ 0.5*theta*x) * dx;
              ;
    t+=h;		
  }
  return h*integral;  
    
}


double quickFfermi1(const double k, const double eta, const double theta)
{
  //tabulated abcissas from Mathematica
  // we can also tabulate tTBL[ii] - exp(-tTBL[ii]) and exp(xTBL[ii])
  #include "tbl_double.h"
  const double minus_exp_eta=exp(-eta);
  const double k_plus_1 = k+1.0;
  double integral = 0.0;
  int ii;
//#pragma omp simd
  for(ii=0;ii<num;ii++)
  {
    integral+= 1.0/(1.0+exp(xTBL[ii])*minus_exp_eta )
	            *exp(   (k_plus_1) * (tTBL[ii] - exp(-tTBL[ii]))   )
                *sqrt(1.0 + 0.5*theta*xTBL[ii]) * dxTBL[ii];
  }
  return h*integral;  
    
}

double quickFfermi2(const double k, const double eta, const double theta)
{
  //tabulated abcissas from Mathematica
  #include "tbl_double.h"
  const double minus_exp_eta=exp(-eta);
  const double k_plus_1 = k+1.0;
  const double half_theta = 0.5*theta;
  double integral = 0.0;
  int ii;
  
  #pragma omp simd
  #pragma vector aligned
  #pragma ivdep
  for(ii=0;ii<num;ii++)
  {
    integral+= 1.0/(1.0+ExpxTBL[ii]*minus_exp_eta )
	            *exp(   k_plus_1 * tExpTBL[ii]   )
                *sqrt(1.0 + half_theta*xTBL[ii]) * dxTBL[ii];
    ;
  }
  return h*integral;  
    
}

double quickFfermi3(const double k_dbl, const double eta_dbl, const double theta_dbl)
{
  //tabulated abcissas from Mathematica
#include "tbl_float.h"

  float k,eta,theta;
  float integral = 0.0f;
  int ii;
  
  k = (float) k_dbl;
  theta = (float) theta_dbl;
  eta = (float) eta_dbl;
  
  #pragma omp simd
  #pragma vector aligned
  #pragma ivdep
  for(ii=0;ii<num;ii++)
  {
    integral+= 1.0f/(1.0f+ExpxTBL[ii]*expf(-eta) )
	            *expf(   (k+1.0f) * tExpTBL[ii]   )
                *sqrtf(1.0f + 0.5f*theta*xTBL[ii]) * dxTBL[ii];
    ;
  }
  return (double) h*integral;  
    
}
