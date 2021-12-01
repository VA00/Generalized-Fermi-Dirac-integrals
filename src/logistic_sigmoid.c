/*
A. Odrzywolek, aodrzywolek@gmail.com, 2021-11-19
*/
#include "../fermidirac.h"
#include <math.h>
#include <float.h>
#include <quadmath.h>

double sigmoid(double z)
{
  return 1.0/(1.0+exp(-z));
  //return 0.5*tanh(0.5*z)+0.5; //the fastest, clean (no conditionals)
}

long double sigmoid_long(long double z)
{
  return 1.0l/(1.0l+expl(-z));
  //return 0.5L*tanhl(0.5L*z)+0.5L; //the fastest, clean (no conditionals)
}

__float128 sigmoid_quad(__float128 z)
{
  //return 1.0q/(1.0q+expq(-z));
  return 0.5q*tanhq(0.5q*z)+0.5q; //the fastest, clean (no conditionals)
}


// Mathematica: (-1)^i Sum[Eulerian[i, j] LogisticSigmoid[z]^(j + 1) (LogisticSigmoid[z] - 1)^(i - j), {j, 0, i-1}]
double sigmoid_derivative(double z, int i)
{
  double sum=0.0, s=sigmoid(z);
  int j;
  

   
  if(i==0) return s; //general formula do not work for i==0

  for(j=0;j<=i-1;j++)
   sum = sum + eulerian(i,j)*power_squaring(s,j + 1)*power_squaring(s - 1.0,i - j);

  return ((i%2) ? -1.0 : 1.0)*sum;

}



// Same as above, but using already computed sigmoid as an argument, without leading sigma !
double sigmoid_derivative_polynomial(double s, int i)
{
  double sum=0.0;
  int j;
   
  if(i==0) return 1.0;

  for(j=0;j<=i-1;j++)
   sum = sum + eulerian(i,j)*power_squaring(s,j)*power_squaring(s - 1.0,i - j);

  return ((i%2) ? -1.0 : 1.0)*sum;

}

long double sigmoid_derivative_polynomial_long(long double s, int i)
{
  long double sum=0.0L;
  int j;
   
  if(i==0) return 1.0L;
  if(i==1) return 1.0L - s;
  if(i==2) return 1.0L + s*(2.0L*s-3.0L);
  if(i==3) return s*((12.0L - 6.0L*s)*s - 7.0L) + 1.0L;

  for(j=0;j<=i-1;j++)
   sum = sum + eulerian(i,j)*power_squaring(s,j)*power_squaring(s - 1.0,i - j);

  return ((i%2) ? -1.0 : 1.0)*sum;

}

__float128 sigmoid_derivative_polynomial_quad(__float128 s, int i)
{
  __float128 sum=0.0q;
  int j;
   
  if(i==0) return 1.0q;
  if(i==1) return 1.0q - s;
  if(i==2) return 1.0q + s*(2.0q*s-3.0q);
  if(i==3) return s*((12.0q - 6.0q*s)*s - 7.0q) + 1.0q;

  for(j=0;j<=i-1;j++)
   sum = sum + eulerian(i,j)*power_squaring(s,j)*power_squaring(s - 1.0,i - j);

  return ((i%2) ? -1.0 : 1.0)*sum;

}


// Same as above, but using already computed sigmoid as an argument, without leading sigma !
void sigmoid_derivative_polynomial_vector(double s, double ds[DERIVATIVE_MATRIX_SIZE])
{
  double sum=0.0;
  int i,j;
   
  ds[0]=1.0;

  for(i=1;i<DERIVATIVE_MATRIX_SIZE;i++)
   {
    sum=0.0;
    for(j=0;j<=i-1;j++)
     sum = sum + eulerian(i,j)*power_squaring(s,j)*power_squaring(s - 1.0,i - j);
    ds[i]=((i%2) ? -1.0 : 1.0)*sum;
   }

  return;

}