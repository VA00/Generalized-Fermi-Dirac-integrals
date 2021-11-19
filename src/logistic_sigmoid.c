/*
A. Odrzywolek, aodrzywolek@gmail.com, 2021-11-19
*/
#include "../fermidirac.h"
#include <math.h>
#include <float.h>

double sigmoid(double z)
{
  //return 1.0/(1.0+exp(-z));
  return 0.5*tanh(0.5*z)+0.5; //the fastest, clean (no conditionals)
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