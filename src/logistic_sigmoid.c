#include <math.h>
#include <float.h>

double sigmoid(double z)
{
  //return 1.0/(1.0+exp(-z));
  return 0.5*tanh(0.5*z)+0.5; //the fastest, clean (no conditionals)
}
