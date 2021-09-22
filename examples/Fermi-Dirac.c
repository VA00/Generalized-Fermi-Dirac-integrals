#include "mathlink.h"
#include <fermidirac.h>



double FFermi( double k, double eta, double theta)
{
  return Ffermi(k,eta,theta);
}

double GFermi( double n, double alpha, double beta)
{
  return Gfermi(n,alpha,beta);
}



/*********************************************************************/

int main(int argc, char **argv)
{
  return MLMain(argc, argv);
}