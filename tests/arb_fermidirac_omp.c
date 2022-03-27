/*
Based on: https://github.com/fredrik-johansson/arb/blob/master/examples/integrals.c

A. Odrzywołek, 2022-02-16, andrzej.odrzywolek@uj.edu.pl

Install Arb https://arblib.org/

sudo apt install libflint-arb-dev
sudo apt install libflint-arb2
sudo apt install libflint-dev


Compile with: gcc arb_fermidirac.c -o arb_fd -lflint -lflint-arb -lm -lfermidirac -lquadmath dfermi200.o fedi_cpc.o -lgfortran


*/


/*

Generalized Fermi–Dirac functions and derivatives: properties and evaluation
Published: 1 June 2001
|
Version 1
|
DOI:
10.17632/57tnc6sby7.1
Contributors: Zhigang Gong, Ladislav Zejda, Werner Däppen, Josep M. Aparicio

DOWNLOAD CODE FROM:

https://elsevier.digitalcommonsdata.com/datasets/57tnc6sby7/1

Compile instructions:

gfortran -c fedi_cpc.for
gfortran -c dfermi200.for

Copy files  fedi_cpc.o, dfermi200.o to test/ subdir

*/
#include <string.h>
#include <acb_calc.h>
#include <arb.h>
#include <acb_hypgeom.h>
//#include "double_interval.h" // Require the most recent Arb, I'm unable to install it A.O. 
#include <fermidirac.h>
#include <quadmath.h>
#include <omp.h>
#include "ULP.c"

double dfermi200_(int *, double *,double *,double *); //Gong GFDI

int idx, gong_idx[4][4] = {{0, 2, 5, 9}, {1, 4, 8, -1}, {3, 7, -1, -1}, {6, -1, -1, -1}}; //Gong et. al, Table 1. p. 299, CPC 136 (2001)



int main(int argc, char *argv[])
{
  double fd_quad,fd_double;
 
  int m,n, i,j, sign, counter=0, underflow=0, overflow=0, failed=0, lg2_max;
  double k, eta, theta, GONG, Arb;
  

 
  sscanf(argv[1],"%lf",&k);
  sscanf(argv[2],"%d",&lg2_max);


for(m=0;m<=3;m++)
 for(n=0;n<=3;n++)
  {
  if(m+n>3) continue; 
  printf("Testing derivative %d%d\n\n",m,n);
#pragma omp parallel for collapse(3) private(eta,theta, fd_double, fd_quad, sign,i,j) shared(lg2_max,m,n) reduction(+ :counter, overflow, underflow, failed)  
  for(sign=-1;sign<=1;sign=sign+2)
   for(i=-lg2_max;i<=lg2_max;i=i+1)
    for(j=-lg2_max;j<=lg2_max;j=j+1)
     {
       counter++; 

       eta = sign*pow(2,i);
       theta = pow(2,j);

       fd_quad   = (double) Ffermi_derivatives_m_n_quad((__float128) k, (__float128) eta, (__float128) theta, m, n);
       fd_double = Ffermi_derivatives_m_n_internal_arb(k, eta, theta, m, n);
       if(fd_quad!=fd_double)
        {
         failed++;
         printf("%d %d k=%lf eta=%e theta=%e, %.17e\t%.17e\tULP=%d\n",m,n,k,eta,theta,fd_double, fd_quad, ULP_distance(fd_double,fd_quad,128));
        }
     }
  }





  printf("Finally tested = %d, overflow=%d, underflow=%d, failed=%d\n",counter, overflow, underflow, failed); 

  return 0;
}
