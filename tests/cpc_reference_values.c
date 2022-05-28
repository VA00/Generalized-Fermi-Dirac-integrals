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

double dfermi200_(int *, double *,double *,double *); //Gong GFDI

int idx, gong_idx[4][4] = {{0, 2, 5, 9}, {1, 4, 8, -1}, {3, 7, -1, -1}, {6, -1, -1, -1}}; //Gong et. al, Table 1. p. 299, CPC 136 (2001)



int main(int argc, char *argv[])
{
  double fd_quad,fd_double;
  acb_t fd_arb, fd;
  arb_t fd_real, rel_err, MachineEpsilon, MaxMachineNumber, MinMachineNumber;

  int m,n, i,j, sign, counter=0, underflow=0, overflow=0, failed=0, lg2_max;
  double k, eta, theta, GONG, Arb;
  //di_t interval; Require the most recent Arb, I'm unable to install it A.O. 
  

  flint_printf("Computed with arb-%s\n", arb_version);

  acb_init(fd_arb);
  acb_init(fd);
  arb_init(fd_real);
  arb_init(rel_err);
  arb_init(MachineEpsilon);
  arb_one(MachineEpsilon);
  arb_mul_2exp_si(MachineEpsilon, MachineEpsilon, -52); 
  arb_init(MaxMachineNumber);
  arb_one(MaxMachineNumber);
  arb_mul_2exp_si(MaxMachineNumber, MaxMachineNumber, 1024); 
  arb_init(MinMachineNumber);
  arb_one(MinMachineNumber);
  arb_mul_2exp_si(MinMachineNumber, MinMachineNumber, -1022); 



#if 0
  sscanf(argv[1],"%lf",&k);
  sscanf(argv[2],"%d",&lg2_max);

for(m=0;m<=3;m++)
 for(n=0;n<=3;n++)
  {
  if(m+n>3) continue; 
  printf("Testing derivative %d%d\n\n",m,n);
//#pragma omp parallel for collapse(3) private(idx, eta,theta, fd, fd_real, fd_arb, fd_quad, rel_err, GONG, sign,i,j) shared(lg2_max,m,n) reduction(+ :counter, overflow, underflow, failed)	  
  for(sign=-1;sign<=1;sign=sign+2)
   for(i=-lg2_max;i<=lg2_max;i=i+1)
    for(j=-lg2_max;j<=lg2_max;j=j+1)
     {
       counter++; 

       eta = sign*pow(2,i);
       theta = pow(2,j);


       if(!(counter%1024))   printf("Total tested = %d, overflow=%d, underflow=%d, failed=%d\n",counter, overflow, underflow, failed); 
       Ffermi_derivatives_m_n_arb(fd_arb, k, eta, theta, m, n);
       acb_get_real(fd_real,fd_arb); 
       if( arb_ge(fd_real, MaxMachineNumber) ){ overflow++;continue;}      
       if( arb_le(fd_real, MinMachineNumber) ){ underflow++;continue;}      

       fd_quad = (double) Ffermi_derivatives_m_n_quad((__float128) k, (__float128) eta, (__float128) theta, m, n);
       //fd_quad = (double) Ffermi_sommerfeld_derivatives_m_n_quad(k, eta, theta, m, n,powq(2.0q,-64.0q),11);
       //fd_quad = (double) Ffermi_dblexp_derivatives_m_n_quad(k,eta,theta, m, n, (__float128) 128*DBL_EPSILON, MAX_REFINE);
       //fd_double = Ffermi_derivatives_m_n_internal_arb(0.5, sign*pow(2,i), pow(2,j), m, n);   
       //fd_quad = Ffermi_derivatives_m_n(k, eta, theta, m, n);

       //printf("%e\n", fd_quad/fd_double-1.0);
     
       acb_set_d(fd, fd_quad);
     
       acb_div(fd, fd, fd_arb, 128);
       acb_add_si(fd, fd, -1, 128);
       acb_abs(rel_err, fd, 128);


       if( arb_ge(rel_err, MachineEpsilon) || (!acb_is_finite(fd_arb)) )
       {
       failed++;
       printf("\nArb=         ");
       //acb_print(fd_arb);printf("\n");
       //acb_printd(fd_arb, 128);printf("\n");
       acb_printn(fd_arb, 128, 0);printf("\n");
       printf("libfermidirac=%.18e\n",  fd_quad);     
      
       idx = gong_idx[m][n];
       GONG = dfermi200_(&idx,&k, &eta, &theta);
       printf("GONG         =%.18e\n",GONG);



     
       printf("%d %d\t%d\t",sign,i,j); 
       printf("k=%.1lf eta=%e theta=%e m=%d n=%d\n", k, eta, theta, m, n);
       printf("Relative error:\t");arb_printn(rel_err, 128, 0);

       acb_set_d(fd, GONG);
     
       acb_div(fd, fd, fd_arb, 128);
       acb_add_si(fd, fd, -1, 128);
       acb_abs(rel_err, fd, 128);
       printf("\nRelative error:\t");arb_printn(rel_err, 128, 0);

       printf("\n0.5*eta*theta=%e\n", eta*theta*0.5);
      
       printf("\n------------------------------------------------------------------------------\n\n");
       }


     }
  }


#endif

#if 1
  double lgeta, fd_dblexp, fd_somm;
  m=3;n=0;  printf("Testing derivative %d%d\n\n",m,n);
  sign=1;
  k=-0.75;
  theta = pow(2.0,1.0);

for(lgeta=1.0;lgeta<=24.0;lgeta=lgeta+1.0)
 {

  eta = pow(2.0,lgeta);
  Arb = Ffermi_derivatives_m_n_internal_arb(k, eta, theta, m, n);
  
  //fd_dblexp = (double) Ffermi_dblexp_derivatives_m_n_quad(k,eta,theta, m, n, (__float128) 8*DBL_EPSILON, MAX_REFINE);
  fd_dblexp = Ffermi_dblexp_derivatives_m_n(k,eta,theta, m, n, 8*DBL_EPSILON, MAX_REFINE);
  fd_somm = (double) Ffermi_sommerfeld_derivatives_m_n_quad(k,eta,theta, m, n, (__float128) 8*DBL_EPSILON, MAX_REFINE);
  //fd_double = (double) Ffermi_derivatives_m_n_quad(k,eta,theta, m, n);
   
  idx = gong_idx[m][n];
  GONG = dfermi200_(&idx,&k, &eta, &theta);

  printf("%lf\t%.18e\t%.18e\t%.18e\t%.18e\n",eta,Arb,GONG,fd_dblexp,fd_somm);
 }

#endif



  acb_clear(fd_arb);
  acb_clear(fd);
  arb_clear(fd_real);
  arb_clear(rel_err);
  arb_clear(MachineEpsilon);
  arb_clear(MaxMachineNumber);
  arb_clear(MinMachineNumber);

  printf("Finally tested = %d, overflow=%d, underflow=%d, failed=%d\n",counter, overflow, underflow, failed); 

  return 0;
}
