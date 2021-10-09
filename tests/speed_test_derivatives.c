/*
Compile with e.g:
gcc speed_test_derivatives.c -o speed_test -lfermidirac fermi_dirac_quadrature.o -lm
Run:
./example
*/
#include <fermidirac.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
/*       subroutine dfermi(dk,eta,theta,fd,fdeta,fdtheta, &
                        fdeta2,fdtheta2,fdetadtheta) 
DOWNLOAD from: http://cococubed.asu.edu/codes/fermi_dirac/fermi_dirac.tbz
*/
void dfermi_(double *, double *,double *,double *,double *,double *,double *,double *,double *);

int main()
{
  
  double k,eta,theta, result;
  double k_START,eta_START,theta_START;
  double k_STOP,eta_STOP,theta_STOP;
  double k_STEP,eta_STEP,theta_STEP;
  double F, dF_deta, d2F_deta2, 
		 dF_dtheta, d2F_dtheta2, d2F_deta_dtheta;
  double eq3_LHS,eq3_RHS;
  
  
  

  k_START     = -0.5;
  k_STOP      =  5.0;
  k_STEP      =  0.125;
			  
  eta_START   = -64.0;
  eta_STOP    = 64.0;
  eta_STEP    =  0.5;
			  
  theta_START =  0.0;
  theta_STOP  =  128.0;
  theta_STEP  =  1.0;
  
  k=k_START;
  do
  {
    eta   = eta_START;
    theta = theta_START;
	do
    {
	  theta = theta_START;
	  do
      {
       //printf("%lf %lf %lf\n", k, eta, theta);
       //fixedFfermi_derivatives   (k,eta,theta, 0.125, -5.0, 5.0, &F, &dF_deta, &d2F_deta2, &dF_dtheta, &d2F_dtheta2, &d2F_deta_dtheta);
       //fixedFfermi_derivatives_v2(k,eta,theta, 0.125, -5.0, 5.0, &F, &dF_deta, &d2F_deta2, &dF_dtheta, &d2F_dtheta2, &d2F_deta_dtheta, NULL, NULL, NULL, NULL, 0);
       //result = fixedFfermi(k,eta,theta, 0.125, -5.0, 5.0);
       //dfermi_(&k,&eta,&theta, &F, &dF_deta, &d2F_deta2, &dF_dtheta, &d2F_dtheta2, &d2F_deta_dtheta);
       result = quickFfermi2(k,eta,theta);

       theta = theta + theta_STEP;
	  }
      while(theta<=theta_STOP);
      eta=eta+eta_STEP;
	}
    while(eta<=eta_STOP);
    k= k+k_STEP;
  }
  while(k<=k_STOP);	  

  

  

  
  

  
  return 0;

}
