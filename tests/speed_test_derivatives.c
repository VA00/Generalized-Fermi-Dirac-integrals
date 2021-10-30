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
		 dF_dtheta, d2F_dtheta2, d2F_deta_dtheta,
		 d3F_dtheta3, d3F_dtheta2_deta, d3F_dtheta_deta2, d3F_deta3;
  double eq3_LHS,eq3_RHS;
  int counter=0;
  
  
  

  k_START     = -0.5;
  k_STOP      =  5.0;
  k_STEP      =  2.5;
			  
  eta_START   = -4.0;
  eta_STOP    =  4.0;
  eta_STEP    =  5.0;
			  
  theta_START =  0.0;
  theta_STOP  =  8.0;
  theta_STEP  =  1;
  
  counter=0;
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
		  counter++;
       //printf("%lf %lf %lf\n", k, eta, theta);
       //fixedFfermi_derivatives   (k,eta,theta, 0.125, -5.0, 5.0, &F, &dF_deta, &d2F_deta2, &dF_dtheta, &d2F_dtheta2, &d2F_deta_dtheta);
       fixedFfermi_derivatives_v2(k,eta,theta, 0.125, -5.0, 5.0, &F, &dF_deta, &d2F_deta2, &dF_dtheta, &d2F_dtheta2, &d2F_deta_dtheta, &d3F_dtheta3, &d3F_dtheta2_deta, &d3F_dtheta_deta2, &d3F_deta3, 3);
       //fixedFfermi_derivatives_v3(k,eta,theta, 0.125, -5.0, 5.0, &F, &dF_deta, &d2F_deta2, &dF_dtheta, &d2F_dtheta2, &d2F_deta_dtheta, &d3F_dtheta3, &d3F_dtheta2_deta, &d3F_dtheta_deta2, &d3F_deta3);
       printf("%.1lf %.1lf %.1lf %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",k,eta,theta,F, dF_deta, d2F_deta2, dF_dtheta, d2F_dtheta2, d2F_deta_dtheta, d3F_dtheta3, d3F_dtheta2_deta, d3F_dtheta_deta2, d3F_deta3);
       //result = fixedFfermi(k,eta,theta, 0.125, -5.0, 5.0);
       //result = Ffermi(k,eta,theta);
       dfermi_(&k,&eta,&theta, &F, &dF_deta, &dF_dtheta, &d2F_deta2, &d2F_dtheta2, &d2F_deta_dtheta);
       printf("%.1lf %.1lf %.1lf %.17e %.17e %.17e %.17e %.17e %.17e\n",k,eta,theta,F, dF_deta, d2F_deta2, dF_dtheta, d2F_dtheta2, d2F_deta_dtheta);
       printf("\n");
       //result = quickFfermi2(k,eta,theta);

       theta = theta + theta_STEP;
	  }
      while(theta<=theta_STOP);
      eta=eta+eta_STEP;
	}
    while(eta<=eta_STOP);
    k= k+k_STEP;
  }
  while(k<=k_STOP);	  

  

  printf("counter = %d\n",counter);

  
  

  
  return 0;

}
