/*
Compile with e.g:
gcc example_derivatives.c -o example -lm -lfermidirac
Run:
./example
*/
#include <fermidirac.h>
#include <stdio.h>
#include <math.h>
#include <float.h>


int main()
{
  
  double k,eta,theta, result, 
         F, dF_deta, d2F_deta2, 
		 dF_dtheta, d2F_dtheta2, d2F_deta_dtheta;
  double eq3_LHS,eq3_RHS;
  
  
  
  result = Ffermi(4.0,1.0,1.0);
  
  printf("Ffermi(4,1,1)=%.16e\n",result);
  
  k     = 4.0;
  eta   = 1.0;
  theta = 1.0;
  
  fixedFfermi_derivatives(k,eta,theta, 0.125, -5.0, 5.0, &F, &dF_deta, &d2F_deta2, &dF_dtheta, &d2F_dtheta2, &d2F_deta_dtheta);
  printf("Ffermi(4,1,1)   = %.16e\n", F);
  printf("dF/deta         = %.16e\n", dF_deta);
  printf("d2F/deta2       = %.16e\n", d2F_deta2);
  printf("dF/dtheta       = %.16e\n", dF_dtheta);
  printf("d2F/dtheta2     = %.16e\n", d2F_dtheta2);
  printf("d2F/dtheta/deta = %.16e\n", d2F_deta_dtheta);
  
  
  eq3_LHS = theta*dF_dtheta+(1.0+k)*F;
  fixedFfermi_derivatives(k+1.0,eta,theta, 0.125, -5.0, 5.0, &F, &dF_deta, &d2F_deta2, &dF_dtheta, &d2F_dtheta2, &d2F_deta_dtheta);
  eq3_RHS = dF_deta;
  printf("\nCox and Giuli test = %.16e\t rel_err/DBL_EPSILON=%lf\n", eq3_LHS-eq3_RHS,(eq3_LHS-eq3_RHS)/eq3_RHS/DBL_EPSILON);
  eq3_LHS = F - 2.0*theta*dF_dtheta;
  fixedFfermi_derivatives(k,eta,theta, 0.125, -5.0, 5.0, &F, &dF_deta, &d2F_deta2, &dF_dtheta, &d2F_dtheta2, &d2F_deta_dtheta);
  eq3_RHS = 4.0*dF_dtheta;
  printf("\nGong test = %.16e\t rel_err/DBL_EPSILON=%lf\n", eq3_LHS-eq3_RHS,(eq3_LHS-eq3_RHS)/eq3_RHS/DBL_EPSILON);
  
  
  //TODO compute ULP difference using nextafter()
  

  
  

  
  return 0;

}
