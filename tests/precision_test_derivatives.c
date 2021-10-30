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
#include "refVALS/20211030/refVALS_F_20211030.h"
/*       subroutine dfermi(dk,eta,theta,fd,fdeta,fdtheta, &
                        fdeta2,fdtheta2,fdetadtheta) 
DOWNLOAD from: http://cococubed.asu.edu/codes/fermi_dirac/fermi_dirac.tbz
*/
void dfermi_(double *, double *,double *,double *,double *,double *,double *,double *,double *);

int main()
{
  
  double k,eta,theta, result;

  double F, dF_deta, d2F_deta2, 
		 dF_dtheta, d2F_dtheta2, d2F_deta_dtheta,
		 d3F_dtheta3, d3F_dtheta2_deta, d3F_dtheta_deta2, d3F_deta3;
  double fdREF, dFdeta, d2Fdeta2, 
		 dFdtheta, d2Fdtheta2, d2Fdetadtheta,
		 d3Fdtheta3, d3Fdtheta2deta, d3Fdthetadeta2, d3Fdeta3;

  double eq3_LHS,eq3_RHS;
  double x[13];
  int ii;
  
    FILE  *datafile;
  
  datafile = fopen("refVALS/20211030/refVALS_F_20211030.bin","r");

  
  if(datafile==NULL){ printf("ERROR: UNABLE TO OPEN FILE\n");return -1;}
  

  for(ii=0;ii<num;ii=ii+1)
   {
	fread(x, 8, 13, datafile);

    k               = x[0];
	eta             = x[1];
	theta           = x[2];
	fdREF           = x[3];
    dFdeta          = x[4];
    d2Fdeta2        = x[5]; 
    dFdtheta        = x[6]; 
    d2Fdtheta2      = x[7]; 
    d2Fdetadtheta   = x[8]; 
    d3Fdtheta3      = x[9]; 
    d3Fdtheta2deta  = x[10]; 
    d3Fdthetadeta2  = x[11]; 
    d3Fdeta3        = x[12];
    
    printf("%.3e %.3e %.3e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",k,eta,theta,fdREF, dFdeta, d2Fdeta2, dFdtheta, d2Fdtheta2, d2Fdetadtheta, d3Fdtheta3, d3Fdtheta2deta, d3Fdthetadeta2, d3Fdeta3);
       

       //printf("%lf %lf %lf\n", k, eta, theta);
       //fixedFfermi_derivatives   (k,eta,theta, 0.125, -5.0, 5.0, &F, &dF_deta, &d2F_deta2, &dF_dtheta, &d2F_dtheta2, &d2F_deta_dtheta);
       //fixedFfermi_derivatives_v2(k,eta,theta, 0.125, -5.0, 5.0, &F, &dF_deta, &d2F_deta2, &dF_dtheta, &d2F_dtheta2, &d2F_deta_dtheta, &d3F_dtheta3, &d3F_dtheta2_deta, &d3F_dtheta_deta2, &d3F_deta3, 3);
       fixedFfermi_derivatives_v3(k,eta,theta, 0.125, -5.0, 5.0, &F, &dF_deta, &d2F_deta2, &dF_dtheta, &d2F_dtheta2, &d2F_deta_dtheta, &d3F_dtheta3, &d3F_dtheta2_deta, &d3F_dtheta_deta2, &d3F_deta3);
       printf("%.3e %.3e %.3e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",k,eta,theta,F, dF_deta, d2F_deta2, dF_dtheta, d2F_dtheta2, d2F_deta_dtheta, d3F_dtheta3, d3F_dtheta2_deta, d3F_dtheta_deta2, d3F_deta3);
       printf("%.3e %.3e %.3e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
          k,eta,theta,     (F/fdREF-1.0)/DBL_EPSILON, 
                    (dF_deta/dFdeta-1.0)/DBL_EPSILON,  
                (d2F_deta2/d2Fdeta2-1.0)/DBL_EPSILON, 
                (dF_dtheta/dFdtheta-1.0)/DBL_EPSILON, 
            (d2F_dtheta2/d2Fdtheta2-1.0)/DBL_EPSILON, 
     (d2F_deta_dtheta/d2Fdetadtheta-1.0)/DBL_EPSILON, 
            (d3F_dtheta3/d3Fdtheta3-1.0)/DBL_EPSILON, 
   (d3F_dtheta2_deta/d3Fdtheta2deta-1.0)/DBL_EPSILON, 
   (d3F_dtheta_deta2/d3Fdthetadeta2-1.0)/DBL_EPSILON, 
                 (d3Fdeta3/d3Fdeta3-1.0)/DBL_EPSILON);
       //result = fixedFfermi(k,eta,theta, 0.125, -5.0, 5.0);
       //result = Ffermi(k,eta,theta);
       //dfermi_(&k,&eta,&theta, &F, &dF_deta, &dF_dtheta, &d2F_deta2, &d2F_dtheta2, &d2F_deta_dtheta);
       //printf("%.1lf %.1lf %.1lf %.17e %.17e %.17e %.17e %.17e %.17e\n",k,eta,theta,F, dF_deta, d2F_deta2, dF_dtheta, d2F_dtheta2, d2F_deta_dtheta);
       printf("\n");
       //result = quickFfermi2(k,eta,theta);
   }

  
  

  
  return 0;

}
