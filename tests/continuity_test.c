#include <fermidirac.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
void dfermi_(double *, double *,double *,double *,double *,double *,double *,double *,double *);

int main()
{
  
  int i,j,m,n;
  double s, lead[DERIVATIVE_MATRIX_SIZE],D[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE];
  double result[10];
  double t=0.0, k=1.0, eta, theta=4.0, z;
  double split_point=8192.0;
  double nearby[256];
    double F, dF_deta, d2F_deta2, 
		 dF_dtheta, d2F_dtheta2, d2F_deta_dtheta,
		 d3F_dtheta3, d3F_dtheta2_deta, d3F_dtheta_deta2, d3F_deta3;


   
  eta=split_point;
  for(i=128;i<256;i++)
   {
     nearby[i]=eta;
     eta=nextafter(eta,DBL_MAX);
   }
   
  eta=split_point; 
  for(i=128;i>=0;i--)
   {
     nearby[i]=eta;
     eta=nextafter(eta,-DBL_MAX);
   }


   //for(i=0;i<256;i++) printf("%.17e\n", nearby[i]);


  
  for(i=0;i<256;i++)
   {   
       eta = nearby[i];
       //dfermi_(&k,&eta,&theta, &F, &dF_deta, &dF_dtheta, &d2F_deta2, &d2F_dtheta2, &d2F_deta_dtheta);
       //D[0][0]=F; D[0][1] = dF_dtheta; D[0][2] = d2F_dtheta2; D[1][0] = dF_deta; D[1][1] = d2F_deta_dtheta; D[2][0] = d2F_deta2;
       Ffermi_derivatives_matrix(k,nearby[i],theta,D);
       for(m=0;m<DERIVATIVE_MATRIX_SIZE;m++)
        for(n=0;n<DERIVATIVE_MATRIX_SIZE;n++)
         if(m+n<=3)
           printf("%.17e\t", D[m][n]);
         else
          continue;
    printf("%d\t", i-128);
    printf("%.17e\t", nearby[i]);
    printf("\n");
   }
  
  

  
 
 
  //printf("\n");  
  return 0;

}
