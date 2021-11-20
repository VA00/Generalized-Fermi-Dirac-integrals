#include <fermidirac.h>
#include <stdio.h>
#include <math.h>
#include <float.h>


int main()
{
  
  int i,j,m,n;
  double s, lead[DERIVATIVE_MATRIX_SIZE],D[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE];
  double result[10];
  double t=0.0, k=1.0, eta=2222.0, theta=4.0, z;
  double split_point=65536.0;
  double nearby[256];

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
   {   Ffermi_derivatives_matrix(k,nearby[i],theta,D);
       for(m=0;m<DERIVATIVE_MATRIX_SIZE;m++)
        for(n=0;n<DERIVATIVE_MATRIX_SIZE;n++)
         if(m+n<=3)
           printf("%.17e\t", D[m][n]);
         else
          continue;
    printf("\n");
   }
  
  

  
 
 
  //printf("\n");  
  return 0;

}
