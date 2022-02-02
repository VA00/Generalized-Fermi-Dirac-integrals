/*
Andrzej Odrzywolek, 2022-02-02, andrzej.odrzywolek@uj.edu.pl
Compile with e.g:
gcc precision_test_derivatives.c -o test -lfermidirac -lquadmath
Run:
./test
*/
#include <fermidirac.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
/*       subroutine dfermi(dk,eta,theta,fd,fdeta,fdtheta, &
                        fdeta2,fdtheta2,fdetadtheta) 
DOWNLOAD from: http://cococubed.asu.edu/codes/fermi_dirac/fermi_dirac.tbz
*/
void dfermi_(double *, double *,double *,double *,double *,double *,double *,double *,double *); //FXT dfermi
double dfermi200_(int *, double *,double *,double *); //Gong GFDI

int main()
{
  
  double k,eta,theta;

  double ref[4][4];
  double val[4][4];
  
  double x[13];
  int m,n;
  
    FILE  *datafile;
  
  datafile = fopen("refVALS/refTBL_double_2022-02-02.bin","r");

  
  if(datafile==NULL){ printf("ERROR: UNABLE TO OPEN FILE\n");return -1;}
  

  while(fread(x, 8, 13, datafile) == 13)
   {  
    k               = x[0];
	eta             = x[1];
	theta           = x[2];
	ref[0][0]       = x[3];
    ref[0][1]       = x[4];
    ref[0][2]       = x[5]; 
    ref[0][3]       = x[6];
    ref[1][0]       = x[7]; 
    ref[1][1]       = x[8]; 
    ref[1][2]       = x[9]; 
    ref[2][0]       = x[10]; 
    ref[2][1]       = x[11]; 
    ref[3][0]       = x[12]; 
    
    

//       printf("%lf %lf %lf\n", k, eta, theta);
    for(m=0;m<=3;m++)
     for(n=0;n<=3;n++)
      {
       if(m+n>3) continue;
       val[m][n] = Ffermi_derivatives_m_n_quad(k,eta,theta,m,n); 
       //printf("%.3e\t%.3e\n", val[m][n], ref[m][n]); 
       printf("%.3e\t", (val[m][n]/ref[m][n]-1.0)/DBL_EPSILON); 
  
      }
    printf("\n");  
    

  }

  
  

  fclose(datafile);
  return 0;

}
