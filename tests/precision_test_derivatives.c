/*
Andrzej Odrzywolek, 2022-02-02, andrzej.odrzywolek@uj.edu.pl
Compile with e.g:
gcc precision_test_derivatives.c -o test -lfermidirac -lquadmath fedi_cpc.o dfermi200.o -lgfortran -lm
Run:
./test refVALS/refTBL_double_2022-02-09.bin
*/
#include <fermidirac.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <quadmath.h>
/*       subroutine dfermi(dk,eta,theta,fd,fdeta,fdtheta, &
                        fdeta2,fdtheta2,fdetadtheta) 
DOWNLOAD from: http://cococubed.asu.edu/codes/fermi_dirac/fermi_dirac.tbz
*/
void dfermi_(double *, double *,double *,double *,double *,double *,double *,double *,double *); //FXT dfermi
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


*/
double dfermi200_(int *, double *,double *,double *); //Gong GFDI

#include "ULP.c"

int main( int argc, char** argv)
{
  
  double k,eta,theta;

  double ref[4][4];
  double val[4][4];
  
  double x[13];
  int m,n,ULP_test, failed_count, count;
  int idx, gong_idx[4][4] = {{0, 2, 5, 9}, {1, 4, 8, -1}, {3, 7, -1, -1}, {6, -1, -1, -1}}; //Gong et. al, Table 1. p. 299, CPC 136 (2001)
  
  FILE  *datafile;
  char  *refDATA;
  

  if(!(argv[1]==NULL)) 
   refDATA = argv[1];
  else
   refDATA = "refVALS/k_slices/0.5/refTBL_double_2022-02-14.bin";
  
  datafile = fopen(refDATA,"r");
  
  if(datafile==NULL){ printf("ERROR: UNABLE TO OPEN FILE\n");return -1;}
  
  printf("Testing started. Perfect results (10 x 0 ULPs) are NOT printed, except the first\n");
  printf("F\tdF dtheta\td2F dtheta2\td3F dtheta3\tdF deta\td2F deta dtheta\t d3F deta dtheta2\td2D deta2\td3F deta^2 dtheta\td^3F deta^3\n");
  count=0;
  failed_count = 0;

  while(fread(x, 8, 13, datafile) == 13)
   {  
    count++;
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
    

    ULP_test=0;

    for(m=0;m<=3;m++)
     for(n=0;n<=3;n++)
      {
       if(m+n>3) continue;
       if(theta>DBL_MAX) continue;
       val[m][n] = Ffermi_derivatives_m_n_quad(k,eta,theta,m,n); 

       ULP_test = ULP_test + ULP_distance(ref[m][n], val[m][n], 1024*1024);
      }
    
    if(ULP_test>100 || count==1)
     {
      failed_count++;
      printf("k=% .1lf eta=% .2e theta=%.2e\t", k, eta, theta);
      for(m=0;m<=3;m++)
       for(n=0;n<=3;n++)
        {
         if(m+n>3) continue;
         if(theta>DBL_MAX) continue;
         //printf("% .1e ", (val[m][n]/ref[m][n]-1.0)/DBL_EPSILON);
         printf("%d\t", ULP_distance(ref[m][n], val[m][n], 1024*1024) );  
        }
      printf("\t(libfermidirac)\n"); 
      printf("k=% .1lf eta=% .2e theta=%.2e\t", k, eta, theta);
      for(m=0;m<=3;m++)
       for(n=0;n<=3;n++)
        {
         if(m+n>3) continue;
         if(theta>DBL_MAX) continue;
         idx = gong_idx[m][n];
         val[m][n] = dfermi200_(&idx,&k, &eta, &theta);
         //printf("% .1e ", (val[m][n]/ref[m][n]-1.0)/DBL_EPSILON);
         printf("%d\t", ULP_distance(ref[m][n], val[m][n], 1024*1024) );  
        }
      printf("\t(GONG)\n\n");   
     }

  }

  printf("Tested:\t%d\tFailed:%d\n", count,failed_count);
  fclose(datafile);

  return 0;

}
