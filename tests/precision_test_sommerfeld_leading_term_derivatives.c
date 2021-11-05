#include <fermidirac.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
//#define TESTDATA "refVALS/20211105/refVALS_S_20211105.h"
#include "refVALS/20211105/refVALS_S_rnd_20211105.h"

int main()
{
  
  double k,z;

  double S,S1,S2,S3;

  double eq3_LHS,eq3_RHS;
  double x[13], results[4] = { [ 0 ... 3 ] = -1.0 };
  int ii;
  
    FILE  *datafile;
  
  datafile = fopen("refVALS/20211105/refVALS_S_rnd_20211105.bin","r");

  
  if(datafile==NULL){ printf("ERROR: UNABLE TO OPEN FILE\n");return -1;}
  /*
   k= 1.0;
   z=-1.0;
   sommerfeld_leading_term_derivatives(k,z,results);
   printf("%.3e %.3e %.17e %.17e %.17e %.17e\n",k,z,results[0],results[1],results[2],results[3]);
*/
  printf("%d\n",num);

  for(ii=0;ii<num;ii=ii+1)
   {
	fread(x, 8, 6, datafile);

    k               = x[0];
	z               = x[1];
	S               = x[2];
	S1              = x[3];
    S2              = x[4];
    S3              = x[5]; 
    
    printf("%.3e %.3e %.17e %.17e %.17e %.17e\n",k,z,S,S1,S2,S3);
    sommerfeld_leading_term_derivatives(k,z,results);
    printf("%.3e %.3e %.17e %.17e %.17e %.17e\n",k,z,results[0],results[1],results[2],results[3]);
    printf("\n");
   }

  
  

  fclose(datafile);
  return 0;

}
