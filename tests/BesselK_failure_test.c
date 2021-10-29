#include <fermidirac.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <float.h>
#include "refVALS/refVALS_K_rnd_20211028.h"

int main()
{
  
  
  double k, z,  fd, fdREF, fdTMP, F1, F2;
  long double kL,etaL, fdL;
  double maxABSerror, minABSerror, maxRELerror, minRELerror;
  double relError, absError;

  
  int ii,ULP,ULP_MAX=0,ULP_AVG=0;

  double h,hMAX;
  double x[3];
  
  FILE  *datafile, *datafile2;
  
  datafile = fopen("refVALS/refVALS_K_rnd_20211028.bin","r");

  
  if(datafile==NULL){ printf("ERROR: UNABLE TO OPEN FILE\n");return -1;}
  

  
  maxABSerror = 0.0;
  minABSerror = 0.0; 
  maxRELerror = 0.0;
  minRELerror = 0.0;
  
  
  for(ii=0;ii<num;ii=ii+1)
  {
	fread(x, 8, 3, datafile);
    //printf("%.16e,%.16e,%.16e\n",x[0],x[1],x[2]);
    k   = x[0];
	z = x[1];
	fdREF = x[2];
    
	//if(fdREF<0.0) continue;
	
    //if( (k!=0.0) || (z<1.0) ) continue;
    if( ( (k!=0.0) && (k!=1.0) ) || (fdREF<1e-300) ) continue;

  //fd     = BesselK_dbl_exp(k, z,0.0,16);
  fd     = BesselK(k, z);
  //fd = BesselK0(z);
  
  fdTMP  = fd;
  
  x[2] = fd;
  //fwrite(x,8, 3, datafile2);
  
   //if(isinf(fd)) continue;
  
  absError = fd-fdREF;
  if( (fdREF==0.0) || isinf(fdREF) ) relError=0.0; else relError = (double) ( ((long double) fd)/((long double) fdREF)-1.0L );
  
  ULP=0;
  
  if( (fdREF!=0.0) && (!isinf(fdREF)) )
  {
    while( (fdREF!=fdTMP) && abs(ULP) <256 ){ ULP++; fdTMP=nextafter(fdTMP,fdREF);}
  }
  
  ULP_AVG += ULP;
  
  
    
  if(absError>maxABSerror) maxABSerror = absError;

  if(absError<minABSerror) minABSerror = absError;

  if(relError>maxRELerror)
   { 
    maxRELerror = relError;printf("%d\tk=%.16lf\t\tz=%.20e\t\t%.20e\t%.20e\t%+.3lf\t%dULP\n", ii+1, k, z, fd, fdREF, relError/DBL_EPSILON, ULP);
   }

  if(relError<minRELerror)
   { 
    minRELerror = relError;printf("%d\tk=%.16lf\t\tz=%.20e\t\t%.20e\t%.20e\t%+.3lf\t%dULP\n", ii+1, k, z, fd, fdREF, relError/DBL_EPSILON, ULP);
   }
  
  if(ULP_MAX<=ULP) 
   {  
       ULP_MAX=ULP;
       //if(ULP>128)  
       //printf("%.18e\t%d\n", log(fdREF), ULP);
   };
	 
	// printf("%d\tk=%.16lf\t\tz=%.20e\t\t%.20e\t%.20e\t%+.3lf\t%dULP\n", ii+1, k, z, fd, fdREF, relError/DBL_EPSILON, ULP);


  }
  
    //printf("DBL_EPSILON = %e\tDBL_MAX=%e\tDBL_MIN=%e\n\n\n",DBL_EPSILON,DBL_MAX,DBL_MIN);
  
  printf("Maximum rel. err. [DBL_EPS]: %.2e\t[FLT_EPS]: %.2e\n", maxRELerror/DBL_EPSILON, maxRELerror/FLT_EPSILON);
  printf("Minimum rel. err. [DBL_EPS]: %.2e\t[FLT_EPS]: %.2e\n", minRELerror/DBL_EPSILON, minRELerror/FLT_EPSILON);
  printf("Maximum ULP's:\t%d\n", ULP_MAX);
  printf("AVG ULP's:\t%lf\t%d\n", ((double) ULP_AVG)/num, ii);
  
  
  
  
 
    fclose(datafile); 
	//fclose(datafile2);
	  
	
	
  return 0;

}
