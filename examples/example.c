#include <fermidirac.h>
#include <stdio.h>
#include <math.h>

int main()
{
  
  double result;
  
  
  result = Ffermi(4.0,1.0,1.0);
  
  printf("%.18lf\n",result);
  
  // Expected reference result
  // Ffermi(4,1,1)=114.066877991379025899925088247759032653609602513093723942663727131064\
03373794999046746440912671766590817078780479429221800302036
  
  return 0;

}
