/*
Compile with e.g:
gcc example.c -o example -lm -lfermidirac -lflint
Run:
./example
with expected result: 114.066877991379016066
*/
#include <fermidirac.h>
#include <stdio.h>
#include <math.h>

int main()
{
  
  double result;
  double k=4.0, eta=1.0, theta=1.0;
  const double expected_result = 114.066877991379025899925088247759032653609602513093723942663727131064;
  
  
  result = Ffermi(k,eta,theta);
  
  printf("For k=%lf, eta=%lf, theta=%lf\n\nFfermi(k,eta,theta)=%.18lf\n",k,eta,theta,result);
  printf("Reference result    %.18lf\n\n",expected_result);
  printf("Absolute error:\t%e\nRelative error:\t%e\n",fabs(result-expected_result), fabs(result-expected_result)/expected_result);
  
  // Expected reference result
  // Ffermi(4,1,1)=114.066877991379025899925088247759032653609602513093723942663727131064\
03373794999046746440912671766590817078780479429221800302036
  
  return 0;

}
