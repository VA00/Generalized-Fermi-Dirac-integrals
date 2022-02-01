/*
Compile with e.g:
gcc example_derivatives.c -o example -lm -lfermidirac -lquadmath
Run:
./example
*/
#include <fermidirac.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <quadmath.h>


int main()
{
  
  __float128 k,eta,theta;
  __float128  F[4][4];
    int m,n;
    
  __float128 refVALS[4][4] = {
1.026144878078815325801646213625781936807394239245833482018455316e14q
,
5.10634949358513041450390379328220989346332028762002609668152987e13q
,
-2.541048159149886847302717894514341747376559125734997854287616874e13q
,
3.793473153118717533976794847545316398357858180206253298491676913e13q
,
1.101765802826048209205534617748197222119091561340770109598235977e12q
,
5.487394820102635628620118605336787675657676248534568650029274609e11q
,
-2.733022012497871135906100926719148671181280002692384068956998357e11q
,
0.0q
,
9.6788771890583213736340480434624177732811575538853375335181888e9q
,
4.8247786292290204316907020535873504020674834635754616875943815e9q
,
0.0q
,
0.0q
,
6.613349281010275954773784815940199822922020428842211626459475e7q
,
0.0q
,
0.0q
,
0.0q
};

/*
Reference Values Mathematica CODE

f = Function[{x, k, \[Eta], \[Theta]}, 
   x^k Sqrt[1 + \[Theta] x/2] LogisticSigmoid[\[Eta] - x]];

F[k_, \[Eta]_, \[Theta]_, m_, n_] := Module[{integrand},
  integrand = (Derivative[0, 0, m, n][f])[x, k, \[Eta], \[Theta]];
  NIntegrate[Evaluate@integrand, {x, 0, \[Eta]}, 
    Method -> "GlobalAdaptive", MinRecursion -> 13, 
    MaxRecursion -> 64, WorkingPrecision -> 64] +
   NIntegrate[Evaluate@integrand, {x, \[Eta], Infinity}, 
    Method -> "GlobalAdaptive", MinRecursion -> 11, 
    MaxRecursion -> 29, WorkingPrecision -> 64]
  ]

For[m = 0, m <= 3, m++,
 For[n = 0, n <= 3, n++,
  If[m + n > 3, Print["0.0q"], Print[ToString[CForm@F[4, 512, 1, m, n]] <> "q"]]]
  ]
 ]
 
*/
  
  
  
  
  k     = 4.0q;
  eta   = 512.0q;
  theta = 1.0q;
  
  for(m=0;m<=3;m++)
    for(n=0;n<=3;n++)
      if(m+n>3) 
       continue;
      else
       F[m][n] = Ffermi_derivatives_m_n_quad(k,eta,theta,m,n);

      


  for(m=0;m<=3;m++)
    for(n=0;n<=3;n++)
      if(m+n>3) 
       continue;
      else
       printf("%.32e\t%.32e\t%e\n",(double) F[m][n], (double) refVALS[m][n], (double)( F[m][n]/refVALS[m][n]-1.0q ));
       

  
  

  
  return 0;

}
