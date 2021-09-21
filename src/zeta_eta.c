/*
A. Odrzywolek, andrzej.odrzywolek@uj.edu.pl, 01-06-2020
*/
#include "../fermidirac.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <stdio.h>
#include <float.h>

/* Binomial coefficient Binomial[n,k] */
double binom(int n, int k)
{
	double prod=1.0;
	int i;
	
	for(i=1;i<=k;i++) prod = prod * (1.0+n-i)/i;
	
	return prod;
	
}

/* Dirichlet eta accelerated series */
double dirichlet_eta(double k, double precision, int SERIES_TERMS_MAX)
{
	double sum_i=0.0,sum=0.0,sum_old=0.0,sum_new=0.0;
	int i,j;
	
	/*
	for(j=0;j<=SERIES_TERMS_MAX;j++)
	{
		sum_i=0.0;
		for(i=0;i<=j;i++)
			sum_i += ( i % 2 == 0 ) ?  binom(j,i)*pow(1.0+i,-k) : -binom(j,i)*pow(1.0+i,-k);
		
		sum = sum + sum_i*pow(2.0,-1.0-j);
	}
	*/
	
	j=0;
	
	do
	{
		sum_old=sum_new;
		sum_i=0.0;
		for(i=0;i<=j;i++)
			sum_i += ( i % 2 == 0 ) ?  binom(j,i)*pow(1.0+i,-k) : -binom(j,i)*pow(1.0+i,-k);
		
		sum_new = sum_old + sum_i*pow(2.0,-1.0-j);
				
		j++;
		
	}
    while ( ( (precision>0.0) ? fabs(sum_old-sum_new)>=precision*sum_new : sum_old!=sum_new )  && j<SERIES_TERMS_MAX );
		
		
	return sum_new;

}



/* Some simple estimates for Riemmann Zeta function */
double zeta1( const double k, const int SERIES_TERM_GOAL)
{
	int i;
	double sum=0.0;
	
	for(i=1;i<=SERIES_TERM_GOAL;i++)
		sum += pow(i,-k);
	
	return sum;
	
}

double zeta2( const double k, const int SERIES_TERM_GOAL)
{
	const double prime[12]={2., 3., 5., 7., 11., 13., 17., 19., 23., 29., 31., 37.};
	double prod=1.0;
	int i;
	
	for(i=1;(i<SERIES_TERM_GOAL)&&(i<12);i++)
		prod = prod/(1.0-pow(prime[i-1],-k));
	
	return prod;
	
}


double zeta3( const double k, const double precision, const int SERIES_TERMS_MAX)
{
	const unsigned char prime[32]={2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, \
                                   67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131};
	double prod_old,prod_new=1.0;
	int i=0;
	
  do
   {
	i++;
	prod_old = prod_new;
    prod_new = prod_new/(1.0-pow(prime[i-1],-k));
   }
  while ( ( (precision>0.0) ? fabs(prod_old-prod_new)>=precision*prod_new : prod_old!=prod_new )  && i<SERIES_TERMS_MAX && i<32);

  return prod_new;
	
}


