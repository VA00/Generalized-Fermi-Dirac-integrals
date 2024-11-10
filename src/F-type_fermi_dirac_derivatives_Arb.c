/*
Based on: https://github.com/fredrik-johansson/arb/blob/master/examples/integrals.c

A. Odrzywołek, 2022-02-16, andrzej.odrzywolek@uj.edu.pl

Install Arb https://arblib.org/

sudo apt install libflint
sudo apt install libflint-dev


Compile with: gcc arb_fermidirac.c -o arb_fd -lflint -lflint-arb -lm

*/

#include <string.h>
#include <math.h>
#include <arb.h>
#include <arb_calc.h>
#include <acb.h>
#include <acb_calc.h>
#include <acb_hypgeom.h>
//#include "double_interval.h" // Require the most recent Arb, I'm unable to install it A.O. 
#include <fermidirac.h>
#include <quadmath.h>
#include "eulerian.h"




/* f(z) = z^k * sqrt(1+0.5*theta*z) / (1 + exp(z-eta)) */
int f_generalized_relativistic_fermi_dirac_integrand(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    
    double k_input, eta_input, theta_input;
    acb_t sigmoid, g;
    acb_t k, eta, theta;

    acb_init(k);
    acb_init(eta);
    acb_init(theta);

    //k_input=0.5;
    //eta_input = pow(2,22);
    //theta_input = pow(2,-50);
    k_input=((double *) param)[0];
    eta_input = ((double *) param)[1];
    theta_input = ((double *) param)[2];
    acb_set_d(k, k_input);
    acb_set_d(eta, eta_input); 
    acb_set_d(theta,theta_input);

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_init(g);
    acb_init(sigmoid);
    
    //Compute logistic sigmoid
    acb_neg(eta, eta); //Zmiana znaku
    acb_add(sigmoid,z,eta,prec);
    acb_exp(sigmoid, sigmoid, prec);
    acb_add_ui(sigmoid, sigmoid, 1, prec);
    acb_inv(sigmoid, sigmoid, prec);

    //Compute z^k
    //acb_pow_analytic(res, z, k, order != 0, prec);
    acb_pow(res, z, k, prec);



    //Compute Sqrt[1+theta*z/2]
    acb_mul(g, theta, z, prec);
    acb_div_ui(g, g, 2, prec);
    acb_add_ui(g, g, 1, prec);
    acb_sqrt(g, g, prec);
    //acb_sqrt_analytic(sqrt_term, sqrt_term, order != 0, prec);
    

    //Multiply three terms x^k*Sqrt*sigmoid
    acb_mul(res, res, sigmoid, prec);
    acb_mul(res, res, g, prec);


    
    acb_clear(g);
    acb_clear(sigmoid);
    acb_clear(k);
    acb_clear(eta);
    acb_clear(theta);
    return 0;
}



int f_generalized_relativistic_fermi_dirac_integrand_m_n(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
  double k_input, eta_input, theta_input, m_input, n_input;
  acb_t s,ds, g, g2, pochhammer, half;
  acb_t k, eta, theta;
  int j;

  acb_init(k);
  acb_init(eta);
  acb_init(theta);

  k_input=((double *) param)[0];
  eta_input = ((double *) param)[1];
  theta_input = ((double *) param)[2];
  m_input =  ((double *) param)[3];
  n_input =  ((double *) param)[4];
  int m = (int) m_input;
  int n = (int) n_input;

  acb_set_d(k, k_input); //acb_one(k);
  acb_set_d(eta, eta_input); 
  acb_set_d(theta,theta_input);


  if (order > 1)
      flint_abort();  /* Would be needed for Taylor method. */

  acb_init(g);
  acb_init(g2);
  acb_init(pochhammer);
  acb_init(half);
  acb_init(s);


  
   //Compute logistic sigmoid
   acb_neg(eta, eta); //Zmiana znaku
   acb_add(s,z,eta,prec);
   acb_exp(s, s, prec);
   acb_add_ui(s, s, 1, prec);
   acb_inv(s, s, prec);

  //Compute logistic sigmoid dervatives
  //ds=sigmoid_derivative_polynomial(s, m);
  if(m>0)
   {
    acb_t sum,term1,term2,s1;
    acb_init(sum);
    acb_init(term1);
    acb_init(term2);
    acb_zero(sum);
    acb_init(s1);
    
    for(j=0;j<=m-1;j++)  
     {
     
     //sum = sum + eulerian(m,j)*power_squaring(s,j)*power_squaring(s - 1.0,m - j);
     acb_pow_ui(term1,s,j,prec); // term1= s^j
     acb_add_si(s1,s,-1,prec);    // s1=s-1
     acb_pow_ui(term2,s1,m-j,prec); // term2= (s-1)^(m-j)
     acb_mul_ui(term1,term1,eulerian_numbers[m][j],prec);
     acb_mul(term2, term1, term2, prec);
     acb_add(sum, sum, term2, prec);
     } 
    if(m%2) acb_neg(sum, sum);
    acb_mul(s,sum,s,prec);
    
    acb_clear(sum);
    acb_clear(term1);
    acb_clear(term2);
    acb_clear(s1);
   }

  //Compute z^k
  acb_pow_analytic(res, z, k, order != 0, prec);
  //acb_pow(res, z, k, prec);
  //acb_one(res);

  //Compute g derivatives g (-1)^n Pochhammer[-1/2,n] (z/2/g^2)^n 
  //dg=g_derivative(x, theta, n);
  //Compute Sqrt[1+theta*z/2]
  acb_mul(g, theta, z, prec);
  acb_div_ui(g, g, 2, prec);
  acb_add_ui(g2, g, 1, prec);
  acb_sqrt(g, g2, prec);
  //compute (z/2/g^2)
  acb_inv(g2,g2,prec);
  acb_div_ui(g2, g2, 4, prec);
  acb_mul(g2,g2,z,prec);
  //compute (z/2/g^2)^n 
  //acb_pow_ui(g2,g2,n,prec);
  //acb_one(half);
  //acb_rising_ui(pochhammer, half, n, prec); //Not implemented in version arb-2.17-0, RENAMED acb_rising
  for(j=1;j<=n;j++)
   {
    acb_mul(g,g,g2,prec);
    acb_mul_si(g,g,3-2*j,prec);
   }

  //Multiply three terms x^k*Sqrt*sigmoid
  acb_mul(res, res, s, prec);
  acb_mul(res, res, g, prec);

  acb_clear(g);
  acb_clear(g2);
  acb_clear(pochhammer);
  acb_clear(half);
  acb_clear(s);
  acb_clear(k);
  acb_clear(eta);
  acb_clear(theta);

  return 0;
}

void Ffermi_derivatives_m_n_arb(acb_t s, const double k, const double eta, const double theta, const int m, const int n)
{
  acb_t a, b;
  mag_t tol;
  slong prec;
  
  double input_parameters[5]={k, eta, theta, (double) m, (double) n};
  int fail;
  double cutoff;

  acb_calc_integrate_opt_t options;
  acb_calc_integrate_opt_init(options);

  options->verbose=0;

  prec = 256;

  acb_init(a);
  acb_init(b);
  mag_init(tol);
  mag_set_ui_2exp_si(tol, 1, -prec); // tol = 1*2^-prec


  if( ( (k-floor(k)==0.5) || (k-floor(k)==0.0) ) && (k>=0.0) )  //half-frac k=1/2,3/2,5/2,... OR integer k=0,1,2,3,4,5,...
    acb_zero(a); // We are able start integration at ZERO
  else
    acb_set_d(a,pow(2.0,-256)); // Unable to integrate near (complex) ZERO, branches?
// TODO: Replace ABOVE with something more serious


  //cutoff = prec + FLINT_BIT_COUNT(prec)+fmax(eta,0.0); // eta added !!!!!
  acb_set_d(b,  fmax(eta,0.0) );
  //acb_set_d(b,  cutoff);
  acb_add_ui(b,b,prec,prec);
  acb_add_ui(b,b,FLINT_BIT_COUNT(prec),prec);

  fail = acb_calc_integrate(s, f_generalized_relativistic_fermi_dirac_integrand_m_n, input_parameters, a, b, prec, tol, options, prec);
  if(fail==ARB_CALC_NO_CONVERGENCE) printf("\nArb integration failed!\n");



  acb_clear(a);
  acb_clear(b);
  mag_clear(tol);
  flint_cleanup();
}

double Ffermi_derivatives_m_n_internal_arb(const double k, const double eta, const double theta, const int m, const int n)
{
  acb_t s;
  double result;
  
  acb_init(s);
  Ffermi_derivatives_m_n_arb(s, k, eta, theta, m, n);
  
  //double output
  arf_t t;
  arf_init(t);
  arb_t res;
  arb_init(res);
 
  acb_get_real(res, s);
  result = arf_get_d( arb_midref(res) , ARF_RND_NEAR);

  arf_clear(t);
  arb_clear(res);
  acb_clear(s);
  
  return result;
}

