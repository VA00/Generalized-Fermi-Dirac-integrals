/*
Based on: https://github.com/fredrik-johansson/arb/blob/master/examples/integrals.c

A. Odrzywołek, 2022-02-16, andrzej.odrzywolek@uj.edu.pl

Install Arb https://arblib.org/

sudo apt install libflint-arb-dev
sudo apt install libflint-arb2
sudo apt install libflint-dev


Compile with: gcc arb_fermidirac.c -o arb_fd -lflint -lflint-arb -lm

*/

#include <string.h>
#include <acb_calc.h>
#include <arb.h>
//#include "double_interval.h" // Require the most recent Arb, I'm unable to install it A.O. 
#include <fermidirac.h>
#include <quadmath.h>



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
    acb_pow_analytic(res, z, k, order != 0, prec);


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

void Ffermi_derivatives_m_n_arb(acb_t s, const double k, const double eta, const double theta, const int m, const int n)
{
  acb_t a, b;
  mag_t tol;
  slong prec;
  //slong N;
  double N;
  double input_parameters[3]={k, eta, theta};

  acb_calc_integrate_opt_t options;
  acb_calc_integrate_opt_init(options);

  prec = 128;

  acb_init(a);
  acb_init(b);
  mag_init(tol);
 
  mag_set_ui_2exp_si(tol, 1, -prec); // tol = 1*2^-prec

  /* error bound (N+1) exp(-N) when truncated at N */
  N = prec + FLINT_BIT_COUNT(prec)+fmax(eta,0.0); // eta added !!!!!
  //flint_printf("N = %wd \n",N);
  //flint_printf("MIN = %wd \n",WORD_MIN);
  //flint_printf("MAX = %wd \n",WORD_MAX);
  //flint_printf("MAX = %wu \n",UWORD_MAX);
  acb_zero(a);
  //acb_set_ui(b, N);
  acb_set_d(b, N);
  acb_calc_integrate(s, f_generalized_relativistic_fermi_dirac_integrand, input_parameters, a, b, prec, tol, options, prec);

  acb_clear(a);
  acb_clear(b);
  mag_clear(tol);
  flint_cleanup();
}




int main(int argc, char *argv[])
{
  double fd_quad;
  acb_t fd_arb, fd;
  arb_t rel_err, MachineEpsilon;

  int i,j, sign;
  //di_t interval; Require the most recent Arb, I'm unable to install it A.O. 
  

  flint_printf("Computed with arb-%s\n", arb_version);

  acb_init(fd_arb);
  acb_init(fd);
  arb_init(rel_err);
  arb_init(MachineEpsilon);
  arb_one(MachineEpsilon);
  arb_mul_2exp_si(MachineEpsilon, MachineEpsilon, -52); 

  for(sign=-1;sign<=1;sign=sign+2)
   for(i=-1022;i<=1024;i=i+66)
    for(j=-1022;j<=1024;j=j+66)
     {

       Ffermi_derivatives_m_n_arb(fd_arb, 0.5, sign*pow(2,i), pow(2,j), 0, 0);    
       
       fd_quad = (double) Ffermi_derivatives_m_n_quad(0.5q, sign*powq(2,i), powq(2,j), 0, 0);
     
       acb_set_d(fd, fd_quad);
     
       acb_div(fd, fd, fd_arb, 128);
       acb_add_si(fd, fd, -1, 128);
     

     
       acb_abs(rel_err, fd, 128);
       if( arb_ge(rel_err, MachineEpsilon) )
       {
       printf("\nArb=\t");
       //acb_print(fd_arb);printf("\n");
       //acb_printd(fd_arb, 128);printf("\n");
       acb_printn(fd_arb, 128, 0);printf("\n");
       printf("libfermidirac=%.18e",  fd_quad);     
       printf("\n\n");
     
       printf("%d %d\t%d\t",sign,i,j); 
       arb_printn(rel_err, 128, 0);
       //arb_printn(MachineEpsilon, 128, 0);
       printf("\n------------------------------------------------------------------------------\n\n");
       }



     }

  acb_clear(fd_arb);
  acb_clear(fd);
  arb_clear(rel_err);
  arb_clear(MachineEpsilon);

  return 0;
}
