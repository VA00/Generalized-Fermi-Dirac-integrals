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
#include "flint/profiler.h"
//#include "arb_hypgeom.h"
//#include "acb_hypgeom.h"
//#include "acb_dirichlet.h"
//#include "acb_modular.h"
#include "acb_calc.h"



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

void Ffermi_derivatives_m_n_arb(acb_ptr result, const double k, const double eta, const double theta, const int m, const int n)
{

    acb_t s, t, a, b;
    mag_t tol;
    slong prec, goal;
    slong N;
    //ulong k;
    int ifrom, ito;
    int i, twice, havegoal, havetol;
    acb_calc_integrate_opt_t options;

    ifrom = ito = -1;


    acb_calc_integrate_opt_init(options);

    prec = 128;
    twice = 0;
    goal = 0;
    havetol = havegoal = 0;
    double input_parameters[3]={k, eta, theta};

    acb_init(a);
    acb_init(b);
    acb_init(s);
    acb_init(t);
    mag_init(tol);

    if (!havegoal)
        goal = prec;

    if (!havetol)
        mag_set_ui_2exp_si(tol, 1, -prec);

 
    if (goal < 0) abort();

    /* error bound (N+1) exp(-N) when truncated at N */
    N = goal + FLINT_BIT_COUNT(goal)+pow(2,22); // eta added !!!!!
    acb_zero(a);
    acb_set_ui(b, N);
    acb_calc_integrate(s, f_generalized_relativistic_fermi_dirac_integrand, input_parameters, a, b, goal, tol, options, prec);
    acb_neg(b, b);
    acb_exp(b, b, prec);
    acb_mul_ui(b, b, N + 1, prec);
    arb_add_error(acb_realref(s), acb_realref(b));
    flint_printf("I = ");
    acb_printn(s, 3.333 * prec, 0);
    flint_printf("\n");

    acb_clear(a);
    acb_clear(b);
    acb_clear(s);
    acb_clear(t);
    mag_clear(tol);

    flint_cleanup();

}




int main(int argc, char *argv[])
{
    acb_t s, t, a, b;
    mag_t tol;
    slong prec, goal;
    slong N;
    //ulong k;
    int ifrom, ito;
    int i, twice, havegoal, havetol;
    acb_calc_integrate_opt_t options;

    ifrom = ito = -1;


    acb_calc_integrate_opt_init(options);

    prec = 128;
    twice = 0;
    goal = 0;
    havetol = havegoal = 0;
    double input_parameters[3]={0.5,pow(2,22),pow(2.0,-50)};

    acb_init(a);
    acb_init(b);
    acb_init(s);
    acb_init(t);
    mag_init(tol);

    if (!havegoal)
        goal = prec;

    if (!havetol)
        mag_set_ui_2exp_si(tol, 1, -prec);

 
    if (goal < 0) abort();

    /* error bound (N+1) exp(-N) when truncated at N */
    N = goal + FLINT_BIT_COUNT(goal)+pow(2,22); // eta added !!!!!
    acb_zero(a);
    acb_set_ui(b, N);
    acb_calc_integrate(s, f_generalized_relativistic_fermi_dirac_integrand, input_parameters, a, b, goal, tol, options, prec);
    acb_neg(b, b);
    acb_exp(b, b, prec);
    acb_mul_ui(b, b, N + 1, prec);
    arb_add_error(acb_realref(s), acb_realref(b));
    flint_printf("I = ");
    acb_printn(s, 3.333 * prec, 0);
    flint_printf("\n");

    printf("\n\n%e\n\n", (double) *s);

    acb_clear(a);
    acb_clear(b);
    acb_clear(s);
    acb_clear(t);
    mag_clear(tol);

    flint_cleanup();

    Ffermi_derivatives_m_n_arb(NULL, 0.5, pow(2,22), pow(2,-50), 0, 0);    

    return 0;
}
