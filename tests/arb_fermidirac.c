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
int
f_generalized_relativistic_fermi_dirac(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t t,sigmoid, sqrt_term;
    acb_t k, eta, theta;

    acb_init(k);
    acb_init(eta);
    acb_init(theta);

    acb_set_d(k,0.5);    // k=0.5
    acb_set_d(eta, 2097152);  // eta = 2^21
    acb_set_d(theta,3.1082702275611665134711390509176302506278509424834232340028998555822468563283335970816e85);  // theta=2^285  ; divided by 2 in advance

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_init(t);
    acb_init(sigmoid);
    
    //Compute logistic sigmoid
    acb_neg(eta, eta); //Zmiana znaku
    acb_add(t,z,eta,prec);
    acb_exp(t, t, prec);
    acb_add_ui(res, t, 1, prec);
    acb_inv(sigmoid, res, prec);

    //Compute z^k
    acb_pow_analytic(res, z, k, order != 0, prec);


    //Compute Sqrt[1+theta*z/2]
    acb_init(sqrt_term);
    acb_mul(sqrt_term, theta, z, prec);
    acb_add_ui(sqrt_term, sqrt_term, 1, prec);
    acb_sqrt(sqrt_term, sqrt_term, prec);
    //acb_sqrt_analytic(sqrt_term, sqrt_term, order != 0, prec);
    

    //Multiply three terms x^k*Sqrt*sigmoid
    acb_mul(res, res, sigmoid, prec);
    acb_mul(res, res, sqrt_term, prec);


    
    acb_clear(t);
    acb_clear(sigmoid);
    acb_clear(k);
    acb_clear(eta);
    acb_clear(theta);
    acb_clear(sqrt_term);

    return 0;
}




int main(int argc, char *argv[])
{
    acb_t s, t, a, b;
    mag_t tol;
    slong prec, goal;
    slong N;
    ulong k;
    int ifrom, ito;
    int i, twice, havegoal, havetol;
    acb_calc_integrate_opt_t options;

    ifrom = ito = -1;


    acb_calc_integrate_opt_init(options);

    prec = 128;
    twice = 0;
    goal = 0;
    havetol = havegoal = 0;

    acb_init(a);
    acb_init(b);
    acb_init(s);
    acb_init(t);
    mag_init(tol);

    if (!havegoal)
        goal = prec;

    if (!havetol)
        mag_set_ui_2exp_si(tol, 1, -prec);

    //TIMEIT_ONCE_START
 
    if (goal < 0) abort();

    /* error bound (N+1) exp(-N) when truncated at N */
    N = goal + FLINT_BIT_COUNT(goal)+2097152; // eta added !!!!!
    acb_zero(a);
    acb_set_ui(b, N);
    acb_calc_integrate(s, f_generalized_relativistic_fermi_dirac, NULL, a, b, goal, tol, options, prec);
    acb_neg(b, b);
    acb_exp(b, b, prec);
    acb_mul_ui(b, b, N + 1, prec);
    arb_add_error(acb_realref(s), acb_realref(b));
    flint_printf("I = ");
    acb_printn(s, 3.333 * prec, 0);
    flint_printf("\n\n");

    acb_clear(a);
    acb_clear(b);
    acb_clear(s);
    acb_clear(t);
    mag_clear(tol);

    flint_cleanup();
    return 0;
}
