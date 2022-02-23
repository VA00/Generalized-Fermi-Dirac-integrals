#define DEBUG 0
/* MAX_REFINE limit recursion depth for FFermi.
Note, that convergence might be slow, and
using very large large MAX_REFINE>16
together with high PRECISION settings
close to DBL_EPSILON results in extremely
slow computations.
*/
#define MAX_REFINE 16 // more than 16 result is significant slow-down
//#define PRECISION sqrt(DBL_EPSILON) // convergence is exponential, so in theory this is enough
#define PRECISION_GOAL 8*DBL_EPSILON   // down to 2*DBL_EPSILON seem harmless, 1*DBL_EPSILON cause problems
#define PRECISION_GOAL_QUAD 8*DBL_EPSILON // 128*FLT128_EPSILON   

#define KAHAN 0 // Enable https://en.wikipedia.org/wiki/Kahan_summation_algorithm ; usually this has no significant effect, but results might be not identical, and computation slow
#define TGAMMA_MAX 170.62437695630272081244437878577 // FindInstance[LogGamma[k + 1] == Log[2^1024] tgamma overflow
// Code is able to compute all derivatives, here we impose some hard-coded limits 
#define DERIVATIVE_MATRIX_SIZE 4 /* Default 4 means derivatives up to $\partial^6 F / \partial^3 \theta \partial^3 \eta$ can be stored */
#define DERIVATIVE_MAX_ORDER   3 /* Default 3 means derivatives up to $\partial^3 F$ are actually computed, remaining matrix entries are unused */


//Standard Fermi-Dirac integrals G-type
double Gfermi(double,double,double);

double Gp(double,double,double);
double Gm(double,double,double);


//Standard Fermi-Dirac integrals F-type
double Ffermi(const double, const double, const double);
long double Ffermi_long(const long double, const long double, const long double);
void   Ffermi_derivatives(const double, const double, const double, double[10]);
void   Ffermi_derivatives_matrix(const double, const double, const double, double[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE]);
    double Ffermi_derivatives_m_n(const double, const double, const double, const int, const int);
__float128 Ffermi_derivatives_m_n_quad(const __float128, const __float128, const __float128, const int, const int);
void   Ffermi_dblexp_derivatives_matrix(const double, const double, const double, const double, const int, double[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE]);
double Ffermi_dblexp_derivatives_m_n(const double, const double, const double, const int, const int, const double, const int);
long double Ffermi_dblexp_derivatives_m_n_long(const long double, const long double, const long double, const int, const int, const long double, const int);
__float128 Ffermi_dblexp_derivatives_m_n_quad(const __float128, const __float128, const __float128, const int, const int, const __float128, const int);


//Fixed-grid version
void fixedFfermi_derivatives(const double, const double, const double,
		   const double, double, double, 
		   double *, double *, double *,
		   double *, double *, double *);
void fixedFfermi_derivatives_v2(const double, const double, const double,
		   const double, double, double, 
		   double *, double *, double *,
		   double *, double *, double *,
		   double *, double *, double *, double *,
		   int
		   );
void fixedFfermi_derivatives_v3(const double, const double, const double,
		   const double, double, double, 
		   double *, double *, double *,
		   double *, double *, double *,
		   double *, double *, double *, double *
		   );
void Ffermi_dblexp_derivatives(const double, const double, const double,
  const double , const int , double[10]);

double fixedFfermi(const double, const double, const double,
		   const double, double, double);
float  fixedFfermif(const float, const float, const float,
		   const float, float, float);
long double  fixedFfermi_long(const long double, const long double, const long double,
		   const long double, const long double, const long double);

double  gaussFfermi(const double, const double, const double);
float  gaussFfermif(const float, const float, const float);


//quick versions - under developement, no accuracy control, vectorized & tabulated
double quickFfermi0(double,double,double);
double quickFfermi1(double,double,double);
double quickFfermi2(double,double,double);
double quickFfermi3(double,double,double);
double quickFfermi4(double,double,double);
//double dfermi(double,double,double);

//internal functions
double r(int, double , double , int , int ); //linear recurrence required for Sommerfeld partial derivatives
double r2(int, double , double , int , int ); //linear recurrence required for Sommerfeld partial derivatives
double sigmoid(double);
long double sigmoid_long(long double);
__float128 sigmoid_quad(__float128);
double sigmoid_derivative(double, int);
double sigmoid_derivative_polynomial(double , int);
long double sigmoid_derivative_polynomial_long(long double , int);
__float128 sigmoid_derivative_polynomial_quad(__float128 , int);
void   sigmoid_derivative_polynomial_vector(double, double[DERIVATIVE_MATRIX_SIZE]);
double g_derivative(double, double, int);
void   g_derivative_vector(double, double, double[DERIVATIVE_MATRIX_SIZE]);
    double fac2(int);
__float128 fac2_quad(int);
    double pochhammer(double, int);
__float128 pochhammer_quad(__float128, int);
double factorial_power(double, int);
double eulerian(int, int);
double power_squaring(double, int);
__float128 power_squaring_quad(__float128, int);
double incomplete_half_fractional_gamma( double, double);
double zeta1( double, int);
double zeta2( double, int);
double zeta3( double, double, int);
    double binom(double, int);
__float128 binom_quad(__float128, int);
double dirichlet_eta(double, double, int);
double hyp1f1(double, double, double);
double hyp2f1(double, double, double, double);
double hyp2f1_series_fd(double, double, double, int);
double hypU_series(double, double, double, double, int);
double hypU_asympt(double, double, double, double, int);
double hypU(double, double, double, double, int);
double U(double, double);
double U_half_frac(double, double);
double BesselK_dbl_exp(const double, const double,  const double, const int);
double BesselK(const double, const double);
double sommerfeld_leading_term(double, double);
__float128 sommerfeld_leading_term_quad(__float128, __float128);
void   sommerfeld_leading_term_derivatives(double, double, double[DERIVATIVE_MATRIX_SIZE]);
void   sommerfeld_leading_term_derivatives_quad(__float128, __float128, __float128[DERIVATIVE_MATRIX_SIZE]);
void   sommerfeld_leading_term_derivatives_matrix(double, double, double, double[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE]);
double sommerfeld_leading_term_derivatives_m_n(const double, const double, const double, const int m, const int n);
__float128 sommerfeld_leading_term_derivatives_m_n_quad(const __float128, const __float128, const __float128, const int m, const int n);
double sommerfeld_leading_term_int(double, double);
double recursion_half_frac_k(double, double);
double recursion_int_k(double, double);
__float128 recursion_half_frac_k_quad(__float128, __float128);
__float128 recursion_int_k_quad(__float128, __float128);
double Ffermi_value         (double, double, double, double, int);
double Ffermi_series_neg    (double, double, double, double, int);
double Ffermi_sommerfeld    (double, double, double, double, int);
void sommerfeld_derivatives(const double, const double, const double, double[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE]);
double sommerfeld_derivatives_m_n(const double, const double, const double, const int, const int);
__float128 sommerfeld_derivatives_m_n_quad(const __float128, const __float128, const __float128, const int, const int);
void   Ffermi_sommerfeld_derivatives(const double, const double, const double, const double, const int, double[10]);
void   Ffermi_sommerfeld_derivatives_matrix(double, double, double, double, int, double[DERIVATIVE_MATRIX_SIZE][DERIVATIVE_MATRIX_SIZE]);
double Ffermi_sommerfeld_derivatives_m_n(const double, const double, const double, const int, const int, const double, const int);
double Ffermi_series_sqrt_a (double, double, double, double, int);
double Ffermi_series_sqrt_b (double, double, double, double, int);


double Ffermi_complete_estimate(double, double, double, double);
long double Ffermi_complete_estimate_long(long double, long double, long double, long double);
double Ffermi_estimate(double, double, double, double, double);
double Gfermi_estimate(double, double, double, double, double);

double Ffermi_complete(const double, const double);
double Ffermi_complete_dbl_exp       (const double, const double, const double, const int);
double Ffermi_complete_series_polylog(const double, const double, const double, const int);
double Ffermi_complete_series_zeta   (const double, const double, const double, const int);
double Ffermi_complete_series_asympt (const double, const double, const double, const int);
long double Ffermi_complete_dblexp_long(const long double, const long double, const long double, const int);
double Ffermi_value(const double, const double, const double, const double, const int);

double integrandF_complete(const double, const double, const double);
double integrandF(const double, const double, const double, const double);
double integrandG(const double, const double, const double, const double);


