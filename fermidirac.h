//Standard Fermi-Dirac integrals G-type
double Gfermi(double,double,double);

double Gp(double,double,double);
double Gm(double,double,double);


//Standard Fermi-Dirac integrals F-type
double Ffermi(const double, const double, const double);
long double Ffermi_long(const long double, const long double, const long double);


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
void Ffermi_value_derivatives(const double, const double, const double,
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
double sigmoid(double);
double incomplete_half_fractional_gamma( double, double);
double zeta1( double, int);
double zeta2( double, int);
double zeta3( double, double, int);
double binom(int, int);
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
void   sommerfeld_leading_term_derivatives(double, double, double[4]);
double sommerfeld_leading_term_int(double, double);
double recursion_half_frac_k(double, double);
double recursion_int_k(double, double);
double Ffermi_value         (double, double, double, double, int);
double Ffermi_series_neg    (double, double, double, double, int);
double Ffermi_sommerfeld    (double, double, double, double, int);
void   Ffermi_sommerfeld_derivatives(const double, const double, const double, const double, const int, double[10]);
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
long double Ffermi_complete_value_long(const long double, const long double, const long double, const int);
double Ffermi_value(const double, const double, const double, const double, const int);

double integrandF_complete(const double, const double, const double);
double integrandF(const double, const double, const double, const double);
double integrandG(const double, const double, const double, const double);


