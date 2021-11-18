#!/usr/bin/env python3
"""
Helper script to generate the compiled FFI interface for usage of C sources for
the generalized Fermi-Dirac functions using the tanh-sinh quadrature.

Can generate either api or abi level package.
"""

import glob
import argparse
from cffi import FFI

CDEF = """
    double Gfermi(double, double, double);
    double Gp(double, double, double);
    double Gm(double, double, double);

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

    double fixedFfermi(const double, const double, const double,
        const double, double, double);
    float fixedFfermif(const float, const float, const float,
        const float, float, float);
    long double fixedFfermi_long(const long double, const long double, const long double,
        const long double, const long double, const long double);

    double gaussFfermi(const double, const double, const double);
    float gaussFfermif(const float, const float, const float);
"""

# definition of functions we want to use
ffibuilder = FFI()
ffibuilder.cdef(CDEF)

def parse_arguments():
    """Argument parsing from command line, using argparse.

    :rtype: Any
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--api",
        help="set compilation in api mode (faster, but requires compiler and sources)",
        action="store_true"
    )

    return parser.parse_args()

if __name__ == "__main__":

    args = parse_arguments()
    if args.api:
        ffibuilder.set_source(
            "_fermidirac_api",  # name of the output C extension
            """
                #include "fermidirac.h"
            """,
            sources=glob.glob("../../src/*.c"),
            libraries=['m']
        )
    else:
        ffibuilder.set_source(
            "_fermidirac_abi",  # name of the output C extension
            """
                #include "fermidirac.h"
            """,
            libraries=['m', 'fermidirac']
        )



    ffibuilder.compile(verbose=True)
