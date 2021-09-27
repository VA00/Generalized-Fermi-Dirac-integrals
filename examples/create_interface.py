import glob
from cffi import FFI
ffibuilder = FFI()

ffibuilder.cdef("""
                    double Ffermi(const double, const double, const double);
                    double Gp(const double, const double, const double);
                    double Gm(const double, const double, const double);
                    double Gfermi(const double, const double, const double);
                """)

ffibuilder.set_source("fermidirac",  # name of the output C extension
"""
    #include "fermidirac.h"
""",
    sources=glob.glob("../src/*.c"),
    libraries=['m'])    # on Unix, link with the math library

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)

