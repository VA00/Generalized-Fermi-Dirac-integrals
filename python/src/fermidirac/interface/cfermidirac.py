#!/usr/bin/env python3
"""
Wrapper file for the api and abi versions of the library.
"""

import sys
from fermidirac.helpers.Derivatives import Derivatives
from cffi import FFI
ffi = FFI()
ffi_stdlib = FFI()

ffi_stdlib.cdef("int snprintf ( char * s, size_t n, const char * format, ... );")
_libc = ffi_stdlib.dlopen("c")
try:
    import _fermidirac_api.lib as _fermi
except ImportError as exc:
    from fermidirac.interface.build_cffi_interface import CDEF
    try:
        ffi.cdef(CDEF)
        _fermi = ffi.dlopen("fermidirac")
    except OSError as exc:
        sys.exit(f"Library not found ({exc})!")

from mpmath import mp
mp.dps = 128

def extract_long_double(ctype_long):
    """Function to extract long double into a high-precision mp.mpf object,
    required as python and ctypes don't support using long double.

    :param Any ctype_long: ctype containing long to extract.
    :return: the extracted long double at 128 digit max precision.
    :rtype: Any
    """
    size = 1024
    buf = ffi_stdlib.new("char[]", size)
    _libc.snprintf(buf, size, b'%.128Le', ctype_long)
    return mp.mpf(ffi_stdlib.string(buf).decode('ascii'))

def Ffermi(k, eta, theta):
    """The generalized F-type Fermi integral quadrature.

    :param float k: The k parameter, indicating the power of x in the numerator of the integral.
    :param float eta: The eta parameter, indicating the change to the exponent in the denominator of the integral.
    :param float theta: The theta parameter, indicating the scaling of the x parameter in the square root.
    :return: The result of the sinh-tanh quadrature.
    :rtype: float
    """
    return _fermi.Ffermi(k, eta, theta)

def fixedFfermi_derivatives(k, eta, theta, h, hmin, hmax):
    """The generalized F-type Fermi integral quadrature, using fixed abscissas and weights.

    :param float k: The k parameter, indicating the power of x in the numerator of the integral.
    :param float eta: The eta parameter, indicating the change to the exponent in the denominator of the integral.
    :param float theta: The theta parameter, indicating the scaling of the x parameter in the square root.
    :param float h: The h parameter, indicating step size.
    :param float hmin: The hmin parameter, indicating lower bound of quadrature to start the stepping on.
    :param float hmax: The hmax parameter, indicating upper bound to end the stepping on.
    :return: The result of the sinh-tanh quadrature.
    :rtype: float
    """
    f = ffi.new("double *")
    df_deta = ffi.new("double *")
    d2f_deta2 = ffi.new("double *")
    df_dtheta = ffi.new("double *")
    d2f_dtheta2 = ffi.new("double *")
    d2f_deta_dtheta = ffi.new("double *")
    _fermi.fixedFfermi_derivatives(k, eta, theta, h, hmin, hmax, f, df_deta, d2f_deta2, df_dtheta, d2f_dtheta2, d2f_deta_dtheta)
    return Derivatives(
        f=f[0],
        df_deta=df_deta[0],
        df_dtheta=df_dtheta[0],
        d2f_deta2=d2f_deta2[0],
        d2f_dtheta2=d2f_dtheta2[0],
        d2f_deta_dtheta=d2f_deta_dtheta[0]
    )

def fixedFfermi(k, eta, theta, h, hmin, hmax):
    """The generalized F-type Fermi integral quadrature, using fixed abscissas and weights.

    :param float k: The k parameter, indicating the power of x in the numerator of the integral.
    :param float eta: The eta parameter, indicating the change to the exponent in the denominator of the integral.
    :param float theta: The theta parameter, indicating the scaling of the x parameter in the square root.
    :param float h: The h parameter, indicating step size.
    :param float hmin: The hmin parameter, indicating lower bound of quadrature to start the stepping on.
    :param float hmax: The hmax parameter, indicating upper bound to end the stepping on.
    :return: The result of the sinh-tanh quadrature.
    :rtype: float
    """
    return _fermi.fixedFfermi(k, eta, theta, h, hmin, hmax)

def Ffermi_long(k, eta, theta):
    """The generalized F-type Fermi integral quadrature, with long precision.

    :param float k: The k parameter, indicating the power of x in the numerator of the integral.
    :param float eta: The eta parameter, indicating the change to the exponent in the denominator of the integral.
    :param float theta: The theta parameter, indicating the scaling of the x parameter in the square root.
    :return: The result of the sinh-tanh quadrature.
    :rtype: float
    """
    return extract_long_double(_fermi.Ffermi_long(k, eta, theta))

def fixedFfermi_long(k, eta, theta, h, hmin, hmax):
    """The generalized F-type Fermi integral quadrature, using fixed abscissas and weights, with long precision.

    :param float k: The k parameter, indicating the power of x in the numerator of the integral.
    :param float eta: The eta parameter, indicating the change to the exponent in the denominator of the integral.
    :param float theta: The theta parameter, indicating the scaling of the x parameter in the square root.
    :param float h: The h parameter, indicating step size.
    :param float hmin: The hmin parameter, indicating lower bound of quadrature to start the stepping on.
    :param float hmax: The hmax parameter, indicating upper bound to end the stepping on.
    :return: The result of the sinh-tanh quadrature.
    :rtype: float
    """
    return extract_long_double(_fermi.fixedFfermi_long(k, eta, theta, h, hmin, hmax))

def gaussFfermi(k, eta, theta):
    """Gauss subdivision F-type Fermi function, for checking results (based on code by F.X.Timmes).

    :param float k: The k parameter, indicating the power of x in the numerator of the integral.
    :param float eta: The eta parameter, indicating the change to the exponent in the denominator of the integral.
    :param float theta: The theta parameter, indicating the scaling of the x parameter in the square root.
    :return: The result of the Gauss quadrature.
    :rtype: float
    """
    return _fermi.gaussFfermi(k, eta, theta)

def Gfermi(n, alpha, beta):
    """The generalized G-type Fermi integral quadrature.

    :param float n: The n parameter, indicating the power of alpha constant in the denominator before integral.
    :param float alpha: The alpha parameter, used before the integral and in the square root.
    :param float beta: The beta parameter, indicating the exponent in the denominator of the integral.
    :return: The result of the sinh-tanh quadrature.
    :rtype: float
    """
    return _fermi.Gfermi(n, alpha, beta)

Gm = Gfermi

def Gp(n, alpha, beta):
    """The generalized G-type Fermi integral quadrature, with negative beta coefficient.

    :param float n: The n parameter, indicating the power of alpha constant in the denominator before integral.
    :param float alpha: The alpha parameter, used before the integral and in the square root.
    :param float beta: The beta parameter, indicating the exponent in the denominator of the integral, taken as negative.
    :return: The result of the sinh-tanh quadrature.
    :rtype: float
    """
    return _fermi.Gp(n, alpha, beta)
