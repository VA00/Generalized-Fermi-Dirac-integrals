#!/usr/bin/env python3
"""
Defines relations between various derivatives, which allows cross-checking of calculations.
"""

from fermidirac.helpers.Derivatives import Derivatives
from fermidirac.interface import cfermidirac

def cox_giuli(k, eta, theta, quadrature=cfermidirac.fixedFfermi_derivatives):
    """Calculates both sides of the following relation:
    θ ∂Fk(η,θ)/∂θ = ∂Fk+1(η,θ)/∂η - (k+1)Fk(η,θ)

    Eq. (3) from 'Generalized Fermi–Dirac functions and derivatives: properties and evaluation'
    by Z. Gong, L. Zejda, W. Däppen, J. M. Aparicio

    :param float k: k parameter
    :param float eta: eta parameter
    :param float theta: theta parameter
    :param func quadrature: function to use for calculations, default fixedFfermi_derivatives.
    :returns: 2-tuple containing result of left and right hand side.
    """
    h, hmin, hmax = 0.125, -4.0, 4.0

    f_k = quadrature(k, eta, theta, h, hmin, hmax)
    f_kplus1 = quadrature(k + 1.0, eta, theta, h, hmin, hmax)

    return (theta * f_k.df_dtheta, f_kplus1.df_deta - (k + 1.0)*f_k.f)


def gong_zejda_k_kplus1(k, eta, theta, quadrature=cfermidirac.fixedFfermi_derivatives):
    """Calculates both sides of the following relation:
    Fk+1(η,θ) = 4 ∂Fk(η,θ)/∂θ + 2θ ∂Fk+1(η,θ)/∂θ

    Eq. (4) from 'Generalized Fermi–Dirac functions and derivatives: properties and evaluation'
    by Z. Gong, L. Zejda, W. Däppen, J. M. Aparicio

    :param float k: k parameter
    :param float eta: eta parameter
    :param float theta: theta parameter
    :param func quadrature: function to use for calculations, default fixedFfermi_derivatives.
    :returns: 2-tuple containing result of left and right hand side.
    """
    h, hmin, hmax = 0.125, -4.0, 4.0

    f_k = quadrature(k, eta, theta, h, hmin, hmax)
    f_kplus1 = quadrature(k + 1.0, eta, theta, h, hmin, hmax)

    return (f_kplus1.f, 4 * f_k.df_dtheta + 2 * theta * f_kplus1.df_dtheta)


def gong_zejda_kminus1_k_kplus1(k, eta, theta, quadrature=cfermidirac.fixedFfermi_derivatives):
    """Calculates both sides of the following relation:
    ∂Fk+1(η,θ)/∂η = (k + 3/2) Fk(η,θ) - 2 ∂Fk-1(η,θ)/∂θ

    Eq. (5) from 'Generalized Fermi–Dirac functions and derivatives: properties and evaluation'
    by Z. Gong, L. Zejda, W. Däppen, J. M. Aparicio

    :param float k: k parameter
    :param float eta: eta parameter
    :param float theta: theta parameter
    :param func quadrature: function to use for calculations, default fixedFfermi_derivatives.
    :returns: 2-tuple containing result of left and right hand side.
    """
    h, hmin, hmax = 0.125, -4.0, 4.0

    f_kminus1 = quadrature(k - 1.0, eta, theta, h, hmin, hmax)
    f_k = quadrature(k, eta, theta, h, hmin, hmax)
    f_kplus1 = quadrature(k + 1.0, eta, theta, h, hmin, hmax)

    return (f_kplus1.df_deta, (k + 1.5) * f_k.f - 2 * f_kminus1.df_dtheta)


def gong_zejda_etaderivatives(k, eta, theta, quadrature=cfermidirac.fixedFfermi_derivatives):
    """Calculates both sides of the following relation:
    k Fk-1(η,θ) + (k/2 + 3/4)θ Fk(η,θ) = ∂Fk(η,θ)/∂η + θ/2 ∂Fk+1(η,θ)/∂η

    Eq. (6) from 'Generalized Fermi–Dirac functions and derivatives: properties and evaluation'
    by Z. Gong, L. Zejda, W. Däppen, J. M. Aparicio

    :param float k: k parameter
    :param float eta: eta parameter
    :param float theta: theta parameter
    :param func quadrature: function to use for calculations, default fixedFfermi_derivatives.
    :returns: 2-tuple containing result of left and right hand side.
    """
    h, hmin, hmax = 0.125, -4.0, 4.0

    f_kminus1 = quadrature(k - 1.0, eta, theta, h, hmin, hmax)
    f_k = quadrature(k, eta, theta, h, hmin, hmax)
    f_kplus1 = quadrature(k + 1.0, eta, theta, h, hmin, hmax)

    return (k * f_kminus1.f + (0.5 * k + 0.75) * theta * f_k.f, f_k.df_deta + 0.5 * theta * f_kplus1.df_deta)
