from fermidirac.helpers.derivative_relations import *
from fermidirac.interface.cfermidirac import Ffermi, gaussFfermi
from itertools import product

import pytest
import math

PARAMETER_TUPLES = list(product(
    [(x - 1) * 0.5 for x in range(10)], # fractional k's from -1/2 to 4.0 with 0.5 step
    [(x + 1) * 0.2 for x in range(25)], # fractional etas from 0 to 5.0 with 0.2 step
    [(x + 1) * 0.2 for x in range(25)], # fractional thetas from 0 to 5.0 with 0.2 step
))

RELATIVE_TOLERANCE = 1e-9 # tolerance relative to larger absolute value in fractional terms
ABSOLUTE_TOLERANCE = 0.0 # useful when one of the comparisons is to 0.0

@pytest.mark.parametrize("k, eta, theta", PARAMETER_TUPLES)
def test_derivatives_using_cox_giuli(k, eta, theta):
    left_hand_side, right_hand_side = cox_giuli(k, eta, theta)
    assert math.isclose(
        left_hand_side,
        right_hand_side,
        rel_tol=RELATIVE_TOLERANCE,
        abs_tol=ABSOLUTE_TOLERANCE
    )


@pytest.mark.parametrize("k, eta, theta", PARAMETER_TUPLES)
def test_derivatives_using_gong_zejda_k_kplus1(k, eta, theta):
    left_hand_side, right_hand_side = gong_zejda_k_kplus1(k, eta, theta)
    assert math.isclose(
        left_hand_side,
        right_hand_side,
        rel_tol=RELATIVE_TOLERANCE,
        abs_tol=ABSOLUTE_TOLERANCE
    )


@pytest.mark.parametrize("k, eta, theta", PARAMETER_TUPLES)
def test_derivatives_using_gong_zejda_kminus1_k_kplus1(k, eta, theta):
    left_hand_side, right_hand_side = gong_zejda_kminus1_k_kplus1(k, eta, theta)
    if k - 1.0 >= 0.5:
        assert math.isclose(
            left_hand_side,
            right_hand_side,
            rel_tol=RELATIVE_TOLERANCE,
            abs_tol=ABSOLUTE_TOLERANCE
        )


@pytest.mark.parametrize("k, eta, theta", PARAMETER_TUPLES)
def test_derivatives_using_gong_zejda_etaderivatives(k, eta, theta):
    left_hand_side, right_hand_side = gong_zejda_etaderivatives(k, eta, theta)
    if k - 1.0 >= 0.5:
        assert math.isclose(
            left_hand_side,
            right_hand_side,
            rel_tol=RELATIVE_TOLERANCE,
            abs_tol=ABSOLUTE_TOLERANCE
        )


@pytest.mark.parametrize("k, eta, theta", PARAMETER_TUPLES)
def test_ffermi_compared_to_reference_gauss(k, eta, theta):
    assert math.isclose(
        Ffermi(k, eta, theta),
        gaussFfermi(k, eta, theta),
        rel_tol=RELATIVE_TOLERANCE,
        abs_tol=ABSOLUTE_TOLERANCE
    )
