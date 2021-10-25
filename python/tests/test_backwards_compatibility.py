import pytest
import csv

from fermidirac.interface.cfermidirac import *

COMPATIBILITY_DATA = []

with open("fermi_functions_reference_data.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        COMPATIBILITY_DATA.append(
            row
        )


@pytest.mark.parametrize("data", COMPATIBILITY_DATA)
def test_backwards_fixedffermi_derivative(data):
    result_now = fixedFfermi_derivatives(
        float(data['k']),
        float(data['eta']),
        float(data['theta']),
        0.125,
        -4.0,
        4.0
    )
    assert result_now.f == float(data['fixedFfermi_derivatives'])
    assert result_now.df_deta == float(data['df_deta'])
    assert result_now.df_dtheta == float(data['df_dtheta'])
    assert result_now.d2f_deta2 == float(data['d2f_deta2'])
    assert result_now.d2f_dtheta2 == float(data['d2f_dtheta2'])
    assert result_now.d2f_deta_dtheta == float(data['d2f_deta_dtheta'])

