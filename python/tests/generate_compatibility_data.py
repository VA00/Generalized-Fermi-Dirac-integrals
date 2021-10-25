from itertools import product
import csv

from fermidirac.interface.cfermidirac import *

PARAMETER_TUPLES = list(product(
    [(x - 1) * 0.5 for x in range(10)], # fractional k's from -1/2 to 4.0 with 0.5 step
    [(x + 1) * 1 for x in range(5)], # etas
    [(x + 1) * 1 for x in range(5)], # thetas
))

def generate_result_csv_fermi_functions():
    fields = [
        "k",
        "eta",
        "theta",
        "Ffermi",
        "fixedFfermi",
        "Ffermi_long",
        "fixedFfermi_long",
        "fixedFfermi_derivatives",
        "df/deta",
        "df/dtheta",
        "d2f/deta2",
        "d2f/dtheta2",
        "d2f/detadtheta"
    ]

    rows = []
    for k, eta, theta in PARAMETER_TUPLES:
        derivatives = fixedFfermi_derivatives(k, eta, theta, 0.125, -4.0, 4.0)
        rows.append([
            k,
            eta,
            theta,
            Ffermi(k, eta, theta),
            fixedFfermi(k, eta, theta, 0.125, -4.0, 4.0),
            Ffermi_long(k, eta, theta),
            fixedFfermi_long(k, eta, theta, 0.125, -4.0, 4.0),
            derivatives.f,
            derivatives.df_deta,
            derivatives.df_dtheta,
            derivatives.d2f_deta2,
            derivatives.d2f_dtheta2,
            derivatives.d2f_deta_dtheta
        ])

    with open("fermi_functions_reference_data.csv", "w") as f:
        write = csv.writer(f)
        write.writerow(fields)
        write.writerows(rows)

if __name__ == '__main__':
    generate_result_csv_fermi_functions()
