import pytest

from fermidirac.interface.cfermidirac import *

def test_benchmark_ffermi(benchmark):
    benchmark(Ffermi, 4.0, 1.0, 1.0)


def test_benchmark_ffermi_long(benchmark):
    benchmark(Ffermi_long, 4.0, 1.0, 1.0)


def test_benchmark_fixed_ffermi(benchmark):
    benchmark(fixedFfermi, 4.0, 1.0, 1.0, 0.125, -4.0, 4.0)


def test_benchmark_fixed_ffermi_long(benchmark):
    benchmark(fixedFfermi_long, 4.0, 1.0, 1.0, 0.125, -4.0, 4.0)


def test_benchmark_fixed_ffermi_derivatives(benchmark):
    benchmark(fixedFfermi_derivatives, 4.0, 1.0, 1.0, 0.125, -4.0, 4.0)


def test_benchmark_gfermi(benchmark):
    benchmark(Gfermi, 1.0, 1.0, 1.0)


def test_benchmark_gm(benchmark):
    benchmark(Gm, 1.0, 1.0, 1.0)


def test_benchmark_gp(benchmark):
    benchmark(Gp, 1.0, 1.0, 1.0)


def test_benchmark_gaussfermi(benchmark):
    benchmark(gaussFfermi, 4.0, 1.0, 1.0)
