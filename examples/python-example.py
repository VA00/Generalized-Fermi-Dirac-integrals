import cfermidirac

k = 4.0
result = cfermidirac.Ffermi(k, 1.0, 1.0)
print(f"Ffermi({k}, 1.0, 1.0): {result}")
result = cfermidirac.fixedFfermi(k, 1.0, 1.0, 0.125, -5.0, 5.0)
print(f"fixedFfermi({k}, 1.0, 1.0, 0.125, -5.0, 5.0): {result}")
result = cfermidirac.fixedFfermi_long(k, 1.0, 1.0, 0.125, -5.0, 5.0)
print(f"fixedFfermi_long({k}, 1.0, 1.0, 0.125, -5.0, 5.0): {result}")
result = cfermidirac.gaussFfermi(k, 1.0, 1.0)
print(f"gaussFfermi({k}, 1.0, 1.0, 0.125, -5.0, 5.0): {result}")
derivatives = cfermidirac.fixedFfermi_derivatives(k, 1.0, 1.0, 0.125, -5.0, 5.0)
print(f"fixedFfermi_derivatives({k}, 1.0, 1.0, 0.125, -5.0, 5.0): {derivatives}")
result = cfermidirac.Ffermi_long(k, 1.0, 1.0)
print(f"Ffermi_long({k}, 1.0, 1.0): {result}")
n = 1.0
result = cfermidirac.Gfermi(n, 1.0, 1.0)
print(f"Gfermi({n}, 1.0, 1.0): {result}")
result = cfermidirac.Gm(n, 1.0, 0.0)
print(f"Gm({n}, 1.0, 1.0): {result}")
result = cfermidirac.Gp(n, 1.0, 0.0)
print(f"Gp({n}, 1.0, 1.0): {result}")

