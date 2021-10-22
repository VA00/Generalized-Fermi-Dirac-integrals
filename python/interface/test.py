import sys
try:
    import _fermidirac_api.lib as _fermi
except ImportError as exc:
    print(f"API not found ({exc}), attempting ABI...")
    try:
        import _fermidirac_abi.lib as _fermi
    except ImportError as exc:
        sys.exit(f"No generated interface found ({exc}). Please generate interface first.")

k = 4.0
result = _fermi.Ffermi(k, 1.0, 1.0)
print(f"{k}: {result}")
