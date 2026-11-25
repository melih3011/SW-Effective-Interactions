#!/usr/bin/env python3
"""Gibbs adsorption of the solute at the liquid-vapor interface.

    python scripts/gibbs_adsorption.py data/runs/profile_q2.00_eps1.3918.dat [...]

For each profile the Gibbs dividing surface z_s of the solvent (Gamma_1 = 0)
is located, the solute adsorption Gamma_2 relative to it is computed, and the
surface tension is estimated from the Gibbs adsorption equation in the dilute
limit, d(beta gamma) = -Gamma_2 d(beta mu_2) = -Gamma_2 d ln rho_2, i.e.

    beta gamma ~ beta gamma_0 - Gamma_2 ,

with gamma_0 the surface tension of the pure solvent (Eq. 26 of the paper).
The estimate is compared with the surface tension from the grand potential.
"""

import re
import sys

import numpy as np

GAMMA0 = 0.20906  # beta gamma sigma_1^2 of the pure SW solvent, beta eps_1 = 1, lambda_1 = 1.5


def adsorption(z: np.ndarray, rho: np.ndarray, zs: float) -> float:
    """Excess of rho over the step profile (rho[0] for z < zs, rho[-1] for z > zs)."""
    return np.trapezoid(rho, z) - rho[0] * (zs - z[0]) - rho[-1] * (z[-1] - zs)


def analyse(path: str) -> None:
    text = open(path).read()
    gamma = float(re.search(r"beta\*gamma\*s1\^2 = ([-0-9.eE+]+)", text).group(1))
    eps2 = float(re.search(r"beta\*eps2 = ([-0-9.eE+]+)", text).group(1))
    q = float(re.search(r"sigma2/sigma1 = ([-0-9.eE+]+)", text).group(1))
    z, rho1, rho2 = np.loadtxt(path, usecols=(0, 1, 2), unpack=True)

    # dividing surface with Gamma_1 = 0 (closed form, Gamma_1 is linear in zs)
    zs = (rho1[-1] * z[-1] - rho1[0] * z[0] - np.trapezoid(rho1, z)) / (rho1[-1] - rho1[0])
    g2 = adsorption(z, rho2, zs)

    print(f"{path}\n  sigma2/sigma1 = {q:g}, beta eps2 = {eps2:.4f}")
    print(f"  Gibbs dividing surface      z_s/sigma1     = {zs:+.4f}")
    print(f"  solute adsorption           Gamma_2 s1^2   = {g2:.4e}")
    print(f"  Gibbs estimate              beta gamma s1^2 = {GAMMA0 - g2:.6f}")
    print(f"  DFT (grand potential)       beta gamma s1^2 = {gamma:.6f}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for p in sys.argv[1:]:
        analyse(p)
