#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Static 2D Mode-I fracture benchmark using a constant
Displacement Discontinuity Method (DDM), a boundary-element formulation.

Problem
-------
Infinite isotropic elastic plate with a central straight crack of length 2a,
subjected to a uniform remote tensile stress sigma_inf normal to the crack.

This is a SINGLE STATIC CASE:
- no time stepping
- no crack propagation
- no dynamic/inertial effects

The crack alone is discretized into constant boundary elements.
The unknown is the normal displacement discontinuity D_n on each element.

Sign convention used by Crouch & Starfield:
    D_n = u_n^- - u_n^+
so a physically opening crack has D_n < 0.
The Crack Opening Displacement (COD) plotted here is therefore
    COD = -D_n > 0

Dependencies
------------
numpy
matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# 1. User parameters
# ============================================================

E = 70.0e9             # Young's modulus [Pa]
nu = 0.33              # Poisson's ratio [-]
sigma_inf = 10.0e6     # remote Mode-I tensile stress [Pa]
a = 10.0e-3            # half crack length [m] -> total crack length = 2a
N = 80                 # number of constant crack elements

# Plane-strain benchmark
G = E / (2.0 * (1.0 + nu))
Eprime = E / (1.0 - nu**2)


# ============================================================
# 2. Crack discretization
# ============================================================

dx = 2.0 * a / N
h = dx / 2.0

# element centers
x = np.linspace(-a + h, a - h, N)


# ============================================================
# 3. DDM / BEM influence matrix
# ============================================================
#
# For a constant normal displacement-discontinuity element centered at x_j,
# extending from x_j-h to x_j+h, the induced sigma_yy on the crack line y=0
# at x_i is
#
# sigma_yy = C * D_n *
#            [ 1/(x_i-x_j-h) - 1/(x_i-x_j+h) ]
#
# where
#
# C = -G / [2*pi*(1-nu)]
#
# The total crack-face traction must vanish:
#
# sigma_inf + sum_j A_ij D_n,j = 0
#
# so
#
# A D_n = -sigma_inf
#

C = -G / (2.0 * np.pi * (1.0 - nu))
A = np.zeros((N, N), dtype=float)

for i, xi in enumerate(x):
    r = xi - x
    A[i, :] = C * (1.0 / (r - h) - 1.0 / (r + h))

rhs = -sigma_inf * np.ones(N)

# Solve for normal displacement discontinuity
Dn = np.linalg.solve(A, rhs)

# Crouch convention: negative Dn corresponds to physical opening
cod = -Dn


# ============================================================
# 4. Analytical benchmark
# ============================================================
#
# Infinite plate, central crack, uniform remote tension:
#
# COD_exact(x) = (4 sigma_inf / E') * sqrt(a^2 - x^2)
#
# K_I,exact = sigma_inf * sqrt(pi a)
#

cod_exact = (4.0 * sigma_inf / Eprime) * np.sqrt(np.maximum(a*a - x*x, 0.0))
KI_exact = sigma_inf * np.sqrt(np.pi * a)

# Fit the numerical COD shape to C_fit*sqrt(a^2-x^2)
# and infer an effective Mode-I SIF from the fitted amplitude.
phi = np.sqrt(np.maximum(a*a - x*x, 0.0))
C_fit = np.dot(phi, cod) / np.dot(phi, phi)
sigma_fit = C_fit * Eprime / 4.0
KI_fit = sigma_fit * np.sqrt(np.pi * a)

rel_l2_cod = np.linalg.norm(cod - cod_exact) / np.linalg.norm(cod_exact)
ki_error = abs(KI_fit - KI_exact) / KI_exact

# Check traction-free residual
traction_total = sigma_inf + A @ Dn


# ============================================================
# 5. Print summary
# ============================================================

print("=" * 68)
print("Static Mode-I Crack: Constant DDM / BEM")
print("=" * 68)
print(f"E                         = {E/1e9:.3f} GPa")
print(f"nu                        = {nu:.4f}")
print(f"Remote tension            = {sigma_inf/1e6:.3f} MPa")
print(f"Half crack length a       = {a*1e3:.3f} mm")
print(f"Total crack length 2a     = {2*a*1e3:.3f} mm")
print(f"Number of crack elements  = {N}")
print("-" * 68)
print(f"Max COD (DDM)             = {np.max(cod)*1e6:.6f} um")
print(f"Max COD (analytical)      = {np.max(cod_exact)*1e6:.6f} um")
print(f"COD relative L2 error     = {100*rel_l2_cod:.4f} %")
print(f"K_I exact                 = {KI_exact/1e6:.6f} MPa*sqrt(m)")
print(f"K_I from COD fit          = {KI_fit/1e6:.6f} MPa*sqrt(m)")
print(f"K_I fit error             = {100*ki_error:.4f} %")
print(f"Max traction residual     = {np.max(np.abs(traction_total))/1e6:.3e} MPa")
print("=" * 68)


# ============================================================
# 6. Plot 1: COD comparison
# ============================================================

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(x * 1e3, cod * 1e6, "o", ms=4, label="DDM / BEM")
ax.plot(x * 1e3, cod_exact * 1e6, "-", lw=2, label="Analytical")
ax.set_xlabel("Crack coordinate x [mm]")
ax.set_ylabel("Crack opening displacement [um]")
ax.set_title("Mode-I Crack Opening Displacement")
ax.grid(True, alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig("bem_modeI_cod.png", dpi=220)


# ============================================================
# 7. Plot 2: opened crack shape
# ============================================================

fig, ax = plt.subplots(figsize=(9, 4.8))

upper = +0.5 * cod * 1e6
lower = -0.5 * cod * 1e6

ax.plot(x * 1e3, upper, "-", lw=2, label="upper crack face")
ax.plot(x * 1e3, lower, "-", lw=2, label="lower crack face")
ax.fill_between(x * 1e3, lower, upper, alpha=0.18)

# remote loading arrows (schematic)
arrow_x = np.linspace(-0.8*a, 0.8*a, 7) * 1e3
y_arrow = 0.62 * np.max(cod) * 1e6

for xx in arrow_x:
    ax.annotate(
        "",
        xy=(xx, y_arrow * 1.45),
        xytext=(xx, y_arrow),
        arrowprops=dict(arrowstyle="->", lw=1.2),
    )
    ax.annotate(
        "",
        xy=(xx, -y_arrow * 1.45),
        xytext=(xx, -y_arrow),
        arrowprops=dict(arrowstyle="->", lw=1.2),
    )

ax.text(
    0.0,
    y_arrow * 1.62,
    f"remote tension = {sigma_inf/1e6:.1f} MPa",
    ha="center",
    va="bottom",
)

ax.set_xlabel("Crack coordinate x [mm]")
ax.set_ylabel("Opening direction [um]")
ax.set_title("Single Static Mode-I Opening Case")
ax.grid(True, alpha=0.25)
ax.legend()
fig.tight_layout()
fig.savefig("bem_modeI_opening_shape.png", dpi=220)


# ============================================================
# 8. Plot 3: traction-free residual
# ============================================================

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.plot(x * 1e3, traction_total / 1e6, "o-", ms=3)
ax.axhline(0.0, lw=1.0)
ax.set_xlabel("Crack coordinate x [mm]")
ax.set_ylabel("Total crack-face traction [MPa]")
ax.set_title("Traction-Free Boundary-Condition Residual")
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig("bem_modeI_traction_residual.png", dpi=220)

plt.show()
