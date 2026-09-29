#2D Mode-I Central Crack: Static BEM/DDM Example
This example solves one static load case only. It does not include time integration or crack propagation.
The problem is a straight central crack in an infinite elastic plate subjected to remote Mode-I tension. The crack is solved with the Displacement Discontinuity Method (DDM), a boundary-element formulation in which the displacement jump across the crack is used as the primary unknown.
##1. Problem definition
2D linear elasticity
Infinite plate
Central straight crack of total length `2a`
Uniform remote tensile stress `sigma_inf`
Plane strain
Crack located on the x-axis
Traction-free crack faces
No crack propagation
No time dependence
Default parameters:
Variable	Value
Young's modulus, `E`	70 GPa
Poisson's ratio, `nu`	0.33
Remote tensile stress, `sigma_inf`	10 MPa
Half crack length, `a`	10 mm
Total crack length, `2a`	20 mm
Crack elements, `N`	80
##2. BEM/DDM formulation
Unlike FEM, the whole domain is not meshed. Only the crack line is discretized into boundary elements.
The unknown on each crack element is the normal displacement discontinuity:
$$ D_n=u_n^- -u_n^+ $$
With the Crouch-Starfield sign convention, a physically opening crack has `Dn < 0`. Therefore the Crack Opening Displacement (COD) used in the plots is
$$ COD=-D_n $$
At each collocation point, the sum of the remote tensile stress and the stresses induced by all crack elements must satisfy the traction-free crack-face condition:
$$ \sigma_\infty+\sum_j A_{ij}D_{n,j}=0 $$
The resulting linear system is
$$ A\mathbf{D}n=-\sigma\infty\mathbf{1} $$
For a constant DDM element on the crack line, the influence coefficient used in the code is
$$ A_{ij}=-{G\over 2\pi(1-\nu)}\left[{1\over x_i-x_j-h}-{1\over x_i-x_j+h}\right] $$
where `h` is the half-length of a crack element.
##3. Analytical benchmark
For a central crack in an infinite plate under uniform remote Mode-I tension, the exact crack opening displacement is
$$ COD_{exact}(x)={4\sigma_\infty\over E'}\sqrt{a^2-x^2} $$
For plane strain,
$$ E'={E\over 1-\nu^2} $$
The exact Mode-I stress intensity factor is
$$ K_I=\sigma_\infty\sqrt{\pi a} $$
The code compares the numerical DDM COD with the analytical COD. It also fits the numerical COD to the `sqrt(a^2-x^2)` shape and uses the fitted amplitude to obtain a verification estimate of `K_I`.
> Note: this `K_I` is a **global COD-fit estimate**, not a direct constant-element crack-tip extrapolation.
##4. Running the code
Required packages:
```bash
pip install numpy matplotlib
```
Run:
```bash
python BEM_ModeI_DDM.py
```
##5. Outputs
The console reports:
maximum numerical COD
maximum analytical COD
relative L2 error of COD
exact `K_I`
COD-fit `K_I`
maximum crack-face traction residual
The following figures are saved:
```text
bem_modeI_cod.png
bem_modeI_opening_shape.png
bem_modeI_traction_residual.png
```
`bem_modeI_cod.png`  
: numerical BEM/DDM COD compared with the analytical solution
`bem_modeI_opening_shape.png`  
: upper and lower crack faces for the single static Mode-I case
`bem_modeI_traction_residual.png`  
: numerical check of the traction-free boundary condition
##6. What is not included
This example does not include:
crack propagation
physical time `t`
inertia
dynamic fracture
a finite rectangular outer boundary
mixed Mode I/II loading
plasticity
a cohesive-zone model
It is intentionally a minimal one-case fracture-mechanics BEM example.
To model a finite rectangular plate and treat both coincident crack faces explicitly, a natural next step is the Dual Boundary Element Method (DBEM).
##7. References
Crouch, S. L. and Starfield, A. M. Boundary Element Methods in Solid Mechanics: With Applications in Rock Mechanics and Geological Engineering. Allen & Unwin, 1983.
Portela, A., Aliabadi, M. H., and Rooke, D. P. “The dual boundary element method: Effective implementation for crack problems.” International Journal for Numerical Methods in Engineering, 33, 1269–1287, 1992.
