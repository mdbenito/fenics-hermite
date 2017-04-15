# Hermite elements for FEniCS #

This is a repo to collect the notes and notebooks on my progress in
the implementation of (Cubic) Hermite and Discrete Kirchhoff Triangle
elements for [FEniCS](https://fenicsproject.org/).

## To do ##

* Fix the problems with the scale of the solutions depending on the
  number of cells in the mesh. This is most surely an issue with the
  quadrature code not properly scaling the derivatives or something
  like that.
* Finish implementing Dirichlet boundary conditions: Check the effect
  that meddling with the stiffness matrix has over the solution of the
  system. **In particular check that symmetry is not broken, cf §6.3 of
  the FEniCS book.**
* Finish essential Neumann BCs for 4th order problems (beyond
  the current hack).
* Implement arbitrary (d > 3) polynomial orders using hierarchic bases
  (cf. [Solin, §6.3.2]).
* Disable `vertex_to_dof_map()` and `dof_to_vertex_map()` for spaces
  using Hermite elements.
* More tests.
* <strike>Implement `evaluate_dof()` for `PointDerivative` to enable
  nodal interpolation.</strike> See 
  [interpolation.ipynb](interpolation.ipynb).
* ???

## Dependencies ##

Most of the code depends on my branches on ffc, fiat and ufl, where
the actual implementation is done.

The code should be compatible with the following revisions of the
other subprojects:

* dolfin: 8d4c7c807451f2301020ec23a2106f0501c632e3
* dijitso: tags/dijitso-2016.2.0
* instant: tags/instant-2016.2.0
* mshr: c4058b5287722fbcc9dd8ec25bebfd31e3e58ea4

## Files ##

Basic implementation:

* `fiat_hermite.ipynb`: Implementation of the FIAT Hermite elements.
* `fiat_dkt.ipynb`: Implementation of the FIAT DKT elements.
* `ffc_tests.ipynb`: Tests and routines used during the implementation
  in FFC.
* `quadrature.ipynb`: Progress during the implementation of
   quadratures.
* `interpolation.ipynb`: Nodal interpolation for elements using
  PointDerivatives.
* `discrete_gradient.ipynb`: An implementation using PETSc matrices of
  the discrete gradient operator from CG1 to DG0 and from DKT to
  Vector<P2>.
* `boundary.ipynb`: Steps towards the implementation of Dirichlet and
  Neumann boundary conditions with Hermite elements. Mostly just a
  hack on the assembled stiffness matrix.

Models:

* `Poisson1D.ipynb`: A simple 1D example to test essential BCs.
* `Poisson2D.ipynb`: Dolfin's included 2D example, but with Hermite
  elements. **TODO:** Find out why the solution is slightly
  skewed.
* `Euler-Bernoulli.ipynb`: The Euler Bernoulli beam model. The
   mathematics are taken from (and possibly out of sync with)
   `hermite.tm`.
* `biharmonic.ipynb`: A test (biharmonic with non conforming interior
   penalty method) using Lagrange or Hermite elements.
* `linear_kirchhoff.ipynb`: A test for the non-conforming method
  using the discrete gradient operator for Kirchhoff elements.

Other stuff:

* `debug.ipynb`: Some random tests.
* `hermite.tm`, `hermite.bib`, `img/`: Document with the math for
  Hermite elements. In progress.
* `instant_test.ipynb`: Tests with instant.
* `test.py`: Use this to run tests with pdb (or realgud:pdb in emacs)
* `utils.py`: Colorizing text, boolean operations on functions and more.
* `README.md`: Suprise!
* `doc/`: duh.
* `misc.ipynb`: double duh.
