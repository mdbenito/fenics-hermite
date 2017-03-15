# Hermite elements for FEniCS #

This is a repo to collect the notes and notebooks on my progress in
the implementation of (Cubic) Hermite elements
for [FEniCS](https://fenicsproject.org/).

## To do ##

* Fix the problems with the scale of the solutions depending on the
  number of cells in the mesh. This is most surely an issue with the
  quadrature code not properly scaling the derivatives or something
  like that.
* Finish implementing boundary conditions: Check the effect that
  meddling with the stiffness matrix has over the solution of the
  system. Finish essential Neumann BCs for 4th order problems (beyond
  the current hack).
* Incorporate the restriction on the space of polynomials for the
  element leading to Kirchhoff triangles.
* Implement arbitrary (d > 3) polynomial orders using hierarchic bases
  (cf. [Solin, §6.3.2]).
* Implement `evaluate_dof()` for `PointDerivative` to enable
  interpolation.
* More tests.
* ???

## Dependencies ##

Most of the code depends on my branches on ffc, fiat and ufl, where
the actual implementation is done.

The code should be compatible with the following revisions of the
other subprojects:

* dolfin: 8d4c7c807451f2301020ec23a2106f0501c632e3
* dijitso: 030e1bb59dbf55ae05b09078597cd01bfe89bdc3
* instant: 4bd92ec7eb114c8dfe944e770ce3c22323a75b0d
* mshr: c4058b5287722fbcc9dd8ec25bebfd31e3e58ea4


## Files ##

* `Poisson1D.ipynb`: A 1D example to test Dirichlet BCs.
* `Poisson2D.ipynb`: The official 2D example in dolfin, with Hermite
  elements. **TODO:** Find out why the solution is slightly
  skewed. Properly apply Neumann.
* `discrete-gradient.ipynb`: An implementation of the discrete
  gradient operator from CG1 to DG0 using PETSc matrices. This should
  be the first step towards the actual operator I need for Kirchhoff
  elements.
* `boundary.ipynb`: Steps towards the implementation of Dirichlet and
  Neumann boundary conditions with Hermite elements. Mostly just a
  hack on the assembled stiffness matrix.
* `Euler-Bernoulli.ipynb`: The Euler Bernoulli beam model in a
   notebook.  The mathematics are taken from (and possibly out of sync
   with) `hermite.tm`.
* `biharmonic.ipynb`: A test (biharmonic with non conforming interior
   penalty method) using Lagrange or Hermite elements.
* `debug.ipynb`: Some random tests.
* `elements.ipynb`: Progress during the implementation in FFC.
* `quadrature.ipynb`: Progress during the implementation of
   quadratures.
* `hermite.tm`, `hermite.bib`, `img/`: Document with the math for
  Hermite elements. In progress.
* `instant-test.ipynb`: Tests with instant.
* `test.py`: Use this to run tests with pdb (or realgud:pdb in emacs)
* `utils.py`: colorizing and other amazing stuff.
* `README.md`: ...
* `doc`: duh.
* `misc.ipynb`: double duh.
* `nbimporter.py`: A module to import notebooks with some custom
  modifications.
