# Hermite elements for FEniCS #

This is a repo to collect the notes and notebooks on my progress in
the implementation of (Cubic) Hermite and Discrete Kirchhoff Triangle
elements for [FEniCS](https://fenicsproject.org/).

## To do ##

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
* Implement Hermite transformation for
  `optimisedquadraturetransformer.py` too.
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

There is a script, `checkout-hermite.sh` to automate the checkouts.

## Contents ##

There is a bunch of Jupyter notebooks with implementation steps,
simple tests and some models for deeper testing. See `index.ipynb`.

Other things

* `checkout_hermite.sh`: A script useful after first pulling the
  docker image with the official sources. It checks out the right
  commits and branches and configures the remotes with my private
  branches.
* `hermite.tm`, `hermite.bib`, `img/`: Document with the math for
  Hermite elements. In progress.
* `discrete_gradient.tm`: A document explaining how these operators
  are built. In progress.
* `test.py`: Use this to run tests with pdb (or realgud:pdb in emacs)
