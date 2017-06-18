# Hermite elements for FEniCS #

This is a repo to collect the notes and notebooks on my progress in
the implementation of Cubic Hermite and Discrete Kirchhoff Triangle
elements for [FEniCS](https://fenicsproject.org/), as well as their application to some models
(Euler-Bernoulli, Poisson, Biharmonic, non-linear Kirchhoff...) .

## Contents ##

* A bunch of Jupyter notebooks with implementation steps, simple tests
  and some models for deeper testing. See `index.ipynb` for details.
* `hermite.tm`, `hermite.bib`, `img/`: Document with the math for
  Hermite elements. In progress.
* `discrete_gradient.tm`: A document explaining how these operators
  are built. Contains descriptions of
  the
  [linear Kirchhoff](https://bitbucket.org/mdbenito/linear-kirchhoff)
  and
  [nonlinear Kirchhoff](https://bitbucket.org/mdbenito/nonlinear-kirchhoff) models. In
  progress.
* `test.py`: Use this to run tests with `pdb` (or `realgud:pdb` in emacs)
* `checkout_hermite.sh`: A script useful after first pulling the
  docker image with the official sources. It checks out the right
  commits and branches and configures the remotes with my private
  branches.

## To do ##

* More tests.
* Finish implementing Dirichlet boundary conditions: Check the effect
  that meddling with the stiffness matrix has over the solution of the
  system. **In particular check that symmetry is not broken, cf §6.3 of
  the FEniCS book.**
* Finish essential Neumann BCs for 4th order problems (beyond
  the current hack).
* Implement arbitrary (d > 3) polynomial orders in Hermite elements
  using hierarchic bases (cf. [Solin, §6.3.2]).
* Disable `vertex_to_dof_map()` and `dof_to_vertex_map()` for spaces
  using Hermite elements.
* <strike>Implement `evaluate_dof()` for `PointDerivative` to enable
  nodal interpolation.</strike> See 
  [interpolation.ipynb](interpolation.ipynb).
* Implement Hermite transformation for
  `optimisedquadraturetransformer.py` too.
* ???

## Dependencies ##

Most of the code in the notebooks and models depends on my branches on
ffc, fiat and ufl, where Hermite and DKT elements are implemented
(although a description of the required modifications is contained in
some of the notebooks).

As to the other subprojects (dolfin, dijitso, instant and mshr), this
code should be compatible with `tags/2017.1.0`.

The script `checkout-hermite.sh` automates the checkouts inside the
container.

## License ##

Licensed under the [GPL v.3](https://www.gnu.org/licenses/gpl-3.0.en.html).
