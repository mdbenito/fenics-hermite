# Implementation of Hermite elements for FEniCS

# Hermite elements for FEniCS #

This is a repo to collect the notes and notebooks on my progress. Most
of the code depends on my branches on ffc, fiat and ufl, where the
actual implementation is done.

## Files ##

* `Euler-Bernoulli.ipynb`: The Euler Bernoulli beam model in a
   notebook. Missing: How to set the boundary conditions for the
   simply supported beam? I'm a bit confused by how to set the bending
   moments to zero (the text in the notebook might be wrong). The
   mathematics are taken from (and possibly out of sync with)
   `hermite.tm`.
* `biharmonic.ipynb`: A test (biharmonic with non conforming interior
   penalty method) using Lagrange or Hermite elements.
* `debug.ipynb`: Some random tests.
* `ffc-elements.ipynb`: Progress during the implementation in FFC.
* `ffc-forms.ipynb`: Progress during the implementation of
   quadratures.
* `hermite.tm`: Math for the Hermite elements. In progress.
* `instant.ipynb`: Tests with instant.
* `test.py`: Use this to run tests with pdb (or realgud:pdb in emacs)
* `utils.py`: colorizing and other amazing stuff.
* `README.md`: ...
* `doc`: duh.
* `misc.ipynb`: double duh.
