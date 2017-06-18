#!/bin/bash
##############################################################################
# Checks those versions of the repos which are known to work and my
# private branches.
#
# NOTE: this requires that user fenics has an ssh key pair configured
# and with access to the bitbucket sources
#

source $HOME/fenics-dev.env.conf && $HOME/bin/fenics-pull

cd $SRC/dolfin
git checkout tags/2017.1.0

cd $SRC/mshr
git checkout tags/2017.1.0

cd $SRC/dijitso
git checkout tags/2017.1.0

cd $SRC/instant
git checkout tags/2017.1.0

cd $SRC/ffc
git remote rename origin fenics
git remote add origin git@bitbucket.org:mdbenito/ffc-fork.git
git fetch origin
git checkout mdbenito/hermite

cd $SRC/fiat
git remote rename origin fenics
git remote add origin git@bitbucket.org:mdbenito/fiat-fork.git
git fetch origin
git checkout mdbenito/hermite

cd $SRC/ufl
git remote rename origin fenics
git remote add origin git@bitbucket.org:mdbenito/ufl-fork.git
git fetch origin
git checkout mdbenito/hermite

