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
git checkout 8d4c7c807451f2301020ec23a2106f0501c632e3

cd $SRC/mshr
git checkout c4058b5287722fbcc9dd8ec25bebfd31e3e58ea4

cd $SRC/dijitso
git checkout tags/dijitso-2016.2.0

cd $SRC/instant
git checkout tags/instant-2016.2.0

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

