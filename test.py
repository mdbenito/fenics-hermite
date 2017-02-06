# -*- fenics/dependencies: ("ffc"); pythonic-environment: "/docker:fenics@fenics-dev:/usr" -*-
# Recall that fenics-dev is 127.0.0.1 in /etc/hosts at the docker host

from __future__ import print_function
from dolfin import *
import ffc
from ffc.log import add_logfile, set_level, DEBUG

parameters["form_compiler"]["representation"] = "quadrature"
set_level(DEBUG)
add_logfile("/tmp/fenics.log")

mesh = UnitSquareMesh(32, 32)
V = FunctionSpace(mesh, "Hermite", 3)
u = TrialFunction(V)
v = TestFunction(V)
a = inner(div(grad(u)), div(grad(v)))*dx

with open("/tmp/biharmonic.h", "wt") as f:
    out = ffc.compile_form(a, prefix="bh")
    #print(out[1])
    f.write(out[0])

 
################################################################################
if False:
    from dolfin import *
    mesh = UnitIntervalMesh(10)
    V = FunctionSpace(mesh, 'Hermite', 3)
    f = Argument(V, 0)
    R = assemble(f*dx)
    print("done")

    def x2(x=0):
        return 2*x

    def x3(x=0):
        return 3*x

    print(x)
