###############################################################################
## For some reason using nbimporter on discrete_gradient.ipynb to
## import the stuff I copied here fails. The error is most weird: some
## dolfin routine complains that 'form_compiler' is not a valid key in
## 'parameters', which makes exactly no sense until one checks the
## object and sees that it is indeed not there...

## Also: evaluation of that notebook kills the kernel !?!?!?!!?

from dolfin import *
from petsc4py import PETSc
from itertools import chain

import numpy as np

class DKTCellGradient(object):
    """ Computes the transformation matrix for $\nabla_h$ over one cell.
    
    This assumes that the dofs of the "source" cell (in DKT) are ordered as
    
         w1, w1_x, w1_y, w2, w2_x, w2_y, w3, w3_x, w3_y,
         
    i.e. point evaluation, partial derivatives for each vertex in turn, and
    those of the "target" cell (in P_2^2) as
    
        t1_0, t1_1, t2_0, t2_1, t3_0, t3_1,
        ts1_0, ts1_1, ts2_0, ts2_1, ts3_0, ts3_1
        
    i.e. values at the three vertices, then values at the the midpoints of
    the opposite sides.
    """
    def __init__(self):
        """ Initialize the local cell gradient matrix.
        The matrix is 12x9  (target's local space_dimension() x source's 
        local space_dimension()."""

        self.M = np.zeros((12, 9), dtype=np.float)
        rows = [[0,1], [2,3], [4,5], [6,7], [8,9], [10,11]]
        # Hack to make fancy indexing with lists work:
        # if one of rows[i], self.cols[j] has shape 2,1, broadcasting will actually
        # compute the outer product of rows[i] and cols[j], instead of a "zip"
        # of sorts.
        self.rows = list(map(lambda x: np.array(x).reshape(2,1), rows))
        self.cols = [0, [1,2], 3, [4,5], 6, [7,8]]
        #self.cols = list(map(lambda x: np.array(x).reshape((-1,1)), self.cols))
        # TODO: Note that storing these copies of the identity is silly
        # When we apply the operator, it would make more sense to just copy
        # the partial derivative dofs (which is what these first 3 multiplications 
        # with the identity do), then multiply with the lower part of M...
        # (But we still need these here if we insist on building a global matrix
        # which acts on all dofs of the space simultaneously)
        self.M[self.rows[0], self.cols[1]] = np.eye(2)
        self.M[self.rows[1], self.cols[3]] = np.eye(2)
        self.M[self.rows[2], self.cols[5]] = np.eye(2)

        # We are only defined for triangles (to do: assert this)
        # tdim = 2, num_vertices = 3

        # Init temporaries
        self.tt = np.empty((3, 2))
        self.TT = np.empty((3, 2, 2))
        self.ss = np.empty(3)

    def update(self, cell:Cell):
        """Update the mutating entries of M for the given cell."""
        cc = cell.get_vertex_coordinates().reshape(-1, 2)
        self.tt[0] = (cc[2] - cc[1])
        self.tt[1] = -(cc[0] - cc[2])  # Magic: why was this minus sign here again?
        self.tt[2] = (cc[1] - cc[0])
        for i in range(3):
            self.ss[i] = np.linalg.norm(self.tt[i], ord=2)
            self.tt[i] /= self.ss[i]
            self.TT[i] = -3./4 * np.outer(self.tt[i], self.tt[i])
            # HERE! if we don't use the hierarchical basis we need to add this, see...
            self.TT[i] += 0.5*np.eye(2)
            self.tt[i] *= -3./(2*self.ss[i])

        self.M[self.rows[3], self.cols[2]] = self.tt[0].copy().reshape(2,1)
        self.M[self.rows[3], self.cols[3]] = self.TT[0].copy()
        self.M[self.rows[3], self.cols[4]] = -self.tt[0].copy().reshape(2,1)
        self.M[self.rows[3], self.cols[5]] = self.TT[0].copy()

        self.M[self.rows[4], self.cols[0]] = self.tt[1].copy().reshape(2,1)
        self.M[self.rows[4], self.cols[1]] = self.TT[1].copy()
        self.M[self.rows[4], self.cols[4]] = -self.tt[1].copy().reshape(2,1)
        self.M[self.rows[4], self.cols[5]] = self.TT[1].copy()

        self.M[self.rows[5], self.cols[0]] = self.tt[2].copy().reshape(2,1)
        self.M[self.rows[5], self.cols[1]] = self.TT[2].copy()
        self.M[self.rows[5], self.cols[2]] = -self.tt[2].copy().reshape(2,1)
        self.M[self.rows[5], self.cols[3]] = self.TT[2].copy()



def dkt_gradient_operator(V, Vp):
    """ Build sparse PETSc matrix of the discrete gradient
    operator between V and Vp.
    
    Arguments
    ---------
        V: Source FunctionSpace of type DKT.
        Vp: Target (Vector)FunctionSpace of type P2
            (WITHOUT a hierarchical basis)
        
    Returns
    -------
        Dh: PETSc sparse matrix such that for u in V,
              Dh * u.vector().array()
            yields the vector of coefficients of the discrete
            gradient of u in Vp.
    """   
    gdim = V.mesh().geometry().dim()
    assert V.ufl_element().shortstr()[:3] == 'DKT',\
           "I need a DKT source space"
    assert Vp.ufl_element().shortstr()[:14] == 'Vector<%d x CG2' % gdim,\
           "I need a CG2 VectorFunctionSpace with %d components" % gdim
    assert Vp.mesh().num_cells() == V.mesh().num_cells(),\
           "Number of cells in the meshes don't match"
    dm = V.dofmap()
    dmp1, dmp2 = Vp.sub(0).dofmap(), Vp.sub(1).dofmap()
    e = V.element()
    ep = Vp.element()
    tdim = e.topological_dimension()
    assert tdim == 2, "I need topological dimension 2 (was %d)" % tdim
    coords = V.tabulate_dof_coordinates().reshape((-1, tdim))

    Dh = PETSc.Mat()
    Dh.create(PETSc.COMM_WORLD)
    Dh.setSizes([Vp.dim(), V.dim()])
    Dh.setType("aij")
    Dh.setUp()
    # In order not to assume the following it would be enough
    # to create the array inside the loop...
    warning("Assuming all cells have the same number of dofs")

    ## Preinit element-wise gradient
    grad = DKTCellGradient()

    # TODO: check this for parallel operation...
    istart, iend = Dh.getOwnershipRange()
    #print ("This process owns the range %d - %d" % (istart, iend))

    for cell_id in range(V.mesh().num_cells()):
        cell = Cell(V.mesh(), cell_id)
        grad.update(cell)
        # Massage dofs into the ordering we require with
        # both coordinates for each vector-valued dof consecutive.
        dest_rows = np.array(list(chain.from_iterable(zip(dmp1.cell_dofs(cell_id),
                                                          dmp2.cell_dofs(cell_id)))),
                             dtype=np.intc)
        dest_rows = dest_rows[(dest_rows>=istart) & (dest_rows < iend)]
        # dofs in V are already in the proper ordering (point eval, gradient eval, ...)
        dest_columns = dm.cell_dofs(cell_id)

        # NOTE: the copy() is necessary!
        #print("Setting:\n\trows=%s\n\tcols=%s\nwith:\n%s" % (dest_rows, dest_columns, M))
        Dh.setValues(dest_rows, dest_columns, grad.M.copy())
        assert len(dest_rows)*len(dest_columns) ==  grad.M.size, "duh"
    Dh.assemble()
    return Dh

def compute_dkt_gradient(V, Vp, u):
    """ Returns the discrete gradient of u.

    Arguments
    ---------
        V: Source FunctionSpace of type DKT. (redundant...)
        Vp: Target VectorFunctionSpace of type CG2.

    Returns
    -------
        A Function in Vp interpolating the gradient of u.
    """
    Dh = dkt_gradient_operator(V, Vp)
    du = Function(Vp)
    petsc_u = as_backend_type(u.vector()).vec()
    petsc_du = as_backend_type(du.vector()).vec()
    
    Dh.mult(petsc_u, petsc_du)
    du.vector().apply('insert')

    return du
