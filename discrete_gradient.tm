<TeXmacs|1.99.4>

<style|<tuple|generic|british|better-amsart>>

<\body>
  <\hide-preamble>
    \;

    <assign|python-code|<\macro|body>
      <\small>
        <\pseudo-code>
          <python|<arg|body>>
        </pseudo-code>
      </small>
    </macro>>

    <assign|verbatim|<macro|body|<small|<with|font-family|tt|language|verbatim|<arg|body>>>>>

    <assign|python|<macro|body|<small|<with|mode|prog|prog-language|python|font-family|rm|<arg|body>>>>>

    <assign|cpp|<macro|body|<small|<with|mode|prog|prog-language|cpp|font-family|rm|<arg|body>>>>>
  </hide-preamble>

  <section|Discrete Kirchhoff Triangles>

  <todo|Description, applications, implementation in <name|fiat>.>

  Developed in the 1970s and 80s <todo|(refs 5,14 de Bar13)>

  <math|H<rsup|2>>-non-conforming continouous elements, based on reduced
  <math|P<rsub|3>> elements, their main ingredient is the interpolation of
  gradients into a <math|P<rsub|2>> space. <strong|Key:> derivatives at nodes
  are decoupled from function values and can be handled in interleaving steps
  of numerical schemes.

  <section|Discrete gradient operators>

  Taking derivatives of functions belonging to finite element spaces is a
  form of interpolation into other finite element spaces. This can always be
  achieved by orthogonal projection onto the target space but in some cases
  it might be desirable to compute them applying some linear operator, the
  simplest example being piecewise linear (CG1) elements, whose derivatives
  are piecewise constant (CG0) and can be computed exactly and faster than
  projecting. After briefly reviewing this case we focus on the discrete
  gradient operator for Discrete Kirchhoff Triangles which does nodal
  interpolation at the vertices while ensuring a linear normal component
  along the edges of the triangles.

  <\remark>
    <label|rem:extending-ufl-discrete-gradient>A major drawback in
    <name|FEniCS> of these interpolation operators for differentiation \ is
    that in order to add an operator to <name|UFL> we need to instruct
    <name|FEniCS> how to generate <name|UFC> code for any form where it
    appears, as well as code for the respective cell or facet integrals, and
    it is not clear that this is possible at all. Indeed, <name|UFL> forms
    are \Pelement-agnostic\Q, meaning that they have no knowledge of the
    basis of shape functions. In particular they are meant to act on <em|one>
    generic basis function at a time, so that it seems impossible to define a
    discrete operator which acts on all basis functions of one cell.
  </remark>

  <subsection|Piecewise linear functions>

  Let <math|u<rsub|h>> be a piecewise linear function over a mesh
  <math|\<cal-T\><rsub|h>>, that is <math|u<rsub|h>\<in\>V<rsub|h>=span<around*|{|\<phi\><rsup|i><rsub|h>|}>>,
  where <math|V<rsub|h>> is a CG1 finite element space and
  <math|\<phi\><rsup|i>> are locally in <math|P<rsub|1><around*|(|T|)>> for
  some simplex <math|T\<in\>\<cal-T\><rsub|h>>. Then
  <math|u<rsub|h><rprime|'>> is piecewise constant, i.e.
  <math|u<rprime|'><rsub|h>\<in\>V<rsub|h><rprime|'>=span<around*|{|\<psi\><rsub|h><rsup|i>|}>>,
  a DG0 space over <math|\<cal-T\><rsub|h>>, where
  <math|\<psi\><rsup|i><rsub|h>> are constant over simplices. Denote by
  <math|U=<around*|(|U<rsub|i>|)>> and <math|U<rprime|'>=<around*|(|u<rprime|'><rsub|i>|)>>
  their vectors of coefficients wrt. the standard bases of nodal functions
  over <math|V<rsub|h>> and <math|V<rsub|h><rprime|'>> respectively. We want
  to compute a matrix <math|M> for the gradient operator

  <\equation*>
    \<nabla\><rsub|h>:V<rsub|h>\<rightarrow\>V<rsub|h><rprime|'>
  </equation*>

  such that <math|<around*|(|M U|)><rsub|j>=u<rprime|'><rsub|j>>.

  Since <math|u<rprime|'><rsub|h>> is not in the same space as
  <math|u<rsub|h>>, one way of computing the coefficients <math|U<rprime|'>>
  in <name|FEniCS> is to <python|df.project(df.grad(u), V)> as mentioned
  above. But here we manually compute the transformation without solving a
  linear system. Note that if we were only interested in evaluation of
  <math|u<rsub|h><rprime|'>>, because <math|u<rsub|h>=<big|sum><rsub|i>U<rsub|i>*\<phi\><rsup|i><rsub|h>>
  and by linearity <math|u<rsub|h><rprime|'>=<big|sum><rsub|i>U<rsub|i>*\<phi\><rsup|i><rsub|h><rprime|'>>.
  We could then use <python|df.evaluate_basis_derivatives()> to evaluate the
  <math|\<phi\><rsub|i><rprime|'>> and we would have finished.

  But again, we are interested in obtaining the new coefficients
  <math|U<rsub|i><rprime|'>> such that

  <\equation*>
    u<rprime|'>=<big|sum><rsub|i>U<rprime|'><rsub|i>*\<psi\><rsub|i>.
  </equation*>

  We assemble <math|M> cell by cell as follows, assuming that both <math|V>
  and <math|V<rprime|'>> are defined over the same mesh
  <math|\<cal-T\><rsub|h>> with <math|N> cells.

  <\itemize-dot>
    <item>First, for any given cell <math|T> compute the mapping
    <math|L<rsub|T>> of the dofs in <math|V<rsub|\|T>> to the dofs in
    <math|V<rprime|'><rsub|\|T>>.

    <item>Then apply this mapping to the components of <math|U<rsub|\|T>> to
    obtain the coefficients <math|U<rprime|'><rsub|\|T>>.
  </itemize-dot>

  To explain the idea, consider first a 1D mesh for simplicity and assume
  that we have an index <math|i\<in\><around*|[|N|]>> over the cells of the
  mesh such that:

  <\itemize-dot>
    <item><math|U<rsub|i>,U<rsub|i+1>> are the coefficients for the dofs in
    <math|V> over cell <math|i>.

    <item><math|U<rsub|i><rprime|'>> is the coefficient for the only dof in
    <math|V<rprime|'>> (the constant 1) over cell <math|i>.
  </itemize-dot>

  This assumption is actually a description of the mapping <math|L<rsub|T>>
  mentioned above for each interval. Then we only need to compute
  <math|\<phi\><rsup|i><rsub|h><rprime|'>> (this is not a typo, we are
  actually computing the derivative) for each <math|i> to obtain the new
  vector of coefficients <math|U<rsub|\|T><rprime|'>>:

  <\equation*>
    M<rsub|T>*U<rsub|\|T>=<around*|(|<tabular*|<tformat|<cwith|1|-1|1|1|cell-halign|c>|<cwith|1|-1|1|1|cell-lborder|0ln>|<cwith|1|-1|2|2|cell-halign|c>|<cwith|1|-1|3|3|cell-halign|c>|<cwith|1|-1|4|4|cell-halign|c>|<cwith|1|-1|5|5|cell-halign|c>|<cwith|1|-1|5|5|cell-rborder|0ln>|<table|<row|<cell|\<phi\><rsub|h><rsup|0><rprime|'>>|<cell|\<phi\><rsub|h><rsup|1><rprime|'>>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|\<phi\><rsub|h><rsup|1><rprime|'>>|<cell|\<phi\><rsub|h><rsup|2><rprime|'>>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|\<phi\><rsub|h><rsup|2><rprime|'>>|<cell|\<phi\><rsub|h><rsup|3><rprime|'>>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|\<phi\><rsub|h><rsup|3><rprime|'>>|<cell|\<phi\><rsub|h><rsup|4><rprime|'>>>>>>|)><around*|(|<tabular*|<tformat|<cwith|1|-1|1|1|cell-halign|c>|<cwith|1|-1|1|1|cell-lborder|0ln>|<cwith|1|-1|1|1|cell-rborder|0ln>|<table|<row|<cell|U<rsub|0>>>|<row|<cell|U<rsub|1>>>|<row|<cell|U<rsub|2>>>|<row|<cell|U<rsub|3>>>|<row|<cell|U<rsub|4>>>>>>|)>=<around*|(|<tabular*|<tformat|<cwith|1|-1|1|1|cell-halign|c>|<cwith|1|-1|1|1|cell-lborder|0ln>|<cwith|1|-1|1|1|cell-rborder|0ln>|<table|<row|<cell|U<rsub|0><rprime|'>>>|<row|<cell|U<rsub|1><rprime|'>>>|<row|<cell|U<rsub|2><rprime|'>>>|<row|<cell|U<rsub|3><rprime|'>>>>>>|)>.
  </equation*>

  Note that in this particular instance and with the assumptions made, we
  obtain a band matrix and computing <math|U<rprime|'>> is actually a matter
  of computing the dot product of a \Psubset\Q of the vector of constants
  <math|<around|(|\<phi\><rsub|h><rsup|0><rprime|'>,...,\<phi\><rsub|h><rsup|3><rprime|'>|)>>
  with <math|U>.

  However in other situations, e.g. in higher dimensions, the ordering of the
  nodes will not be as convenient. We proceed then more generally as follows,
  now for the global matrix. Let <math|i> run over all cell indices and let
  <math|\<iota\><rsub|i>> be its local-to-global mapping for <math|V> and
  <math|\<iota\><rprime|'><rsub|i>> for <math|V<rprime|'>>. Initialise the
  global matrix <math|M> of the discrete gradient to zero. At each step of
  the construction of <math|M> we edit row
  <math|r=\<iota\><rsub|i><rprime|'><around|(|0|)>>, i.e. the row whose
  product by <math|U> will return the coefficient for the dof in
  <math|V<rprime|'>> corresponding to cell <math|i>.

  Now, for every vertex <math|z<rsub|j>,j\<in\><around*|[|d|]>> in cell
  <math|i> compute the value of the derivative of the global dof
  <math|\<phi\><rsub|h><rsup|\<iota\><rsub|r><around|(|j|)>>> and set
  <todo|check this, I'm not sure after so many changes in notation>

  <\equation*>
    M<rsub|r*j>=\<phi\><rsub|h><rsup|\<iota\><rsub|r><around|(|j|)>><rprime|'><around*|(|z<rsub|j>|)>,j\<in\><around|[|d|]>.
  </equation*>

  Here is a straightforward implementation of this procedure using <name|AIJ>
  sparse matrices in <name|PETSc>:<\footnote>
    Note that in order to use this we need to instruct FEniCS to create the
    derivatives of the shape functions of elements with
    <python|parameters["form_compiler"]["no-evaluate_basis_derivatives"] =
    False>.
  </footnote>

  <\python-code>
    def discrete_gradient_operator(V:FunctionSpace, Vp:VectorFunctionSpace):

    \ \ \ \ """ (...) """ \ \ 

    \ \ \ \ gdim = V.mesh().geometry().dim()

    \ \ \ \ # Input validation skipped here...

    \ \ \ \ 

    \ \ \ \ dm, dmp = V.dofmap(), Vp.dofmap()

    \ \ \ \ e = V.element()

    \ \ \ \ coords = V.tabulate_dof_coordinates().reshape((-1, gdim))

    \ \ \ \ coordsp = Vp.tabulate_dof_coordinates().reshape((-1, gdim))

    \;

    \ \ \ \ Dh = PETSc.Mat()

    \ \ \ \ Dh.create(PETSc.COMM_WORLD) \ # For parallel assembly

    \ \ \ \ Dh.setSizes([Vp.dim(), V.dim()])

    \ \ \ \ Dh.setType("aij")

    \ \ \ \ Dh.setUp()

    \ \ \ \ warning("Assuming all cells have the same number of dofs")

    \ \ \ \ vals = np.zeros(dm.num_element_dofs(0)*gdim)

    \ \ \ \ # TODO: check this for parallel operation...

    \ \ \ \ istart, iend = Dh.getOwnershipRange()

    \ \ \ \ for cell_id in range(Vp.mesh().num_cells()):

    \ \ \ \ \ \ \ \ rows = dmp.cell_dofs(cell_id)

    \ \ \ \ \ \ \ \ rows = rows[(rows\<gtr\>=istart) & (rows \<less\> iend)]

    \ \ \ \ \ \ \ \ columns = dm.cell_dofs(cell_id)

    \ \ \ \ \ \ \ \ dof_coords = coords[columns,:]

    \ \ \ \ \ \ \ \ point = dof_coords[0] # Any point in the (closure of the)
    cell will do.

    \ \ \ \ \ \ \ \ e.evaluate_basis_derivatives_all(1, vals, point,
    dof_coords, 1)

    \ \ \ \ \ \ \ \ # NOTE: the copy() is necessary!

    \ \ \ \ \ \ \ \ Dh.setValues(rows, columns, vals.reshape((-1,
    gdim)).T.copy())

    \ \ \ \ Dh.assemble()

    \ \ \ \ return Dh
  </python-code>

  <subsection|Discrete Kirchhoff Triangles>

  Recall that <math|\<nabla\><rsub|h>:W<rsub|h>\<rightarrow\>\<Theta\><rsub|h>>
  is <todo|defined as...>

  The operator <math|\<nabla\><rsub|h>> can be seen as an interpolation
  operator computing the gradients of functions in <math|W<rsub|h>>.

  Interpolation into <math|P<rsub|2><rsup|2>>:

  Fix some cell <math|T\<in\>\<cal-T\><rsub|h>> in the triangulation.
  Assuming the dofs in <math|W<rsub|h>> local to <math|T> are ordered as
  <math|\<b-w\><rsub|T>=<around*|(|w<around*|(|z<rsub|1>|)>,w<rsub|,1><around*|(|z<rsub|1>|)>,w<rsub|,2><around*|(|z<rsub|1>|)>,w<around*|(|z<rsub|2>|)>,w<rsub|,1><around*|(|z<rsub|2>|)>,w<rsub|,2><around*|(|z<rsub|2>|)>,w<around*|(|z<rsub|3>|)>,w<rsub|,1><around*|(|z<rsub|3>|)>,w<rsub|,2><around*|(|z<rsub|3>|)>|)>>
  and those in <math|\<Theta\><rsub|h>> as
  <math|\<b-theta\><rsub|T>=<around*|(|\<theta\><rsub|z<rsub|1>>,\<theta\><rsub|z<rsub|2>>,\<theta\><rsub|z<rsub|3>>,\<theta\><rsub|S<rsub|1>>,\<theta\><rsub|S<rsub|2>>,\<theta\><rsub|S<rsub|3>>|)>>,
  the action of <math|\<nabla\><rsub|h>> on <math|T> is given by (see
  <cite-detail|bartels_numerical_2015|Ÿ8.2.2> for the details):

  <\equation*>
    \<b-theta\><rsub|T>=M<rsub|T>*\<b-w\><rsub|T>
  </equation*>

  with

  <\equation>
    <label|eq:local-dkt-gradient>M<rsub|T>=<matrix|<tformat|<table|<row|<cell|0>|<cell|I<rsub|2>>|<cell|0>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|I<rsub|2>>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|I<rsub|2>>>|<row|<cell|0>|<cell|0>|<cell|<wide|t|~><rsub|S<rsub|1>>>|<cell|<wide|T|~><rsub|S<rsub|1>>>|<cell|-<wide|t|~><rsub|S<rsub|1>>>|<cell|<wide|T|~><rsub|S<rsub|1>>>>|<row|<cell|<wide|t|~><rsub|S<rsub|2>>>|<cell|<wide|T|~><rsub|S<rsub|2>>>|<cell|0>|<cell|0>|<cell|-<wide|t|~><rsub|S<rsub|2>>>|<cell|<wide|T|~><rsub|S<rsub|2>>>>|<row|<cell|<wide|t|~><rsub|S<rsub|3>>>|<cell|<wide|T|~><rsub|S<rsub|3>>>|<cell|-<wide|t|~><rsub|S<rsub|3>>>|<cell|<wide|T|~><rsub|S<rsub|3>>>|<cell|0>|<cell|0>>>>>\<in\>\<bbb-R\><rsup|12\<times\>9>,
  </equation>

  where

  <\equation*>
    <wide|T|~><rsub|S<rsub|l>>=<tfrac|-3|4>*t<rsub|S<rsub|l>>*t<rsub|S<rsub|l>><rsup|\<top\>>+<with|color|red|<tfrac|1|2>*I<rsub|2>>,<space|1em><wide|t|~><rsub|S<rsub|l>>=<tfrac|-3|2*<around*|\||S<rsub|l>|\|>>*t<rsub|S<rsub|l>>,
  </equation*>

  and <math|t<rsub|S<rsub|l>>> is the vector tangent to face
  <math|S<rsub|l>>, which is the one containing the two vertices other than
  <math|z<rsub|l>> (i.e. opposite to <math|z<rsub|l>>). In the notation of
  <cite-detail|bartels_numerical_2015|Ÿ8.2>, choosing the hierarchical basis
  described there for the space <math|\<Theta\><rsub|h>>, one has for the
  coefficient at each midpoint that:

  <\equation>
    <label|eq:coeff-midpoint>\<theta\><rsub|s>=\<theta\><rsub|h><around*|(|z<rsub|s>|)>-<frac|1|2>*<around*|(|\<theta\><rsub|h><around*|(|z<rsub|s><rsup|1>|)>+\<theta\><rsub|h><around*|(|z<rsub|s><rsup|2>|)>|)>,
  </equation>

  but because we are not using a hierarchical basis for the
  <math|P<rsub|2><rsup|2>> element, the computations don't apply verbatim and
  our <math|<wide|T|~><rsub|S<rsub|l>>> differs by a constant
  <math|+<tfrac|1|2>I<rsub|2>>.

  Application of <math|\<nabla\><rsub|h>> is best done locally, during
  assembly, but it is straightforward to implement this globally: it is
  enough to construct once a global (sparse) matrix by traversing all the
  cells in the mesh and then multiply the vector of coefficients of any
  <math|w<rsub|h>\<in\>W<rsub|h>> by it.<\footnote>
    Actually this would perform many computations twice so we would rather
    only traverse the facets, carefully keeping track of the dofs involved.
  </footnote>

  <subsubsection|Computation of the local tensor <math|M>>

  As mentioned before (Remark <reference|rem:extending-ufl-discrete-gradient>),
  extending UFL is not an (easy?) option. The alternative is to assemble the
  system matrix ourselves for each problem where we need
  <math|\<nabla\><rsub|h>>. Here is one way of doing this which reuses much
  of what <name|FFC> already provides. We pick a specific very simple
  bilinear form for clarity but this can (and will) be extended to others.

  We consider the discrete problem which is detailed below in Section
  <reference|sec:linear-kirchhoff>:

  <\equation*>
    <around|(|\<nabla\>\<nabla\><rsub|h>*w<rsub|h>,\<nabla\>\<nabla\><rsub|h>*v<rsub|h>|)>=<around|(|f<rsub|3>,v<rsub|h>|)>,
  </equation*>

  where <math|w<rsub|h>\|<rsub|T>\<in\>P<rsub|3><rsup|red><around*|(|T|)>>
  for every cell <math|T> in the triangulation <math|\<cal-T\><rsub|h>>. We
  want to assemble the <strong|local tensor> for the bilinear form over the
  DKT element, <math|P<rsup|red><rsub|3><around*|(|T|)>>:

  <\equation>
    <label|eq:local-tensor-lki><around|(|\<nabla\>\<nabla\><rsub|h>*\<phi\><rsup|j><rsub|h>,\<nabla\>\<nabla\><rsub|h>*\<phi\><rsup|i><rsub|h>|)>=<around|(|f<rsub|3>,\<phi\><rsup|i><rsub|h>|)>,<text|
    for all shape functions >\<phi\><rsup|i><rsub|h>,\<phi\><rsup|j><rsub|h>\<in\>P<rsup|red><rsub|3><around*|(|T|)>.
  </equation>

  Recal that the image of <math|\<nabla\><rsub|h>\|<rsub|P<rsub|3><rsup|red>>>
  is in <math|P<rsub|2><rsup|2><around*|(|T|)>>. Let
  <math|\<theta\><rsup|l><rsub|h>,\<theta\><rsup|k><rsub|h>\<in\>P<rsub|2><rsup|2><around*|(|T|)>>
  denote shape functions and let

  <\equation*>
    A<rsup|<around*|(|2|)>><rsub|l\<nocomma\>k>=<around*|(|\<nabla\>\<theta\><rsup|k><rsub|h>,\<nabla\>\<theta\><rsup|l><rsub|h>|)>,<application-space|1em>k,l\<in\><around*|[|12|]>
  </equation*>

  be the <strong|local tensor> for the <math|P<rsub|2><rsup|2><around*|(|T|)>>
  element. Notice that this matrix is provided by the
  <cpp|ufc::cell_integral> compiled by <name|FFC>. Let
  <math|M<rsub|T>=<around*|(|m<rsub|p\<nocomma\>q>|)>\<in\>\<bbb-R\><rsup|12\<times\>9>>
  be the matrix for the local discrete gradient <eqref|eq:local-dkt-gradient>
  and take one shape function <math|\<phi\><rsup|i><rsub|h>\<in\>P<rsup|3><rsub|red><around*|(|T|)>,i\<in\><around*|[|9|]>>.
  With some abuse of notation where we confuse functions with their
  coefficient vectors in their respective finite elements we can write:

  <\equation*>
    \<nabla\><rsub|h> \<phi\><rsup|i><rsub|h>=M<rsub|T>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|\<vdots\>>>|<row|<cell|1>>|<row|<cell|\<vdots\>>>|<row|<cell|0>>>>>=m<rsub|l\<nocomma\>i>*\<theta\><rsup|l><rsub|h>,
  </equation*>

  and analogously <math|\<nabla\><rsub|h>\<phi\><rsup|j><rsub|h>=m<rsub|k\<nocomma\>j>*\<theta\><rsup|k><rsub|h>>.
  Substituting back in <eqref|eq:local-tensor-lki> we obtain

  <\equation*>
    A<rsub|i\<nocomma\>j>=<around|(|\<nabla\>\<nabla\><rsub|h>*\<phi\><rsup|j><rsub|h>,\<nabla\>\<nabla\><rsub|h>*\<phi\><rsup|i><rsub|h>|)>=m<rsub|k\<nocomma\>j>*m<rsub|l\<nocomma\>i>*<around*|(|\<nabla\>\<theta\><rsup|k><rsub|h>,\<nabla\>\<theta\><rsup|l><rsub|h>|)>=m<rsub|k\<nocomma\>j>*m<rsub|l\<nocomma\>i>*A<rsub|l\<nocomma\>k><rsup|<around*|(|2|)>>,
  </equation*>

  or in matrix form <math|A=M<rsup|\<top\>>*A<rsup|<around*|(|2|)>>*M>. The
  concrete steps required are:

  <\algorithm>
    <\enumerate>
      <item>Compile with <name|FFC> the form
      <python|inner(\<nabla\>u,\<nabla\>v)*dx> for a
      <python|VectorFunctionSpace> of type <verbatim|P2>. This provides us
      with the local cell tensor <math|A<rsup|<around|(|2|)>>> in
      <cpp|cell_integral::tabulate_tensor()>.

      <item>Compile with <name|FFC> the form
      <verbatim|inner<around|(|u,v|)>*dx> for
      <verbatim|FunctionSpace\<less\>DKT\<gtr\>>. The dofmaps computed by
      <name|FEniCS> will be the ones required to preinitialise the sparsity
      pattern of the global tensor,<\footnote>
        This is done in <cpp|AssemblerBase::init_global_tensor()>, using
        information of the form, e.g. whether it has cell, interior or
        exterior facet integrals. Since this won't change in the real form we
        are considering <todo|things should work out of the box>.
      </footnote> but this provides us with the <strong|wrong> local cell
      tensor. In order to compute the right one, denoted <math|A>:

      <\enumerate>
        <item>Compute the local matrix for <math|\<nabla\><rsub|h>> over
        simplex <math|T>: <math|M<rsub|T>=<around|(|m<rsub|p*q>|)>>.

        <item>Compute

        <\equation*>
          A<rsub|i\<nocomma\>j>=m<rsub|k*j>*m<rsub|l*i>*A<rsub|l*k><rsup|<around|(|2|)>>=<around|(|M<rsub|:i>|)><rsup|\<top\>>*A<rsup|<around|(|2|)>>*M<rsub|:j>
        </equation*>
      </enumerate>
    </enumerate>
  </algorithm>

  Here one might be tempted (or at least the author was) to compute
  <math|A<rsub|i\<nocomma\>j>> inside <verbatim|cell_integral::tabulate_tensor()>
  for the form, compile and link the new form and make it available to
  <name|dolfin> to let the it assemble the full system matrix. This proves
  technically tricky and provides little advantage over assembling the matrix
  ourselves using the tensor from 1. and the dofmaps from 2. This method
  might have the advantage that all issues related to the preinitialisation
  of the system matrices are left to <name|FEniCS>.

  <section|Linear Kirchhoff model><label|sec:linear-kirchhoff>

  Let <math|\<omega\>\<subseteq\>\<bbb-R\><rsup|2>> be a bounded Lipschitz
  domain and <math|f:\<omega\>\<rightarrow\>\<bbb-R\><rsup|3>>. As a first
  test for DKT elements and the operator <math|\<nabla\><rsub|h>>, we look
  for a vertical <strong|displacement> <math|u:\<omega\>\<rightarrow\>\<bbb-R\>>
  minimizing the energy

  <\equation*>
    I<rsub|lKi><around|(|w|)>=<frac|1|2>*<big|int><rsub|\<omega\>><around|\||D<rsup|2>*w|\|><rsup|2><math-up|d>x-<big|int><rsub|\<omega\>>f<rsub|3>*w<math-up|d>x
  </equation*>

  under homogeneous Dirichlet boundary conditions. This means solving the
  linear problem

  <\equation*>
    <around|(|D<rsup|2>*w,D<rsup|2>*v|)>=<around|(|f<rsub|3>,v|)>,
  </equation*>

  for all <math|v\<in\>H<rsup|1><rsub|0><around*|(|\<omega\>|)>>. Again
  following <cite|bartels_numerical_2015> we implement the non-conforming
  discretisation

  <\equation*>
    <around|(|\<nabla\>\<nabla\><rsub|h>*w<rsub|h>,\<nabla\>\<nabla\><rsub|h>*v<rsub|h>|)>=<around|(|f<rsub|3>,v<rsub|h>|)>,
  </equation*>

  for all <math|v<rsub|h>\<in\>W<rsub|h,D>\<assign\><around|{|w<rsub|h>\<in\>C<around|(|<wide|\<omega\>|\<bar\>>|)>:w<rsub|h>\|<rsub|T>\<in\>P<rsub|3><rsup|red><around|(|T|)>*<text|
  for all >T\<in\><with|math-font|cal|T><rsub|h><text| with continuous
  gradients at all nodes and >v<rsub|h><around*|(|z|)>=0<text| for all
  >z<text| at the boundary >\<gamma\><rsub|D>|}>> where the crux of the
  matter is the discrete gradient operator <math|\<nabla\><rsub|h>>. An
  alternative non-conforming discretisation with an interior penalty method
  is <todo|detailed in ...> for completeness (see also
  <cite|brenner_c0_2005>).

  <big-figure|<clipped|<image|img/linear-kirchhoff-cpp.eps|0.8par|||>||||4cm>|Solution
  of the linear Kirchhoff model for a clamped plate under a constant force.>

  We provide two implementations in <name|FEniCS>, a first one in
  <name|Python> which served as test-bed and a second one in <c++>. The two
  main issues where the computation of the local discrete gradient operator
  and its application during manual assembly of the system matrix and careful
  assignment of essential boundary conditions. See Appendix
  <inactive|<reference|>> for the details.

  <section|Non linear Kirchhoff model>

  As before, let <math|\<omega\>\<subseteq\>\<bbb-R\><rsup|2>> be a bounded
  Lipschitz domain and <math|f\<in\>L<rsup|1><around*|(|\<omega\>;\<bbb-R\><rsup|3>|)>>.
  We look for a minimizing <strong|deformation>
  <math|y\<in\>H<rsup|2><around*|(|\<omega\>;\<bbb-R\><rsup|3>|)>> for the
  energy<\footnote>
    Recall that in the linear case we were looking for a displacement
    <math|u> whereas here we are interested in the deformation <math|y=id+u>.
  </footnote>

  <\equation*>
    I<rsub|Ki><around|(|y|)>=<frac|1|2>*<big|int><rsub|\<omega\>><around|\||D<rsup|2>*y|\|><rsup|2><math-up|d>x-<big|int><rsub|\<omega\>>f<rsub|3>*y*<math-up|d>x
  </equation*>

  subject to the constraint that <math|y> be an <dfn|isometry>

  <\equation*>
    \<nabla\><rsup|\<top\>>y*\<nabla\>y=I<rsub|2><text| a.e. in >\<omega\>,
  </equation*>

  and under homogeneous Dirichlet boundary conditions for <math|y> and
  <math|\<nabla\>y> over a region <math|\<gamma\><rsub|D>\<subseteq\>\<omega\>>
  with positive measure. The (non-conforming) discretization we choose
  (<cite|bartels_approximation_2013>, <cite-detail|bartels_numerical_2015|Ÿ8.3>)
  uses the previously introduced spaces <math|W<rsub|h>> and
  <math|\<Theta\><rsub|h>> and makes clever use of the decoupling of function
  values and gradients in <math|W<rsub|h>>. The discrete problem is:

  <\problem*>
    <dueto|discrete Kirchhoff model>Minimize

    <\equation*>
      I<rsup|h><rsub|Ki><around*|(|y<rsub|h>|)>=<frac|1|2>*<big|int><rsub|\<omega\>><around*|\||\<nabla\>\<nabla\><rsub|h>y<rsub|h>|\|><rsup|2>*\<mathd\>x-<big|int><rsub|\<omega\>>f*\<cdot\>y<rsub|h>*\<mathd\>x
    </equation*>

    subject to

    <\equation*>
      y<rsub|h>\<in\>\<cal-A\><rsub|h>\<assign\><around*|{|v<rsub|h>\<in\>W<rsub|h,D><rsup|3>:\<nabla\><rsup|\<top\>>v<rsub|h><around*|(|z|)>*\<nabla\><rsub|h><around*|(|z|)>=I<rsub|2><text|
      for all >z\<in\>\<cal-N\><rsub|h>|}>
    </equation*>

    where

    <\equation*>
      W<rsub|h,D><rsup|3>=<around*|(|W<rsub|h,D>|)><rsup|3>=<around*|{|v<rsub|h>\<in\>W<rsub|h><rsup|3>:v<rsub|h><around*|(|z|)>=y<rsub|D><around*|(|z|)><text|
      and >\<nabla\>v<rsub|h><around*|(|z|)>=\<Phi\><rsub|D><text| for all
      >z\<in\>\<cal-N\><rsub|h>\<cap\>\<gamma\><rsub|D>|}>.
    </equation*>
  </problem*>

  For the approximation properties of this discretization, see
  <cite-detail|bartels_approximation_2013|Theorem 3.1>. It is solved using a
  discrete gradient flow of the energy <cite|bartels_approximation_2013>,<cite-detail|bartels_numerical_2015|Algorithm
  8.1>:

  <\named-algorithm|(discrete <math|H<rsup|2>>-isometry-flow)>
    Let <math|\<tau\>\<gtr\>0> and <math|y<rsub|h><rsup|0>\<in\>\<cal-A\><rsub|h>>.
    For <math|k=1,2,\<ldots\>>, define

    <\equation*>
      \<cal-F\><rsub|h><around*|[|y<rsup|k-1><rsub|h>|]>\<assign\><around*|{|w<rsub|h>\<in\>W<rsub|h,D><rsup|3>:\<nabla\><rsup|\<top\>>w<rsub|h><around*|(|z|)>*\<nabla\>y<rsup|k-1><rsub|h><around*|(|z|)>+\<nabla\><rsup|\<top\>>y<rsup|k-1><rsub|h><around*|(|z|)>*\<nabla\>w<rsub|h><around*|(|z|)>=0<text|
      f.a. >z\<in\>\<cal-N\><rsub|h>|}>
    </equation*>

    and set <math|y<rsup|k><rsub|h>=y<rsup|k-1><rsub|h>+\<tau\>*d<rsub|t>y<rsub|h><rsup|k>>
    where the <em|update> <math|d<rsub|t>y<rsup|k><rsub|h>> satisfies the
    constraint <math|d<rsub|t> y<rsup|k><rsub|h>\<in\>\<cal-F\><rsub|h><around*|[|y<rsup|k-1><rsub|h>|]>>
    and solves the system

    <\equation>
      <label|eq:h2-flow-update-system><around*|(|\<nabla\>\<nabla\><rsub|h>
      d<rsub|t>y<rsup|k><rsub|h>,\<nabla\>\<nabla\><rsub|h>w<rsub|h>|)>+\<alpha\>*<around*|(|\<nabla\>\<nabla\><rsub|h>
      <around*|(|y<rsup|k-1><rsub|h>+\<tau\>*d<rsub|t>y<rsub|h>|)>,\<nabla\>\<nabla\><rsub|h>
      w<rsub|h>|)>=<around*|(|f,w<rsub|h>|)>
    </equation>

    for all <math|w<rsub|h>\<in\>\<cal-F\><rsub|h><around*|[|y<rsup|k-1><rsub|h>|]>>.
    Stop the iteration if <math|<around*|\<\|\|\>|\<nabla\>\<nabla\><rsub|h>
    d<rsub|t> y<rsup|k><rsub|h>|\<\|\|\>>\<leqslant\>\<varepsilon\><rsub|stop>>.
  </named-algorithm>

  <cite-detail|bartels_approximation_2013|Theorem 3.2> proves that the
  iterates <math|<around*|(|y<rsup|k><rsub|h>|)><rsub|k=0,1,\<ldots\>>>
  produced by this algorithm reduce the energy and approximate a solution
  even if they might not always fulfil the isometry constraint.

  If we express the constraint <math|d<rsub|t>
  y<rsup|k><rsub|h>\<in\>\<cal-F\><rsub|h><around*|[|y<rsup|k-1><rsub|h>|]>>
  component-wise we have (we drop the subindex <math|h> for notational
  simplicity)

  <\equation>
    <label|eq:nodal-isometry-constraint-components>d<rsub|t>
    y<rsup|k><rsub|l,i>*y<rsup|k-1><rsub|l,j>+y<rsub|l,i><rsup|k-1>*d<rsub|t>
    y<rsup|k><rsub|l,j>=0<text| for each >i,j\<in\><around*|{|1,2|}><text| at
    every node >z\<in\>\<cal-N\><rsub|h>
  </equation>

  and\ 

  <\equation>
    <label|eq:nodal-dirichlet-constraint>d<rsub|t> y<rsup|k><rsub|l>=0<text|
    at every node >z\<in\>\<cal-N\><rsub|h>\<cap\>\<Gamma\><rsub|D>.<todo|fix
    notation>
  </equation>

  Using the local dofmaps we translate <eqref|eq:nodal-isometry-constraint-components>
  into matrix <math|\<cdot\>> vector form, for one cell. First, assume that
  the <verbatim|VectorFunctionSpace \<less\>DKT\<gtr\>> sorts the local dofs
  of its subspaces simply by concatenating them, that is: for a fixed cell,
  the local indices 0-8 are for the first subspace, 9-17 for the second and
  18-26 for the third.<\footnote>
    This is actually the case, at least in <name|FEniCS> 2017.1
  </footnote> Then we have that the local components
  <math|d<rsub|t>Y<rsup|k>> of <math|d<rsub|t> y<rsup|k><rsub|h>> are
  <math|d<rsub|t>Y<rsup|k><rsub|3*l+n*i>=d<rsub|t>y<rsub|l,i><rsup|k><around*|(|z<rsub|n>|)>>,
  <math|l\<in\><around*|{|0,1,2|}>>, explicitly:

  <\equation*>
    d<rsub|t> Y<rsup|k>=<around*|(|d<rsub|t>y<rsub|1><rsup|k><around*|(|z<rsub|1>|)>,d<rsub|t>y<rsub|1,1><rsup|k><around*|(|z<rsub|1>|)>,d<rsub|t>
    y<rsub|1,2><rsup|k><around*|(|z<rsub|1>|)>,\<ldots\>,d<rsub|t>y<rsub|3><rsup|k><around*|(|z<rsub|3>|)>,d<rsub|t>y<rsub|3,1><rsup|k><around*|(|z<rsub|3>|)>,d<rsub|t>y<rsub|3,2><rsup|k><around*|(|z<rsub|3>|)>|)>.
  </equation*>

  We construct a matrix <math|B<rsup|k-1>> whose rows will realize each one
  of the components of the constraint. By inspecting the products in
  <eqref|eq:nodal-isometry-constraint-components> and
  <eqref|eq:nodal-dirichlet-constraint> we can construct the matrix for one
  cell

  <\equation*>
    B<rsup|k-1>\<assign\><matrix|<tformat|<table|<row|<cell|B<rsub|1><rsup|k-1>>|<cell|B<rsub|2><rsup|k-1>>|<cell|B<rsub|3><rsup|k-1>>>>>>\<in\>\<bbb-R\><rsup|7\<times\>27>
  </equation*>

  whose components are given in Table <reference|tab:matrix-bk>.<\footnote>
    The full matrix is

    <\equation*>
      <tabular|<tformat|<cwith|1|1|12|12|cell-row-span|1>|<cwith|1|1|12|12|cell-col-span|9>|<cwith|1|1|21|21|cell-row-span|1>|<cwith|1|1|21|21|cell-col-span|9>|<cwith|2|2|3|3|cell-row-span|1>|<cwith|2|2|3|3|cell-col-span|3>|<cwith|1|1|3|3|cell-row-span|1>|<cwith|1|1|3|3|cell-col-span|9>|<cwith|2|2|6|6|cell-row-span|1>|<cwith|2|2|6|6|cell-col-span|3>|<cwith|2|2|9|9|cell-row-span|1>|<cwith|2|2|9|9|cell-col-span|3>|<cwith|2|2|12|12|cell-row-span|1>|<cwith|2|2|12|12|cell-col-span|3>|<cwith|2|2|15|15|cell-row-span|1>|<cwith|2|2|15|15|cell-col-span|3>|<cwith|2|2|18|18|cell-row-span|1>|<cwith|2|2|18|18|cell-col-span|3>|<cwith|2|2|21|21|cell-row-span|1>|<cwith|2|2|21|21|cell-col-span|3>|<cwith|2|2|24|24|cell-row-span|1>|<cwith|2|2|24|24|cell-col-span|3>|<cwith|2|2|27|27|cell-row-span|1>|<cwith|2|2|27|27|cell-col-span|3>|<cwith|3|3|6|6|cell-row-span|1>|<cwith|3|3|6|6|cell-col-span|3>|<cwith|3|3|9|9|cell-row-span|1>|<cwith|3|3|9|9|cell-col-span|3>|<cwith|3|3|15|15|cell-row-span|1>|<cwith|3|3|15|15|cell-col-span|3>|<cwith|3|3|18|18|cell-row-span|1>|<cwith|3|3|18|18|cell-col-span|3>|<cwith|3|3|24|24|cell-row-span|1>|<cwith|3|3|24|24|cell-col-span|3>|<cwith|3|3|27|27|cell-row-span|1>|<cwith|3|3|27|27|cell-col-span|3>|<cwith|3|8|4|4|cell-background|pastel
      green>|<cwith|4|8|7|7|cell-background|pastel
      green>|<cwith|4|8|10|10|cell-background|pastel
      green>|<cwith|4|8|8|8|cell-background|pastel
      orange>|<cwith|4|8|5|5|cell-background|pastel
      orange>|<cwith|4|8|11|11|cell-background|pastel
      orange>|<cwith|6|6|4|4|cell-background|pastel
      green>|<cwith|6|6|7|7|cell-background|pastel
      green>|<cwith|6|6|10|10|cell-background|pastel
      green>|<cwith|6|6|8|8|cell-background|pastel
      orange>|<cwith|6|6|5|5|cell-background|pastel
      orange>|<cwith|6|6|11|11|cell-background|pastel
      orange>|<cwith|4|4|3|29|cell-tborder|0ln>|<cwith|3|3|3|29|cell-bborder|0ln>|<cwith|4|4|3|29|cell-bborder|1ln>|<cwith|5|5|3|29|cell-tborder|1ln>|<cwith|4|4|3|3|cell-lborder|0ln>|<cwith|4|4|29|29|cell-rborder|0ln>|<cwith|5|5|3|29|cell-tsep|0.3em>|<cwith|3|3|4|4|cell-background|pastel
      green>|<cwith|3|3|5|5|cell-background|pastel
      orange>|<cwith|5|8|29|29|cell-lborder|0ln>|<cwith|5|8|28|28|cell-rborder|0ln>|<cwith|5|8|29|29|cell-rborder|1ln>|<cwith|1|-1|1|-1|cell-halign|c>|<cwith|5|5|2|2|cell-row-span|4>|<cwith|5|5|2|2|cell-col-span|1>|<cwith|4|4|3|3|cell-bborder|0ln>|<cwith|5|8|2|5|cell-bborder|0ln>|<cwith|5|8|2|5|cell-rborder|0ln>|<cwith|5|8|2|5|cell-lborder|0ln>|<cwith|5|5|29|29|cell-tborder|0ln>|<cwith|4|4|29|29|cell-bborder|0ln>|<cwith|5|5|29|29|cell-lborder|0ln>|<cwith|5|5|28|28|cell-rborder|0ln>|<cwith|5|5|29|29|cell-rborder|0ln>|<cwith|8|8|29|29|cell-bborder|0ln>|<cwith|8|8|29|29|cell-lborder|0ln>|<cwith|8|8|28|28|cell-rborder|0ln>|<cwith|8|8|29|29|cell-rborder|0ln>|<cwith|6|6|29|29|cell-tborder|0ln>|<cwith|5|5|29|29|cell-bborder|0ln>|<cwith|6|6|29|29|cell-lborder|0ln>|<cwith|6|6|28|28|cell-rborder|0ln>|<cwith|6|6|29|29|cell-rborder|0ln>|<cwith|7|7|29|29|cell-tborder|0ln>|<cwith|6|6|29|29|cell-bborder|0ln>|<cwith|7|7|29|29|cell-bborder|0ln>|<cwith|8|8|29|29|cell-tborder|0ln>|<cwith|7|7|29|29|cell-lborder|0ln>|<cwith|7|7|28|28|cell-rborder|0ln>|<cwith|7|7|29|29|cell-rborder|0ln>|<cwith|5|5|30|30|cell-row-span|4>|<cwith|5|5|30|30|cell-col-span|1>|<cwith|4|4|4|28|cell-tborder|0ln>|<cwith|3|3|4|28|cell-bborder|0ln>|<cwith|4|4|4|28|cell-bborder|0ln>|<cwith|5|8|2|28|cell-tborder|0ln>|<cwith|4|4|4|4|cell-lborder|0ln>|<cwith|4|4|3|3|cell-rborder|0ln>|<cwith|4|4|28|28|cell-rborder|0ln>|<cwith|4|4|29|29|cell-lborder|0ln>|<cwith|4|4|1|-1|cell-bsep|0px>|<cwith|4|4|1|-1|cell-tsep|0px>|<cwith|4|4|1|-1|cell-valign|c>|<cwith|4|4|1|-1|cell-height|0.8em>|<cwith|4|4|1|-1|cell-vmode|exact>|<cwith|5|5|2|2|cell-bsep|0px>|<cwith|5|5|2|2|cell-tsep|0px>|<cwith|5|5|30|30|cell-bsep|0px>|<cwith|5|5|30|30|cell-tsep|0px>|<cwith|8|8|1|-1|cell-bsep|0px>|<cwith|3|8|13|13|cell-background|pastel
      green>|<cwith|3|8|22|22|cell-background|pastel
      green>|<cwith|3|8|14|14|cell-background|pastel
      orange>|<cwith|3|8|23|23|cell-background|pastel
      orange>|<cwith|5|8|16|16|cell-background|pastel
      green>|<cwith|5|8|19|19|cell-background|pastel
      green>|<cwith|5|8|25|25|cell-background|pastel
      green>|<cwith|5|8|28|28|cell-background|pastel
      green>|<cwith|5|8|17|17|cell-background|pastel
      orange>|<cwith|5|8|20|20|cell-background|pastel
      orange>|<cwith|5|8|26|26|cell-background|pastel
      orange>|<cwith|5|8|29|29|cell-background|pastel
      orange>|<cwith|5|5|2|2|cell-valign|b>|<cwith|5|5|30|30|cell-valign|b>|<cwith|3|3|12|12|cell-tborder|0ln>|<cwith|2|2|12|12|cell-bborder|0ln>|<cwith|8|8|12|12|cell-bborder|0ln>|<cwith|3|8|12|12|cell-lborder|0.5ln>|<cwith|3|8|11|11|cell-rborder|0.5ln>|<cwith|3|8|12|12|cell-rborder|0ln>|<cwith|3|8|13|13|cell-lborder|0ln>|<cwith|3|3|21|21|cell-tborder|0ln>|<cwith|2|2|21|21|cell-bborder|0ln>|<cwith|8|8|21|21|cell-bborder|0ln>|<cwith|3|8|21|21|cell-lborder|0.5ln>|<cwith|3|8|20|20|cell-rborder|0.5ln>|<cwith|3|8|21|21|cell-rborder|0ln>|<cwith|3|8|22|22|cell-lborder|0ln>|<cwith|4|4|30|30|cell-background|>|<cwith|4|4|3|29|cell-background|>|<cwith|5|5|1|1|cell-row-span|4>|<cwith|5|5|1|1|cell-col-span|1>|<cwith|5|5|1|1|cell-valign|c>|<twith|table-lsep|0px>|<cwith|1|-1|1|-1|cell-lsep|0.1em>|<cwith|1|-1|1|-1|cell-rsep|0.1em>|<cwith|5|8|6|11|cell-lsep|0.2em>|<cwith|5|8|6|11|cell-rsep|0.2em>|<cwith|5|8|15|19|cell-lsep|0.2em>|<cwith|5|8|15|19|cell-rsep|0.2em>|<cwith|5|8|24|29|cell-lsep|0.2em>|<cwith|5|8|24|29|cell-rsep|0.2em>|<cwith|1|-1|20|20|cell-lsep|0.2em>|<cwith|1|-1|20|20|cell-rsep|0.2em>|<cwith|1|2|3|27|color|darker
      grey>|<cwith|4|4|3|23|color|darker grey>|<cwith|3|3|6|9|color|darker
      grey>|<cwith|3|3|15|18|color|darker grey>|<cwith|3|3|24|27|color|darker
      grey>|<cwith|3|3|3|27|math-level|1>|<cwith|2|2|3|27|math-level|1>|<cwith|1|1|3|21|math-level|1>|<cwith|1|1|1|-1|cell-bsep|0px>|<cwith|2|2|1|-1|cell-tsep|0px>|<cwith|1|-1|2|2|cell-rsep|0px>|<cwith|1|-1|30|30|cell-lsep|0px>|<cwith|6|6|7|7|cell-valign|c>|<cwith|6|6|8|8|cell-row-span|2>|<cwith|6|6|8|8|cell-col-span|1>|<cwith|6|6|7|7|cell-row-span|2>|<cwith|6|6|7|7|cell-col-span|1>|<cwith|6|6|8|8|cell-valign|c>|<cwith|6|6|11|11|cell-row-span|2>|<cwith|6|6|11|11|cell-col-span|1>|<cwith|6|6|10|10|cell-row-span|2>|<cwith|6|6|10|10|cell-col-span|1>|<cwith|6|6|10|10|cell-valign|c>|<cwith|6|6|11|11|cell-valign|c>|<table|<row|<cell|>|<cell|>|<cell|<text|subspace
      0>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|<text|subspace
      1>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|<text|subspace
      2>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|<text|eval.
      at >z<rsub|1>>|<cell|>|<cell|>|<cell|z<rsub|2>>|<cell|>|<cell|>|<cell|z<rsub|3>>|<cell|>|<cell|>|<cell|<text|eval.
      at >z<rsub|1>>|<cell|>|<cell|>|<cell|z<rsub|2>>|<cell|>|<cell|>|<cell|z<rsub|3>>|<cell|>|<cell|>|<cell|<text|eval.
      at >z<rsub|1>>|<cell|>|<cell|>|<cell|z<rsub|2>>|<cell|>|<cell|>|<cell|z<rsub|3>>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|y<rsup|k><rsub|1>>|<cell|y<rsup|k><rsub|1,1>>|<cell|y<rsup|k><rsub|1,2>>|<cell|\<cdots\>>|<cell|>|<cell|>|<cell|\<cdots\>>|<cell|>|<cell|>|<cell|y<rsup|k><rsub|2>>|<cell|y<rsup|k><rsub|2,1>>|<cell|y<rsup|k><rsub|2,2>>|<cell|\<cdots\>>|<cell|>|<cell|>|<cell|\<cdots\>>|<cell|>|<cell|>|<cell|y<rsup|k><rsub|3>>|<cell|y<rsup|k><rsub|3,1>>|<cell|y<rsup|k><rsub|3,2>>|<cell|\<cdots\>>|<cell|>|<cell|>|<cell|\<cdots\>>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|<very-small|0>>|<cell|<very-small|1>>|<cell|<very-small|2>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|<very-small|9>>|<cell|<very-small|10>>|<cell|<very-small|11>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|<very-small|18>>|<cell|<very-small|19>>|<cell|<very-small|20>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|B<rsub|k-1>\<assign\><rsub|>>|<cell|<smash|<shift|<around*|(|<resize|||||5.2em>|\<nobracket\>>||0.3em>>>|<cell|>|<cell|2*y<rsup|k-1><rsub|1,1>>|<cell|0>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|2*y<rsup|k-1><rsub|2,1>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|2*y<rsup|k-1><rsub|3,1>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|<shift|||-0.2em><smash|<shift|<around*|\<nobracket\>|<resize|||||5.2em>|)>||0.3em>>>>|<row|<cell|>|<cell|>|<cell|>|<cell|y<rsup|k-1><rsub|1,2>>|<cell|y<rsup|k-1><rsub|1,1>>|<cell|>|<cell|\<vdots\>>|<cell|\<vdots\>>|<cell|>|<cell|\<vdots\>>|<cell|\<vdots\>>|<cell|>|<cell|y<rsup|k-1><rsub|2,2>>|<cell|y<rsup|k-1><rsub|2,1>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|y<rsup|k-1><rsub|3,2>>|<cell|y<rsup|k-1><rsub|3,1>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|y<rsup|k-1><rsub|1,2>>|<cell|y<rsup|k-1><rsub|1,1>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|y<rsup|k-1><rsub|2,2>>|<cell|y<rsup|k-1><rsub|2,1>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|y<rsup|k-1><rsub|3,2>>|<cell|y<rsup|k-1><rsub|3,1>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|0>|<cell|2*y<rsup|k-1><rsub|1,2>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|2*y<rsup|k-1><rsub|2,2>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|2*y<rsup|k-1><rsub|3,2>>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>>>>
    </equation*>
  </footnote> Row 1 encodes <eqref|eq:nodal-isometry-constraint-components>
  for <math|i=j=1>, rows 2,3 encode <math|i=1,2,j=2,1> and row 4 encodes
  <math|i=j=2> Note that we have three times all the products in
  <eqref|eq:nodal-isometry-constraint-components>, once per node in a cell.

  <\big-table>
    <\equation*>
      <matrix|<tformat|<cwith|1|4|2|2|cell-background|pastel
      green>|<cwith|1|4|5|5|cell-background|pastel
      green>|<cwith|1|4|8|8|cell-background|pastel
      green>|<cwith|1|4|6|6|cell-background|pastel
      orange>|<cwith|1|4|3|3|cell-background|pastel
      orange>|<cwith|1|4|9|9|cell-background|pastel
      orange>|<cwith|2|2|2|2|cell-background|pastel
      green>|<cwith|1|1|1|9|cell-tborder|1ln>|<cwith|1|1|1|9|cell-tsep|0.3em>|<cwith|1|4|1|9|cell-halign|c>|<cwith|1|4|1|3|cell-bborder|0ln>|<cwith|1|4|1|3|cell-rborder|0ln>|<cwith|1|4|1|3|cell-lborder|0ln>|<cwith|1|4|1|9|cell-tborder|0ln>|<cwith|4|4|1|9|cell-bsep|0px>|<cwith|1|4|9|9|cell-rborder|0.5ln>|<cwith|1|4|1|9|cell-lsep|0.1em>|<cwith|1|4|1|9|cell-rsep|0.1em>|<cwith|1|4|4|9|cell-lsep|0.2em>|<cwith|1|4|4|9|cell-rsep|0.2em>|<cwith|2|2|5|5|cell-valign|c>|<cwith|2|2|6|6|cell-valign|c>|<cwith|2|2|8|8|cell-valign|c>|<cwith|2|2|9|9|cell-valign|c>|<cwith|4|4|8|9|cell-bborder|0ln>|<cwith|1|4|8|9|cell-lborder|0ln>|<cwith|1|4|7|7|cell-rborder|0ln>|<cwith|1|4|8|9|cell-rborder|0ln>|<cwith|2|2|6|6|cell-row-span|1>|<cwith|2|2|6|6|cell-col-span|1>|<cwith|2|2|9|9|cell-row-span|1>|<cwith|2|2|9|9|cell-col-span|1>|<cwith|2|2|8|8|cell-row-span|1>|<cwith|2|2|8|8|cell-col-span|1>|<cwith|2|2|5|5|cell-row-span|1>|<cwith|2|2|5|5|cell-col-span|1>|<cwith|1|1|6|6|cell-tborder|1ln>|<cwith|1|1|6|6|cell-tsep|0.3em>|<cwith|1|4|6|6|cell-halign|c>|<cwith|1|4|6|6|cell-bborder|0ln>|<cwith|1|4|6|6|cell-rborder|0ln>|<cwith|1|4|6|6|cell-lborder|0ln>|<cwith|1|4|6|6|cell-tborder|0ln>|<cwith|4|4|6|6|cell-bsep|0px>|<cwith|1|4|6|6|cell-lsep|0.1em>|<cwith|1|4|6|6|cell-rsep|0.1em>|<cwith|1|1|9|9|cell-tborder|1ln>|<cwith|1|1|9|9|cell-tsep|0.3em>|<cwith|1|4|9|9|cell-halign|c>|<cwith|1|4|9|9|cell-bborder|0ln>|<cwith|1|4|9|9|cell-rborder|0ln>|<cwith|1|4|9|9|cell-lborder|0ln>|<cwith|1|4|9|9|cell-tborder|0ln>|<cwith|4|4|9|9|cell-bsep|0px>|<cwith|1|4|9|9|cell-lsep|0.1em>|<cwith|1|4|9|9|cell-rsep|0.1em>|<cwith|1|4|3|3|cell-background|pastel
      green>|<cwith|1|4|5|6|cell-background|pastel
      orange>|<cwith|1|4|8|9|cell-background|pastel
      yellow>|<cwith|1|4|7|7|cell-background|pastel
      yellow>|<cwith|4|4|4|4|cell-background|pastel
      orange>|<cwith|1|3|4|4|cell-background|pastel
      orange>|<cwith|4|4|1|1|cell-background|pastel
      green>|<cwith|3|3|1|1|cell-background|pastel
      green>|<cwith|2|2|1|1|cell-background|pastel
      green>|<cwith|1|1|1|1|cell-background|pastel
      green>|<cwith|5|7|1|3|cell-background|pastel
      green>|<cwith|5|7|4|6|cell-background|pastel
      orange>|<cwith|5|7|7|9|cell-background|pastel
      yellow>|<table|<row|<cell|<application-space|1em>>|<cell|2*y<rsup|k-1><rsub|i,1>>|<cell|>|<cell|<application-space|1em>>|<cell|2*y<rsup|k-1><rsub|i,1>>|<cell|>|<cell|<application-space|1em>>|<cell|2*y<rsup|k-1><rsub|i,1>>|<cell|>>|<row|<cell|>|<cell|y<rsup|k-1><rsub|i,2>>|<cell|y<rsup|k-1><rsub|i,1>>|<cell|>|<cell|y<rsup|k-1><rsub|i,2>>|<cell|y<rsup|k-1><rsub|i,1>>|<cell|>|<cell|y<rsup|k-1><rsub|i,2>>|<cell|y<rsup|k-1><rsub|i,1>>>|<row|<cell|>|<cell|y<rsup|k-1><rsub|i,2>>|<cell|y<rsup|k-1><rsub|i,1>>|<cell|>|<cell|y<rsup|k-1><rsub|i,2>>|<cell|y<rsup|k-1><rsub|i,1>>|<cell|>|<cell|y<rsup|k-1><rsub|i,2>>|<cell|y<rsup|k-1><rsub|i,1>>>|<row|<cell|>|<cell|>|<cell|2*y<rsup|k-1><rsub|i,2>>|<cell|>|<cell|>|<cell|2*y<rsup|k-1><rsub|i,2>>|<cell|>|<cell|>|<cell|2*y<rsup|k-1><rsub|i,2>>>|<row|<cell|1>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|1>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|>|<cell|1>|<cell|>|<cell|>>>>>*<matrix|<tformat|<cwith|1|1|1|2|cell-bborder|0ln>|<cwith|1|1|1|2|cell-halign|c>|<cwith|1|1|1|2|cell-lsep|0.1em>|<cwith|1|1|1|2|cell-rsep|0.1em>|<cwith|1|1|1|2|math-level|1>|<cwith|1|1|1|2|color|darker
      grey>|<cwith|2|2|1|2|cell-bborder|0ln>|<cwith|2|2|1|2|cell-halign|c>|<cwith|2|2|1|2|cell-bborder|0ln>|<cwith|2|2|1|2|cell-lsep|0.1em>|<cwith|2|2|1|2|cell-rsep|0.1em>|<cwith|2|2|1|2|math-level|1>|<cwith|2|2|1|2|color|darker
      grey>|<cwith|3|3|1|2|cell-bborder|0ln>|<cwith|3|3|1|2|cell-halign|c>|<cwith|3|3|1|2|cell-bborder|0ln>|<cwith|3|3|1|2|cell-lsep|0.1em>|<cwith|3|3|1|2|cell-rsep|0.1em>|<cwith|3|3|1|2|math-level|1>|<cwith|3|3|1|2|color|darker
      grey>|<cwith|4|4|1|2|cell-bborder|0ln>|<cwith|4|4|1|2|cell-halign|c>|<cwith|4|4|1|2|cell-bborder|0ln>|<cwith|4|4|1|2|cell-lsep|0.1em>|<cwith|4|4|1|2|cell-rsep|0.1em>|<cwith|4|4|1|2|math-level|1>|<cwith|4|4|1|1|cell-row-span|1>|<cwith|4|4|1|1|cell-col-span|1>|<cwith|4|4|1|2|cell-bborder|0ln>|<cwith|4|4|1|2|cell-halign|c>|<cwith|4|4|1|2|cell-lsep|0.1em>|<cwith|4|4|1|2|cell-rsep|0.1em>|<cwith|4|4|1|2|math-level|1>|<cwith|4|4|1|2|color|darker
      grey>|<cwith|5|5|1|2|cell-bborder|0ln>|<cwith|5|5|1|2|cell-halign|c>|<cwith|5|5|1|2|cell-bborder|0ln>|<cwith|5|5|1|2|cell-lsep|0.1em>|<cwith|5|5|1|2|cell-rsep|0.1em>|<cwith|5|5|1|2|math-level|1>|<cwith|5|5|1|2|cell-bborder|0ln>|<cwith|5|5|1|2|cell-halign|c>|<cwith|5|5|1|2|cell-bborder|0ln>|<cwith|5|5|1|2|cell-lsep|0.1em>|<cwith|5|5|1|2|cell-rsep|0.1em>|<cwith|5|5|1|2|math-level|1>|<cwith|5|5|1|2|color|darker
      grey>|<cwith|6|6|1|2|cell-bborder|0ln>|<cwith|6|6|1|2|cell-halign|c>|<cwith|6|6|1|2|cell-bborder|0ln>|<cwith|6|6|1|2|cell-lsep|0.1em>|<cwith|6|6|1|2|cell-rsep|0.1em>|<cwith|6|6|1|2|math-level|1>|<cwith|6|6|1|1|cell-bborder|0ln>|<cwith|6|6|1|1|cell-halign|c>|<cwith|6|6|1|1|cell-bborder|0ln>|<cwith|6|6|1|1|cell-lsep|0.1em>|<cwith|6|6|1|1|cell-rsep|0.1em>|<cwith|6|6|1|1|math-level|1>|<cwith|6|6|1|2|color|darker
      grey>|<cwith|7|7|1|2|cell-bborder|0ln>|<cwith|7|7|1|2|cell-halign|c>|<cwith|7|7|1|2|cell-bborder|0ln>|<cwith|7|7|1|2|cell-lsep|0.1em>|<cwith|7|7|1|2|cell-rsep|0.1em>|<cwith|7|7|1|2|math-level|1>|<cwith|7|7|1|1|cell-row-span|1>|<cwith|7|7|1|1|cell-col-span|1>|<cwith|7|7|1|2|cell-bborder|0ln>|<cwith|7|7|1|2|cell-halign|c>|<cwith|7|7|1|2|cell-lsep|0.1em>|<cwith|7|7|1|2|cell-rsep|0.1em>|<cwith|7|7|1|2|math-level|1>|<cwith|7|7|1|2|color|darker
      grey>|<cwith|8|8|1|2|cell-bborder|0ln>|<cwith|8|8|1|2|cell-halign|c>|<cwith|8|8|1|2|cell-bborder|0ln>|<cwith|8|8|1|2|cell-lsep|0.1em>|<cwith|8|8|1|2|cell-rsep|0.1em>|<cwith|8|8|1|2|math-level|1>|<cwith|8|8|1|2|cell-bborder|0ln>|<cwith|8|8|1|1|cell-rborder|0ln>|<cwith|8|8|1|2|cell-bborder|0ln>|<cwith|8|8|1|2|cell-halign|c>|<cwith|8|8|1|2|cell-bborder|0ln>|<cwith|8|8|1|2|cell-lsep|0.1em>|<cwith|8|8|1|2|cell-rsep|0.1em>|<cwith|8|8|1|2|math-level|1>|<cwith|8|8|1|2|color|darker
      grey>|<cwith|9|9|1|1|cell-bborder|0ln>|<cwith|9|9|1|1|cell-halign|c>|<cwith|9|9|1|1|cell-bborder|0ln>|<cwith|9|9|1|1|cell-rborder|0.5ln>|<cwith|9|9|1|1|cell-lsep|0.1em>|<cwith|9|9|1|1|cell-rsep|0.1em>|<cwith|9|9|1|1|math-level|1>|<cwith|9|9|1|1|cell-bborder|0ln>|<cwith|9|9|1|1|cell-tborder|0ln>|<cwith|9|9|1|1|cell-bborder|0ln>|<cwith|9|9|1|1|cell-lborder|0ln>|<cwith|9|9|1|1|cell-rborder|0ln>|<cwith|9|9|1|1|cell-bborder|0ln>|<cwith|9|9|1|1|cell-halign|c>|<cwith|9|9|1|1|cell-bborder|0ln>|<cwith|9|9|1|1|cell-lsep|0.1em>|<cwith|9|9|1|1|cell-rsep|0.1em>|<cwith|9|9|1|1|math-level|1>|<cwith|9|9|1|1|color|darker
      grey>|<cwith|1|3|1|1|cell-background|pastel
      green>|<cwith|4|6|1|1|cell-background|pastel
      orange>|<cwith|7|9|1|1|cell-background|pastel
      yellow>|<table|<row|<cell|d<rsub|t>
      y<rsup|k><rsub|i>>>|<row|<cell|d<rsub|t>
      y<rsup|k><rsub|i,1>>>|<row|<cell|d<rsub|t>
      y<rsup|k><rsub|i,2>>>|<row|<cell|d<rsub|t>
      y<rsup|k><rsub|i>>>|<row|<cell|d<rsub|t>
      y<rsup|k><rsub|i,1>>>|<row|<cell|d<rsub|t>
      y<rsup|k><rsub|i,2>>>|<row|<cell|d<rsub|t>
      y<rsup|k><rsub|i>>>|<row|<cell|d<rsub|t>
      y<rsup|k><rsub|i,1>>>|<row|<cell|d<rsub|t> y<rsup|k><rsub|i,2>>>>>>
    </equation*>
  </big-table|<label|tab:matrix-bk>One of the (local) products
  <math|B<rsub|i><rsup|k-1>*d<rsub|t><rsub|>Y<rsup|k><rsub|i>> (for subspace
  <math|W<rsup|<around*|(|i|)>><rsub|h>>) for a single element. Colours
  denote evaluation at each of the three nodes of a triangle. Empty entries
  are zeroes. Column indices correspond to local dof indices:
  <math|0,1,2,\<ldots\>>>

  For instance, the products <eqref|eq:nodal-isometry-constraint-components>
  for <math|i=1,j=2> are given by row 2 of all 3 matrices
  <math|B<rsup|k-1><rsub|i>>:

  <\eqnarray*>
    <tformat|<table|<row|<cell|<around*|(|B<rsup|k-1>|)><rsub|2l>*d<rsub|t
    >Y<rsup|k><rsub|l>>|<cell|=>|<cell|y<rsup|k-1><rsub|1,2><around*|(|z<rsub|1>|)>*d<rsub|t
    >y<rsup|k><rsub|1,1><around*|(|z<rsub|1>|)>+y<rsup|k-1><rsub|1,1><around*|(|z<rsub|1>|)>*d<rsub|t
    >y<rsup|k><rsub|1,2><around*|(|z<rsub|1>|)>>>|<row|<cell|>|<cell|>|<cell|+y<rsup|k-1><rsub|1,2><around*|(|z<rsub|2>|)>*d<rsub|t
    >y<rsup|k><rsub|1,1><around*|(|z<rsub|2>|)>+y<rsup|k-1><rsub|1,1><around*|(|z<rsub|2>|)>*d<rsub|t
    >y<rsup|k><rsub|1,2><around*|(|z<rsub|2>|)>>>|<row|<cell|>|<cell|>|<cell|+y<rsup|k-1><rsub|1,2><around*|(|z<rsub|3>|)>*d<rsub|t
    >y<rsup|k><rsub|1,1><around*|(|z<rsub|3>|)>+y<rsup|k-1><rsub|1,1><around*|(|z<rsub|3>|)>*d<rsub|t
    >y<rsup|k><rsub|1,2><around*|(|z<rsub|3>|)>>>|<row|<cell|>|<cell|>|<cell|+\<ldots\>>>>>
  </eqnarray*>

  Using <math|B<rsup|k-1>> the (local) linear system to solve
  <eqref|eq:h2-flow-update-system> at time step <math|k> is:

  <\equation>
    <label|eq:kirchhoff-local-timestep-system><matrix|<tformat|<table|<row|<cell|<around*|(|1+\<alpha\>*\<tau\>|)>*M<rsup|\<top\>>*A<rsup|<around*|(|2|)>>*M>|<cell|<around*|(|B<rsup|k-1>|)><rsup|\<top\>>>>|<row|<cell|B<rsup|k-1>>|<cell|0>>>>>*<matrix|<tformat|<table|<row|<cell|d<rsub|t>
    Y<rsup|k>>>|<row|<cell|\<Lambda\>>>>>>=<matrix|<tformat|<table|<row|<cell|-\<alpha\>*M<rsup|\<top\>>*A<rsup|<around*|(|2|)>>*M*Y<rsup|k-1>+\<tau\>*F>>|<row|<cell|0>>>>>.
  </equation>

  This is a square <math|<around*|(|27+4|)>\<times\><around*|(|27+4|)>>
  matrix and a <math|27+4> right hand side. <todo|The unknown
  <math|\<Lambda\>> can be discarded>

  The matrix <math|M> realizes the (local) operator
  <math|\<nabla\><rsub|h>:W<rsub|h><rsup|3>\<rightarrow\>\<Theta\><rsub|h><rsup|2>>.
  Recall from the scalar case that <math|M> uses geometric information of the
  cell to combine the partial derivatives of its argument. Given the
  definition of <math|\<nabla\><rsub|h>>, the full matrix
  <math|M\<in\>\<bbb-R\><rsup|36\<times\>27>> is block diagonal with each of
  its three diagonal blocks given by a copy of the scalar version.

  <subsection|Implementation>

  Assembly of the full system from all the pieces given in
  <eqref|eq:kirchhoff-local-timestep-system> requires some dof-juggling.
  Also, it is essential to precompute a valid sparsity pattern or assembly at
  each timestep will be too costly.

  <\itemize-dot>
    <item>The <strong|local> <math|M<rsup|\<top\>>*A<rsup|<around*|(|2|)>>*M>
    is dense, since <cpp|ufc::cell_integral::tabulate_tensor()> returns a
    dense matrix. However:

    <item>Is the number of zeroes for the global tensor for DKT elements
    correctly precomputed by <name|FEniCS>? How is it done? See
    <cpp|AssemblerBase::init_global_tensor()>

    <item>Using a <cpp|dolfin::BlockMatrix> to store the system matrix seemed
    convenient but no solver can use it since it doesn't implement
    <cpp|dolfin::GenericLinearOperator>. A quick fix was to write a
    <cpp|BlockMatrixAdapter> to conctruct a \Pflattened\Q version of the
    block matrix which kept track of the positions of the blocks so as to
    read from and write to them when updating was needed (in particular,
    updating of the isometry constraint). This doubles memory usage and
    completely broke parallel operation but \Pquickly\Q got the job
    done.<\footnote>
      It still took several days of coding... :(
    </footnote>
  </itemize-dot>

  <\bibliography|bib|tm-alpha|hermite.bib>
    <\bib-list|3>
      <bibitem*|Bar13><label|bib-bartels_approximation_2013>S.<nbsp>Bartels.<newblock>
      Approximation of Large Bending Isometries with Discrete Kirchhoff
      Triangles.<newblock> <with|font-shape|italic|SIAM Journal on Numerical
      Analysis>, 51(1):516\U525, jan 2013.<newblock>

      <bibitem*|Bar15><label|bib-bartels_numerical_2015>Sören
      Bartels.<newblock> <with|font-shape|italic|Numerical Methods for
      Nonlinear Partial Differential Equations>,
      <localize|volume><nbsp>47<localize| of
      ><with|font-shape|italic|Springer Series in Computational
      Mathematics>.<newblock> Springer International Publishing, Cham,
      2015.<newblock>

      <bibitem*|BS05><label|bib-brenner_c0_2005>Susanne<nbsp>C.<nbsp>Brenner<localize|
      and >Li-Yeng Sung.<newblock> <math|C<rsup|0>> Interior Penalty Methods
      for Fourth Order Elliptic Boundary Value Problems on Polygonal
      Domains.<newblock> <with|font-shape|italic|Journal of Scientific
      Computing>, 22-23(1-3):83\U118, jun 2005.<newblock>
    </bib-list>
  </bibliography>
</body>

<\initial>
  <\collection>
    <associate|font|stix>
    <associate|font-base-size|11>
    <associate|indent-indentation|1.5tab>
    <associate|info-flag|detailed>
    <associate|math-font|math-stix>
    <associate|page-medium|papyrus>
    <associate|preamble|false>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-10|<tuple|4.1|?>>
    <associate|auto-11|<tuple|7|?>>
    <associate|auto-2|<tuple|2|1>>
    <associate|auto-3|<tuple|2.1|1>>
    <associate|auto-4|<tuple|2.2|3>>
    <associate|auto-5|<tuple|2.2.1|4>>
    <associate|auto-6|<tuple|3|5>>
    <associate|auto-7|<tuple|1|5>>
    <associate|auto-8|<tuple|4|7>>
    <associate|auto-9|<tuple|1|7>>
    <associate|bib-bartels_approximation_2013|<tuple|Bar13|7>>
    <associate|bib-bartels_numerical_2015|<tuple|Bar15|7>>
    <associate|bib-brenner_c0_2005|<tuple|BS05|8>>
    <associate|eq:coeff-midpoint|<tuple|2|?>>
    <associate|eq:h2-flow-update-system|<tuple|4|6>>
    <associate|eq:kirchhoff-local-timestep-system|<tuple|7|?>>
    <associate|eq:local-dkt-gradient|<tuple|1|3>>
    <associate|eq:local-tensor-lki|<tuple|3|4>>
    <associate|eq:nodal-dirichlet-constraint|<tuple|6|?>>
    <associate|eq:nodal-isometry-constraint-components|<tuple|5|6>>
    <associate|footnote-1|<tuple|1|2>>
    <associate|footnote-2|<tuple|2|4>>
    <associate|footnote-3|<tuple|3|5>>
    <associate|footnote-4|<tuple|4|6>>
    <associate|footnote-5|<tuple|5|?>>
    <associate|footnote-6|<tuple|6|?>>
    <associate|footnote-7|<tuple|7|?>>
    <associate|footnr-1|<tuple|1|2>>
    <associate|footnr-2|<tuple|2|4>>
    <associate|footnr-3|<tuple|3|5>>
    <associate|footnr-4|<tuple|4|6>>
    <associate|footnr-5|<tuple|5|?>>
    <associate|footnr-6|<tuple|6|?>>
    <associate|footnr-7|<tuple|7|?>>
    <associate|rem:extending-ufl-discrete-gradient|<tuple|1|1>>
    <associate|sec:linear-kirchhoff|<tuple|3|5>>
    <associate|tab:matrix-bk|<tuple|1|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      bartels_numerical_2015

      bartels_numerical_2015

      bartels_numerical_2015

      brenner_c0_2005

      bartels_approximation_2013

      bartels_numerical_2015

      bartels_approximation_2013

      bartels_approximation_2013

      bartels_numerical_2015

      bartels_approximation_2013
    </associate>
    <\associate|figure>
      <tuple|normal|Solution of the linear Kirchhoff model for a clamped
      plate under a constant force.|<pageref|auto-7>>
    </associate>
    <\associate|table>
      <tuple|normal|One of the (local) products
      <with|mode|<quote|math>|B<rsub|i><rsup|k-1>*d<rsub|t><rsub|>Y<rsup|k><rsub|i>>
      (for subspace <with|mode|<quote|math>|W<rsup|<around*|(|i|)>><rsub|h>>)
      for a single element. Colours denote evaluation at each of the three
      nodes of a triangle. Empty entries are zeroes. Column indices
      correspond to local dof indices: <with|mode|<quote|math>|0,1,2,\<ldots\>>|<pageref|auto-9>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Discrete
      Kirchhoff Triangles> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Discrete
      gradient operators> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>

      <with|par-left|<quote|1tab>|2.1<space|2spc>Piecewise linear functions
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|1tab>|2.2<space|2spc>Discrete Kirchhoff Triangles
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|2tab>|2.2.1<space|2spc>Computation of the local
      tensor <with|mode|<quote|math>|M> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Linear
      Kirchhoff model> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|4<space|2spc>Non
      linear Kirchhoff model> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8><vspace|0.5fn>

      <with|par-left|<quote|1tab>|4.1<space|2spc>Implementation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>