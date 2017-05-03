<TeXmacs|1.99.4>

<style|generic>

<\body>
  <section|Piecewise linear functions>

  Let <math|u> be a piecewise linear function over a mesh <math|\<cal-M\>>,
  that is <math|u\<in\>V=span<around*|{|\<phi\><rsub|i>|}>>, where <math|V>
  is a CG1 finite element space and <math|\<phi\><rsub|i>> are locally in
  <math|P<rsub|1><around*|(|T|)>> for some triangle <math|T\<in\>\<cal-M\>>.
  Then <math|u<rprime|'>> is piecewise constant, i.e.
  <math|u<rprime|'>\<in\>V<rprime|'>=span<around*|{|\<psi\><rsub|i>|}>>, a
  DG0 space over <math|\<cal-M\>>, where <math|\<psi\><rsub|i>> are constant
  over triangles. Denote by <math|\<b-u\>> and <math|\<b-u\><rprime|'>> their
  vectors of coefficients wrt. the standard bases of nodal functions over
  <math|V> and <math|V<rprime|'>> respectively. \ We want to compute a matrix
  <math|D<rsub|h>> for the gradient operator

  <\equation*>
    \<nabla\><rsub|h>:V\<rightarrow\>V<rprime|'>
  </equation*>

  such that <math|<around*|(|D<rsub|h> \<b-u\>|)><rsub|j>=\<b-u\><rprime|'><rsub|j>>.

  Since <math|u<rprime|'>> is not in the same space as <math|u>, one way of
  computing <math|\<b-u\><rprime|'>> in FEniCS is to `project(grad(u), V')`,
  but what we are doing here is to manually compute the transformation
  without solving a linear system. Note that if we were only interested in
  evaluation of <math|u<rprime|'>>, because
  <math|u=<big|sum><rsub|i>u<rsub|i>*\<phi\><rsub|i>> and by linearity
  <math|u<rprime|'>=<big|sum><rsub|i>u<rsub|i>*\<phi\><rsub|i><rprime|'>>. We
  could then use `evaluate_basis_derivatives()` to evaluate the
  <math|\<phi\><rsub|i><rprime|'>> and we would have finished.

  But again, we are interested in obtaining the new coefficients
  <math|\<b-u\><rsub|i><rprime|'>> such that

  <\equation*>
    u<rprime|'>=<big|sum><rsub|i>\<b-u\><rprime|'><rsub|i>*\<psi\><rsub|i>.
  </equation*>

  We assemble the transformation matrix cell by cell as follows. We need to
  assume that both <math|V> and <math|V<rprime|'>> are defined over the same
  mesh <math|\<cal-M\>>.

  <\itemize-dot>
    <item>First, for any given cell <math|T> compute the mapping
    <math|L<rsub|T>> of the dofs in <math|V<rsub|\|T>> to the dofs in
    <math|V<rprime|'><rsub|\|T>>.

    <item>Then apply this mapping to the components of
    <math|\<b-u\><rsub|\|T>> to obtain the coefficients
    <math|\<b-u\><rprime|'><rsub|\|T>>.
  </itemize-dot>

  To explain the idea, consider first a 1D mesh for simplicity and assume
  that we have an index <math|i\<in\><around|{|0,1,...,M|}>> over the cells
  of the mesh such that:

  <\itemize-dot>
    <item><math|\<b-u\><rsub|i>,\<b-u\><rsub|i+1>> are the coefficients for
    the dofs in <math|V> over cell <math|i>.

    <item><math|\<b-u\><rsub|i><rprime|'>> is the coefficient for the only
    dof in <math|V<rprime|'>> (the constant 1) over cell <math|i>.
  </itemize-dot>

  This assumption is actually a description of the mapping <math|L<rsub|T>>
  mentioned above for each interval. Then we only need to compute
  <math|\<phi\><rsub|i><rprime|'>> (not a typo, we are actually computing the
  derivative!) for each <math|i> to obtain the new vector of coefficients
  <math|\<b-u\><rsub|\|T><rprime|'>>:

  <\equation*>
    <around*|(|<tabular*|<tformat|<cwith|1|-1|1|1|cell-halign|c>|<cwith|1|-1|1|1|cell-lborder|0ln>|<cwith|1|-1|2|2|cell-halign|c>|<cwith|1|-1|3|3|cell-halign|c>|<cwith|1|-1|4|4|cell-halign|c>|<cwith|1|-1|5|5|cell-halign|c>|<cwith|1|-1|5|5|cell-rborder|0ln>|<table|<row|<cell|\<phi\><rsub|0><rprime|'>>|<cell|\<phi\><rsub|1><rprime|'>>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|\<phi\><rsub|1><rprime|'>>|<cell|\<phi\><rsub|2><rprime|'>>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|\<phi\><rsub|2><rprime|'>>|<cell|\<phi\><rsub|3><rprime|'>>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|\<phi\><rsub|3><rprime|'>>|<cell|\<phi\><rsub|4><rprime|'>>>>>>|)><around*|(|<tabular*|<tformat|<cwith|1|-1|1|1|cell-halign|c>|<cwith|1|-1|1|1|cell-lborder|0ln>|<cwith|1|-1|1|1|cell-rborder|0ln>|<table|<row|<cell|u<rsub|0>>>|<row|<cell|u<rsub|1>>>|<row|<cell|u<rsub|2>>>|<row|<cell|u<rsub|3>>>|<row|<cell|u<rsub|4>>>>>>|)>=<around*|(|<tabular*|<tformat|<cwith|1|-1|1|1|cell-halign|c>|<cwith|1|-1|1|1|cell-lborder|0ln>|<cwith|1|-1|1|1|cell-rborder|0ln>|<table|<row|<cell|u<rsub|0><rprime|'>>>|<row|<cell|u<rsub|1><rprime|'>>>|<row|<cell|u<rsub|2><rprime|'>>>|<row|<cell|u<rsub|3><rprime|'>>>>>>|)>.
  </equation*>

  Note that in this particular instance and with the assumptions made, we
  obtain a band matrix and computing <math|u<rprime|'>> is actually a matter
  of computing the dot product of a "subset" of the vector of constants
  <math|<around|(|\<phi\><rsub|0><rprime|'>,...,\<phi\><rsub|N<rsub|V<rprime|'>>><rprime|'>|)>>
  with <math|u>.

  However in other situations, e.g. in higher dimensions, the ordering of the
  nodes won't be as convenient. We proceed then more generally as follows:
  Let <math|i> run over all cell indices and let <math|\<iota\><rsub|i>> be
  the local-to-global mappings for <math|V> and
  <math|\<iota\><rsub|i><rprime|'>> for <math|V<rprime|'>>. Initialise the
  matrix <math|M> of the discrete gradient to zero. At each step of the
  construction of <math|M> we edit row <math|r=\<iota\><rsub|i><rprime|'><around|(|0|)>>,
  i.e. the row whose product by <math|<around|(|u|)><rsub|j>> will return the
  coefficient for the dof in <math|V<rprime|'>> corresponding to cell
  <math|i>.

  Now, for every vertex <math|v<rsub|j>,j\<in\><around|{|0,...,d|}>> in cell
  <math|i> compute the value of the derivative of the global dof
  <math|\<phi\><rsub|\<iota\><rsub|r><around|(|j|)>>> and set

  <\equation*>
    M<rsub|r*j>=\<phi\><rsub|\<iota\><rsub|r><around|(|j|)>>,j\<in\><around|{|0,...,d|}>.
  </equation*>

  <todo|...>

  <section|Discrete Kirchoff Triangles>

  As a first test for DKT elements, we find optimal configurations of the
  vertical displacements <math|u> for the energy

  <\equation*>
    I<rsub|K*i><around|(|w|)>=<frac|1|2>*<big|int><rsub|\<omega\>><around|\||D<rsup|2>*w|\|><rsup|2><math-up|d>x-<big|int><rsub|\<omega\>>f<rsub|3>*w<math-up|d>x
  </equation*>

  under homogeneous Dirichlet boundary conditions. This means solving the
  linear problem

  <\equation*>
    <around|(|D<rsup|2>*w,D<rsup|2>*v|)>=<around|(|f<rsub|3>,v|)>,
  </equation*>

  for all <math|v\<in\>H<rsup|1><rsub|0><around*|(|\<Omega\>|)>>. Following
  Bartels2016... we implement the non-conforming discretisation

  <\equation*>
    <around|(|\<nabla\>*\<nabla\><rsub|h>*w<rsub|h>,\<nabla\>*\<nabla\><rsub|h>*v<rsub|h>|)>=<around|(|f<rsub|3>,v<rsub|h>|)>,
  </equation*>

  for all <math|v<rsub|h>\<in\>W<rsub|h,D>\<assign\>>..., where the crux of
  the matter is the discrete gradient operator <math|\<nabla\><rsub|h>>.
  Recall that <math|\<nabla\><rsub|h>:W<rsub|h>\<rightarrow\>\<Theta\><rsub|h>>
  is defined as...

  <subsection|Interpolation into <math|P<rsub|2><rsup|2>>>

  The operator <math|\<nabla\><rsub|h>> can be easily implemented as an
  interpolation operator from <math|W<rsub|h>> into <math|\<Theta\><rsub|h>>.
  It is enough to construct a global (sparse) matrix by traversing all the
  cells in the mesh (actually we would rather only traverse the facets, as
  long as we can then assign the dofs properly) and multiply the vector of
  coefficients by it. <todo|See ... for the details>.

  But if we wanted to use this operator in UFL we would need to instruct
  FEniCS how to generate UFC code for any form using it and for their
  respective cell or facet integrals. However, it is not clear that this is
  possible at all: UFL forms are element-agnostic, meaning that they have no
  knowledge of the basis of shape functions. In particular, they are meant to
  act on <em|one> generic basis function at a time, so that it seems
  impossible to define a discrete operator which acts on all basis functions
  of one cell.

  <subsection|Computation of the local tensor>

  The alternative to extending UFL is assembling the system matrix ourselves.
  Here is one possible way of doing this which reuses what FFC already
  provides. Let <math|T\<in\>\<cal-T\><rsub|h>> be one cell in the
  triangulation. We want to assemble the <strong|local tensor> for the
  bilinear form over the DKT element, <math|P<rsup|red><rsub|3><around*|(|T|)>>:

  <\equation>
    <label|eq:local-tensor-lki><around|(|\<nabla\>\<nabla\><rsub|h>*\<phi\><rsub|j>,\<nabla\>\<nabla\><rsub|h>*\<phi\><rsub|i>|)>=<around|(|f<rsub|3>,\<phi\><rsub|i>|)>,<text|
    for all >\<phi\><rsub|i>,\<phi\><rsub|j>\<in\>P<rsup|red><rsub|3><around*|(|T|)>.
  </equation>

  Let <math|\<theta\><rsub|l>,\<theta\><rsub|k>\<in\>P<rsub|2><rsup|2><around*|(|T|)>>
  and let

  <\equation*>
    A<rsup|<around*|(|2|)>><rsub|l\<nocomma\>k>=<around*|(|\<nabla\>\<theta\><rsub|k>,\<nabla\>\<theta\><rsub|l>|)>,<application-space|1em>k,l\<in\><around*|[|12|]>
  </equation*>

  be the <strong|local tensor> for the <math|P<rsub|2><rsup|2>> element over
  <math|T>. Recall that <math|M<rsub|T>=<around*|(|m<rsub|p\<nocomma\>q>|)>\<in\>\<bbb-R\><rsup|12\<times\>9>>
  is the matrix for the local discrete gradient <todo|defined...>. Then, for
  every <math|\<phi\><rsub|i>,i\<in\><around*|[|9|]>>:

  <\equation*>
    \<nabla\><rsub|h> \<phi\><rsub|i>=M<rsub|T>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|\<vdots\>>>|<row|<cell|1>>|<row|<cell|\<vdots\>>>|<row|<cell|0>>>>>=<around*|(|M<rsub|T>|)><rsub|l\<nocomma\>i>\<theta\><rsub|l>=m<rsub|l\<nocomma\>i>*\<theta\><rsub|l>,
  </equation*>

  and analogously

  <\equation*>
    \<nabla\><rsub|h>\<phi\><rsub|j>=m<rsub|k\<nocomma\>j>*\<theta\><rsub|k>.
  </equation*>

  Substituting back in <eqref|eq:local-tensor-lki> we obtain

  <\equation*>
    A<rsub|i\<nocomma\>j>=<around|(|\<nabla\>\<nabla\><rsub|h>*\<phi\><rsub|j>,\<nabla\>\<nabla\><rsub|h>*\<phi\><rsub|i>|)>=m<rsub|k\<nocomma\>j>*m<rsub|l\<nocomma\>i>*<around*|(|\<nabla\>\<theta\><rsub|k>,\<nabla\>\<theta\><rsub|l>|)>=m<rsub|k\<nocomma\>j>*m<rsub|l\<nocomma\>i>*A<rsub|l\<nocomma\>k><rsup|<around*|(|2|)>>.
  </equation*>

  To summarize:

  <\enumerate>
    <item>Compile with FFC the form <verbatim|inner<around|(|\<nabla\>u,\<nabla\>v|)>\<ast\>dx>
    for <verbatim|VectorFunctionSpace\<less\>P2\<gtr\>>. This provides us
    with the local cell tensor <math|A<rsup|<around|(|2|)>>> in
    <verbatim|cell_integral::tabulate_tensor()>.

    <item>Compile with FFC the form <verbatim|inner<around|(|u,v|)>\<ast\>dx>
    for <verbatim|FunctionSpace\<less\>DKT\<gtr\>>. The dofmaps, etc. should
    be ok. This provides us with the <strong|wrong> local cell tensor. In
    order to compute the right one, denoted <math|A>:

    <\enumerate>
      <item>Compute the local matrix for <math|\<nabla\><rsub|h>> over
      simplex <math|t>: <math|M<rsub|t>=<around|(|m<rsub|p*q>|)>>.

      <item>Compute

      <\equation*>
        A<rsub|i\<nocomma\>j>=m<rsub|k*j>*m<rsub|l*i>*A<rsub|l*k><rsup|<around|(|2|)>>=<around|(|M<rsub|:i>|)><rsup|T>*A<rsup|<around|(|2|)>>*M<rsub|:j>
      </equation*>

      <item>Do this inside <verbatim|cell_integral::tabulate_tensor()> for
      the form.
    </enumerate>

    <item>Let dolfin assemble the full system matrix.
  </enumerate>
</body>

<initial|<\collection>
</collection>>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-2|<tuple|2|?>>
    <associate|auto-3|<tuple|2.1|?>>
    <associate|auto-4|<tuple|2.2|?>>
    <associate|eq:local-tensor-dkt|<tuple|2.1|?>>
    <associate|eq:local-tensor-lki|<tuple|1|?>>
  </collection>
</references>