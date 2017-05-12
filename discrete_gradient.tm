<TeXmacs|1.99.4>

<style|generic>

<\body>
  Taking derivatives of functions belonging to finite element spaces is a
  form of interpolation into other finite element spaces. This can be always
  achieved by orthogonal projection onto the target space but in some cases
  it can be computed exactly (and faster) with some linear operator, the
  simplest example being CG1 elements (piecewise linear functions), whose
  derivatives are CG0 (piecewise constant). After briefly reviewing this we
  focus in a slightly different situation. For Discrete Kirchhoff Triangles,
  <inactive|<cite-detail|bartels|Ÿ8.2>> defines a discrete gradient operator
  which does nodal interpolation at the vertices while ensuring a linear
  normal component along the edges of the triangles. The actual
  implementation is quite more involved than before since we cannot easily
  extend <name|UFL> with this operator and hence we need to manually assemble
  the system matrix for each problem.

  <section|Piecewise linear functions>

  Let <math|u> be a piecewise linear function over a mesh
  <math|\<cal-T\><rsub|h>>, that is <math|u\<in\>V=span<around*|{|\<phi\><rsub|i>|}>>,
  where <math|V> is a CG1 finite element space and <math|\<phi\><rsub|i>> are
  locally in <math|P<rsub|1><around*|(|T|)>> for some simplex
  <math|T\<in\>\<cal-T\><rsub|h>>. Then <math|u<rprime|'>> is piecewise
  constant, i.e. <math|u<rprime|'>\<in\>V<rprime|'>=span<around*|{|\<psi\><rsub|i>|}>>,
  a DG0 space over <math|\<cal-T\><rsub|h>>, where <math|\<psi\><rsub|i>> are
  constant over simplices. Denote by <math|\<b-u\>> and
  <math|\<b-u\><rprime|'>> their vectors of coefficients wrt. the standard
  bases of nodal functions over <math|V> and <math|V<rprime|'>> respectively.
  We want to compute a matrix <math|D<rsub|h>> for the gradient operator

  <\equation*>
    \<nabla\><rsub|h>:V\<rightarrow\>V<rprime|'>
  </equation*>

  such that <math|<around*|(|D<rsub|h> \<b-u\>|)><rsub|j>=\<b-u\><rprime|'><rsub|j>>.

  Since <math|u<rprime|'>> is not in the same space as <math|u>, one way of
  computing <math|\<b-u\><rprime|'>> in <name|FEniCS> is to
  <python|project(grad(u), V)>, but what we are doing here is to manually
  compute the transformation without solving a linear system. Note that if we
  were only interested in evaluation of <math|u<rprime|'>>, because
  <math|u=<big|sum><rsub|i>u<rsub|i>*\<phi\><rsub|i>> and by linearity
  <math|u<rprime|'>=<big|sum><rsub|i>u<rsub|i>*\<phi\><rsub|i><rprime|'>>. We
  could then use <python|evaluate_basis_derivatives()> to evaluate the
  <math|\<phi\><rsub|i><rprime|'>> and we would have finished.

  But again, we are interested in obtaining the new coefficients
  <math|\<b-u\><rsub|i><rprime|'>> such that

  <\equation*>
    u<rprime|'>=<big|sum><rsub|i>\<b-u\><rprime|'><rsub|i>*\<psi\><rsub|i>.
  </equation*>

  We assemble the transformation matrix cell by cell as follows, assuming
  that both <math|V> and <math|V<rprime|'>> are defined over the same mesh
  <math|\<cal-T\><rsub|h>>.

  <\itemize-dot>
    <item>First, for any given cell <math|T> compute the mapping
    <math|L<rsub|T>> of the dofs in <math|V<rsub|\|T>> to the dofs in
    <math|V<rprime|'><rsub|\|T>>.

    <item>Then apply this mapping to the components of
    <math|\<b-u\><rsub|\|T>> to obtain the coefficients
    <math|\<b-u\><rprime|'><rsub|\|T>>.
  </itemize-dot>

  To explain the idea, consider first a 1D mesh for simplicity and assume
  that we have an index <math|i\<in\><around*|[|N|]>> over the cells of the
  mesh such that:

  <\itemize-dot>
    <item><math|\<b-u\><rsub|i>,\<b-u\><rsub|i+1>> are the coefficients for
    the dofs in <math|V> over cell <math|i>.

    <item><math|\<b-u\><rsub|i><rprime|'>> is the coefficient for the only
    dof in <math|V<rprime|'>> (the constant 1) over cell <math|i>.
  </itemize-dot>

  This assumption is actually a description of the mapping <math|L<rsub|T>>
  mentioned above for each interval. Then we only need to compute
  <math|\<phi\><rsub|i><rprime|'>> (not a typo, we are actually computing the
  derivative) for each <math|i> to obtain the new vector of coefficients
  <math|\<b-u\><rsub|\|T><rprime|'>>:

  <\equation*>
    <around*|(|<tabular*|<tformat|<cwith|1|-1|1|1|cell-halign|c>|<cwith|1|-1|1|1|cell-lborder|0ln>|<cwith|1|-1|2|2|cell-halign|c>|<cwith|1|-1|3|3|cell-halign|c>|<cwith|1|-1|4|4|cell-halign|c>|<cwith|1|-1|5|5|cell-halign|c>|<cwith|1|-1|5|5|cell-rborder|0ln>|<table|<row|<cell|\<phi\><rsub|0><rprime|'>>|<cell|\<phi\><rsub|1><rprime|'>>|<cell|>|<cell|>|<cell|>>|<row|<cell|>|<cell|\<phi\><rsub|1><rprime|'>>|<cell|\<phi\><rsub|2><rprime|'>>|<cell|>|<cell|>>|<row|<cell|>|<cell|>|<cell|\<phi\><rsub|2><rprime|'>>|<cell|\<phi\><rsub|3><rprime|'>>|<cell|>>|<row|<cell|>|<cell|>|<cell|>|<cell|\<phi\><rsub|3><rprime|'>>|<cell|\<phi\><rsub|4><rprime|'>>>>>>|)><around*|(|<tabular*|<tformat|<cwith|1|-1|1|1|cell-halign|c>|<cwith|1|-1|1|1|cell-lborder|0ln>|<cwith|1|-1|1|1|cell-rborder|0ln>|<table|<row|<cell|u<rsub|0>>>|<row|<cell|u<rsub|1>>>|<row|<cell|u<rsub|2>>>|<row|<cell|u<rsub|3>>>|<row|<cell|u<rsub|4>>>>>>|)>=<around*|(|<tabular*|<tformat|<cwith|1|-1|1|1|cell-halign|c>|<cwith|1|-1|1|1|cell-lborder|0ln>|<cwith|1|-1|1|1|cell-rborder|0ln>|<table|<row|<cell|u<rsub|0><rprime|'>>>|<row|<cell|u<rsub|1><rprime|'>>>|<row|<cell|u<rsub|2><rprime|'>>>|<row|<cell|u<rsub|3><rprime|'>>>>>>|)>.
  </equation*>

  Note that in this particular instance and with the assumptions made, we
  obtain a band matrix and computing <math|u<rprime|'>> is actually a matter
  of computing the dot product of a \Psubset\Q of the vector of constants
  <math|<around|(|\<phi\><rsub|0><rprime|'>,...,\<phi\><rsub|N<rsub|V<rprime|'>>><rprime|'>|)>>
  with <math|u>.

  However in other situations, e.g. in higher dimensions, the ordering of the
  nodes will not be as convenient. We proceed then more generally as follows.
  Let <math|i> run over all cell indices and let <math|\<iota\><rsub|i>> be
  the local-to-global mappings for <math|V> and
  <math|\<iota\><rsub|i><rprime|'>> for <math|V<rprime|'>>. Initialise the
  matrix <math|M> of the discrete gradient to zero. At each step of the
  construction of <math|M> we edit row <math|r=\<iota\><rsub|i><rprime|'><around|(|0|)>>,
  i.e. the row whose product by <math|<around|(|u|)><rsub|j>> will return the
  coefficient for the dof in <math|V<rprime|'>> corresponding to cell
  <math|i>.

  Now, for every vertex <math|v<rsub|j>,j\<in\><around*|[|d|]>> in cell
  <math|i> compute the value of the derivative of the global dof
  <math|\<phi\><rsub|\<iota\><rsub|r><around|(|j|)>>> and set

  <\equation*>
    M<rsub|r*j>=\<phi\><rsub|\<iota\><rsub|r><around|(|j|)>>,j\<in\><around|[|d|]>.
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
  <todo|Bartels2016...> we implement the non-conforming discretisation

  <\equation*>
    <around|(|\<nabla\>*\<nabla\><rsub|h>*w<rsub|h>,\<nabla\>*\<nabla\><rsub|h>*v<rsub|h>|)>=<around|(|f<rsub|3>,v<rsub|h>|)>,
  </equation*>

  for all <math|v<rsub|h>\<in\>W<rsub|h,D>\<assign\>>..., where the crux of
  the matter is the discrete gradient operator <math|\<nabla\><rsub|h>>.
  Recall that <math|\<nabla\><rsub|h>:W<rsub|h>\<rightarrow\>\<Theta\><rsub|h>>
  is <todo|defined as...>

  The operator <math|\<nabla\><rsub|h>> can be seen as an interpolation
  operator from <math|W<rsub|h>> into <math|\<Theta\><rsub|h>>.

  <subsection|Interpolation into <math|P<rsub|2><rsup|2>>>

  Fix some cell <math|T\<in\>\<cal-T\><rsub|h>> in the triangulation.
  Assuming the dofs in <math|W<rsub|h>> local to <math|T> are ordered as
  <math|\<b-w\><rsub|T>=<around*|(|w<rsub|z<rsub|1>>,\<partial\><rsub|1>
  w<rsub|z<rsub|1>>,\<partial\><rsub|2> w<rsub|z<rsub|1>>,w<rsub|z<rsub|2>>,\<partial\><rsub|1>
  w<rsub|z<rsub|2>>,\<partial\><rsub|2> w<rsub|z<rsub|2>>,w<rsub|z<rsub|3>>,\<partial\><rsub|1>
  w<rsub|z<rsub|3>>,\<partial\><rsub|2> w<rsub|z<rsub|3>>|)>> and those in
  <math|\<Theta\><rsub|h>> as <math|\<b-theta\><rsub|T>=<around*|(|\<theta\><rsub|z<rsub|1>>,\<theta\><rsub|z<rsub|2>>,\<theta\><rsub|z<rsub|3>>,\<theta\><rsub|S<rsub|1>>,\<theta\><rsub|S<rsub|2>>,\<theta\><rsub|S<rsub|3>>|)>>,
  the action of <math|\<nabla\><rsub|h>> on <math|T> is given by <todo|See
  <inactive|<cite-detail|bartels|Ÿ8.2.2>> for the details>

  <\equation*>
    \<b-theta\><rsub|T>=M<rsub|T>*\<b-w\><rsub|T>
  </equation*>

  where

  <\equation>
    <label|eq:local-dkt-gradient>M<rsub|T>=<matrix|<tformat|<table|<row|<cell|0>|<cell|I<rsub|2>>|<cell|0>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|I<rsub|2>>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|0>|<cell|I<rsub|2>>>|<row|<cell|0>|<cell|0>|<cell|<wide|t|~><rsub|S<rsub|1>>>|<cell|<wide|T|~><rsub|S<rsub|1>>>|<cell|-<wide|t|~><rsub|S<rsub|1>>>|<cell|<wide|T|~><rsub|S<rsub|1>>>>|<row|<cell|<wide|t|~><rsub|S<rsub|2>>>|<cell|<wide|T|~><rsub|S<rsub|2>>>|<cell|0>|<cell|0>|<cell|-<wide|t|~><rsub|S<rsub|2>>>|<cell|<wide|T|~><rsub|S<rsub|2>>>>|<row|<cell|<wide|t|~><rsub|S<rsub|3>>>|<cell|<wide|T|~><rsub|S<rsub|3>>>|<cell|-<wide|t|~><rsub|S<rsub|3>>>|<cell|<wide|T|~><rsub|S<rsub|3>>>|<cell|0>|<cell|0>>>>>\<in\>\<bbb-R\><rsup|12\<times\>9>,<space|1em><wide|T|~><rsub|S<rsub|l>>=<tfrac|-3|4>*t<rsub|S<rsub|l>>*t<rsub|S<rsub|l>><rsup|\<top\>>,<space|1em><wide|t|~><rsub|S<rsub|l>>=<tfrac|-3|2*<around*|\||S<rsub|l>|\|>>*t<rsub|S<rsub|l>>
  </equation>

  and <math|t<rsub|S<rsub|l>>> is the vector tangent to face
  <math|S<rsub|l>>, which is the one containing the two vertices other than
  <math|z<rsub|l>> (i.e. opposite to <math|z<rsub|l>>).

  In order to implement this globally, its enough to construct once a global
  (sparse) matrix by traversing all the cells in the mesh and then multiply
  the vector of coefficients of any <math|w<rsub|h>\<in\>W<rsub|h>> by
  it.<\footnote>
    Actually this would perform many computations twice so we would rather
    only traverse the facets, carefully keeping track of the dofs involved.
  </footnote>

  However, if we want to use this operator in <name|UFL> we need to instruct
  <name|FEniCS> how to generate <name|UFC> code for any form where it
  appears, as well as code for the respective cell or facet integrals, and it
  is not clear that this is possible at all. Indeed, <name|UFL> forms are
  \Pelement-agnostic\Q, meaning that they have no knowledge of the basis of
  shape functions. In particular they are meant to act on <em|one> generic
  basis function at a time, so that it seems impossible to define a discrete
  operator which acts on all basis functions of one cell.

  <subsection|Computation of the local tensor>

  The alternative to extending <name|UFL> is assembling the system matrix
  ourselves for each problem where we need <math|\<nabla\><rsub|h>>. Here is
  one way of doing this which reuses what <name|FFC> already provides. Let
  <math|T\<in\>\<cal-T\><rsub|h>> be one cell in the triangulation. We first
  want to assemble the <strong|local tensor> for the bilinear form over the
  DKT element, <math|P<rsup|red><rsub|3><around*|(|T|)>>:

  <\equation>
    <label|eq:local-tensor-lki><around|(|\<nabla\>\<nabla\><rsub|h>*\<phi\><rsub|j>,\<nabla\>\<nabla\><rsub|h>*\<phi\><rsub|i>|)>=<around|(|f<rsub|3>,\<phi\><rsub|i>|)>,<text|
    for all shape functions >\<phi\><rsub|i>,\<phi\><rsub|j>\<in\>P<rsup|red><rsub|3><around*|(|T|)>.
  </equation>

  Let <math|\<theta\><rsub|l>,\<theta\><rsub|k>\<in\>P<rsub|2><rsup|2><around*|(|T|)>>
  be two shape functions and let

  <\equation*>
    A<rsup|<around*|(|2|)>><rsub|l\<nocomma\>k>=<around*|(|\<nabla\>\<theta\><rsub|k>,\<nabla\>\<theta\><rsub|l>|)>,<application-space|1em>k,l\<in\><around*|[|12|]>
  </equation*>

  be the <strong|local tensor> for the <math|P<rsub|2><rsup|2><around*|(|T|)>>
  element. Let <math|M<rsub|T>=<around*|(|m<rsub|p\<nocomma\>q>|)>\<in\>\<bbb-R\><rsup|12\<times\>9>>
  be the matrix for the local discrete gradient <eqref|eq:local-dkt-gradient>
  and take one shape function <math|\<phi\><rsub|i>\<in\>P<rsup|3><rsub|red><around*|(|T|)>,i\<in\><around*|[|9|]>>.
  With some abuse of notation where we confuse functions with their
  coefficient vectors in their respective finite elements we can write:

  <\equation*>
    \<nabla\><rsub|h> \<phi\><rsub|i>=M<rsub|T>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|\<vdots\>>>|<row|<cell|1>>|<row|<cell|\<vdots\>>>|<row|<cell|0>>>>>=m<rsub|l\<nocomma\>i>*\<theta\><rsub|l>,
  </equation*>

  and analogously <math|\<nabla\><rsub|h>\<phi\><rsub|j>=m<rsub|k\<nocomma\>j>*\<theta\><rsub|k>>.
  Substituting back in <eqref|eq:local-tensor-lki> we obtain

  <\equation*>
    A<rsub|i\<nocomma\>j>=<around|(|\<nabla\>\<nabla\><rsub|h>*\<phi\><rsub|j>,\<nabla\>\<nabla\><rsub|h>*\<phi\><rsub|i>|)>=m<rsub|k\<nocomma\>j>*m<rsub|l\<nocomma\>i>*<around*|(|\<nabla\>\<theta\><rsub|k>,\<nabla\>\<theta\><rsub|l>|)>=m<rsub|k\<nocomma\>j>*m<rsub|l\<nocomma\>i>*A<rsub|l\<nocomma\>k><rsup|<around*|(|2|)>>.
  </equation*>

  Now, the steps required are:

  <\enumerate>
    <item>Compile with <name|FFC> the form
    <verbatim|inner<around|(|\<nabla\>u,\<nabla\>v|)>\<ast\>dx> for a
    <verbatim|VectorFunctionSpace\<less\>P2\<gtr\>>. This provides us with
    the local cell tensor <math|A<rsup|<around|(|2|)>>> in
    <verbatim|cell_integral::tabulate_tensor()>.

    <item>Compile with <name|FFC> the form
    <verbatim|inner<around|(|u,v|)>\<ast\>dx> for
    <verbatim|FunctionSpace\<less\>DKT\<gtr\>>. The dofmaps, etc. should be
    ok. This provides us with the <strong|wrong> local cell tensor. In order
    to compute the right one, denoted <math|A>:

    <\enumerate>
      <item>Compute the local matrix for <math|\<nabla\><rsub|h>> over
      simplex <math|T>: <math|M<rsub|T>=<around|(|m<rsub|p*q>|)>>.

      <item>Compute

      <\equation*>
        A<rsub|i\<nocomma\>j>=m<rsub|k*j>*m<rsub|l*i>*A<rsub|l*k><rsup|<around|(|2|)>>=<around|(|M<rsub|:i>|)><rsup|\<top\>>*A<rsup|<around|(|2|)>>*M<rsub|:j>
      </equation*>

      <item>Do this inside <verbatim|cell_integral::tabulate_tensor()> for
      the form.
    </enumerate>

    <item>Compile and link the new form and make it available to
    <name|dolfin>.

    <item>Let <name|dolfin> assemble the full system matrix.
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
    <associate|eq:local-dkt-gradient|<tuple|1|?>>
    <associate|eq:local-tensor-lki|<tuple|2|?>>
    <associate|footnote-1|<tuple|1|?>>
    <associate|footnr-1|<tuple|1|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Piecewise
      linear functions> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Discrete
      Kirchoff Triangles> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>

      <with|par-left|<quote|1tab>|2.1<space|2spc>Interpolation into
      <with|mode|<quote|math>|P<rsub|2><rsup|2>>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|1tab>|2.2<space|2spc>Computation of the local
      tensor <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>
    </associate>
  </collection>
</auxiliary>