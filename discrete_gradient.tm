<TeXmacs|1.99.4>

<style|generic>

<\body>
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
</body>

<initial|<\collection>
</collection>>