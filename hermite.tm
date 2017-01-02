<TeXmacs|1.99.4>

<style|<tuple|generic|british>>

<\body>
  <\hide-preamble>
    <assign|vector|<macro|x|<shrink-inline|<left|(><tformat|<cwith|1|-1|1|-1|cell-halign|c><cwith|1|-1|1|1|cell-lsep|0><cwith|1|-1|1|1|cell-rsep|0><cwith|1|-1|1|1|cell-bsep|0.1em><cwith|1|-1|1|1|cell-tsep|0.1em>|<arg|x>><right|)>>>>

    <assign|quotation|<\macro|body>
      <\padded>
        <\indent-both|<value|quote-left-indentation>|<value|quote-right-indentation>>
          <surround|<yes-indent>||<em|<arg|body>>>
        </indent-both>
      </padded>
    </macro>>
  </hide-preamble>

  <doc-data|<doc-title|Hermite elements: introduction<line-break> and
  implementation in FEniCS>|<doc-author|<author-data|<author-name|Miguel de
  Benito>|<\author-affiliation>
    Universität Augsburg
  </author-affiliation>>>|<doc-date|December 2016>>

  <subsection|Setting>

  We fix a polygonal domain <math|\<Omega\>\<subset\>\<bbb-R\><rsup|2>> (the
  <dfn|physical domain>) and a triangulation<\footnote>
    A <em|triangulation> is...
  </footnote> <math|\<cal-T\>=<around*|{|K<rsub|1>,\<ldots\>,K<rsub|N>|}>>.
  Let <math|V> be a (not necessarily proper) subspace of
  <math|H<rsup|2><around*|(|\<Omega\>|)>> and
  <math|L:V\<rightarrow\>\<bbb-R\>> and <math|a:V\<times\>V\<rightarrow\>\<bbb-R\>>
  respectively linear and bilinear forms. We consider the linear variational
  problem:

  <\em>
    Find <math|u\<in\>V> such that

    <\equation>
      <label|eq:linear-problem>a<around*|(|v,u|)>=L<around*|(|v|)><text|<space|1em>for
      every >v\<in\>V.
    </equation>
  </em>

  In broad terms, the standard construction of a discrete polynomial space
  <math|V<rsub|h>> for this problem involves:\ 

  <\itemize-dot>
    <item>The triangulation <math|\<cal-T\>>.

    <item>A reference simplex <math|<wide|K|^>\<subset\>\<bbb-R\><rsup|2>>
    with an associated set of linear forms (<dfn|degrees of freedom>) over
    <math|P<rsub|p><around*|(|<wide|K|^>|)>>, the space of polynomials of
    degree up to <math|p>, defining a (dual) local polynomial basis
    <math|<wide|\<cal-B\>|^>> of <math|P<rsub|p><around*|(|<wide|K|^>|)>>.

    <item>A set of invertible affine maps
    <math|<around*|{|\<b-x\><rsup|K>|}><rsub|K\<in\>\<cal-T\>>> between
    <math|<wide|K|^>> and each triangle in <math|\<cal-T\>>.

    <item>A transformation <math|\<cal-F\>:<wide|\<cal-B\>|^><around*|(|P<rsub|p><around*|(|<wide|K|^>|)>|)>\<rightarrow\>\<cal-B\><around*|(|P<rsub|p><around*|(|K|)>|)>>
    preserving certain properties of <math|<wide|\<cal-B\>|^>> and
    <math|<wide|\<cal-B\>|^><rprime|'>>(e.g. duality).

    <item>A local-to-global dof mapping <math|\<iota\><rsup|K>>, assigning
    local dofs in each <math|K\<in\>\<cal-T\>> to global dofs in
    <math|V<rsub|h>>.

    <item><todo|...>
  </itemize-dot>

  <subsection|Conforming and non-conforming discretizations>

  If <math|V<rsub|h>> is a subspace of <math|H<rsup|2><around*|(|\<Omega\>|)>>,
  we say that the discretization is <dfn|conforming>. This is not the case we
  will be handling. For Hermite elements in 2D:

  <\quotation>
    The cubic Hermite elements have a continuous normal derivative but not
    full <math|C<rsup|1>> continuity. In particular, the normal derivatives
    may not match at the boundary of two elements, away from the vertices. If
    you want full <math|C<rsup|1>> continuity you will have to use the
    Argyris element or Hsieh-Clough-Tucker or something. I recommend the
    discussion in chapter 6 of Ciarlet's finite element book.<htab|5mm>
    (<hlink|CompSci SE|http://scicomp.stackexchange.com/questions/2012/construction-of-c1-h2-conforming-finite-element-basis-for-triangular-or-te>)
  </quotation>

  Convergence results for non-conforming discretizations...

  <subsection|Cubic Hermite elements in 2D>

  <\notation*>
    In the following, we use the indices <math|\<alpha\>\<in\><around*|{|1,4,7|}>>
    and <math|i\<in\><around*|{|1,2|}>>. We also set
    <math|\<alpha\><rsub|i>=\<alpha\>+i>. <todo|ugly!>
  </notation*>

  Let <math|<wide|K|^>> be a reference simplex in 2D with vertices
  <math|<around*|{|<wide|v|^><rsub|1>,<wide|v|^><rsub|2>,<wide|v|^><rsub|3>|}>>
  and <math|\<b-x\>:<wide|K|^>\<rightarrow\>K> be a non-degenerate affine
  transformation into one of the cells of the physical domain
  <math|\<Omega\>>. Denote the inverse mapping by
  <math|<wide|\<b-x\>|^>:K\<rightarrow\><wide|K|^>> (Figure
  <reference|fig:reference-triangle>).

  Let <math|P<rsub|3>=P<rsub|3><around*|(|<wide|K|^>|)>> be the space of
  polynomials of order up to 3 defined over <math|<wide|K|^>>. Define 10
  linear forms (the <dfn|degrees of freedom> or <dfn|dofs>) in the dual
  <math|P<rsub|3><rprime|'>\<assign\>P<rsub|3><around*|(|<wide|K|^>|)><rprime|'>>
  as follows: for each vertex <math|v<rsub|\<alpha\>>> in <math|<wide|K|^>>
  add one function evaluation <math|\<psi\><rsub|\<alpha\>><around*|(|f|)>\<assign\>f<around*|(|<wide|v|^><rsub|<around*|(|\<alpha\>-1|)>/3>|)>>
  and evaluation of both partial derivatives
  <math|\<psi\><rsub|\<alpha\><rsub|i>><around*|(|f|)>=f<rsub|,i><around*|(|<wide|v|^><rsub|<around*|(|\<alpha\>-1|)>/3>|)>>
  for <math|i=1,2>. Finally at the barycenter <math|<wide|v|^><rsub|0>> of
  <math|<wide|K|^>> add a new evaluation <math|\<psi\><rsub|0><around*|(|f|)>=f<around*|(|<wide|v|^><rsub|0>|)>>.

  <big-figure|<with|gr-mode|<tuple|edit|math-at>|gr-frame|<tuple|scale|1cm|<tuple|0.270004gw|0.180009gh>>|gr-geometry|<tuple|geometry|0.280007par|0.200003par|center>|gr-grid|<tuple|cartesian|<point|0|0>|1>|gr-grid-old|<tuple|cartesian|<point|0|0>|1>|gr-edit-grid-aspect|<tuple|<tuple|axes|none>|<tuple|1|none>|<tuple|10|none>>|gr-edit-grid|<tuple|cartesian|<point|0|0>|1>|gr-edit-grid-old|<tuple|cartesian|<point|0|0>|1>|<graphics||<line|<point|0|0>|<point|2.0|0.0>|<point|0.0|2.0>|<point|0.0|0.0>>>>|<todo|<label|fig:reference-triangle>The
  reference triangle and the mapping <math|\<b-x\>>...>>

  The forms <math|\<cal-B\><rprime|'>\<assign\><around*|{|<wide|\<psi\>|^><rsub|\<alpha\>>,<wide|\<psi\>|^><rsub|\<alpha\><rsub|1>>,<wide|\<psi\>|^><rsub|\<alpha\><rsub|2>>:\<alpha\>=1,4,7|}>\<cup\><around*|{|<wide|\<psi\>|^><rsub|0>|}>>
  are a basis of <math|P<rsub|3><rprime|'>> with dual
  <math|\<cal-B\>=<around*|{|<wide|\<varphi\>|^><rsub|\<alpha\>>,<wide|\<varphi\>|^><rsub|\<alpha\><rsub|1>>,<wide|\<varphi\>|^><rsub|\<alpha\><rsub|2>>:\<alpha\>=1,4,7|}>\<cup\><around*|{|<wide|\<varphi\>|^><rsub|0>|}>>.
  This means that the latter consists of 10 polynomials
  <math|<wide|\<varphi\>|^><rsub|k>\<in\>P<rsub|3><around*|(|<wide|K|^>|)>,k\<in\><around*|{|0,\<ldots\>,9|}>>
  such that

  <\equation>
    <label|eq:delta-property><wide|\<psi\>|^><rsub|j><around*|(|<wide|\<varphi\>|^><rsub|k>|)>=\<delta\><rsub|j\<nocomma\>k>.
  </equation>

  In particular <math|<wide|\<psi\>|^><rsub|1><around*|(|<wide|\<varphi\>|^><rsub|1>|)>=<wide|\<varphi\>|^><rsub|1><around*|(|<wide|v|^><rsub|1>|)>=1>,
  <math|<wide|\<psi\>|^><rsub|2><around*|(|<wide|\<varphi\>|^><rsub|2>|)>=\<partial\><rsub|1><wide|\<varphi\>|^><rsub|2><around*|(|<wide|v|^><rsub|1>|)>=1>,
  but <math|<wide|\<psi\>|^><rsub|1><around*|(|<wide|\<varphi\>|^><rsub|2>|)>=<wide|\<varphi\>|^><rsub|2><around*|(|<wide|v|^><rsub|1>|)>=0>
  and so on.

  This <dfn|reference polynomial basis> is then mapped onto a <dfn|local
  basis over the physical cell> <math|K> by means of a mapping

  <\equation*>
    <wide|\<varphi\>|^>\<mapsto\>\<varphi\>=\<cal-F\><around*|(|<wide|\<varphi\>|^>|)>.
  </equation*>

  For standard Lagrange elements (i.e. when <math|\<cal-B\><rprime|'>>
  consists only of point evaluation forms), the immediate choice
  <math|\<cal-F\><around*|(|<wide|\<varphi\>|^>|)>\<assign\><wide|\<varphi\>|^>\<circ\><wide|\<b-x\>|^>>
  maps the basis <math|\<cal-B\>> into a basis of
  <math|P<rsub|3><around*|(|K|)>> which fulfills the delta property
  <eqref|eq:delta-property>. However, in our case, in order for
  <eqref|eq:delta-property> to hold in <math|K> for the images of the
  \Phermite\Q pairs of dofs <math|<wide|\<Phi\>|^><rsub|\<alpha\>>\<assign\><around*|(|<wide|\<varphi\>|^><rsub|\<alpha\><rsub|1>>,<wide|\<varphi\>|^><rsub|\<alpha\><rsub|2>>|)>>
  we require a different transformation, namely:

  <\equation>
    <label|eq:hermite-transform>\<Phi\><rsub|\<alpha\>>=\<cal-F\><around*|(|<wide|\<Phi\>|^><rsub|\<alpha\>>|)>\<assign\><around*|(|D
    \<b-x\>*<wide|\<Phi\>|^><rsub|\<alpha\>>|)>\<circ\><wide|\<b-x\>|^>,
  </equation>

  or, component-wise:

  <\equation*>
    \<varphi\><rsub|\<alpha\><rsub|i>><around*|(|x|)>=\<nabla\>\<b-x\><rsub|i><around*|(|<wide|\<b-x\>|^><around*|(|x|)>|)>\<cdot\><wide|\<Phi\>|^><rsub|\<alpha\>><around*|(|<wide|\<b-x\>|^><around*|(|x|)>|)>=<around*|(|\<b-x\><rsub|i,j>*<wide|\<varphi\>|^><rsub|\<alpha\><rsub|j>>|)><around*|(|<wide|\<b-x\>|^><around*|(|x|)>|)>.
  </equation*>

  We can check that <eqref|eq:delta-property> holds: for the point
  evaluations this is clear and for the others:

  <math|\<psi\><rsub|2><around*|(|\<varphi\><rsub|2>|)>=\<partial\><rsub|1>\<varphi\><rsub|2><around*|(|v<rsub|1>|)>=\<ldots\>>,

  Note that at vertex <math|v<rsub|\<alpha\>>> of <math|K>, each new basis
  function is a linear combination of all the old ones at
  <math|<wide|v|^><rsub|\<alpha\>>> (the coefficients for the point
  evaluations happen to be zero), so we can write
  <math|<wide|\<Phi\>|^><rsub|\<alpha\>>\<assign\><around*|(|<wide|\<varphi\>|^><rsub|\<alpha\>>,<wide|\<varphi\>|^><rsub|\<alpha\><rsub|1>>,<wide|\<varphi\>|^><rsub|\<alpha\><rsub|2>>|)>>
  and define

  <\equation*>
    <vector|<tformat|<cwith|1|3|1|1|cell-lsep|0>|<cwith|1|3|1|1|cell-rsep|0>|<cwith|1|3|1|1|cell-bsep|0>|<cwith|1|3|1|1|cell-tsep|0>|<table|<row|<cell|\<varphi\><rsub|><rsub|\<alpha\>>>>|<row|<cell|\<varphi\><rsub|><rsub|\<alpha\><rsub|1>>>>|<row|<cell|\<varphi\><rsub|><rsub|\<alpha\><rsub|2>>>>>>><around*|(|x|)>=\<cal-F\><around*|(|<wide|\<Phi\>|^><rsub|\<alpha\>>|)>\<assign\><wide*|<matrix|<tformat|<table|<row|<cell|1>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|\<b-x\><rsub|1,1>>|<cell|\<b-x\><rsub|1,2>>>|<row|<cell|0>|<cell|\<b-x\><rsub|2,1>>|<cell|\<b-x\><rsub|2,2>>>>>>|\<wide-squnderbrace\>><rsub|H><around*|(|<wide|\<b-x\>|^><around*|(|x|)>|)>*<vector|<tformat|<cwith|1|3|1|1|cell-lsep|0>|<cwith|1|3|1|1|cell-rsep|0>|<cwith|1|3|1|1|cell-bsep|0>|<cwith|1|3|1|1|cell-tsep|0>|<table|<row|<cell|<wide|\<varphi\>|^><rsub|><rsub|\<alpha\>>>>|<row|<cell|<wide|\<varphi\>|^><rsub|><rsub|\<alpha\><rsub|1>>>>|<row|<cell|<wide|\<varphi\>|^><rsub|><rsub|\<alpha\><rsub|2>>>>>>><rsub|\|<wide|\<b-x\>|^><around*|(|x|)>>
  </equation*>

  Using this <math|\<cal-F\>> means that we have to take the factor <math|D
  \<b-x\>> into account when differentiating the basis functions. Tedious
  application of the chain rule yields that each partial derivative
  <math|\<varphi\><rsub|\<alpha\><rsub|i>,j>> is a linear combination of the
  derivatives of both reference basis functions
  <math|\<nabla\><wide|\<varphi\>|^><rsub|\<alpha\><rsub|k><rsub|>>> with
  weights given by <math|\<nabla\> \<b-x\><rsub|i>>:

  <\equation>
    <label|eq:hermite-first-derivatives>\<varphi\><rsub|\<alpha\><rsub|i>,j><around*|(|x|)>=<around*|(|\<b-x\><rsub|i,k>*<wide|\<varphi\>|^><rsub|\<alpha\><rsub|k>,l>|)><rsub|\|<wide|\<b-x\>|^><around*|(|x|)>>*<wide|\<b-x\>|^><rsub|l,j><around*|(|x|)>
  </equation>

  or, in matrix form for the three dofs over <math|v<rsub|\<alpha\>>> as
  above:

  <\equation*>
    D \<Phi\><rsub|\<alpha\>><around*|(|x|)>=<around*|(|H*D<wide|\<Phi\>|^><rsub|\<alpha\>>|)><around*|(|<wide|\<b-x\>|^><around*|(|x|)>|)>*D<wide|\<b-x\>|^><around*|(|x|)>.
  </equation*>

  The last term <math|D<wide|\<b-x\>|^>> comes from the affine change of
  coordinates and is also present in the Lagrange case.

  For the second derivatives we have, thanks to <math|D<rsup|2> \<b-x\>=0>:

  <\equation>
    <label|eq:hermite-second-derivatives>\<varphi\><rsub|\<alpha\><rsub|i>,j\<nocomma\>k>=\<b-x\><rsub|i,p>*<wide|\<varphi\>|^><rsub|\<alpha\><rsub|p>,q\<nocomma\>r>*<wide|\<b-x\>|^><rsub|r,k>*<wide|\<b-x\>|^><rsub|q,j><text|,
    i.e. >D<rsup|2> \<varphi\><rsub|\<alpha\><rsub|i>>=D<rsup|\<top\>><wide|\<b-x\>|^>*<around*|(|\<b-x\><rsub|i,1>*D<rsup|2><wide|\<varphi\>|^><rsub|\<alpha\><rsub|1>>+\<b-x\><rsub|i,2>*D<rsup|2><wide|\<varphi\>|^><rsub|\<alpha\><rsub|2>>|)>*D
    <wide|\<b-x\>|^>.
  </equation>

  <todo|Write in terms of <math|H>>

  Here as before, we obtain a linear combination of the second derivatives of
  both reference basis functions <math|D<rsup|2>
  <wide|\<varphi\>|^><rsub|\<alpha\><rsub|j>>> with weights given by
  <math|\<nabla\> \<b-x\><rsub|i>>. The terms
  <math|D<rsup|\<top\>><wide|\<b-x\>|^>> and <math|D<wide|\<b-x\>|^>> are
  again the same as if we only had used the inverse mapping
  <math|<wide|\<b-x\>|^>> to change variables.

  <subsection|Assembly of the stiffness matrix>

  In order to discretize <eqref|eq:linear-problem> we require the action of
  <math|a> on the discrete space <math|V<rsub|h>=span<around*|{|\<varphi\><rsub|1>,\<ldots\>,\<varphi\><rsub|<todo|M*N>>|}>>:

  <\equation*>
    a:V<rsup|1><rsub|h>\<times\>V<rsup|2><rsub|h>\<rightarrow\>\<bbb-R\>,
  </equation*>

  where we used superindices to differentiate two copies of <math|V<rsub|h>>.
  This can be split into the contribution of each element, meaning that we
  compute the values

  <\equation*>
    A<rsub|i>=a<around*|(|\<varphi\><rsub|i<rsub|1>><rsup|1>,\<varphi\><rsup|2><rsub|i<rsub|2>>|)>=<big|sum><rsub|K\<in\>\<cal-T\>>a<rsup|K><around*|(|\<varphi\><rsub|i<rsub|1>><rsup|1>,\<varphi\><rsup|2><rsub|i<rsub|2>>|)>
  </equation*>

  with <math|i> a multindex with <math|0\<leqslant\><around*|\||i|\|>\<less\>2*<around*|\||\<cal-T\>|\|>>
  and <math|a<rsup|K>> the contribution from element <math|K>. Because of the
  compact supports of the basis functions, most terms in the sum will be
  zero. <todo|This results in a global element matrix <math|A> and a system
  <math|A*u=b>...>

  <subsection|Implementation in FEniCS >

  This required:

  <\enumerate>
    <item>Fixing the FIAT Hermite element definition (e.g. it was returning
    an inconsistent format of <tt|entity_ids>) in dimensions 1 to 3.\ 

    <item>Implementing in FFC the \PHermite transformation\Q described above
    for the element basis functions which are associated to the evaluation of
    partial derivatives at the nodes of the simplex. This means:

    <\itemize-dot>
      <item><tt|representation.py>: for every dof
      <math|<wide|\<psi\>|^><rsub|\<alpha\><rsub|i>>> include the
      coefficients of the other dofs <math|<wide|\<psi\>|^><rsub|\<alpha\><rsub|j>>>
      for <math|j\<neq\>i> in the <tt|dof_data> structure passed to the next
      compiler stage.

      <item><tt|evaluatebasis.py>: using the Jacobian of the geometric
      transformation, compute the new transformed values as per
      <eqref|eq:hermite-first-derivatives>.

      <item><tt|quadratures/quadraturetransformer.py>:\ 
    </itemize-dot>

    <item><todo|Implementing the same transformation for the quadratures>.
    See <tt|ffc/quadratures/quadraturetransformerbase.py>.
  </enumerate>

  \;
</body>

<\initial>
  <\collection>
    <associate|font|stix>
    <associate|font-base-size|11>
    <associate|info-flag|detailed>
    <associate|math-font|math-stix>
    <associate|preamble|false>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-2|<tuple|2|?>>
    <associate|auto-3|<tuple|3|?>>
    <associate|auto-4|<tuple|1|?>>
    <associate|auto-5|<tuple|4|?>>
    <associate|auto-6|<tuple|5|?>>
    <associate|eq:delta-property|<tuple|2|?>>
    <associate|eq:hermite-first-derivatives|<tuple|4|?>>
    <associate|eq:hermite-second-derivatives|<tuple|5|?>>
    <associate|eq:hermite-transform|<tuple|3|?>>
    <associate|eq:linear-problem|<tuple|1|?>>
    <associate|fig:reference-triangle|<tuple|1|?>>
    <associate|footnote-1|<tuple|1|?>>
    <associate|footnr-1|<tuple|1|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|figure>
      <tuple|normal|<with|color|<quote|dark
      red>|<datoms|<macro|x|<resize|<tabular|<tformat|<cwith|1|1|1|1|cell-background|pastel
      red>|<cwith|1|1|1|1|cell-lsep|0fn>|<cwith|1|1|1|1|cell-rsep|0fn>|<cwith|1|1|1|1|cell-bsep|<value|marked-padding>>|<cwith|1|1|1|1|cell-tsep|<value|marked-padding>>|<table|<row|<cell|<arg|x>>>>>>|<plus|1l|0fn>|<plus|1b|<value|marked-padding>>|<minus|1r|0fn>|<minus|1t|<value|marked-padding>>>>|[The
      reference triangle and the mapping <with|mode|<quote|math>|\<b-x\>>...]>>|<pageref|auto-3>>
    </associate>
    <\associate|toc>
      <with|par-left|<quote|1tab>|1<space|2spc>Setting
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1>>

      <with|par-left|<quote|1tab>|2<space|2spc>Cubic Hermite elements in 2D
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1tab>|3<space|2spc>Assembly of the stiffness
      matrix <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|1tab>|4<space|2spc>Implementation in FEniCS
      \ <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>
    </associate>
  </collection>
</auxiliary>