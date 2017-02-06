<TeXmacs|1.99.4>

<style|<tuple|amsart|british>>

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

  <doc-data|<doc-title|Hermite elements:<new-line>introduction and
  implementation>|<doc-author|<author-data|<author-name|Miguel de
  Benito>|<\author-affiliation>
    Universität Augsburg
  </author-affiliation>>>|<doc-date|December 2016>>

  <abstract-data|<abstract|After quickly reviewing the setting and notation
  for finite element approximations, we explain (cubic) Hermite elements and
  the coordinate transformation they require
  <cite-detail|solin_partial_2005|Chapter 6>. We then detail some aspects of
  our implementation for the FEniCS library
  <cite|alnaes_fenics_2015|logg_automated_2012> and finish with some
  examples.>>

  <section|Setting>

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
  <math|V<rsub|h>> for this problem involves:

  <\itemize-dot>
    <item>The triangulation <math|\<cal-T\>>.

    <item>A reference simplex <math|<wide|K|^>\<subset\>\<bbb-R\><rsup|2>>
    with an associated set of linear forms (<dfn|degrees of freedom> or
    <em|dofs>) over <math|P<rsub|p><around*|(|<wide|K|^>|)>>, the space of
    polynomials of degree up to <math|p>, defining a (dual) local polynomial
    basis <math|<wide|\<cal-B\>|^>> of <math|P<rsub|p><around*|(|<wide|K|^>|)>>.<\footnote>
      In <name|FIAT>, the polynomials of <math|<wide|\<cal-B\>|^>> (the shape
      functions) are themselves expressed in a fixed polynomial basis with
      good numerical properties. In particular the one picked
      (<todo|Dubiner?>) yields much better condition numbers for the mass
      matrices than e.g. the basis of monomials
      <cite-detail|logg_automated_2012|Ÿ4.3.2>.
    </footnote>

    <item>A set of invertible affine maps
    <math|<around*|{|\<b-x\><rsup|K>|}><rsub|K\<in\>\<cal-T\>>> between
    <math|<wide|K|^>> and each triangle in <math|\<cal-T\>>.

    <item>A transformation <math|\<cal-F\>:<wide|\<cal-B\>|^><around*|(|P<rsub|p><around*|(|<wide|K|^>|)>|)>\<rightarrow\>\<cal-B\><around*|(|P<rsub|p><around*|(|K|)>|)>>
    preserving certain properties of <math|<wide|\<cal-B\>|^>> and
    <math|<wide|\<cal-B\>|^><rprime|'>>(e.g. duality).

    <item>A <em|local-to-global dof mapping> <math|\<iota\><rsup|K>>,
    assigning local dofs in each <math|K\<in\>\<cal-T\>> to global dofs in
    <math|V<rsub|h>>. This mapping is used in the definition of the global
    space <math|V<rsub|h>> to impose continuity properties of functions in
    <math|V<rsub|h>>: by requiring local dofs to agree\ 

    <item>A method of interpolating functions
    <math|g\<in\>V\<subset\>H<rsup|2><around*|(|\<Omega\>|)>> into
    <math|V<rsub|h>>.

    <item><todo|...>
  </itemize-dot>

  In this note we work in 2 dimensions, with <math|p=3> and a particular
  choice of of degrees of freedom giving rise to the classic <em|Hermite
  finite element> <inactive|<cite|>>. This element is a generalization of
  cubic Hermite interpolating polynomials on 1 dimension and have been around
  in the literature since at least <inactive|<cite|ciarlet_raviart>>.

  <subsection|Conforming and non-conforming discretizations>

  If <math|V<rsub|h>> is a subspace of <math|H<rsup|2><around*|(|\<Omega\>|)>>,
  we say that the discretization is <dfn|conforming>. Hermite elements are
  conforming in 1D but not in 2D <cite-detail|brenner_mathematical_2008|Proposition
  3.3.17>:

  <\quotation>
    The cubic Hermite elements have a continuous normal derivative but not
    full <math|C<rsup|1>> continuity. In particular, the normal derivatives
    may not match at the boundary of two elements, away from the vertices. If
    you want full <math|C<rsup|1>> continuity you will have to use the
    Argyris element or Hsieh-Clough-Tucker or something. I recommend the
    discussion in chapter 6 of Ciarlet's finite element book.<htab|5mm>
    (<hlink|CompSci SE|http://scicomp.stackexchange.com/questions/2012/construction-of-c1-h2-conforming-finite-element-basis-for-triangular-or-te>)
  </quotation>

  <todo|Approximations and errors for non-conforming discretizations...>

  <section|Cubic Hermite elements in 2D>

  <\notation*>
    In the following, we use the indices <math|\<alpha\>\<in\><around*|{|1,4,7|}>>
    and <math|i\<in\><around*|{|0,1,2|}>>. We also set
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
  add one function evaluation <math|\<psi\><rsub|\<alpha\><rsub|0>><around*|(|f|)>\<assign\>f<around*|(|<wide|v|^><rsub|<around*|(|\<alpha\>-1|)>/3>|)>>
  and evaluation of both partial derivatives
  <math|\<psi\><rsub|\<alpha\><rsub|i>><around*|(|f|)>=f<rsub|,i><around*|(|<wide|v|^><rsub|<around*|(|\<alpha\>-1|)>/3>|)>>
  for <math|i=1,2>. Finally at the barycenter <math|<wide|v|^><rsub|0>> of
  <math|<wide|K|^>> add a new evaluation <math|\<psi\><rsub|0><around*|(|f|)>=f<around*|(|<wide|v|^><rsub|0>|)>>.

  <big-figure|<with|gr-mode|<tuple|edit|math-at>|gr-frame|<tuple|scale|1cm|<tuple|0.270004gw|0.180009gh>>|gr-geometry|<tuple|geometry|0.280007par|0.200003par|center>|gr-grid|<tuple|cartesian|<point|0|0>|1>|gr-grid-old|<tuple|cartesian|<point|0|0>|1>|gr-edit-grid-aspect|<tuple|<tuple|axes|none>|<tuple|1|none>|<tuple|10|none>>|gr-edit-grid|<tuple|cartesian|<point|0|0>|1>|gr-edit-grid-old|<tuple|cartesian|<point|0|0>|1>|<graphics||<line|<point|0|0>|<point|2.0|0.0>|<point|0.0|2.0>|<point|0.0|0.0>>>>|<todo|<label|fig:reference-triangle>The
  reference triangle and the mapping <math|\<b-x\>>...>>

  The forms <math|\<cal-B\><rprime|'>\<assign\><around*|{|<wide|\<psi\>|^><rsub|\<alpha\><rsub|0>>,<wide|\<psi\>|^><rsub|\<alpha\><rsub|1>>,<wide|\<psi\>|^><rsub|\<alpha\><rsub|2>>:\<alpha\>=1,4,7|}>\<cup\><around*|{|<wide|\<psi\>|^><rsub|0>|}>>
  are a basis of <math|P<rsub|3><rprime|'>> with dual
  <math|\<cal-B\>=<around*|{|<wide|\<varphi\>|^><rsub|\<alpha\><rsub|0>>,<wide|\<varphi\>|^><rsub|\<alpha\><rsub|1>>,<wide|\<varphi\>|^><rsub|\<alpha\><rsub|2>>:\<alpha\>=1,4,7|}>\<cup\><around*|{|<wide|\<varphi\>|^><rsub|0>|}>>.
  Reindexing, this means that the latter consists of 10 polynomials, called
  <dfn|shape functions>, <math|<wide|\<varphi\>|^><rsub|k>\<in\>P<rsub|3><around*|(|<wide|K|^>|)>,k\<in\><around*|{|0,\<ldots\>,9|}>>
  such that

  <\equation>
    <label|eq:delta-property><wide|\<psi\>|^><rsub|j><around*|(|<wide|\<varphi\>|^><rsub|k>|)>=\<delta\><rsub|j\<nocomma\>k>.
  </equation>

  In particular <math|<wide|\<psi\>|^><rsub|1><around*|(|<wide|\<varphi\>|^><rsub|1>|)>=<wide|\<varphi\>|^><rsub|1><around*|(|<wide|v|^><rsub|1>|)>=1>,
  <math|<wide|\<psi\>|^><rsub|2><around*|(|<wide|\<varphi\>|^><rsub|2>|)>=\<partial\><rsub|1><wide|\<varphi\>|^><rsub|2><around*|(|<wide|v|^><rsub|1>|)>=1>,
  but <math|<wide|\<psi\>|^><rsub|1><around*|(|<wide|\<varphi\>|^><rsub|2>|)>=<wide|\<varphi\>|^><rsub|2><around*|(|<wide|v|^><rsub|1>|)>=0>
  and so on.

  This reference shape functions are then mapped onto a <dfn|local basis over
  the physical cell> <math|K> by means of a mapping

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
  \Phermite\Q pairs of dofs <math|<around*|(|<wide|\<varphi\>|^><rsub|\<alpha\><rsub|1>>,<wide|\<varphi\>|^><rsub|\<alpha\><rsub|2>>|)>>
  we require a different transformation, namely:

  <\equation>
    <label|eq:hermite-transform><around*|(|\<varphi\><rsub|\<alpha\><rsub|1>>,\<varphi\><rsub|\<alpha\><rsub|2>>|)><rsup|\<top\>>=\<cal-F\><around*|(|<wide|\<varphi\>|^><rsub|\<alpha\><rsub|1>>,<wide|\<varphi\>|^><rsub|\<alpha\><rsub|2>>|)>\<assign\><around*|(|D
    \<b-x\>\<cdot\><around*|(|<wide|\<varphi\>|^><rsub|\<alpha\><rsub|1>>,<wide|\<varphi\>|^><rsub|\<alpha\><rsub|2>>|)><rsup|\<top\>>|)>\<circ\><wide|\<b-x\>|^>,
  </equation>

  or, component-wise:

  <\equation*>
    \<varphi\><rsub|\<alpha\><rsub|i>><around*|(|x|)>=\<nabla\>\<b-x\><rsub|i><around*|(|<wide|\<b-x\>|^><around*|(|x|)>|)>\<cdot\><around*|(|<wide|\<varphi\>|^><rsub|\<alpha\><rsub|1>>,<wide|\<varphi\>|^><rsub|\<alpha\><rsub|2>>|)><around*|(|<wide|\<b-x\>|^><around*|(|x|)>|)>=<around*|(|\<b-x\><rsub|i,j>*<wide|\<varphi\>|^><rsub|\<alpha\><rsub|j>>|)><around*|(|<wide|\<b-x\>|^><around*|(|x|)>|)>.
  </equation*>

  <todo|We can check that <eqref|eq:delta-property> holds>: for the point
  evaluations this is clear and for the others:

  <math|\<psi\><rsub|2><around*|(|\<varphi\><rsub|2>|)>=\<partial\><rsub|1>\<varphi\><rsub|2><around*|(|v<rsub|1>|)>=\<ldots\>>,

  Notice that at vertex <math|v<rsub|\<alpha\>>> of <math|K>, each new basis
  function is a linear combination of all the old ones at
  <math|<wide|v|^><rsub|\<alpha\>>> with some coefficients set to zero: we
  can write

  <\equation*>
    \<Phi\><rsub|\<alpha\>>\<assign\><around*|(|\<varphi\><rsub|\<alpha\><rsub|0>>,\<varphi\><rsub|\<alpha\><rsub|1>>,\<varphi\><rsub|\<alpha\><rsub|2>>|)><rsup|\<top\>>,<space|2em><wide|\<Phi\>|^><rsub|\<alpha\>>\<assign\><around*|(|<wide|\<varphi\>|^><rsub|\<alpha\><rsub|0>>,<wide|\<varphi\>|^><rsub|\<alpha\><rsub|1>>,<wide|\<varphi\>|^><rsub|\<alpha\><rsub|2>>|)><rsup|\<top\>>,
  </equation*>

  and define

  <\equation*>
    \<Phi\><rsub|\<alpha\>><around*|(|x|)>=\<cal-F\><around*|(|<wide|\<Phi\>|^><rsub|\<alpha\>>|)>\<assign\><wide*|<matrix|<tformat|<table|<row|<cell|1>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|\<b-x\><rsub|1,1>>|<cell|\<b-x\><rsub|1,2>>>|<row|<cell|0>|<cell|\<b-x\><rsub|2,1>>|<cell|\<b-x\><rsub|2,2>>>>>>|\<wide-squnderbrace\>><rsub|H><around*|(|<wide|\<b-x\>|^><around*|(|x|)>|)>*<vector|<tformat|<cwith|1|3|1|1|cell-lsep|0>|<cwith|1|3|1|1|cell-rsep|0>|<cwith|1|3|1|1|cell-bsep|0>|<cwith|1|3|1|1|cell-tsep|0>|<table|<row|<cell|<wide|\<varphi\>|^><rsub|><rsub|\<alpha\><rsub|0>>>>|<row|<cell|<wide|\<varphi\>|^><rsub|><rsub|\<alpha\><rsub|1>>>>|<row|<cell|<wide|\<varphi\>|^><rsub|><rsub|\<alpha\><rsub|2>>>>>>><rsub|\|<wide|\<b-x\>|^><around*|(|x|)>>.
  </equation*>

  Using this <math|\<cal-F\>> means that we have to take the factor <math|D
  \<b-x\>> into account when differentiating the basis functions. Using the
  chain rule and the fact that <math|D H=0>, we obtain

  <\equation*>
    D \<Phi\><rsub|\<alpha\>><around*|(|x|)>=<around*|(|H*D
    <wide|\<Phi\>|^><rsub|\<alpha\>>|)><around*|(|<wide|\<b-x\>|^><around*|(|x|)>|)>*D<wide|\<b-x\>|^><around*|(|x|)>,
  </equation*>

  where the last term <math|D<wide|\<b-x\>|^>> comes from the affine change
  of coordinates and is also present in the Lagrange case. Component-wise:

  <\equation*>
    \<varphi\><rsub|\<alpha\><rsub|i>,j><around*|(|x|)>=<around*|(|H<rsub|<around*|(|i+1|)>\<nocomma\>k>*<wide|\<varphi\>|^><rsub|\<alpha\><rsub|k>,l>|)><rsub|\|<wide|\<b-x\>|^><around*|(|x|)>>*<wide|\<b-x\>|^><rsub|l,j><around*|(|x|)>,<space|2em>i\<in\><around*|{|0,1,2|}>.
  </equation*>

  Each partial derivative <math|\<varphi\><rsub|\<alpha\><rsub|i>,j>> for
  <math|i\<in\><around*|{|1,2|}>> is a linear combination of the derivatives
  of both reference basis functions <math|\<nabla\><wide|\<varphi\>|^><rsub|\<alpha\><rsub|k><rsub|>>>
  with weights given by <math|\<nabla\> \<b-x\><rsub|i>>:

  <\equation>
    <label|eq:hermite-first-derivatives>\<varphi\><rsub|\<alpha\><rsub|i>,j><around*|(|x|)>=<around*|(|\<b-x\><rsub|i,k>*<wide|\<varphi\>|^><rsub|\<alpha\><rsub|k>,l>|)><rsub|\|<wide|\<b-x\>|^><around*|(|x|)>>*<wide|\<b-x\>|^><rsub|l,j><around*|(|x|)>,<space|2em>i\<in\><around*|{|1,2|}>.
  </equation>

  For the second derivatives we have a similar expression, again thanks to
  <math|D<rsup|2> \<b-x\>=0>:

  <\equation>
    <label|eq:hermite-second-derivatives>\<varphi\><rsub|\<alpha\><rsub|i>,j\<nocomma\>k>=\<b-x\><rsub|i,p>*<wide|\<varphi\>|^><rsub|\<alpha\><rsub|p>,q\<nocomma\>r>*<wide|\<b-x\>|^><rsub|q,j>*<wide|\<b-x\>|^><rsub|r,k>,<space|2em>i\<in\><around*|{|1,2|}>,
  </equation>

  i.e.

  <\equation*>
    D<rsup|2> \<varphi\><rsub|\<alpha\><rsub|i>>=D<rsup|\<top\>><wide|\<b-x\>|^>*<around*|(|\<b-x\><rsub|i,1>*D<rsup|2><wide|\<varphi\>|^><rsub|\<alpha\><rsub|1>>+\<b-x\><rsub|i,2>*D<rsup|2><wide|\<varphi\>|^><rsub|\<alpha\><rsub|2>>|)>*D
    <wide|\<b-x\>|^>,<space|2em>i\<in\><around*|{|1,2|}>.
  </equation*>

  And if we use <math|H> (the sum is over <math|p\<in\><around*|{|0,1,2|}>>):

  <\equation*>
    \<varphi\><rsub|\<alpha\><rsub|i>,j\<nocomma\>k>=H<rsub|i+1,p+1>*<wide|\<varphi\>|^><rsub|\<alpha\><rsub|p>,q\<nocomma\>r>*<wide|\<b-x\>|^><rsub|q,j>*<wide|\<b-x\>|^><rsub|r,k>,<space|2em>i\<in\><around*|{|0,1,2|}>.
  </equation*>

  Here as before, we obtain a linear combination of the second derivatives of
  both reference basis functions <math|D<rsup|2>
  <wide|\<varphi\>|^><rsub|\<alpha\><rsub|j>>> with weights given by
  <math|\<nabla\> \<b-x\><rsub|i>>. The terms
  <math|D<rsup|\<top\>><wide|\<b-x\>|^>> and <math|D<wide|\<b-x\>|^>> are
  once more the same as if we had only used the inverse mapping
  <math|<wide|\<b-x\>|^>> to change variables.

  With <math|3> derivatives <todo|generalise to <math|k> derivatives with a
  multiindex>:

  <\equation*>
    \<varphi\><rsub|\<alpha\><rsub|i>,j\<nocomma\>k\<nocomma\>l>=H<rsub|i+1,p+1>*<wide|\<varphi\>|^><rsub|\<alpha\><rsub|p>,q\<nocomma\>r\<nocomma\>s>*<wide|\<b-x\>|^><rsub|q,j>*<wide|\<b-x\>|^><rsub|r,k>*<wide|\<b-x\>|^><rsub|s,l>.
  </equation*>

  <section|Assembly of the stiffness matrix>

  In order to discretise <eqref|eq:linear-problem> we require the action of
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

  <subsection|Computing integrals>

  The contribution <math|a<rsup|K>> from element <math|K> to the stiffness
  matrix will typically be of integral form. For example, if we are
  discretising the biharmonic equation <math|\<Delta\><rsup|2> u=f>, assuming
  homogeneous boundary conditions, the corresponding bilinear form is
  <math|a<around*|(|u,v|)>=<big|int>\<Delta\>u*\<Delta\>v>, so we need to
  compute

  <\equation*>
    a<rsup|K><around*|(|\<varphi\><rsub|i<rsub|1>><rsup|1>,\<varphi\><rsup|2><rsub|i<rsub|2>>|)>=<big|int><rsub|K>\<Delta\>\<varphi\><rsub|i<rsub|1>><rsup|1>*\<mathLaplace\>\<varphi\><rsup|2><rsub|i<rsub|2>>*\<mathd\>x=<big|int><rsub|<wide|K|^>><todo|\<ldots\>>.
  </equation*>

  <todo|As can be seen, the transformation introduced above plays here a role
  as well...> Notice that even though <math|V<rsub|h>\<nsubset\>H<rsup|2>> we
  can still compute the Laplacian in the interior of each element <math|K>.
  However, the normal derivative might present discontinuities, <todo|so...>

  \;

  <section|Interpolation>

  <todo|See <cite-detail|brenner_mathematical_2008|Ÿ4.8, p. 121> for info on
  interpolation of non-smooth functions.>

  We have at least the following options for interpolating functions
  <math|g\<in\>H<rsup|2><around*|(|\<Omega\>|)>> into <todo|notation Hermite
  space>, from slowest and most accurate to fastest and less accurate
  <cite-detail|solin_partial_2005|Ÿ6.3.8>.

  <\enumerate>
    <item>Use a global orthogonal projection to compute the best interpolant
    in the <math|H<rsup|2>>-sense. This means solving the system: <todo|...>.

    <item>Use nodal interpolation of vertex and derivatives plus local
    orthogonal projections in the element interiors. <todo|...>

    <item>Compute the nodal interpolant.

    <\question*>
      Does the nodal interpolant include the dofs for derivative evaluation?
      If yes, how should one approximate the derivatives? Recall that the
      nodal interpolant is given by

      <\equation*>
        I<rsub|K><around*|(|g|)>=<big|sum><rsub|i=1><rsup|n<rsub|K>>L<rsub|i><around*|(|g|)>*\<varphi\><rsub|i>
      </equation*>
    </question*>
  </enumerate>

  <section|Implementation in FEniCS >

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
      transformation, compute the basis evaluations as per
      <eqref|eq:hermite-first-derivatives>. The intermediate variable
      <python|dof_data['hermite_node_offset']> is used in
      <python|_compute_values()> to compute the right linear combination of
      reference basis functions.

      <item><tt|quadratures/quadraturetransformer.py>:\ 
    </itemize-dot>

    <item><todo|Implementing the same transformation for the quadratures>.
    See <tt|ffc/quadratures/quadraturetransformerbase.py>.

    <item><todo|Implementing a new interpolation method>. FIXME: why does
    this work out of the box, even though I didn't code it yet? See
    <tt|ffc/interpolatevertexvalues.py>.
  </enumerate>

  <section|The Euler-Bernoulli beam model>

  Model: 2D beam reduced to 1D problem. We study the deformation of the
  midplane under the assumptions that after the deformation the normals to
  the midplane:

  <\itemize-dot>
    <item>do not bend,

    <item>do not stretch,

    <item>remain orthogonal to the midplane.
  </itemize-dot>

  This theory is adequate for thin beams and is intended for small strains
  even with large global deformations: it is a physically linear but
  geometrically non-linear theory.

  For thicker beams, Timoshenko's theory, which accounts for internal shear
  forces, yields more accurate predictions.

  Fix <math|\<omega\>=<around*|(|a,b|)>> to be the midplane of the beam
  <math|\<Omega\>=\<omega\>\<times\><around*|(|-h/2,h/2|)>\<subset\>\<bbb-R\><rsup|2>>.

  <subsection|Derivation>

  Conservation of momentum and linear constitutive relations yield
  <cite-detail|solin_partial_2005|Ÿ6.1.1>:

  <\equation>
    <label|eq:euler-bernoulli><frac|\<mathd\><rsup|2>|\<mathd\>x<rsup|2>>
    <around*|(|b<around|(|x|)>*<frac|\<mathd\><rsup|2>|\<mathd\>x<rsup|2>>
    u<around|(|x|)>|)>=f<around*|(|x|)>
  </equation>

  where <math|b<around|(|x|)>=E<around|(|x|)>*I<around|(|x|)>> is the product
  of Young's modulus <math|E> and the area moment of inertia of the beam
  <math|I>.

  For steel (with \<less\> 0.3% carbon) at 21\<degree\>C,
  <math|E=203.4\<cdot\>10<rsup|9>> Pa and since our beam has a constant
  square cross section of side 0.01m, we have <math|I=8\<cdot\>10<rsup|-10>>.
  <todo|compute...>

  <subsection|Weak formulation>

  Write <math|\<nabla\>u=<frac|\<mathd\>|\<mathd\>x> u> and
  <math|\<Delta\>u=<frac|\<mathd\><rsup|2>|\<mathd\>x<rsup|2>>*u>. Let
  <math|V> be a subspace of <math|H<rsup|2><around*|(|\<omega\>|)>> to be
  specified later. Multiplying <eqref|eq:euler-bernoulli> by a test function
  <math|v\<in\>V> and integrating by parts we arrive at

  <\equation*>
    <big|int><rsub|\<omega\>>b*\<Delta\>u*\<Delta\>v*\<mathd\>x-<around|[|b*\<Delta\>u*\<nabla\>v|]><rsub|a><rsup|b>+<around|[|\<nabla\><around|(|b*\<Delta\>u|)>*v|]><rsub|a><rsup|b>=<big|int><rsub|\<omega\>>f*v\<mathd\>x.
  </equation*>

  In order for these integrals to make sense we may take
  <math|b\<in\>L<rsup|\<infty\>><around*|(|\<omega\>|)>> and
  <math|f\<in\>L<rsup|2><around|(|\<omega\>|)>> (or even
  <math|H<rsup|-2><around*|(|\<omega\>|)>>). As always, the definition of
  <math|V> and the final form of the equation are determined by our choices
  for the four boundary conditions that have to be specified:

  <paragraph|Essential boundary conditions>Split the boundary
  <math|\<gamma\>=<around*|{|a,b|}>> into (possibly empty) sets
  <math|\<gamma\><rsub|u>,\<gamma\><rsub|\<theta\>>\<subset\><around*|{|a,b|}>>.
  We fix either the <strong|deflections> <math|u<around|(|\<alpha\>|)>> at
  <math|\<alpha\>\<in\>\<gamma\><rsub|u>> or the <strong|slopes>
  <math|u<rprime|'><around|(|\<alpha\>|)>> at
  <math|\<alpha\>\<in\>\<gamma\><rsub|\<theta\>>> or both. These conditions
  are incorporated into the definition of <math|V>. For example if we
  <em|clamp> the beam at an horizontal position we have

  <\equation*>
    V=V<rsub|<text|clamped>>=<around*|{|v\<in\>H<rsup|2><around*|(|\<omega\>|)>:v<around*|(|\<alpha\>|)>=v<rsub|\<alpha\>>,v<rprime|'><around*|(|\<alpha\>|)>=v<rprime|'><rsub|\<alpha\>>,\<alpha\>=1,2|}>.
  </equation*>

  <paragraph|Natural boundary conditions>For
  <math|\<alpha\>\<in\>\<gamma\><rsub|M>\<subset\><around|{|a,b|}>>, we can
  fix the <strong|bending moment>:

  <\equation*>
    M<around*|(|\<alpha\>|)>=<around|(|b*\<Delta\>u|)><around|(|\<alpha\>|)>,
  </equation*>

  which for general <math|x\<in\>\<omega\>> is the torque exerted by forces
  surrounding <math|x>. If, for example we find solutions such that
  <math|M<around*|(|a|)>=0>, then we are assuming that the left end of the
  beam is free to rotate, i.e. that it undergoes no bending due to torque.
  Alternatively we can set the <strong|shear force> at
  <math|\<alpha\>\<in\>\<gamma\><rsub|F>\<subset\><around*|{|a,b|}>>:

  <\equation*>
    F<around*|(|\<alpha\>|)>=<around*|[|\<nabla\>*<around|(|b*\<Delta\>*u|)>|]><around|(|\<alpha\>|)>,
  </equation*>

  which is the resultant of transversal forces at <math|x\<in\>\<omega\>>.
  This will in general be zero at the ends.

  After choosing some combination of the conditions the problem is: Find
  <math|u\<in\>V> such that for all <math|v\<in\>V>:

  <\equation*>
    <big|int><rsub|\<omega\>>b*\<Delta\>u*\<Delta\>v*\<mathd\>x=<big|int><rsub|\<omega\>>f*v*\<mathd\>x+<big|int><rsub|\<gamma\><rsub|M>>M*\<nabla\>v*\<mathd\>s-<big|int><rsub|\<gamma\><rsub|F>>F*v*\<mathd\>s
  </equation*>

  <math|V,\<gamma\><rsub|M>,\<gamma\><rsub|F>> to be specified and such that
  <math|\<gamma\><rsub|u>\<uplus\>\<gamma\><rsub|F>=<around*|{|a,b|}>> and
  <math|\<gamma\><rsub|\<theta\>>\<uplus\>\<gamma\><rsub|M>=<around*|{|a,b|}>>.

  <subsection|Existence and uniqueness>

  A simple application of Hölder and Poincaré is required to check the
  conditions for existence in Lax-Milgram. For ellipticity, assume that
  <math|b> is a.e. bounded away from and above zero.

  <subsection|Discretization>

  Recall that <math|H<rsup|2><around|(|a,b|)>\<in\>C<rsup|1><around|(|a,b|)>>
  by the Sobolev embeddings.

  Even though cuadratic polynomials might be enough, we want to construct a
  Ciarlet finite element, i.e. with a unisolvent set of degrees of freedom,
  which requires at least cubic polynomials ... [elaborate, see
  <cite-detail|solin_partial_2005|p. 218>]

  So we use cubic Hermite elements, <todo|which are <math|H<rsup|2>>
  conforming in <math|<with|math-font|Bbb|\<bbb-R\>>>>.

  \;

  \;

  <\bibliography|bib|tm-alpha|hermite.bib>
    <\bib-list|4>
      <bibitem*|ABH+15><label|bib-alnaes_fenics_2015>Martin<nbsp>S.<nbsp>Alnaes,
      Jan Blechta, Johan Hake, August Johansson, Benjamin Kehlet, Anders
      Logg, Chris Richardson, Johannes Ring,
      Marie<nbsp>E.<nbsp>Rognes<localize|, and
      >Garth<nbsp>N.<nbsp>Wells.<newblock> The FEniCS Project Version
      1.5.<newblock> <with|font-shape|italic|Archive of Numerical Software>,
      3(100), 2015.<newblock>

      <bibitem*|BS08><label|bib-brenner_mathematical_2008>Susanne<nbsp>C.<nbsp>Brenner<localize|
      and >L.<nbsp>Ridgway Scott.<newblock> <with|font-shape|italic|The
      Mathematical Theory of Finite Element Methods>.<newblock>
      <localize|Number><nbsp>15<localize| in >Texts in Applied Mathematics.
      Springer New York, New York, NY, 3<localize| edition>, 2008.<newblock>

      <bibitem*|LMW12><label|bib-logg_automated_2012>Anders Logg, Kent-Andre
      Mardal<localize|, and >Garth<nbsp>N.<nbsp>Wells<localize|,
      editors>.<newblock> <with|font-shape|italic|Automated solution of
      differential equations by the Finite Element method>.<newblock>
      <localize|Number><nbsp>84<localize| in >Lecture notes in computational
      science and engineering. Springer, 2012.<newblock> DOI
      10.1007/978-3-642-23099-8.<newblock>

      <bibitem*|Sol05><label|bib-solin_partial_2005>Pavel Solin.<newblock>
      <with|font-shape|italic|Partial differential equations and the finite
      element method>.<newblock> Pure and applied Mathematics. Dec
      2005.<newblock>
    </bib-list>
  </bibliography>

  \;
</body>

<\initial>
  <\collection>
    <associate|font|stix>
    <associate|font-base-size|11>
    <associate|info-flag|detailed>
    <associate|math-font|math-stix>
    <associate|page-medium|papyrus>
    <associate|preamble|false>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-10|<tuple|6.1|?>>
    <associate|auto-11|<tuple|6.2|?>>
    <associate|auto-12|<tuple|1|?>>
    <associate|auto-13|<tuple|2|?>>
    <associate|auto-14|<tuple|6.3|?>>
    <associate|auto-15|<tuple|6.4|?>>
    <associate|auto-16|<tuple|6.4|?>>
    <associate|auto-2|<tuple|1.1|1>>
    <associate|auto-3|<tuple|2|2>>
    <associate|auto-4|<tuple|1|2>>
    <associate|auto-5|<tuple|3|3>>
    <associate|auto-6|<tuple|3.1|3>>
    <associate|auto-7|<tuple|4|3>>
    <associate|auto-8|<tuple|5|4>>
    <associate|auto-9|<tuple|6|4>>
    <associate|bib-alnaes_fenics_2015|<tuple|ABH+15|4>>
    <associate|bib-brenner_mathematical_2008|<tuple|BS08|?>>
    <associate|bib-logg_automated_2012|<tuple|LMW12|4>>
    <associate|bib-solin_partial_2005|<tuple|Sol05|4>>
    <associate|eq:delta-property|<tuple|2|2>>
    <associate|eq:euler-bernoulli|<tuple|6|?>>
    <associate|eq:hermite-first-derivatives|<tuple|4|3>>
    <associate|eq:hermite-second-derivatives|<tuple|5|3>>
    <associate|eq:hermite-transform|<tuple|3|2>>
    <associate|eq:linear-problem|<tuple|1|1>>
    <associate|fig:reference-triangle|<tuple|1|2>>
    <associate|footnote-1|<tuple|1|1>>
    <associate|footnote-2|<tuple|2|?>>
    <associate|footnr-1|<tuple|1|1>>
    <associate|footnr-2|<tuple|2|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      solin_partial_2005

      alnaes_fenics_2015

      logg_automated_2012

      brenner_mathematical_2008

      brenner_mathematical_2008

      solin_partial_2005

      solin_partial_2005

      solin_partial_2005
    </associate>
    <\associate|figure>
      <tuple|normal|<with|color|<quote|dark
      red>|<datoms|<macro|x|<resize|<tabular|<tformat|<cwith|1|1|1|1|cell-background|pastel
      red>|<cwith|1|1|1|1|cell-lsep|0fn>|<cwith|1|1|1|1|cell-rsep|0fn>|<cwith|1|1|1|1|cell-bsep|<value|marked-padding>>|<cwith|1|1|1|1|cell-tsep|<value|marked-padding>>|<table|<row|<cell|<arg|x>>>>>>|<plus|1l|0fn>|<plus|1b|<value|marked-padding>>|<minus|1r|0fn>|<minus|1t|<value|marked-padding>>>>|[The
      reference triangle and the mapping <with|mode|<quote|math>|\<b-x\>>...]>>|<pageref|auto-4>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1.<space|2spc>Setting>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1tab>|1.1.<space|2spc>Conforming and
      non-conforming discretizations <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2.<space|2spc>Cubic
      Hermite elements in 2D> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3.<space|2spc>Assembly
      of the stiffness matrix> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5><vspace|0.5fn>

      <with|par-left|<quote|1tab>|3.1.<space|2spc>Computing integrals
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|4.<space|2spc>Interpolation>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|5.<space|2spc>Implementation
      in FEniCS > <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|6.<space|2spc>The
      Euler-Bernoulli beam model> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9><vspace|0.5fn>

      <with|par-left|<quote|1tab>|6.1.<space|2spc>Derivation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10>>

      <with|par-left|<quote|1tab>|6.2.<space|2spc>Weak formulation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <with|par-left|<quote|4tab>|Essential boundary conditions:
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Natural boundary conditions:
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13><vspace|0.15fn>>

      <with|par-left|<quote|1tab>|6.3.<space|2spc>Existence and uniqueness
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-14>>

      <with|par-left|<quote|1tab>|6.4.<space|2spc>Discretization
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-16><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>