<TeXmacs|1.99.4>

<style|<tuple|amsart|british|better-amsart>>

<\body>
  <\hide-preamble>
    <assign|vector|<macro|x|<shrink-inline|<left|(><tformat|<cwith|1|-1|1|-1|cell-halign|c><cwith|1|-1|1|1|cell-lsep|0><cwith|1|-1|1|1|cell-rsep|0><cwith|1|-1|1|1|cell-bsep|0.2em><cwith|1|-1|1|1|cell-tsep|0.2em>|<arg|x>><right|)>>>>

    <assign|quotation|<\macro|body>
      <\padded>
        <\indent-both|<value|quote-left-indentation>|<value|quote-right-indentation>>
          <surround|<yes-indent>||<em|<arg|body>>>
        </indent-both>
      </padded>
    </macro>>

    <assign|dfn|<macro|body|<strong|<arg|body>>>>
  </hide-preamble>

  <doc-data|<doc-title|Hermite elements:<new-line>introduction and
  implementation>|<doc-author|<author-data|<author-name|Miguel de
  Benito>|<\author-affiliation>
    Universität Augsburg
  </author-affiliation>>>|<doc-date|December 2016>>

  <abstract-data|<abstract|We briefly review the setting and notation for
  finite element approximations, using cubic Hermite elements and the
  coordinate transformation they require. We then detail some aspects of our
  implementation for the FEniCS library <cite|alnaes_fenics_2015|logg_automated_2012>
  and finish with some applications.>>

  <\small>
    <\table-of-contents|toc>
      <vspace*|1fn><with|font-series|bold|math-font-series|bold|1.<space|2spc>Setting>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|bold|math-font-series|bold|2.<space|2spc>Cubic
      Hermite elements in 2D> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>

      <with|par-left|1tab|2.1.<space|2spc>The reference simplex and shape
      functions <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|1tab|2.2.<space|2spc>The local basis
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <vspace*|1fn><with|font-series|bold|math-font-series|bold|3.<space|2spc>Assembly
      of the stiffness matrix> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6><vspace|0.5fn>

      <with|par-left|1tab|3.1.<space|2spc>The global basis
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>

      <with|par-left|1tab|3.2.<space|2spc>The stiffness matrix
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8>>

      <with|par-left|1tab|3.3.<space|2spc>Computing integrals
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9>>

      <vspace*|1fn><with|font-series|bold|math-font-series|bold|4.<space|2spc>Interpolation>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10><vspace|0.5fn>

      <with|par-left|1tab|4.1.<space|2spc>The local nodal interpolant
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <with|par-left|1tab|4.2.<space|2spc><with|mode|math|W<rsup|m,2>>-projection
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12>>

      <vspace*|1fn><with|font-series|bold|math-font-series|bold|5.<space|2spc>Implementation
      in <with|font-shape|small-caps|FEniCS> >
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13><vspace|0.5fn>

      <vspace*|1fn><with|font-series|bold|math-font-series|bold|6.<space|2spc>The
      Euler-Bernoulli beam model> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-14><vspace|0.5fn>

      <with|par-left|1tab|6.1.<space|2spc>Derivation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15>>

      <with|par-left|1tab|6.2.<space|2spc>Weak formulation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-16>>

      <with|par-left|4tab|Essential boundary conditions
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-17><vspace|0.15fn>>

      <with|par-left|4tab|Natural boundary conditions
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-18><vspace|0.15fn>>

      <with|par-left|1tab|6.3.<space|2spc>Existence and uniqueness
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-19>>

      <with|par-left|1tab|6.4.<space|2spc>Discretization
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-20>>

      <vspace*|1fn><with|font-series|bold|math-font-series|bold|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-21><vspace|0.5fn>
    </table-of-contents>
  </small>

  <section|Setting>

  We fix a polygonal domain <math|\<Omega\>\<subset\>\<bbb-R\><rsup|2>> (the
  <dfn|physical domain>) and a sufficiently regular triangulation<\footnote>
    <em|\SA triangulation of a polygonal domain \<Omega\> is a subdivision
    consisting of triangles having the property that no vertex of any
    triangle lies in the interior of an edge of another triangle.\T>
    <cite-detail|brenner_mathematical_2008|Ÿ3.3.11>.
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

  The finite element method builds a suitable discretisation <math|V<rsub|h>>
  of the space <math|V> and projects Problem <eqref|eq:linear-problem> onto
  it for its numerical solution. In broad terms, the standard construction is
  a discrete polynomial space <math|V<rsub|h>> involving:

  <\itemize-dot>
    <item>The triangulation <math|\<cal-T\>>.

    <item>A reference simplex <math|<wide|K|^>\<subset\>\<bbb-R\><rsup|2>>
    with an associated set of linear forms (<em|degrees of freedom> or
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
    <math|V<rsub|h>> by requiring local dofs to agree on matching vertices,
    faces, etc.

    <item>A method of interpolating functions
    <math|g\<in\>V\<subset\>H<rsup|2><around*|(|\<Omega\>|)>> into
    <math|V<rsub|h>>, e.g. <math|H<rsup|2>>-projection.

    <item>A method of computing integrals of linear forms defined over
    <math|V<rsub|h>>.

    <item><todo|...>
  </itemize-dot>

  In this note we work in 2 dimensions, with polynomials of degree <math|p=3>
  and a particular choice of degrees of freedom giving rise to the classic
  <em|Hermite finite element>. This element is a generalization of cubic
  Hermite interpolating polynomials on 1 dimension and have been around in
  the literature since at least <inactive|<cite|ciarlet_raviart>>. Hermite
  elements are <math|H<rsup|2>>-conforming in 1D but not in 2D
  <cite-detail|brenner_mathematical_2008|Proposition 3.3.17>.<\footnote>
    From <hlink|CompSci SE|http://scicomp.stackexchange.com/questions/2012/construction-of-c1-h2-conforming-finite-element-basis-for-triangular-or-te>:

    <\quotation>
      \SThe cubic Hermite elements have a continuous normal derivative but
      not full <math|C<rsup|1>> continuity. In particular, the normal
      derivatives may not match at the boundary of two elements, away from
      the vertices. If you want full <math|C<rsup|1>> continuity you will
      have to use the Argyris element or Hsieh-Clough-Tucker or something. I
      recommend the discussion in chapter 6 of Ciarlet's finite element
      book.\T
    </quotation>
  </footnote>

  <todo|Approximations and errors for non-conforming discretizations...>

  <section|Cubic Hermite elements in 2D>

  <\notation*>
    In this section <math|\<alpha\>,i,p\<in\><around*|{|1,2,3|}>> and
    <math|j,k,l,q,r,s\<in\><around*|{|1,2|}>>.
  </notation*>

  <subsection|The reference simplex and shape functions>

  Let <math|<wide|K|^>> be a <dfn|reference simplex> in 2D with vertices
  <math|<wide|v|^><rsup|\<alpha\>>> and <math|\<b-x\>:<wide|K|^>\<rightarrow\>K>
  be a non-degenerate affine transformation into one of the cells of the
  physical domain <math|\<Omega\>>. Denote the inverse mapping by
  <math|<wide|\<b-x\>|^>:K\<rightarrow\><wide|K|^>> (Figure
  <reference|fig:reference-triangle>).

  Let <math|P<rsub|3>=P<rsub|3><around*|(|<wide|K|^>|)>> be the space of
  polynomials of order up to 3 defined over <math|<wide|K|^>>. We define 10
  linear forms (the <dfn|degrees of freedom> or <dfn|dofs>) in the dual
  <math|P<rsub|3><rprime|'>\<assign\>P<rsub|3><around*|(|<wide|K|^>|)><rprime|'>>
  as follows: for each vertex <math|<wide|v|^><rsup|\<alpha\>>> in
  <math|<wide|K|^>> add one point evaluation

  <\equation*>
    \<psi\><rsup|\<alpha\>><rsub|1><around*|(|f|)>\<assign\>f<around*|(|<wide|v|^><rsup|\<alpha\>>|)>,
  </equation*>

  and evaluation of both partial derivatives

  <\equation*>
    \<psi\><rsub|i+1><rsup|\<alpha\>><around*|(|f|)>\<assign\>f<rsub|,i><around*|(|<wide|v|^><rsup|\<alpha\>>|)>,<space|2em>i\<in\><around*|{|1,2|}>.
  </equation*>

  Finally at the barycenter <math|<wide|v|^><rsup|0>> of <math|<wide|K|^>>
  add a last evaluation

  <\equation*>
    \<psi\><rsub|0><around*|(|f|)>\<assign\>f<around*|(|<wide|v|^><rsup|0>|)>.
  </equation*>

  <big-figure|<with|gr-mode|<tuple|edit|math-at>|gr-frame|<tuple|scale|1cm|<tuple|0.270004gw|0.180009gh>>|gr-geometry|<tuple|geometry|0.280007par|0.200003par|center>|gr-grid|<tuple|empty>|gr-grid-old|<tuple|cartesian|<point|0|0>|1>|gr-edit-grid-aspect|<tuple|<tuple|axes|none>|<tuple|1|none>|<tuple|10|none>>|gr-edit-grid|<tuple|empty>|gr-edit-grid-old|<tuple|cartesian|<point|0|0>|1>|<graphics||<line|<point|0|0>|<point|2.0|0.0>|<point|0.0|2.0>|<point|0.0|0.0>>>>|<todo|<label|fig:reference-triangle>The
  reference triangle and the mapping <math|\<b-x\>>...>>

  The linear forms <math|\<cal-B\><rprime|'>\<assign\><around*|{|<wide|\<psi\>|^><rsup|\<alpha\>><rsub|i>|}>\<cup\><around*|{|<wide|\<psi\>|^><rsub|0>|}>>
  are a basis of <math|P<rsub|3><rprime|'>> with dual
  <math|\<cal-B\>=<around*|{|<wide|\<varphi\>|^><rsub|i><rsup|\<alpha\>>|}>\<cup\><around*|{|<wide|\<varphi\>|^><rsub|0>|}>>
  <cite-detail|brenner_mathematical_2008|3.x.11>. Reindexing with
  <math|k=3*<around*|(|\<alpha\>-1|)>+i>, this means that <math|\<cal-B\>>
  consists of 10 polynomials <math|<wide|\<varphi\>|^><rsub|k>\<in\>P<rsub|3><around*|(|<wide|K|^>|)>>
  called <dfn|shape functions>, such that the following duality property
  holds:

  <\equation>
    <label|eq:delta-property><wide|\<psi\>|^><rsub|j><around*|(|<wide|\<varphi\>|^><rsub|k>|)>=\<delta\><rsub|j\<nocomma\>k>,<space|2em>j,k\<in\><around*|{|0,\<ldots\>,9|}>.
  </equation>

  In particular <math|<wide|\<psi\>|^><rsub|1><around*|(|<wide|\<varphi\>|^><rsub|1>|)>=<wide|\<varphi\>|^><rsub|1><around*|(|<wide|v|^><rsub|1>|)>=1>,
  <math|<wide|\<psi\>|^><rsub|2><around*|(|<wide|\<varphi\>|^><rsub|2>|)>=\<partial\><rsub|1><wide|\<varphi\>|^><rsub|2><around*|(|<wide|v|^><rsub|1>|)>=1>
  and <math|<wide|\<psi\>|^><rsub|3><around*|(|<wide|\<varphi\>|^><rsub|3>|)>=\<partial\><rsub|2><wide|\<varphi\>|^><rsub|3><around*|(|<wide|v|^><rsub|1>|)>=1>
  but <math|<wide|\<psi\>|^><rsub|1><around*|(|<wide|\<varphi\>|^><rsub|2>|)>=<wide|\<varphi\>|^><rsub|2><around*|(|<wide|v|^><rsub|1>|)>=0>
  and so on.

  <subsection|The local basis>

  The reference shape functions are then mapped onto a <dfn|local basis over
  the physical cell> <math|K> by means of a mapping

  <\equation*>
    <wide|\<varphi\>|^>\<mapsto\>\<varphi\>=\<cal-F\><around*|(|<wide|\<varphi\>|^>|)>.
  </equation*>

  For standard Lagrange elements (i.e. when <math|\<cal-B\><rprime|'>>
  consists only of point evaluation forms), the immediate choice
  <math|\<cal-F\><around*|(|<wide|\<varphi\>|^>|)>\<assign\><wide|\<varphi\>|^>\<circ\><wide|\<b-x\>|^>>
  maps the basis <math|\<cal-B\>> into a basis of
  <math|P<rsub|3><around*|(|K|)>> which fulfills the <dfn|delta property>
  <eqref|eq:delta-property>. However, in order for <eqref|eq:delta-property>
  to hold in <math|K> (the linear functionals remaining unchanged) for the
  images of the Hermite pairs of dofs <math|<around*|(|<wide|\<varphi\>|^><rsub|2><rsup|\<alpha\>>,<wide|\<varphi\>|^><rsub|3><rsup|\<alpha\>>|)>>
  we require a different transformation which at vertex
  <math|v<rsup|\<alpha\>>> of <math|K>, sets the shape functions to a linear
  combination of all the shape functions at
  <math|<wide|v|^><rsup|\<alpha\>>>.

  <\definition*>
    Let

    <\equation*>
      \<Phi\><rsub|\<alpha\>>\<assign\><around*|(|\<varphi\><rsub|1><rsup|\<alpha\>>,\<varphi\><rsub|2><rsup|\<alpha\>>,\<varphi\><rsub|3><rsup|\<alpha\>>|)><rsup|\<top\>>,<space|2em><wide|\<Phi\>|^><rsub|\<alpha\>>\<assign\><around*|(|<wide|\<varphi\>|^><rsub|1><rsup|\<alpha\>>,<wide|\<varphi\>|^><rsub|2><rsup|\<alpha\>>,<wide|\<varphi\>|^><rsub|3><rsup|\<alpha\>>|)><rsup|\<top\>>,
    </equation*>

    and define

    <\equation*>
      H\<assign\><matrix|<tformat|<table|<row|<cell|1>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|\<b-x\><rsub|1,1>>|<cell|\<b-x\><rsub|1,2>>>|<row|<cell|0>|<cell|\<b-x\><rsub|2,1>>|<cell|\<b-x\><rsub|2,2>>>>>>.
    </equation*>

    The <dfn|Hermite transformation> <math|\<cal-F\>> is given by the matrix
    <math|H> as:

    <\equation>
      <label|eq:hermite-transform>\<Phi\><rsub|\<alpha\>><around*|(|x|)>=\<cal-F\><around*|(|<wide|\<Phi\>|^><rsub|\<alpha\>><around*|(|<wide|\<b-x\>|^><around*|(|x|)>|)>|)>\<assign\>H<around*|(|<wide|\<b-x\>|^><around*|(|x|)>|)>*<wide|\<Phi\>|^><rsub|\<alpha\>><around*|(|<wide|\<b-x\>|^><around*|(|x|)>|)>.
    </equation>

    This transformation respects property
    <eqref|eq:delta-property>.<\footnote>
      This can be seen by solving the linear system resulting from imposing
      condition <eqref|eq:delta-property> for a linear combination of the
      relevant shape functions <cite-detail|solin_partial_2005|Ÿ6.4.3>. We
      can easily check that it holds in the physical domain: for the point
      evaluations this is clear and for the others it is a matter of applying
      the chain rule: <math|\<psi\><rsub|2><around*|(|\<varphi\><rsub|2>|)>=\<partial\><rsub|1>\<varphi\><rsub|2><around*|(|v<rsub|1>|)>=\<b-x\><rsub|1,k><around*|(|<wide|v|^><rsub|1>|)>*<wide|\<varphi\>|^><rsup|\<alpha\>><rsub|k,l><around*|(|<wide|v|^><rsub|1>|)>*<wide|\<b-x\>|^><rsub|l,1><around*|(|v<rsub|1>|)>=\<b-x\><rsub|1,k><around*|(|<wide|v|^><rsub|1>|)>*\<delta\><rsub|k\<nocomma\>l>*<wide|\<b-x\>|^><rsub|l,1><around*|(|v<rsub|1>|)>=1>,
      where <math|k,l\<in\><around*|{|1,2|}>>.
    </footnote>
  </definition*>

  If we look at the entries related to partial derivatives we have:

  <\equation*>
    \<varphi\><rsub|i+1><rsup|\<alpha\>><around*|(|x|)>=<around*|(|\<b-x\><rsub|i,j>*<wide|\<varphi\>|^><rsup|\<alpha\>><rsub|j>|)><around*|(|<wide|\<b-x\>|^><around*|(|x|)>|)>.
  </equation*>

  Using the chain rule and the fact that <math|D H=0>, we obtain

  <\equation*>
    D \<Phi\><rsub|\<alpha\>><around*|(|x|)>=<around*|(|H*D
    <wide|\<Phi\>|^><rsub|\<alpha\>>|)><around*|(|<wide|\<b-x\>|^><around*|(|x|)>|)>*D<wide|\<b-x\>|^><around*|(|x|)>,
  </equation*>

  where the last term <math|D<wide|\<b-x\>|^>> comes from the affine change
  of coordinates and is also present in the Lagrange case. Component-wise:

  <\equation>
    <label|eq:hermite-first-derivatives>\<varphi\><rsub|i,j><rsup|\<alpha\>><around*|(|x|)>=<around*|(|H<rsub|i\<nocomma\>k>*<wide|\<varphi\>|^><rsup|\<alpha\>><rsub|k,l>|)><rsub|\|<wide|\<b-x\>|^><around*|(|x|)>>*<wide|\<b-x\>|^><rsub|l,j><around*|(|x|)>,<space|2em>i\<in\><around*|{|1,2,3|}>.
  </equation>

  If we look again at the matrix product, each partial derivative
  <math|\<varphi\><rsub|i,j><rsup|\<alpha\>>>,
  <math|j\<in\><around*|{|1,2|}>> is a linear combination of the derivatives
  of both reference basis functions <math|\<nabla\><wide|\<varphi\>|^><rsub|k><rsup|\<alpha\>>>
  with weights given by <math|\<nabla\> \<b-x\><rsub|i>>:

  <\equation*>
    \<varphi\><rsub|i+1,j><rsup|\<alpha\>><around*|(|x|)>=<around*|(|\<b-x\><rsub|i,k>*<wide|\<varphi\>|^><rsup|\<alpha\>><rsub|k,l>|)><rsub|\|<wide|\<b-x\>|^><around*|(|x|)>>*<wide|\<b-x\>|^><rsub|l,j><around*|(|x|)>,<space|2em>i\<in\><around*|{|1,2|}>.
  </equation*>

  For the second derivatives we have a similar expression, again thanks to
  <math|D<rsup|2> \<b-x\>=0>:

  <\equation*>
    \<varphi\><rsub|i+1,j\<nocomma\>k><rsup|\<alpha\>>=\<b-x\><rsub|i,p>*<wide|\<varphi\>|^><rsup|\<alpha\>><rsub|p,q\<nocomma\>r>*<wide|\<b-x\>|^><rsub|q,j>*<wide|\<b-x\>|^><rsub|r,k>,<space|2em>i\<in\><around*|{|1,2|}>,
  </equation*>

  i.e.

  <\equation*>
    D<rsup|2> \<varphi\><rsup|\<alpha\>><rsub|i+1>=D<rsup|\<top\>><wide|\<b-x\>|^>*<around*|(|\<b-x\><rsub|i,1>*D<rsup|2><wide|\<varphi\>|^><rsub|2><rsup|\<alpha\>>+\<b-x\><rsub|i,2>*D<rsup|2><wide|\<varphi\>|^><rsub|3><rsup|\<alpha\>>|)>*D
    <wide|\<b-x\>|^>,<space|2em>i\<in\><around*|{|1,2|}>.
  </equation*>

  And if we use <math|H> (the sums are over
  <math|p\<in\><around*|{|1,2,3|}>,q,r\<in\><around*|{|1,2|}>>):

  <\equation>
    <label|eq:hermite-second-derivatives>\<varphi\><rsub|i,j\<nocomma\>k><rsup|\<alpha\>>=H<rsub|i\<nocomma\>p>*<wide|\<varphi\>|^><rsup|\<alpha\>><rsub|p,q\<nocomma\>r>*<wide|\<b-x\>|^><rsub|q,j>*<wide|\<b-x\>|^><rsub|r,k>.
  </equation>

  Here as before, we obtain a linear combination of the second derivatives of
  both reference basis functions <math|D<rsup|2>
  <wide|\<varphi\>|^><rsup|\<alpha\>><rsub|j>> with weights given by
  <math|\<nabla\> \<b-x\><rsub|i>>. The terms
  <math|D<rsup|\<top\>><wide|\<b-x\>|^>> and <math|D<wide|\<b-x\>|^>> are
  once more the same as if we had only used the inverse mapping
  <math|<wide|\<b-x\>|^>> to change variables.

  With <math|3> derivatives <todo|generalise to <math|k> derivatives with a
  multiindex>:

  <\equation*>
    \<varphi\><rsub|i,j\<nocomma\>k\<nocomma\>l><rsup|\<alpha\>>=H<rsub|i\<nocomma\>p>*<wide|\<varphi\>|^><rsup|\<alpha\>><rsub|p,q\<nocomma\>r\<nocomma\>s>*<wide|\<b-x\>|^><rsub|q,j>*<wide|\<b-x\>|^><rsub|r,k>*<wide|\<b-x\>|^><rsub|s,l>.
  </equation*>

  <section|Assembly of the stiffness matrix>

  We now review the steps required to piece together the linear system
  discretising <eqref|eq:linear-problem>.

  <subsection|The global basis>

  After defining the local shape functions
  <math|\<varphi\><rsup|\<alpha\>><rsub|i>> for each simplex <math|K> of the
  triangulation, we gather all of them into a basis for a polynomial space
  <math|V<rsub|h>> with the help of a <dfn|local-to-global map>
  <math|\<iota\>>. This assigns local indices to a global set, matching
  vertices, etc. <todo|elaborate?>

  <subsection|The stiffness matrix>

  In order to discretise <eqref|eq:linear-problem> we require the action of
  <math|a> on the polynomial space <math|V<rsub|h>=span<around*|{|\<varphi\><rsub|1>,\<ldots\>,\<varphi\><rsub|10*N>|}>>:

  <\equation*>
    a:V<rsup|<around*|(|1|)>><rsub|h>\<times\>V<rsup|<around*|(|2|)>><rsub|h>\<rightarrow\>\<bbb-R\>,
  </equation*>

  where we used superindices to tell the two copies of <math|V<rsub|h>>
  apart. This can be split into the contribution of each element, meaning
  that we compute the values

  <\equation*>
    A<rsub|i>=a<around*|(|\<varphi\><rsub|i><rsup|<around*|(|1|)>>,\<varphi\><rsup|<around*|(|2|)>><rsub|j>|)>=<big|sum><rsub|K\<in\>\<cal-T\>>a<rsup|K><around*|(|\<varphi\><rsub|i><rsup|<around*|(|1|)>>,\<varphi\><rsup|<around*|(|2|)>><rsub|j>|)>
  </equation*>

  with \ <math|0\<leqslant\>i,j\<less\>10*N> and <math|a<rsup|K>> the
  contribution from element <math|K>. Because of the compact supports of the
  basis functions, most terms in the sum will be zero. <todo|This results in
  a global element matrix <math|A> and a system <math|A*u=b>...>

  <subsection|Computing integrals>

  The contribution <math|a<rsup|K>> from element <math|K> to the stiffness
  matrix will typically be of integral form. For example, if we are
  discretising the biharmonic equation <math|\<Delta\><rsup|2> u=f>, assuming
  homogeneous boundary conditions, the corresponding bilinear form is
  <math|a<around*|(|u,v|)>=<big|int>\<Delta\>u*\<Delta\>v>, so we need to
  compute

  <\equation*>
    a<rsup|K><around*|(|\<varphi\><rsub|i><rsup|<around*|(|1|)>>,\<varphi\><rsup|<around*|(|2|)>><rsub|j>|)>=<big|int><rsub|K>\<Delta\>\<varphi\><rsub|i><rsup|<around*|(|1|)>>*\<mathLaplace\>\<varphi\><rsup|<around*|(|2|)>><rsub|j>*\<mathd\>x=<big|int><rsub|<wide|K|^>><todo|\<ldots\>>.
  </equation*>

  <todo|As can be seen, the transformation introduced above plays here a role
  as well...> Notice that even though <math|V<rsub|h>\<nsubset\>H<rsup|2>> we
  can still compute the Laplacian in the interior of each element <math|K>.
  However, the normal derivative might present discontinuities, <todo|so...>

  We implement a quadrature representation (as opposed to a tensor
  representation), as in <cite|olgaard_optimizations_2010>.

  \;

  <section|Interpolation>

  <\note*>
    See <cite-detail|brenner_mathematical_2008|Ÿ4.8> and
    <cite|girault_hermite_2002> for info on interpolation of non-smooth
    functions.
  </note*>

  We have at least the following options for interpolating functions
  <math|g\<in\>H<rsup|2><around*|(|\<Omega\>|)>> into <math|V<rsub|h>>, from
  slowest and most accurate to fastest and less accurate
  <cite-detail|solin_partial_2005|Ÿ6.3.8>.

  <\enumerate>
    <item>Use a global orthogonal projection to compute the best interpolant
    in the <math|H<rsup|2>>-sense. This means solving the system: <todo|...>.
    Does FEniCS do this by default with <python|project()>?

    <item>Use nodal interpolation of vertex and derivatives plus local
    orthogonal projections in the element interiors. <todo|...>

    <item>Compute the nodal interpolant.
  </enumerate>

  <subsection|The local nodal interpolant>

  The local nodal interpolant is given by

  <\equation*>
    I<rsub|K><around*|(|g|)>=<big|sum><rsub|\<alpha\>=1><rsup|3>\<psi\><rsup|\<alpha\>><rsub|i><around*|(|g|)>*\<varphi\><rsub|i>
  </equation*>

  where <math|\<psi\><rsup|\<alpha\>><rsub|i><around*|(|g|)>=\<partial\><rsub|i-1>
  g<around*|(|v<rsup|\<alpha\>>|)>> and <math|\<partial\><rsub|0>f=f>.

  <\question*>
    How should one compute the derivatives? Should we extend the interface of
    <cpp|ufc::function> to include differentiation? Can we use AD?
  </question*>

  <subsection|<math|W<rsup|m,2>>-projection>

  See <cite-detail|solin_partial_2005|Ÿ6.3.8>,
  <cite|brenner_mathematical_2008>...

  <section|Implementation in <name|FEniCS> >

  An incomplete list of things that needed doing:

  <\enumerate>
    <item>Fixing the FIAT Hermite element definition (e.g. it was returning
    an inconsistent format of <tt|entity_ids>) in dimensions 1 to 3.\ 

    <item>Implementing in FFC the \PHermite transformation\Q described above
    for the element basis functions which are associated to the evaluation of
    partial derivatives at the nodes of the simplex. This means:

    <\itemize-dot>
      <item><tt|representation.py>: for every dof
      <math|<wide|\<psi\>|^><rsub|i><rsup|\<alpha\>>> include the
      coefficients of the other dofs <math|<wide|\<psi\>|^><rsub|j><rsup|\<alpha\>>>
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

    <item>Implementing the same transformation for the quadratures. See
    <python|create_argument()> and <python|create_function()> in
    <tt|ffc/quadratures/quadraturetransformerbase.py>.

    <item><todo|Implementing the nodal interpolant>.
  </enumerate>

  <section|The Euler-Bernoulli beam model>

  This classical model reduces a 2D beam to a 1D problem. We study the
  deformation of the midplane of a beam subject to transveral load, under the
  assumptions that after the deformation the normals to the midplane:

  <\itemize-dot>
    <item>do not bend,

    <item>do not stretch,

    <item>remain orthogonal to the midplane.
  </itemize-dot>

  This theory is adequate for thin beams and is intended for small strains
  even with large global deformations: it is a <em|physically linear> but
  <em|geometrically non-linear> theory. For thicker beams, Timoshenko's
  theory, which accounts for internal shear forces, yields more accurate
  predictions <inactive|<cite|>>.

  <subsection|Derivation>

  Fix <math|\<omega\>=<around*|(|a,b|)>> to be the midplane of the beam
  represented by the domain <math|\<Omega\>=\<omega\>\<times\><around*|(|-h/2,h/2|)>\<subset\>\<bbb-R\><rsup|2>>.
  Conservation of momentum and certain linear constitutive relations yield
  <cite-detail|solin_partial_2005|Ÿ6.1.1> the equation

  <\equation>
    <label|eq:euler-bernoulli><frac|\<mathd\><rsup|2>|\<mathd\>x<rsup|2>>
    <around*|(|b<around|(|x|)>*<frac|\<mathd\><rsup|2>|\<mathd\>x<rsup|2>>
    u<around|(|x|)>|)>=f<around*|(|x|)>,
  </equation>

  where <math|b<around|(|x|)>=E<around|(|x|)>*I<around|(|x|)>> is the product
  of <dfn|Young's modulus> <math|E> and the <dfn|area moment of inertia> of
  the beam <math|I>.

  For steel (with \<less\> 0.3% carbon) at 21\<degree\>C, the modulus is
  <math|E=203.4\<cdot\>10<rsup|9>> Pa and since our beam has a constant
  square cross section of side 0.01m, we have <math|I=8\<cdot\>10<rsup|-10>>.
  <todo|compute...>

  <subsection|Weak formulation>

  For consistency with notation for higher dimensions, write
  <math|\<nabla\>u=<frac|\<mathd\>|\<mathd\>x> u> and
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
  We fix either the <dfn|deflections> <math|u<around|(|\<alpha\>|)>> at
  <math|\<alpha\>\<in\>\<gamma\><rsub|u>> or the <dfn|slopes>
  <math|u<rprime|'><around|(|\<alpha\>|)>> at
  <math|\<alpha\>\<in\>\<gamma\><rsub|\<theta\>>> or both. These conditions
  are incorporated into the definition of <math|V>. For example if we
  <em|clamp> the beam at an horizontal position we have

  <\equation*>
    V=V<rsub|<text|clamped>>=<around*|{|u\<in\>H<rsup|2><around*|(|\<omega\>|)>:u<around*|(|\<alpha\>|)>=u<rsub|\<alpha\>>,u<rprime|'><around*|(|\<alpha\>|)>=u<rprime|'><rsub|\<alpha\>>,\<alpha\>=1,2|}>.
  </equation*>

  <paragraph|Natural boundary conditions>For
  <math|\<alpha\>\<in\>\<gamma\><rsub|M>\<subset\><around|{|a,b|}>>, we can
  fix the <dfn|bending moment>:

  <\equation*>
    M<around*|(|\<alpha\>|)>=<around|(|b*\<Delta\>u|)><around|(|\<alpha\>|)>,
  </equation*>

  which for general <math|x\<in\>\<omega\>> is the torque exerted by forces
  surrounding <math|x>. If for example we find solutions such that
  <math|M<around*|(|a|)>=0>, then we are assuming that the left end of the
  beam is free to rotate, i.e. that it undergoes no bending due to torque.
  Alternatively we can set the <dfn|shear force> at
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
    <\bib-list|6>
      <bibitem*|ABH+15><label|bib-alnaes_fenics_2015>Martin<nbsp>S.<nbsp>Alnaes,
      Jan Blechta, Johan Hake, August Johansson, Benjamin Kehlet, Anders
      Logg, Chris Richardson, Johannes Ring,
      Marie<nbsp>E.<nbsp>Rognes<localize|, and
      >Garth<nbsp>N.<nbsp>Wells.<newblock> The FEniCS Project Version
      1.5.<newblock> <with|font-shape|italic|Archive of Numerical Software>,
      3(100), 2015.<newblock>

      <bibitem*|BS08><label|bib-brenner_mathematical_2008>Susanne<nbsp>C.<nbsp>Brenner<localize|
      and >L.<nbsp>Ridgway Scott.<newblock> <with|font-shape|italic|The
      mathematical theory of finite element methods>.<newblock>
      <localize|Number><nbsp>15<localize| in >Texts in Applied Mathematics.
      Springer New York, New York, NY, 3<localize| edition>, 2008.<newblock>

      <bibitem*|GS02><label|bib-girault_hermite_2002>V.<nbsp>Girault<localize|
      and >L.<nbsp>Scott.<newblock> Hermite interpolation of nonsmooth
      functions preserving boundary conditions.<newblock>
      <with|font-shape|italic|Mathematics of Computation>,
      71(239):1043\U1074, 2002.<newblock>

      <bibitem*|LMW12><label|bib-logg_automated_2012>Anders Logg, Kent-Andre
      Mardal<localize|, and >Garth<nbsp>N.<nbsp>Wells<localize|,
      editors>.<newblock> <with|font-shape|italic|Automated solution of
      differential equations by the finite element method>.<newblock>
      <localize|Number><nbsp>84<localize| in >Lecture notes in computational
      science and engineering. Springer, 2012.<newblock> DOI
      10.1007/978-3-642-23099-8.<newblock>

      <bibitem*|Sol05><label|bib-solin_partial_2005>Pavel Solin.<newblock>
      <with|font-shape|italic|Partial differential equations and the finite
      element method>.<newblock> Pure and applied Mathematics. Dec
      2005.<newblock>

      <bibitem*|ØW10><label|bib-olgaard_optimizations_2010>Kristian<nbsp>B.<nbsp>Ølgaard<localize|
      and >Garth<nbsp>N.<nbsp>Wells.<newblock> Optimizations for Quadrature
      Representations of Finite Element Tensors Through Automated Code
      Generation.<newblock> <with|font-shape|italic|ACM Trans. Math. Softw.>,
      37(1):8\U1, jan 2010.<newblock>
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
    <associate|auto-10|<tuple|4|5>>
    <associate|auto-11|<tuple|4.1|5>>
    <associate|auto-12|<tuple|4.2|5>>
    <associate|auto-13|<tuple|5|5>>
    <associate|auto-14|<tuple|6|5>>
    <associate|auto-15|<tuple|6.1|6>>
    <associate|auto-16|<tuple|6.2|6>>
    <associate|auto-17|<tuple|1|6>>
    <associate|auto-18|<tuple|2|7>>
    <associate|auto-19|<tuple|6.3|7>>
    <associate|auto-2|<tuple|2|2>>
    <associate|auto-20|<tuple|6.4|7>>
    <associate|auto-21|<tuple|6.4|7>>
    <associate|auto-3|<tuple|2.1|2>>
    <associate|auto-4|<tuple|1|2>>
    <associate|auto-5|<tuple|2.2|3>>
    <associate|auto-6|<tuple|3|3>>
    <associate|auto-7|<tuple|3.1|4>>
    <associate|auto-8|<tuple|3.2|4>>
    <associate|auto-9|<tuple|3.3|4>>
    <associate|bib-alnaes_fenics_2015|<tuple|ABH+15|7>>
    <associate|bib-brenner_mathematical_2008|<tuple|BS08|7>>
    <associate|bib-girault_hermite_2002|<tuple|GS02|7>>
    <associate|bib-logg_automated_2012|<tuple|LMW12|7>>
    <associate|bib-olgaard_optimizations_2010|<tuple|ØW10|?>>
    <associate|bib-solin_partial_2005|<tuple|Sol05|7>>
    <associate|eq:delta-property|<tuple|2|3>>
    <associate|eq:euler-bernoulli|<tuple|6|6>>
    <associate|eq:hermite-first-derivatives|<tuple|4|4>>
    <associate|eq:hermite-second-derivatives|<tuple|5|4>>
    <associate|eq:hermite-transform|<tuple|3|3>>
    <associate|eq:linear-problem|<tuple|1|1>>
    <associate|fig:reference-triangle|<tuple|1|3>>
    <associate|footnote-1|<tuple|1|1>>
    <associate|footnote-2|<tuple|2|2>>
    <associate|footnote-3|<tuple|3|2>>
    <associate|footnote-4|<tuple|4|3>>
    <associate|footnr-1|<tuple|1|1>>
    <associate|footnr-2|<tuple|2|2>>
    <associate|footnr-3|<tuple|3|2>>
    <associate|footnr-4|<tuple|4|3>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      alnaes_fenics_2015

      logg_automated_2012

      brenner_mathematical_2008

      logg_automated_2012

      olgaard_optimizations_2010

      brenner_mathematical_2008

      brenner_mathematical_2008

      solin_partial_2005

      brenner_mathematical_2008

      girault_hermite_2002

      solin_partial_2005

      solin_partial_2005

      brenner_mathematical_2008

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

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2.<space|2spc>Cubic
      Hermite elements in 2D> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>

      <with|par-left|<quote|1tab>|2.1.<space|2spc>The reference simplex and
      shape functions <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|1tab>|2.2.<space|2spc>The local basis
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3.<space|2spc>Assembly
      of the stiffness matrix> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6><vspace|0.5fn>

      <with|par-left|<quote|1tab>|3.1.<space|2spc>The global basis
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>

      <with|par-left|<quote|1tab>|3.2.<space|2spc>The stiffness matrix
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8>>

      <with|par-left|<quote|1tab>|3.3.<space|2spc>Computing integrals
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|4.<space|2spc>Interpolation>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10><vspace|0.5fn>

      <with|par-left|<quote|1tab>|4.1.<space|2spc>The local nodal interpolant
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <with|par-left|<quote|1tab>|4.2.<space|2spc><with|mode|<quote|math>|W<rsup|m,2>>-projection
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|5.<space|2spc>Implementation
      in <with|font-shape|<quote|small-caps>|FEniCS> >
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|6.<space|2spc>The
      Euler-Bernoulli beam model> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-14><vspace|0.5fn>

      <with|par-left|<quote|1tab>|6.1.<space|2spc>Derivation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15>>

      <with|par-left|<quote|1tab>|6.2.<space|2spc>Weak formulation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-16>>

      <with|par-left|<quote|4tab>|Essential boundary conditions
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-17><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Natural boundary conditions
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-18><vspace|0.15fn>>

      <with|par-left|<quote|1tab>|6.3.<space|2spc>Existence and uniqueness
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-19>>

      <with|par-left|<quote|1tab>|6.4.<space|2spc>Discretization
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-20>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-21><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>