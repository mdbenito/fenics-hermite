<TeXmacs|1.99.5>

<style|<tuple|amsart|british|better-amsart>>

<\body>
  <\hide-preamble>
    \;

    <assign|python-code|<macro|body|<\small>
      <\pseudo-code>
        <python|<arg|body>>
      </pseudo-code>
    </small>>>
  </hide-preamble>

  <doc-data|<doc-title|Hermite elements:<new-line>introduction and
  implementation>|<doc-running-title|Hermite
  elements>|<doc-author|<author-data|<author-name|Miguel de Benito
  Delgado>|<\author-affiliation>
    Universitšt Augsburg
  </author-affiliation>>>|<doc-date|March 2017>>

  <abstract-data|<abstract|We briefly review the setting and notation for
  finite element approximations using cubic Hermite elements and the
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

      <with|par-left|1tab|2.2.<space|2spc>The local basis over a physical
      element <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|1tab|2.3.<space|2spc>The global basis for
      <with|mode|math|V<rsub|h>> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>

      <vspace*|1fn><with|font-series|bold|math-font-series|bold|3.<space|2spc>Construction
      of the linear system> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8><vspace|0.5fn>

      <with|par-left|1tab|3.1.<space|2spc>Assembly of the stiffness matrix
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9>>

      <with|par-left|1tab|3.2.<space|2spc>Imposing boundary conditions
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10>>

      <with|par-left|1tab|3.3.<space|2spc>Computing integrals
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <vspace*|1fn><with|font-series|bold|math-font-series|bold|4.<space|2spc>Interpolation>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12><vspace|0.5fn>

      <with|par-left|1tab|4.1.<space|2spc>The local nodal interpolant
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13>>

      <with|par-left|1tab|4.2.<space|2spc><with|mode|math|H<rsup|2>>-projection
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-14>>

      <vspace*|1fn><with|font-series|bold|math-font-series|bold|5.<space|2spc>The
      Euler-Bernoulli beam model> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15><vspace|0.5fn>

      <with|par-left|1tab|5.1.<space|2spc>Derivation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-16>>

      <with|par-left|1tab|5.2.<space|2spc>Weak formulation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-17>>

      <with|par-left|4tab|Essential boundary conditions
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-18><vspace|0.15fn>>

      <with|par-left|4tab|Natural boundary conditions
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-19><vspace|0.15fn>>

      <with|par-left|1tab|5.3.<space|2spc>Existence and uniqueness
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-20>>

      <with|par-left|1tab|5.4.<space|2spc>Discretization
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-21>>

      <with|par-left|1tab|5.5.<space|2spc>Clamped beam
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-22>>

      <with|par-left|1tab|5.6.<space|2spc>Cantilevered beam
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-24>>

      <with|par-left|1tab|5.7.<space|2spc>Simply supported beam
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-26>>

      <vspace*|1fn><with|font-series|bold|math-font-series|bold|6.<space|2spc>The
      biharmonic equation> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-28><vspace|0.5fn>

      <vspace*|1fn><with|font-series|bold|math-font-series|bold|Appendix
      A.<space|2spc>Implementation details>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-32><vspace|0.5fn>

      <with|par-left|1tab|A.1.<space|2spc><with|font-shape|small-caps|FIAT>
      implementation <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-33>>

      <with|par-left|1tab|A.2.<space|2spc>Computing the Hermite shape
      functions <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-34>>

      <with|par-left|1tab|A.3.<space|2spc>Essential boundary conditions for
      Hermite elements <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-35>>

      <with|par-left|1tab|A.4.<space|2spc>Automatic differentiation for nodal
      interpolation <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-36>>

      <vspace*|1fn><with|font-series|bold|math-font-series|bold|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-37><vspace|0.5fn>
    </table-of-contents>
  </small>

  <section|Setting>

  We fix a polygonal domain <math|\<Omega\>\<subset\>\<bbb-R\><rsup|2>> (the
  <dfn|physical domain>) and a sufficiently regular triangulation<\footnote>
    <em|\SA triangulation of a polygonal domain \<Omega\> is a subdivision
    consisting of triangles having the property that no vertex of any
    triangle lies in the interior of an edge of another triangle.\T>
    <cite-detail|brenner_mathematical_2008|ü3.3.11>.
  </footnote> <math|\<cal-T\>=<around*|{|K<rsub|1>,\<ldots\>,K<rsub|N>|}>>.
  Let <math|V> be a (not necessarily proper) subspace of
  <math|H<rsup|2><around*|(|\<Omega\>|)>> and
  <math|L:V\<rightarrow\>\<bbb-R\>> and <math|a:V\<times\>V\<rightarrow\>\<bbb-R\>>
  respectively linear and bilinear forms. We consider the linear variational
  problem:

  <\em>
    <\with|par-par-sep|1fn|par-sep|0fn>
      Find <math|u\<in\>V> such that

      <\equation>
        <label|eq:linear-problem>a<around*|(|v,u|)>=L<around*|(|v|)><text|<space|1em>for
        every >v\<in\>V.
      </equation>
    </with>
  </em>

  Assuming coerciveness of <math|a>, this problem has a solution by the Lemma
  of Lax-Milgram. The finite element method builds a suitable discretisation
  <math|V<rsub|h>> of the space <math|V> and projects Problem
  <eqref|eq:linear-problem> onto it for its numerical solution as a linear
  system: More precisely, letting <math|<around*|{|\<varphi\><rsub|i>|}>> be
  a basis for <math|V<rsub|h>>, a function
  <math|u<rsub|h>=u<rsub|i>*\<varphi\><rsub|i>\<in\>V<rsub|h>> solves
  <eqref|eq:linear-problem><\footnote>
    In this note we use implicit summation over repeated indices.
  </footnote> iff

  <\equation*>
    a<around*|(|u<rsub|h>,\<varphi\><rsub|j>|)>=u<rsub|i>*a<around*|(|\<varphi\><rsub|i>,\<varphi\><rsub|j>|)>=L<around*|(|\<varphi\><rsub|j>|)><text|,
    for every >\<varphi\><rsub|j>,
  </equation*>

  or, more succintly iff

  <\equation>
    <label|eq:linear-system>A*u<rsub|h>=b,
  </equation>

  where <math|A<rsub|j\<nocomma\>i>=a<around*|(|\<varphi\><rsub|i>,\<varphi\><rsub|j>|)>>
  and <math|b<rsub|j>=L<around*|(|\<varphi\><rsub|j>|)>>.<\footnote>
    Note the transposition: this can be a source of confusion when
    implementing the discrete system.
  </footnote> If one proves that <math|<around*|\<\|\|\>|u<rsub|h>-u|\<\|\|\>>\<rightarrow\>0>
  as <math|h\<rightarrow\>0> in some norm, then \Pall\Q that is left is to
  actually construct a suitable <math|V<rsub|h>> and solve the linear system
  <eqref|eq:linear-system>.<\footnote>
    The error between <math|u<rsub|h>> and <math|u> will be bounded above by
    some constant multiple of the \Pbest approximation error\Q (see e.g.
    Cea's lemma). One typically bounds the latter by the error made by an
    interpolation operator, which in turn is bounded above by the norm of the
    solution times some power of the grid size (cf.
    <cite-detail|grossmann_numerical_2007|ü4.4 and ü3.3> for the conforming
    case).
  </footnote> In broad terms, the standard procedure to build
  <math|V<rsub|h>> is to construct a polynomial space involving:

  <\itemize-dot>
    <item>The triangulation <math|\<cal-T\>>.

    <item>A reference simplex <math|<wide|K|^>\<subset\>\<bbb-R\><rsup|d>>
    with an associated set <math|<wide|\<cal-L\>|^>> of linear forms
    (<em|degrees of freedom> or <em|dofs>) over
    <math|P<rsub|p><around*|(|<wide|K|^>|)>>, the space of polynomials of
    degree up to <math|p>, defining a (dual) local polynomial basis
    <math|<wide|\<cal-V\>|^>> of <math|P<rsub|p><around*|(|<wide|K|^>|)>>
    (<em|shape functions>).<\footnote>
      In <name|FIAT>, the polynomials of <math|<wide|V|^>> are themselves
      expressed in a fixed polynomial basis with good numerical properties.
      In particular the one picked (<todo|Dubiner?>) yields much better
      condition numbers for the mass matrices than e.g. the basis of
      monomials <cite-detail|logg_automated_2012|ü4.3.2>.
    </footnote>

    <item>A set of invertible affine maps
    <math|<around*|{|\<b-x\><rsup|K>|}><rsub|K\<in\>\<cal-T\>>> between
    <math|<wide|K|^>> and each simplex <math|K\<in\>\<cal-T\>>.

    <item>A transformation <math|\<cal-F\><rsub|K>:<wide|\<cal-V\>|^><around*|(|P<rsub|p><around*|(|<wide|K|^>|)>|)>\<rightarrow\>\<cal-V\><around*|(|P<rsub|p><around*|(|K|)>|)>>
    preserving the duality property of <math|<wide|\<cal-V\>|^>> and
    <math|<wide|\<cal-L\>|^>>.

    <item>A <em|local-to-global dof mapping> <math|\<iota\><rsub|K>>,
    assigning local dofs in each <math|K\<in\>\<cal-T\>> to global dofs in
    <math|V<rsub|h>>.
  </itemize-dot>

  Furthermore one requires methods to:

  <\itemize-dot>
    <item>Interpolate functions <math|g\<in\>V> into <math|V<rsub|h>>, e.g.
    <math|H<rsup|1>>-projection.

    <item>Compute integrals of linear forms defined over <math|V<rsub|h>>.

    <item>Set boundary conditions.
  </itemize-dot>

  In this note we work in <math|d\<in\><around*|{|1,2|}>> dimensions, with
  polynomials of degree <math|p=3> and a particular choice of degrees of
  freedom giving rise to the classic <em|Hermite finite element>. This
  element is a generalization of cubic Hermite interpolating polynomials in 1
  dimension which has been in the literature since at least
  <inactive|<cite|ciarlet_raviart>>. The dofs are given by function and
  partial derivative evaluation at the nodes of the triangulation (plus
  evaluation at the barycenters in 2D), which results in <math|C<rsup|1>>
  function spaces in 1D. In 2D, since continuity of partial derivatives is
  only ensured at the nodes but not at the edges, the resulting function
  spaces are not in <math|C<rsup|1>> <cite-detail|brenner_mathematical_2008|Proposition
  3.3.17>, hence not in <math|H<rsup|2><around*|(|\<Omega\>|)>>. This will be
  relevant for some <math|4<rsup|th>> order problems of interest below, since
  the discretisation won't be conforming as described above.

  <todo|Approximations and errors for non-conforming discretizations...>

  <section|Cubic Hermite elements in 2D>

  <\notation*>
    In this section <math|\<alpha\>,i,p\<in\><around*|{|1,2,3|}>> and
    <math|j,k,m,q,r,s\<in\><around*|{|1,2|}>>.
  </notation*>

  <subsection|The reference simplex and shape functions>

  Let <math|<wide|K|^>> be a <dfn|reference simplex> in 2D with vertices
  <math|<wide|v|^><rsup|\<alpha\>>> and <math|\<b-x\>:<wide|K|^>\<rightarrow\>K>
  be a non-degenerate affine transformation into one of the cells of the
  physical domain <math|\<Omega\>> with inverse mapping
  <math|<wide|\<b-x\>|^>:K\<rightarrow\><wide|K|^>> (Figure
  <reference|fig:reference-triangle>).

  <big-figure|<with|gr-mode|<tuple|edit|point>|gr-frame|<tuple|scale|1cm|<tuple|0.270004gw|0.180009gh>>|gr-geometry|<tuple|geometry|0.280007par|0.200003par|center>|gr-grid|<tuple|empty>|gr-grid-old|<tuple|cartesian|<point|0|0>|1>|gr-edit-grid-aspect|<tuple|<tuple|axes|none>|<tuple|1|none>|<tuple|10|none>>|gr-edit-grid|<tuple|empty>|gr-edit-grid-old|<tuple|cartesian|<point|0|0>|1>|<graphics||<line|<point|0|0>|<point|2.0|0.0>|<point|0.0|2.0>|<point|0.0|0.0>>|<point|0|0>|<point|0|0>|<point|2|0>|<point|0|2>|<point|0.7|0.6>>>|<todo|<label|fig:reference-triangle>The
  reference triangle and the mapping <math|\<b-x\>>...>>

  Let <math|P<rsub|3>=P<rsub|3><around*|(|<wide|K|^>|)>> be the space of
  polynomials of order up to 3 defined over <math|<wide|K|^>>. We define 10
  linear forms (the <dfn|degrees of freedom> or <dfn|dofs>) in the dual
  <math|P<rsub|3><around*|(|<wide|K|^>|)><rprime|'>> as follows: for each
  vertex <math|<wide|v|^><rsup|\<alpha\>>> in <math|<wide|K|^>> add one point
  evaluation

  <\equation*>
    l<rsup|\<alpha\>><rsub|1><around*|(|f|)>\<assign\>f<around*|(|<wide|v|^><rsup|\<alpha\>>|)>,
  </equation*>

  and evaluation of both partial derivatives

  <\equation*>
    l<rsub|i+1><rsup|\<alpha\>><around*|(|f|)>\<assign\>f<rsub|,i><around*|(|<wide|v|^><rsup|\<alpha\>>|)>,<space|2em>i\<in\><around*|{|1,2|}>.
  </equation*>

  Finally at the barycenter <math|<wide|v|^><rsup|0>> of <math|<wide|K|^>>
  add a last evaluation

  <\equation*>
    l<rsub|0><around*|(|f|)>\<assign\>f<around*|(|<wide|v|^><rsup|0>|)>.
  </equation*>

  The linear forms <math|<wide|\<cal-L\>|^>\<assign\><around*|{|<wide|l|^><rsup|\<alpha\>><rsub|i>|}>\<cup\><around*|{|<wide|l|^><rsub|0>|}>>
  are a basis of <math|P<rsub|3><around*|(|<wide|K|^>|)><rprime|'>> with dual
  <math|<wide|\<cal-V\>|^>=<around*|{|<wide|\<varphi\>|^><rsub|0>|}>\<cup\><around*|{|<wide|\<varphi\>|^><rsub|i><rsup|\<alpha\>>|}>>
  <cite-detail|brenner_mathematical_2008|3.x.11>. Reindexing the last 9
  functions with <math|k=3*<around*|(|\<alpha\>-1|)>+i>, this means that
  <math|<wide|\<cal-V\>|^>> consists of 10 polynomials
  <math|<wide|\<varphi\>|^><rsub|k>\<in\>P<rsub|3><around*|(|<wide|K|^>|)>>
  called <dfn|shape functions>, such that the following <dfn|duality
  property> holds:

  <\equation>
    <label|eq:delta-property><wide|l|^><rsub|j><around*|(|<wide|\<varphi\>|^><rsub|k>|)>=\<delta\><rsub|j\<nocomma\>k>,<space|2em>j,k\<in\><around*|{|0,\<ldots\>,9|}>.
  </equation>

  In particular <math|<wide|l|^><rsub|1><around*|(|<wide|\<varphi\>|^><rsub|1>|)>=<wide|\<varphi\>|^><rsub|1><around*|(|<wide|v|^><rsub|1>|)>=1>,
  <math|<wide|l|^><rsub|2><around*|(|<wide|\<varphi\>|^><rsub|2>|)>=\<partial\><rsub|1><wide|\<varphi\>|^><rsub|2><around*|(|<wide|v|^><rsub|1>|)>=1>
  and <math|<wide|l|^><rsub|3><around*|(|<wide|\<varphi\>|^><rsub|3>|)>=\<partial\><rsub|2><wide|\<varphi\>|^><rsub|3><around*|(|<wide|v|^><rsub|1>|)>=1>
  but <math|<wide|l|^><rsub|1><around*|(|<wide|\<varphi\>|^><rsub|2>|)>=<wide|\<varphi\>|^><rsub|2><around*|(|<wide|v|^><rsub|1>|)>=0>
  and so on.<\footnote>
    To avoid any possible confusions, we note here that in our
    implementation, <math|<wide|\<varphi\>|^><rsub|0>> is actually assigned
    the <math|10<rsup|th>> position in the array of basis functions.
  </footnote>

  <\big-figure>
    <image|img/hermite2d-w-0.eps|0.32par|||><image|img/hermite2d-w-1.eps|0.32par|||><image|img/hermite2d-w-2.eps|0.32par|||>

    <image|img/hermite2d-w-3.eps|0.32par|||><image|img/hermite2d-w-4.eps|0.32par|||><image|img/hermite2d-w-5.eps|0.32par|||>

    <image|img/hermite2d-w-6.eps|0.32par|||><image|img/hermite2d-w-7.eps|0.32par|||><image|img/hermite2d-w-8.eps|0.32par|||>

    <htab|48mm><image|img/hermite2d-w-9.eps|0.32par|||>
  </big-figure|The ten Hermite shape functions over the reference triangle.
  Here <math|<wide|\<varphi\>|^><rsub|0>> is depicted last.>

  <subsection|The local basis over a physical element>

  The degrees of freedom <math|<wide|\<cal-L\>|^>> are mapped into another
  basis <math|\<cal-L\>=<around*|{|l<rsub|i>:P<rsub|3><around*|(|K|)>\<rightarrow\>\<bbb-R\>|}>>
  of <math|P<rsub|3><around*|(|K|)><rprime|'>> by simple identification of
  the vertices: <math|l<rsup|\<alpha\>><rsub|1><around*|(|f|)>\<assign\>f<around*|(|v<rsup|\<alpha\>>|)>,l<rsup|\<alpha\>><rsub|2><around*|(|f|)>\<assign\>\<partial\><rsub|1>f<around*|(|v<rsup|\<alpha\>>|)>,l<rsup|\<alpha\>><rsub|3><around*|(|f|)>\<assign\>\<partial\><rsub|2>f<around*|(|v<rsup|\<alpha\>>|)>>,
  whereas the reference shape functions are mapped onto a <dfn|local basis of
  shape functions> over the physical cell <math|K> by means of a mapping

  <\equation*>
    <wide|\<varphi\>|^>\<mapsto\>\<varphi\>\<assign\>\<cal-F\><rsub|K><around*|(|<wide|\<varphi\>|^>|)>.
  </equation*>

  For standard Lagrange elements (i.e. when <math|\<cal-L\>> consists only of
  point evaluation forms), the immediate choice
  <math|\<cal-F\><rsub|K><around*|(|<wide|\<varphi\>|^>|)>\<assign\><wide|\<varphi\>|^>\<circ\><wide|\<b-x\>|^>>
  would map <math|<wide|\<cal-V\>|^>> into a basis of
  <math|P<rsub|3><around*|(|K|)>> which fulfills the duality property
  <eqref|eq:delta-property> wrt. <math|\<cal-L\>>. However, note that in
  general:

  <\equation*>
    <wide|\<varphi\>|^><rsub|,i><around*|(|<wide|x|^>|)>\<neq\><around*|(|<wide|\<varphi\>|^>\<circ\><wide|\<b-x\>|^>|)><rsub|,i><around*|(|\<b-x\><around*|(|<wide|x|^>|)>|)>,
  </equation*>

  so that <eqref|eq:delta-property> will not hold in <math|K> for the partial
  derivatives. For the images of the Hermite pairs of dofs
  <math|<around*|(|<wide|\<varphi\>|^><rsub|2><rsup|\<alpha\>>,<wide|\<varphi\>|^><rsub|3><rsup|\<alpha\>>|)>>
  we require a different transformation which at vertex
  <math|v<rsup|\<alpha\>>> of <math|K> sets the shape functions to a linear
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

    The <dfn|Hermite transformation> <math|\<cal-F\><rsub|K>> is given by the
    matrix <math|H> as:

    <\equation>
      <label|eq:hermite-transform>\<Phi\><rsub|\<alpha\>><around*|(|x|)>=\<cal-F\><rsub|K><around*|(|<wide|\<Phi\>|^><rsub|\<alpha\>><around*|(|<wide|\<b-x\>|^><around*|(|x|)>|)>|)>\<assign\>H<around*|(|<wide|\<b-x\>|^><around*|(|x|)>|)>*<wide|\<Phi\>|^><rsub|\<alpha\>><around*|(|<wide|\<b-x\>|^><around*|(|x|)>|)>.
    </equation>
  </definition*>

  <\proposition>
    The Hermite transformation <math|\<cal-F\><rsub|K>>:

    <\enumerate>
      <item>maps <math|<wide|\<cal-V\>|^>> into a basis <math|\<cal-V\>> of
      <math|P<rsub|3><around*|(|K|)>>,

      <item>respects property <eqref|eq:delta-property>.
    </enumerate>
  </proposition>

  <\proof>
    The matrix <math|H> is obtained by solving the linear system resulting
    from imposing condition <eqref|eq:delta-property> for a linear
    combination of the relevant shape functions
    <cite-detail|solin_partial_2006|ü6.4.3>. Property 1 follows from the fact
    that <math|\<b-x\>> is non-singular, hence <math|det H=det
    D\<b-x\>\<neq\>0>. Property 2 for the point evaluations is clear and for
    the others it is a matter of applying the chain rule:
    <math|l<rsub|2><around*|(|\<varphi\><rsub|2>|)>=\<partial\><rsub|1>\<varphi\><rsub|2><around*|(|v<rsub|1>|)>=\<b-x\><rsub|1,k><around*|(|<wide|v|^><rsub|1>|)>*<wide|\<varphi\>|^><rsup|\<alpha\>><rsub|k,l><around*|(|<wide|v|^><rsub|1>|)>*<wide|\<b-x\>|^><rsub|l,1><around*|(|v<rsub|1>|)>=\<b-x\><rsub|1,k><around*|(|<wide|v|^><rsub|1>|)>*\<delta\><rsub|k\<nocomma\>l>*<wide|\<b-x\>|^><rsub|l,1><around*|(|v<rsub|1>|)>=1>,
    where <math|k,l\<in\><around*|{|1,2|}>>.
  </proof>

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
    <label|eq:hermite-first-derivatives>\<varphi\><rsub|i,j><rsup|\<alpha\>><around*|(|x|)>=<around*|(|H<rsub|i\<nocomma\>k>*<wide|\<varphi\>|^><rsup|\<alpha\>><rsub|k,l>|)><around*|(|<wide|\<b-x\>|^><around*|(|x|)>|)>*<wide|\<b-x\>|^><rsub|l,j><around*|(|x|)>,<space|2em>i\<in\><around*|{|1,2,3|}>.
  </equation>

  Each partial derivative <math|\<varphi\><rsub|i,j><rsup|\<alpha\>>>,
  <math|j\<in\><around*|{|1,2|}>> is a linear combination of the derivatives
  of both reference basis functions <math|<wide|\<varphi\>|^><rsub|k,l><rsup|\<alpha\>>>
  with weights given by <math|\<nabla\> \<b-x\><rsub|l>>:

  <\equation*>
    \<varphi\><rsub|i+1,j><rsup|\<alpha\>><around*|(|x|)>=<around*|(|\<b-x\><rsub|i,k>*<wide|\<varphi\>|^><rsup|\<alpha\>><rsub|k,l>|)><around*|(|<wide|\<b-x\>|^><around*|(|x|)>|)>*<wide|\<b-x\>|^><rsub|l,j><around*|(|x|)>,<space|2em>i\<in\><around*|{|1,2|}>.
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

  With <math|3> derivatives:

  <\equation*>
    \<varphi\><rsub|i,j\<nocomma\>k\<nocomma\>m><rsup|\<alpha\>>=H<rsub|i\<nocomma\>p>*<wide|\<varphi\>|^><rsup|\<alpha\>><rsub|p,q\<nocomma\>r\<nocomma\>s>*<wide|\<b-x\>|^><rsub|q,j>*<wide|\<b-x\>|^><rsub|r,k>*<wide|\<b-x\>|^><rsub|s,m>.
  </equation*>

  <\proposition>
    <todo|generalise to <math|k> derivatives with a multiindex>.
  </proposition>

  <subsection|The global basis for <math|V<rsub|h>>>

  With the previous construction we obtain a collection of triples
  <math|<around*|{|<around*|(|K,\<cal-V\><rsub|K>,\<cal-L\><rsub|K>|)>|}><rsub|K\<in\>\<cal-T\>>>.
  We now gather all of the shape functions into a basis for a polynomial
  space <math|V<rsub|h>> with the help of a collection of
  <dfn|local-to-global maps> <math|\<iota\><rsub|K>>. These map indices local
  to a simplex into a global set and ensure any continuity properties by
  requiring local dofs to agree. Specifically, for each
  <math|K\<in\>\<cal-T\>> we define

  <\equation*>
    \<iota\><rsub|K>:I<rsub|K>=<around*|{|1,\<ldots\>,n<rsub|K>|}>\<rightarrow\>I=<around*|{|1,\<ldots\>,N|}>.
  </equation*>

  This defines a family of linear functionals
  <math|\<cal-L\>=<around*|{|l<rsub|i>|}><rsub|i=1><rsup|N>> via

  <\equation*>
    l<rsub|\<iota\><rsub|K><around*|(|i|)>><around*|(|v|)>\<assign\>l<rsub|i><rsup|K><around*|(|v<rsub|\|K>|)>,<space|2em>i\<in\>I<rsub|K>,
  </equation*>

  for every <math|v\<in\>V<rsub|h>>, where <math|V<rsub|h>> consists of the
  functions <math|v> satisfying:

  <\enumerate>
    <item><math|v<rsub|\|K>\<in\>\<cal-V\><rsub|K>>.

    <item>For all <math|K,K<rprime|'>> and
    <math|<around*|(|i,i<rprime|'>|)>\<in\>I<rsub|K>\<times\>I<rsub|K<rprime|'>>>
    such that <math|\<iota\><rsub|K><around*|(|i|)>=\<iota\><rsub|K<rprime|'>><around*|(|i<rprime|'>|)>>
    it holds that

    <\equation*>
      l<rsub|i><rsup|K><around*|(|v<rsub|\|K>|)>=l<rsup|K<rprime|'>><rsub|i<rprime|'>><around*|(|v<rsub|\|K<rprime|'>>|)>.
    </equation*>
  </enumerate>

  <section|Construction of the linear system>

  We now review the steps required to piece together the linear system
  discretising <eqref|eq:linear-problem>.

  <subsection|Assembly of the stiffness matrix>

  In order to discretise <eqref|eq:linear-problem> we require the action of
  <math|a> on the polynomial space <math|V<rsub|h>=span<around*|{|\<varphi\><rsub|1>,\<ldots\>,\<varphi\><rsub|M>|}>>,
  where <math|M=10*N> in the case of 2D Cubic Hermite elements. Using
  superindices to tell the two copies of <math|V<rsub|h>> apart we want to
  compute the values of

  <\equation*>
    a:V<rsup|<around*|(|1|)>><rsub|h>\<times\>V<rsup|<around*|(|2|)>><rsub|h>\<rightarrow\>\<bbb-R\>,
  </equation*>

  over all basis functions. This can be split into the contribution of each
  element, meaning that we compute the values

  <\equation*>
    A<rsub|i\<nocomma\>j>=a<around*|(|\<varphi\><rsub|j><rsup|<around*|(|1|)>>,\<varphi\><rsup|<around*|(|2|)>><rsub|i>|)>=<big|sum><rsub|K\<in\>\<cal-T\>>a<rsup|K><around*|(|\<varphi\><rsub|j><rsup|<around*|(|1|)>>,\<varphi\><rsup|<around*|(|2|)>><rsub|i>|)>
  </equation*>

  with \ <math|0\<leqslant\>i,j\<less\>M> and <math|a<rsup|K>> the
  contribution from element <math|K>. Applying the same idea to the linear
  form at the right hand side of the equation yields a global linear system
  <math|A*u=b> as in <eqref|eq:linear-system>. For performance reasons, this
  process of <em|assembly> of the matrix <math|A> is actually conducted in
  reverse fashion, by iterating over all cells <math|K> in the triangulation
  and adding up their contributions (the relevant integrals) to each entry
  <math|A<rsub|i\<nocomma\>j>>.

  <subsection|Imposing boundary conditions>

  A standard way of imposing Dirichlet boundary conditions (and the one used
  in <name|FEniCS>) is to force the coefficients of the discrete solution
  corresponding to nodal Lagrange shape functions at the boundary to have the
  desired values. If the solution has to equal some value <math|g<rsub|0>> at
  a point whose corresponding global dof has index <math|i<rsub|0>>, then
  this is easily achieved by setting <math|A<rsub|i<rsub|0>\<nocomma\>j>=\<delta\><rsub|i<rsub|0>\<nocomma\>j>>
  for all <math|j> and <math|b<rsub|i<rsub|0>>=g<rsub|0>> before solving the
  system, see <cite-detail|quarteroni_numerical_2009|ü8.4.5>.<\footnote>
    Actually, one would like to preserve any symmetries that the system may
    have without the boundary conditions and the operation described would
    break them. This can be avoided if the boundary condition is applied
    element-wise during the assembly, as is done in <name|FEniCS>.
  </footnote>

  One must be careful however with Hermite elements, since not all dofs at a
  boundary node in the mesh are of Lagrange type: setting the Hermite ones
  effectively fixes the first derivatives of the solution at that point, not
  its value.

  This turns out to be an issue with the current version of
  <name|FEniCS>,<\footnote>
    Version <tt|2017.1.0.dev0>.
  </footnote> where <cpp|DirichletBC> assumes that all dofs at a mesh node
  stem from point evaluations (as does most of the library). See the
  implementation details in Section <reference|sec:bcs-implementation> below.

  <subsection|Computing integrals>

  The contribution <math|a<rsup|K>> from element <math|K> to the stiffness
  requires computing some integral. For example, if we are discretising the
  biharmonic equation <math|\<Delta\><rsup|2>u=f>, assuming homogeneous
  boundary conditions, the corresponding bilinear form is
  <math|a<around*|(|u,v|)>=<big|int>\<Delta\>u*\<Delta\>v>, so we need to
  compute

  <\equation*>
    a<rsup|K><around*|(|\<varphi\><rsub|j><rsup|<around*|(|1|)>>,\<varphi\><rsup|<around*|(|2|)>><rsub|i>|)>=<big|int><rsub|K>\<Delta\>\<varphi\><rsub|j><rsup|<around*|(|1|)>>*\<mathLaplace\>\<varphi\><rsup|<around*|(|2|)>><rsub|i>*\<mathd\>x=<big|int><rsub|<wide|K|^>><todo|\<ldots\>>.
  </equation*>

  <todo|As can be seen, the transformation introduced above plays here a role
  as well...> Notice that even if <math|V<rsub|h>\<nsubset\>H<rsup|2>>
  globally, we can still compute the Laplacian in the interior of each
  element <math|K>. However, the normal derivatives might present
  discontinuities, so we must modify <math|a<rsup|K>> to enforce some sort of
  weak continuity (see below).

  In order to actually compute the integrals, we implement a quadrature
  representation of the bilinear form as described in
  <cite|olgaard_optimizations_2010> (instead of a <em|tensor representation>,
  cf. <inactive|<cite-detail|logg...|Chapter about tensor repr.>>).

  \;

  <section|Interpolation>

  <\note*>
    See <cite-detail|brenner_mathematical_2008|ü4.8> and
    <cite|girault_hermite_2002> for info on interpolation of non-smooth
    functions.
  </note*>

  As already mentioned, we must interpolate the data of the problem into the
  relevant function space, which in our examples above and below is
  <math|H<rsup|2><around*|(|\<Omega\>|)>>. We have at least the following
  options for interpolating functions <math|g\<in\>H<rsup|2><around*|(|\<Omega\>|)>>
  into <math|V<rsub|h>>, from fastest and less accurate to slowest and most
  accurate <cite-detail|solin_partial_2006|ü6.3.8>.\ 

  <\enumerate>
    <item>Compute the nodal interpolant, cf.
    ü<reference|sec:nodal-interpolant>.

    <item>Use nodal interpolation of vertex and derivatives plus local
    orthogonal projections in the element interiors. This improves the nodal
    interpolant in the interior for dimensions 2 or greater or higher order
    Hermite elements.

    <item>Use a global orthogonal projection to compute the best interpolant
    in the <math|H<rsup|2>>-sense. If the boundary conditions allow it (e.g.
    if we have Dirichlet and Poincarť's inequality holds) then one can use
    the <math|H<rsup|2><rsub|0>> seminorm, cf.
    ü<reference|sec:H2-projection>.<\footnote>
      <todo|Does FEniCS do this by default with <python|project()>?>.
    </footnote>
  </enumerate>

  <subsection|The local nodal interpolant><label|sec:nodal-interpolant>

  The local nodal interpolant over a triangle <math|K> of the mesh is given
  by

  <\equation*>
    I<rsub|K><around*|(|g|)>=l<rsub|0><around*|(|g|)>*\<varphi\><rsub|0>+<big|sum><rsub|\<alpha\>=1><rsup|3><big|sum><rsub|i=1><rsup|3>l<rsup|\<alpha\>><rsub|i><around*|(|g|)>*\<varphi\><rsub|i>
  </equation*>

  where <math|l<rsup|\<alpha\>><rsub|i><around*|(|g|)>=\<partial\><rsub|i-1>
  g<around*|(|v<rsup|\<alpha\>>|)>> and <math|\<partial\><rsub|0>
  f\<assign\>f>, and <math|\<varphi\><rsub|i>> are the shape functions.

  Here a difficulty particular to Hermite elements presents itself:

  <\question*>
    How should one compute the derivatives? Should we extend the interface of
    <cpp|ufc::function> to include differentiation?
  </question*>

  A simple first solution is to only allow interpolation of expressions whose
  derivative can be evaluated, e.g. via automatic differentiation. In
  <name|Python>, using <name|autograd> it is straightforward to implemement a
  decorator extending any <python|Callable> for the evaluation of partial
  derivatives and this is the path we first take.

  <subsection|<math|H<rsup|2>>-projection><label|sec:H2-projection>

  This means solving the system <todo|...>. See
  <cite-detail|solin_partial_2006|ü6.3.8>,
  <cite|brenner_mathematical_2008>...

  Implementation is straightforward after defining the relevant forms in the
  UFL file and writing a small <c++> wrapper for the solution of the system.

  <section|The Euler-Bernoulli beam model>

  As a first test of our implementation of Hermite elements we use a one
  dimensional model for which these provide a conforming approximation. The
  classical Euler-Bernoulli model studies the deformation of the midplane of
  a beam subject to a transversal load, under the assumptions that after the
  deformation the normals to the midplane:

  <\itemize-dot>
    <item>do not bend,

    <item>do not stretch,

    <item>remain orthogonal to the midplane.
  </itemize-dot>

  This effectively reduces a 2D beam problem to a 1D problem. This theory is
  adequate for thin beams and is intended for small strains even with large
  global deformations: it is a <em|physically linear> but <em|geometrically
  non-linear> theory. For thicker beams, Timoshenko's theory, which accounts
  for internal shear forces, yields more accurate predictions
  <todo|<inactive|<cite|>>>.

  <subsection|Derivation>

  Fix <math|\<omega\>=<around*|(|a,b|)>> to be the midplane of the beam
  represented by the domain <math|\<Omega\>=\<omega\>\<times\><around*|(|-h/2,h/2|)>\<subset\>\<bbb-R\><rsup|2>>.
  Conservation of momentum and certain linear constitutive relations yield
  <cite-detail|solin_partial_2006|ü6.1.1> the equation

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
  for the four boundary conditions that in this case have to be specified:

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

  A simple application of HŲlder and Poincarť is required to check the
  conditions for existence in Lax-Milgram's Lemma. For ellipticity, assume
  that <math|b> is a.e. positive and bounded away from zero.

  <subsection|Discretization>

  Recall that <math|H<rsup|2><around|(|a,b|)>\<subset\>C<rsup|1,\<gamma\>><around|(|a,b|)>,\<gamma\>=1/2>
  by the Sobolev embeddings. Therefore we use cubic Hermite elements, which
  are <math|C<rsup|1>> in <math|\<bbb-R\>> and thus provide a conforming
  discretization (see <cite-detail|solin_partial_2006|p. 218>). The
  specification of the linear form and the assembly of the stiffness matrix
  are done using <name|dolfin> compiled to use the extensions to <name|FIAT>
  and <name|FFC> detailed below. Boundary conditions need special handling,
  see also the Appendix.

  <subsection|Clamped beam>

  We look for a solution with clamped boundaries, i.e.
  <math|u<around|(|a|)>=u<around|(|b|)>=u<rprime|'><around|(|a|)>=u<rprime|'><around|(|b|)>=0>.
  Because of the essential constraints, the terms corresponding to the
  natural ones in the weak formulation vanish and we are left with:

  <\equation*>
    V=V<rsub|clamped>=<around|{|v\<in\>H<rsup|2><around|(|\<Omega\>|)>:v<around|(|\<alpha\>|)>=v<rsub|\<alpha\>>,v<rprime|'><around|(|\<alpha\>|)>=v<rprime|'><rsub|\<alpha\>>,\<alpha\>=1,2|}>.
  </equation*>

  and <math|\<gamma\><rsub|M>=\<gamma\><rsub|F>=\<emptyset\>>.

  <big-figure|<image|img/clamped-beam.eps|0.8par|||>|Clamped steel beam under
  the effect of gravity alone.>

  <subsection|Cantilevered beam>

  One boundary is clamped, e.g. <math|u<around|(|a|)>=u<rprime|'><around|(|a|)>=0>,
  and the other is left hanging freely, so that no bending moment (i.e. no
  torque) and no shear force appear, e.g. assuming constant <math|E> and
  <math|I>: <math|u<rprime|''><around|(|b|)>=u<rprime|'''><around|(|b|)>=0>.

  We can generalise this to the case where we hang a mass from the free end
  of the beam. The downward force will introduce shear, e.g.
  <math|u<rprime|'''><around|(|b|)>=-m*g>, with <math|m> the mass of the
  object and <math|g> gravity. Just set <math|m=0> for the first situation.

  <big-figure|<image|img/cantilevered-beam.eps|0.8par|||>|Cantilevered steel
  beam with a weight hanging from the right endpoint.>

  <subsection|Simply supported beam>

  The beam rests on two point supports at its ends. Displacements are fixed
  at both ends with <math|u<around|(|a|)>=u<around|(|b|)>=0>, but since the
  beam is free to rotate, it experiences no torque at these points, i.e. the
  bending moments are zero: <math|u<rprime|''><around|(|a|)>=u<rprime|''><around|(|b|)>=0>,
  so we can plug in the natural boundary condition <math|M=0> and the problem
  is: Find <math|u\<in\>V=H<rsup|2><rsub|0><around|(|\<omega\>|)>> such that
  for all <math|v\<in\>V>:

  <\equation*>
    <big|int><rsub|\<omega\>>b*\<Delta\>u*\<Delta\>v*\<mathd\>x+<big|int><rsub|\<gamma\>>\<nabla\>*<around|(|b*\<Delta\>u|)>*v\<mathd\>s=<big|int><rsub|\<omega\>>f*v*\<mathd\>x.
  </equation*>

  <big-figure|<image|img/supported-beam.eps|0.8par|||>|Simply supported steel
  beam under the effect of gravity alone.>

  <section|The biharmonic equation>

  We now use Hermite elements for a non-conforming discretisation of a fourth
  order problem in <math|\<bbb-R\><rsup|2>>. Following the example in
  <name|dolfin>'s repository, we consider the Biharmonic equation:

  <\equation*>
    \<nabla\><rsup|4>u=f<text| over >\<Omega\>=<around|[|0,\<pi\>|]>\<times\>[-\<mathpi\>/2,\<mathpi\>/2]
  </equation*>

  with source <math|f> given by

  <\equation*>
    f<around|(|x,y|)>=4*sin <around|(|x|)>*cos <around|(|y|)>
  </equation*>

  and homogeneous boundary conditions

  <\equation*>
    u=0 <text| and >\<nabla\><rsup|2>u=0
  </equation*>

  on the whole boundary. The analytic solution is

  <\equation*>
    u<around*|(|x,y|)>=sin<around*|(|x|)>*cos<around*|(|y|)>.
  </equation*>

  We use a discontinuous Galerkin formulation (interior penalty method) which
  imposes weak continuity of the normal derivatives across facets and enables
  us to implement a non-conforming discretisation (see e.g.
  <cite|brenner_c0_2005>). <todo|Why does this make any sense at all for
  Hermite elements?> This results in the following weak formulation. Let
  <math|V<rsub|h>> be a finite element space of Lagrange or Hermite type,
  with polynomial degree 3. Find <math|u\<in\>V<rsub|h>> such that

  <\equation*>
    a<around*|(|u,v|)>=<big|int><rsub|\<Omega\>>v*f*\<mathd\>x
  </equation*>

  holds for all <math|v\<in\>V<rsub|h>>, where\ 

  <\equation*>
    a<around*|(|u,v|)>\<assign\><big|sum><rsub|K\<in\><with|math-font|cal|\<cal-T\><rsub|h>>><big|int><rsub|K>\<nabla\><rsup|2>u*\<nabla\><rsup|2>v*\<mathd\>x-<big|int><rsub|\<Gamma\><rsub|K>><around*|\<langle\>|\<nabla\><rsup|2>*u|\<rangle\>>*<around*|\<llbracket\>|\<nabla\>v|\<rrbracket\>><space|0.17em>\<mathd\>s-<big|int><rsub|\<Gamma\><rsub|K>><around*|\<llbracket\>|\<nabla\>u|\<rrbracket\>><around*|\<langle\>|\<nabla\><rsup|2>*v|\<rangle\>>*<space|0.17em>\<mathd\>s+<big|int><rsub|\<Gamma\><rsub|K>><frac|\<alpha\>|h>*<around*|\<llbracket\>|\<nabla\>u|\<rrbracket\>>*<around*|\<llbracket\>|\<nabla\>v|\<rrbracket\>>*\<mathd\>s,
  </equation*>

  <math|<around*|\<langle\>|u|\<rangle\>>\<assign\><frac|1|2>*<around|(|u<rsub|+>+u<rsub|->|)>>
  is the average of <math|u> across some edge (<math|u<rsub|\<pm\>>> being
  its value at each side of the edge), <math|<around*|\<llbracket\>|w|\<rrbracket\>>=w<rsub|+>\<cdot\>n<rsub|+>+w<rsub|->\<cdot\>n<rsub|->>
  is the jump in normal direction of <math|w>, <math|\<alpha\>\<geqslant\>0>
  is a penalty term and <math|h> is some measure of the size of <math|K>.

  The implementation in FEniCS is straightforward and follows closely the
  standard demo in the package, so that we can safely omit it here. Suffice
  to note that the only delicate point is the enforcement of the essential
  (Dirichlet) boundary conditions, where one must only fix degrees of freedom
  related to Lagrange evaluations and not Hermite.

  <big-figure|<image|img/biharmonic-solutions.eps|0.9par|||>|A visual
  inspection of the solutions isn't very revealing...>

  <big-figure|<image|img/biharmonic-diffs.eps|0.9par|||>|The relative
  difference between the Lagrange and Hermite solutions decreases [...] as
  the grid resolution increases.>

  <big-figure|<image|img/biharmonic-times.eps|1par|||>|Enforcing essential
  boundary conditions after assembly requires manipulation of the system
  matrix, invalidating its sparsity pattern. The cost of a naive <em|a
  posteriori> modification in <name|Python> is exponential in the number of
  degrees of freedom.>

  \;

  <appendix|Implementation details>

  Here is an incomplete list of things that needed doing for the
  implementation in <name|FEniCS>:

  <\enumerate>
    <item>Fixing the <name|FIAT> Hermite element definition to instantiate
    the right dual basis and build the right <python|entity_ids>. See Section
    <reference|sec:fiat-implementation> below.

    <item>Implementing in <name|FFC> the Hermite transformation described
    above for the element basis functions which are associated to the
    evaluation of partial derivatives at the nodes of the simplex. This meant
    modifying:

    <\itemize-dot>
      <item><tt|ffc/representation.py>: for every dof
      <math|<wide|l|^><rsub|i><rsup|\<alpha\>>> include the coefficients of
      the other dofs <math|<wide|l|^><rsub|j><rsup|\<alpha\>>> for
      <math|j\<neq\>i> in the <tt|dof_data> structure passed to the next
      compiler stage.

      <item><tt|ffc/evaluatebasis.py>: using the Jacobian of the geometric
      transformation, compute the basis evaluations as per
      <eqref|eq:hermite-first-derivatives>. The intermediate variable
      <python|dof_data['hermite_node_offset']> is used in
      <python|_compute_values()> to compute the right linear combination of
      reference basis functions.
    </itemize-dot>

    <item>Implementing the same transformation for the quadratures. See
    <python|create_argument()> and <python|create_function()> in
    <tt|ffc/quadratures/quadraturetransformerbase.py>.

    <item>Amending the assembled stiffness matrix to account for the presence
    of Hermite dofs in boundary conditions. In particular undoing the changes
    done by <cpp|DirichletBC.apply()> to the corresponding rows.

    <item>Implementing the nodal interpolant. See Section
    <reference|sec:ad-nodal-interpolation> below.
  </enumerate>

  <subsection|<name|FIAT> implementation><label|sec:fiat-implementation>

  <\python-code>
    from FIAT import finite_element, polynomial_set, dual_set, functional,
    reference_element

    \;

    class CubicHermiteDualSet(dual_set.DualSet):

    \ \ \ \ """The dual basis for Cubic Hermite elements.

    \ \ \ \ This class works for simplices of any dimension.

    \ \ \ \ Nodes are point evaluation and derivatives at vertices,

    \ \ \ \ plus point evaluation at the barycenter of each face if

    \ \ \ \ the spatial dimension is \<gtr\>= 2. """

    \;

    \ \ \ \ def __init__(self, ref_el):

    \ \ \ \ \ \ \ \ entity_ids = {}

    \ \ \ \ \ \ \ \ nodes = []

    \ \ \ \ \ \ \ \ cur = 0

    \;

    \ \ \ \ \ \ \ \ top = ref_el.get_topology()

    \ \ \ \ \ \ \ \ verts = ref_el.get_vertices()

    \ \ \ \ \ \ \ \ sd = ref_el.get_spatial_dimension()

    \;

    \ \ \ \ \ \ \ \ # Vertex dofs: one point evaluation and sd partial
    derivatives

    \ \ \ \ \ \ \ \ entity_ids[0] = {}

    \ \ \ \ \ \ \ \ for v in sorted(top[0]):

    \ \ \ \ \ \ \ \ \ \ \ \ nodes.append(functional.PointEvaluation(ref_el,
    verts[v]))

    \ \ \ \ \ \ \ \ \ \ \ \ pd = functional.PointDerivative

    \ \ \ \ \ \ \ \ \ \ \ \ for i in range(sd):

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ alpha = [0] * sd

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ alpha[i] = 1

    \;

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ nodes.append(pd(ref_el, verts[v], alpha))

    \;

    \ \ \ \ \ \ \ \ \ \ \ \ entity_ids[0][v] = list(range(cur, cur + 1 + sd))

    \ \ \ \ \ \ \ \ \ \ \ \ cur += sd + 1

    \;

    \ \ \ \ \ \ \ \ # no edge dofs

    \ \ \ \ \ \ \ \ entity_ids[1] = {}

    \ \ \ \ \ \ \ \ for e in sorted(top[1]):

    \ \ \ \ \ \ \ \ \ \ \ \ entity_ids[1][e] = []

    \;

    \ \ \ \ \ \ \ \ if sd \<gtr\>= 2:

    \ \ \ \ \ \ \ \ \ \ \ \ # face dof

    \ \ \ \ \ \ \ \ \ \ \ \ # point evaluation at barycenter

    \ \ \ \ \ \ \ \ \ \ \ \ entity_ids[2] = {}

    \ \ \ \ \ \ \ \ \ \ \ \ for f in sorted(top[2]):

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ pt = ref_el.make_points(2, f, 3)[0]

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ n = functional.PointEvaluation(ref_el,
    pt)

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ nodes.append(n)

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ entity_ids[2][f] = list(range(cur, cur +
    1))

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ cur += 1

    \;

    \ \ \ \ \ \ \ \ \ \ \ \ # No more dof in higher dimensions

    \ \ \ \ \ \ \ \ \ \ \ \ for dim in range(3, sd + 1):

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ entity_ids[dim] = {x: [] for x in
    sorted(top[dim])}

    \;

    \ \ \ \ \ \ \ \ super(CubicHermiteDualSet, self).__init__(nodes, ref_el,
    entity_ids)
  </python-code>

  The finite element is defined as:

  <\python-code>
    class CubicHermite(finite_element.CiarletElement):

    \ \ \ \ """ The cubic Hermite finite element.

    \ \ \ \ There are three degrees of freedom per vertex and one additional

    \ \ \ \ for each barycenter of a 2d face. """

    \;

    \ \ \ \ def __init__(self, ref_el, degree=3):

    \ \ \ \ \ \ \ \ # Degree is fixed to 3 for cubic Hermite elements.

    \ \ \ \ \ \ \ \ poly_set = polynomial_set.ONPolynomialSet(ref_el, 3)

    \ \ \ \ \ \ \ \ dual = CubicHermiteDualSet(ref_el)

    \ \ \ \ \ \ \ \ super(CubicHermite, self).__init__(poly_set, dual, 3)

    \;

    \ \ \ \ def mapping(self):

    \ \ \ \ \ \ \ \ # MBD FIXME: is this correct?

    \ \ \ \ \ \ \ \ mappings = []

    \ \ \ \ \ \ \ \ for f in self.dual_basis():

    \ \ \ \ \ \ \ \ \ \ \ \ if isinstance(f, functional.PointEvaluation):

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ mappings.append("affine")

    \ \ \ \ \ \ \ \ \ \ \ \ elif isinstance(f, functional.PointDerivative):

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ mappings.append("hermite")

    \ \ \ \ \ \ \ \ return mappings
  </python-code>

  <subsection|Computing the Hermite shape
  functions><label|sec:shape-implementation>

  In order to test the <name|FIAT> implementation of Hermite elements we
  would like to \Smanually\T compute the shape functions. Let
  <math|<wide|K|^>> be the reference triangle with vertices at
  <math|<around*|(|0,0|)>,<around*|(|1,0|)>> and <math|<around*|(|0,1|)>>. We
  choose the monomial basis <math|G=<around|{|1,x,y,x<rsup|2>,x*y,y<rsup|2>,x<rsup|3>,x<rsup|2>*y,x*y<rsup|2>,y<rsup|3>|}>>
  to express the shape functions in. The 10 degrees of freedom
  <math|L=<around|{|l<rsub|i>|}><rsub|i=0><rsup|9>\<subset\>P<rsub|3><around|(|<wide|K|^>|)><rprime|'>>
  are point evaluation and partial differentiation at the vertices, plus
  evaluation at the barycenter, as described above. We build the Vandermonde
  matrix <math|V<rsub|i*j>=L<rsub|i><around|(|G<rsub|j>|)>> and invert it:
  the columns of <math|V<rsup|-1>> are the coefficients of the shape
  functions expressed in the basis <math|G>.

  This computation can be easily done manually or using Automatic
  Differentitation with <name|autograd> and the following bit of code:

  <\python-code>
    import autograd as ad

    import autograd.numpy as np

    \;

    # Partial derivatives

    def dx(f):

    \ \ \ \ return ad.grad(f, 0)

    def dy(f):

    \ \ \ \ return ad.grad(f, 1)

    \;

    # Monomial basis for $P_3(R^2)$

    G = [lambda x,y: 1, lambda x,y: x, lambda x,y: y, lambda x,y: x**2,\ 

    \ \ \ \ \ lambda x,y: x*y, lambda x,y: y**2, lambda x,y: x**3, lambda
    x,y: x**2*y,

    \ \ \ \ \ lambda x,y: x*y**2, lambda x,y: y**3]

    \;

    # Hermite degrees of freedom on the reference triangle
    [(0,0),(1,0),(0,1)]:

    # point evaluation and partial derivatives at each vertex, plus
    evaluation at

    # the barycenter.

    L = [lambda f: f(0., 0.), lambda f: dx(f)(0., 0.), lambda f: dy(f)(0.,
    0.),

    \ \ \ \ \ lambda f: f(1., 0.), lambda f: dx(f)(1., 0.), lambda f:
    dy(f)(1., 0.),

    \ \ \ \ \ lambda f: f(0., 1.), lambda f: dx(f)(0., 1.), lambda f:
    dy(f)(0., 1.),

    \ \ \ \ \ lambda f: f(1./3., 1./3.)]

    \;

    # Build Vandermonde matrix

    V = np.zeros((10,10))

    for i, j in np.ndindex(V.shape):

    \ \ \ \ V[i,j] = L[i](G[j])

    \;

    # Invert and read linear combinations of basis functions from G from the
    columns

    Vinv = np.linalg.inv(V)

    \;

    # Write strings to be pasted into tests/unit/test_elements.py

    S = ["1", "x[0]", "x[1]", "x[0]**2", "x[0]*x[1]", "x[1]**2",

    \ \ \ \ \ "x[0]**3", "x[0]**2*x[1]", "x[0]*x[1]**2", "x[1]**3"]

    s = "sh = [None] * 10\\n"

    for j in range(10):

    \ \ \ \ s += "sh[%d] = lambda x: " % j

    \ \ \ \ l = []

    \ \ \ \ for i in range(10):

    \ \ \ \ \ \ \ \ coeff = Vinv[i,j]

    \ \ \ \ \ \ \ \ # Try to simplify coefficients a bit

    \ \ \ \ \ \ \ \ # We have alreay checked that they are integers

    \ \ \ \ \ \ \ \ if not np.isclose(coeff, 0):

    \ \ \ \ \ \ \ \ \ \ \ \ sign = "+" if coeff \<gtr\> 0 else "-"

    \ \ \ \ \ \ \ \ \ \ \ \ coeff = np.abs(coeff)

    \ \ \ \ \ \ \ \ \ \ \ \ prefix = "%s %d *" % (sign, coeff) if not
    np.isclose(coeff, 1)\\

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ else
    sign

    \ \ \ \ \ \ \ \ \ \ \ \ l.append("%s %s" % (prefix, S[i]))

    \ \ \ \ s += " ".join(l) + "\\n"

    print(s)
  </python-code>

  <subsection|Essential boundary conditions for Hermite
  elements><label|sec:bcs-implementation>

  A simple implementation of Dirichlet BCs (without enforcing symmetry) is as
  follows. The case of Neumann conditions is a bit more involved because we
  need to single out Hermite dofs at the boundary, but the idea is similar.

  <\python-code>
    def apply_dirichlet_hermite(A:Matrix, b:Vector, bc:DirichletBC):

    \ \ \ \ warning("This won't work in parallel!")

    \ \ \ \ rows = list_hermite_boundary_dofs(bc)

    \ \ \ \ cols = np.arange(A.size(1), dtype=np.intc)

    \ \ \ \ saved_block = np.empty((rows.size, cols.size), dtype=np.float)

    \ \ \ \ A.get(saved_block, rows, cols)

    \ \ \ \ saved_vec = b[rows].copy() \ \ 

    \ \ \ \ bc.apply(A, b)

    \ \ \ \ b[rows] = saved_vec

    \ \ \ \ M = as_backend_type(A)

    \ \ \ \ if isinstance(M, dolfin.cpp.la.PETScMatrix):

    \ \ \ \ \ \ \ \ M.mat().setOption(PETSc.MAT_NEW_NONZERO_ALLOCATION_ERR,
    PETSc.PETSc_FALSE)

    \ \ \ \ M.set(saved_block, rows, cols) \ # Slow!!!

    \ \ \ \ M.apply('insert')
  </python-code>

  <subsection|Automatic differentiation for nodal
  interpolation><label|sec:ad-nodal-interpolation>

  For consistency with the rest of <name|FEniCS> we extend the interface of
  <python|dolfin.Expression> to take partial derivatives. Of course we can
  only do this for expressions where we can actually compute these and this
  is most easily done with autograd applied to callables:

  <\python-code>
    import autograd as ad

    from typing import Callable

    from FIAT import functional

    from ffc import fiatinterface

    ##### Docstrings, comments and checks omitted #####

    \;

    def make_derivatives(fun:Callable) -\<gtr\> Callable:

    \ \ \ \ def diff_fun(x, derivatives=()):

    \ \ \ \ \ \ \ \ ret = fun

    \ \ \ \ \ \ \ \ for i, n in enumerate(derivatives):

    \ \ \ \ \ \ \ \ \ \ \ \ for j in range(n):

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ ret = ad.grad(ret, i)

    \ \ \ \ \ \ \ \ return ret(*x)

    \ \ \ \ # Some boilerplate omitted...

    \ \ \ \ return diff_fun

    \;

    class ExpressionAD(Expression):

    \ \ \ \ def __init__(self, **kwargs):

    \ \ \ \ \ \ \ \ self.fun = make_derivatives(kwargs['fun'])

    \ \ \ \ def eval(self, value, x):

    \ \ \ \ \ \ \ \ value[0] = self.fun(x)

    \ \ \ \ def partial(self, value, x, derivatives):

    \ \ \ \ \ \ \ \ value[0] = self.fun(x, derivatives)

    \ \ \ \ def __call__(self, x, derivatives=()):

    \ \ \ \ \ \ \ \ return self.fun(x, derivatives)

    \;

    def interpolate_hermite(f:ExpressionAD, V:FunctionSpace) -\<gtr\>
    Function:

    \ \ \ \ tdim = V.element().topological_dimension()

    \ \ \ \ e = fiatinterface.create_element(V.ufl_element())

    \ \ \ \ hmask = np.array(list(map(lambda f: isinstance(f,
    functional.PointDerivative),

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ e.dual_basis())))

    \ \ \ \ lmask = np.array(list(map(lambda f: isinstance(f,
    functional.PointEvaluation),

    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ e.dual_basis())))

    \ \ \ \ uv = np.zeros(V.dim())

    \ \ \ \ dm = V.dofmap()

    \ \ \ \ dof_coordinates = V.tabulate_dof_coordinates().reshape((-1,
    tdim))

    \ \ \ \ derivatives = partial_derivatives(tdim) # e.g.
    [(1,0,0),(0,1,0),(0,0,1)]*4

    \ \ \ \ for i in range(V.mesh().num_cells()):

    \ \ \ \ \ \ \ \ dofs = dm.cell_dofs(i)

    \ \ \ \ \ \ \ \ hermite_dofs = dofs[hmask]

    \ \ \ \ \ \ \ \ lagrange_dofs = dofs[lmask]

    \ \ \ \ \ \ \ \ for dd, hdof in zip(derivatives, hermite_dofs):

    \ \ \ \ \ \ \ \ \ \ \ \ uv[hdof] = f(dof_coordinates[hdof],
    derivatives=dd)

    \ \ \ \ \ \ \ \ for ldof in lagrange_dofs:

    \ \ \ \ \ \ \ \ \ \ \ \ uv[ldof] = f(dof_coordinates[ldof],
    derivatives=())

    \ \ \ \ u = Function(V)

    \ \ \ \ u.vector().set_local(uv)

    \ \ \ \ return u
  </python-code>

  Usage is simple and homogeneous with <name|dolfin>'s <python|project()>:

  <\python-code>
    V = FunctionSpace(UnitSquareMesh(20, 20), 'Hermite', 3)

    f = ExpressionAD(fun=lambda x,y: x*(y+2)+x**2, degree=3)

    u = interpolate_hermite(f, V) \ # Nodal interpolation

    v = project(f, V)

    w = Function(V, u.vector() - v.vector()) \ # difference
  </python-code>

  <\bibliography|bib|tm-alpha|hermite.bib>
    <\bib-list|9>
      <bibitem*|ABH+15><label|bib-alnaes_fenics_2015>Martin<nbsp>S.<nbsp>Alnaes,
      Jan Blechta, Johan Hake, August Johansson, Benjamin Kehlet, Anders
      Logg, Chris Richardson, Johannes Ring,
      Marie<nbsp>E.<nbsp>Rognes<localize|, and
      >Garth<nbsp>N.<nbsp>Wells.<newblock> The FEniCS Project Version
      1.5.<newblock> <with|font-shape|italic|Archive of Numerical Software>,
      3(100), 2015.<newblock>

      <bibitem*|BS05><label|bib-brenner_c0_2005>Susanne<nbsp>C.<nbsp>Brenner<localize|
      and >Li-Yeng Sung.<newblock> C0 Interior Penalty Methods for Fourth
      Order Elliptic Boundary Value Problems on Polygonal Domains.<newblock>
      <with|font-shape|italic|Journal of Scientific Computing>,
      22-23(1-3):83\U118, jun 2005.<newblock>

      <bibitem*|BS08><label|bib-brenner_mathematical_2008>Susanne<nbsp>C.<nbsp>Brenner<localize|
      and >L.<nbsp>Ridgway Scott.<newblock> <with|font-shape|italic|The
      mathematical theory of finite element methods>.<newblock>
      <localize|Number><nbsp>15<localize| in >Texts in Applied Mathematics.
      Springer New York, New York, NY, 3<localize| edition>, 2008.<newblock>

      <bibitem*|GRS07><label|bib-grossmann_numerical_2007>Christian
      Grossmann, Hans-GŲrg Roos<localize|, and >Martin Stynes.<newblock>
      <with|font-shape|italic|Numerical treatment of partial differential
      equations>.<newblock> Universitext. Springer Berlin Heidelberg, Berlin,
      Heidelberg, 2007.<newblock>

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

      <bibitem*|Qua09><label|bib-quarteroni_numerical_2009>Alfio
      Quarteroni.<newblock> <with|font-shape|italic|Numerical Models for
      Differential Problems>.<newblock> <localize|Number><nbsp>8<localize| in
      >Modeling, Simulation and Applications. Springer, 1<localize| edition>,
      2009.<newblock>

      <bibitem*|Sol06><label|bib-solin_partial_2006>Pavel Solin.<newblock>
      <with|font-shape|italic|Partial differential equations and the finite
      element method>.<newblock> Pure and applied Mathematics. John Wiley &
      Sons, 2006.<newblock> DOI: 10.1002/0471764108.<newblock>

      <bibitem*|ōW10><label|bib-olgaard_optimizations_2010>Kristian<nbsp>B.<nbsp>ōlgaard<localize|
      and >Garth<nbsp>N.<nbsp>Wells.<newblock> Optimizations for Quadrature
      Representations of Finite Element Tensors Through Automated Code
      Generation.<newblock> <with|font-shape|italic|ACM Trans. Math. Softw.>,
      37(1):8\U1, jan 2010.<newblock>
    </bib-list>
  </bibliography>
</body>

<\initial>
  <\collection>
    <associate|font|stix>
    <associate|font-base-size|11>
    <associate|info-flag|detailed>
    <associate|math-font|math-stix>
    <associate|page-medium|paper>
    <associate|page-screen-margin|false>
    <associate|preamble|false>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|2>>
    <associate|auto-10|<tuple|3.2|6>>
    <associate|auto-11|<tuple|3.3|7>>
    <associate|auto-12|<tuple|4|7>>
    <associate|auto-13|<tuple|4.1|7>>
    <associate|auto-14|<tuple|4.2|8>>
    <associate|auto-15|<tuple|5|8>>
    <associate|auto-16|<tuple|5.1|8>>
    <associate|auto-17|<tuple|5.2|8>>
    <associate|auto-18|<tuple|1|9>>
    <associate|auto-19|<tuple|2|9>>
    <associate|auto-2|<tuple|2|3>>
    <associate|auto-20|<tuple|5.3|9>>
    <associate|auto-21|<tuple|5.4|9>>
    <associate|auto-22|<tuple|5.5|9>>
    <associate|auto-23|<tuple|3|10>>
    <associate|auto-24|<tuple|5.6|10>>
    <associate|auto-25|<tuple|4|10>>
    <associate|auto-26|<tuple|5.7|10>>
    <associate|auto-27|<tuple|5|11>>
    <associate|auto-28|<tuple|6|11>>
    <associate|auto-29|<tuple|6|12>>
    <associate|auto-3|<tuple|2.1|3>>
    <associate|auto-30|<tuple|7|12>>
    <associate|auto-31|<tuple|8|12>>
    <associate|auto-32|<tuple|A|12>>
    <associate|auto-33|<tuple|A.1|13>>
    <associate|auto-34|<tuple|A.2|14>>
    <associate|auto-35|<tuple|A.3|15>>
    <associate|auto-36|<tuple|A.4|16>>
    <associate|auto-37|<tuple|A.4|17>>
    <associate|auto-4|<tuple|1|3>>
    <associate|auto-5|<tuple|2|4>>
    <associate|auto-6|<tuple|2.2|4>>
    <associate|auto-7|<tuple|2.3|6>>
    <associate|auto-8|<tuple|3|6>>
    <associate|auto-9|<tuple|3.1|6>>
    <associate|bib-alnaes_fenics_2015|<tuple|ABH+15|17>>
    <associate|bib-brenner_c0_2005|<tuple|BS05|17>>
    <associate|bib-brenner_mathematical_2008|<tuple|BS08|17>>
    <associate|bib-girault_hermite_2002|<tuple|GS02|17>>
    <associate|bib-grossmann_numerical_2007|<tuple|GRS07|17>>
    <associate|bib-logg_automated_2012|<tuple|LMW12|17>>
    <associate|bib-olgaard_optimizations_2010|<tuple|ōW10|17>>
    <associate|bib-quarteroni_numerical_2009|<tuple|Qua09|17>>
    <associate|bib-solin_partial_2006|<tuple|Sol06|17>>
    <associate|eq:delta-property|<tuple|3|3>>
    <associate|eq:euler-bernoulli|<tuple|7|8>>
    <associate|eq:hermite-first-derivatives|<tuple|5|5>>
    <associate|eq:hermite-second-derivatives|<tuple|6|5>>
    <associate|eq:hermite-transform|<tuple|4|5>>
    <associate|eq:linear-problem|<tuple|1|2>>
    <associate|eq:linear-system|<tuple|2|2>>
    <associate|fig:reference-triangle|<tuple|1|3>>
    <associate|footnote-1|<tuple|1|2>>
    <associate|footnote-2|<tuple|2|2>>
    <associate|footnote-3|<tuple|3|2>>
    <associate|footnote-4|<tuple|4|2>>
    <associate|footnote-5|<tuple|5|2>>
    <associate|footnote-6|<tuple|6|3>>
    <associate|footnote-7|<tuple|7|6>>
    <associate|footnote-8|<tuple|8|7>>
    <associate|footnote-9|<tuple|9|7>>
    <associate|footnr-1|<tuple|1|2>>
    <associate|footnr-2|<tuple|2|2>>
    <associate|footnr-3|<tuple|3|2>>
    <associate|footnr-4|<tuple|4|2>>
    <associate|footnr-5|<tuple|5|2>>
    <associate|footnr-6|<tuple|6|3>>
    <associate|footnr-7|<tuple|7|6>>
    <associate|footnr-8|<tuple|8|7>>
    <associate|footnr-9|<tuple|9|7>>
    <associate|sec:H2-projection|<tuple|4.2|8>>
    <associate|sec:ad-nodal-interpolation|<tuple|A.4|16>>
    <associate|sec:bcs-implementation|<tuple|A.3|15>>
    <associate|sec:fiat-implementation|<tuple|A.1|13>>
    <associate|sec:nodal-interpolant|<tuple|4.1|7>>
    <associate|sec:shape-implementation|<tuple|A.2|14>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      alnaes_fenics_2015

      logg_automated_2012

      brenner_mathematical_2008

      grossmann_numerical_2007

      logg_automated_2012

      brenner_mathematical_2008

      brenner_mathematical_2008

      solin_partial_2006

      quarteroni_numerical_2009

      olgaard_optimizations_2010

      brenner_mathematical_2008

      girault_hermite_2002

      solin_partial_2006

      solin_partial_2006

      brenner_mathematical_2008

      solin_partial_2006

      solin_partial_2006

      brenner_c0_2005
    </associate>
    <\associate|figure>
      <tuple|normal|<with|color|<quote|dark
      red>|<datoms|<macro|x|<resize|<tabular|<tformat|<cwith|1|1|1|1|cell-background|pastel
      red>|<cwith|1|1|1|1|cell-lsep|0fn>|<cwith|1|1|1|1|cell-rsep|0fn>|<cwith|1|1|1|1|cell-bsep|<value|marked-padding>>|<cwith|1|1|1|1|cell-tsep|<value|marked-padding>>|<table|<row|<cell|<arg|x>>>>>>|<plus|1l|0fn>|<plus|1b|<value|marked-padding>>|<minus|1r|0fn>|<minus|1t|<value|marked-padding>>>>|[The
      reference triangle and the mapping <with|mode|<quote|math>|\<b-x\>>...]>>|<pageref|auto-4>>

      <tuple|normal|The ten Hermite shape functions over the reference
      triangle. Here <with|mode|<quote|math>|<wide|\<varphi\>|^><rsub|0>> is
      depicted last.|<pageref|auto-5>>

      <tuple|normal|Clamped steel beam under the effect of gravity
      alone.|<pageref|auto-23>>

      <tuple|normal|Cantilevered steel beam with a weight hanging from the
      right endpoint.|<pageref|auto-25>>

      <tuple|normal|Simply supported steel beam under the effect of gravity
      alone.|<pageref|auto-27>>

      <tuple|normal|A visual inspection of the solutions isn't very
      revealing...|<pageref|auto-29>>

      <tuple|normal|The relative difference between the Lagrange and Hermite
      solutions decreases [...] as the grid resolution
      increases.|<pageref|auto-30>>

      <tuple|normal|Enforcing essential boundary conditions after assembly
      requires manipulation of the system matrix, invalidating its sparsity
      pattern. The cost of a naive <with|font-shape|<quote|italic>|a
      posteriori> modification in <with|font-shape|<quote|small-caps>|Python>
      is exponential in the number of degrees of freedom.|<pageref|auto-31>>
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

      <with|par-left|<quote|1tab>|2.2.<space|2spc>The local basis over a
      physical element <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|<quote|1tab>|2.3.<space|2spc>The global basis for
      <with|mode|<quote|math>|V<rsub|h>> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3.<space|2spc>Construction
      of the linear system> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8><vspace|0.5fn>

      <with|par-left|<quote|1tab>|3.1.<space|2spc>Assembly of the stiffness
      matrix <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9>>

      <with|par-left|<quote|1tab>|3.2.<space|2spc>Imposing boundary
      conditions <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10>>

      <with|par-left|<quote|1tab>|3.3.<space|2spc>Computing integrals
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|4.<space|2spc>Interpolation>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12><vspace|0.5fn>

      <with|par-left|<quote|1tab>|4.1.<space|2spc>The local nodal interpolant
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13>>

      <with|par-left|<quote|1tab>|4.2.<space|2spc><with|mode|<quote|math>|H<rsup|2>>-projection
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-14>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|5.<space|2spc>The
      Euler-Bernoulli beam model> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15><vspace|0.5fn>

      <with|par-left|<quote|1tab>|5.1.<space|2spc>Derivation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-16>>

      <with|par-left|<quote|1tab>|5.2.<space|2spc>Weak formulation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-17>>

      <with|par-left|<quote|4tab>|Essential boundary conditions
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-18><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Natural boundary conditions
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-19><vspace|0.15fn>>

      <with|par-left|<quote|1tab>|5.3.<space|2spc>Existence and uniqueness
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-20>>

      <with|par-left|<quote|1tab>|5.4.<space|2spc>Discretization
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-21>>

      <with|par-left|<quote|1tab>|5.5.<space|2spc>Clamped beam
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-22>>

      <with|par-left|<quote|1tab>|5.6.<space|2spc>Cantilevered beam
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-24>>

      <with|par-left|<quote|1tab>|5.7.<space|2spc>Simply supported beam
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-26>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|6.<space|2spc>The
      biharmonic equation> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-28><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Appendix
      A.<space|2spc>Implementation details>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-32><vspace|0.5fn>

      <with|par-left|<quote|1tab>|A.1.<space|2spc><with|font-shape|<quote|small-caps>|FIAT>
      implementation <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-33>>

      <with|par-left|<quote|1tab>|A.2.<space|2spc>Computing the Hermite shape
      functions <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-34>>

      <with|par-left|<quote|1tab>|A.3.<space|2spc>Essential boundary
      conditions for Hermite elements <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-35>>

      <with|par-left|<quote|1tab>|A.4.<space|2spc>Automatic differentiation
      for nodal interpolation <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-36>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-37><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>