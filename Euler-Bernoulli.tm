<TeXmacs|1.99.4>

<style|generic>

<\body>
  <chapter|The Euler-Bernoulli beam model>

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

  Conservation of momentum and linear constitutive relations yield (see e.g.
  <inactive|<cite-detail|>>):

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

  <paragraph|Essential boundary conditions:>Split the boundary
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

  <paragraph|Natural boundary conditions:>For
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
  which requires at least cubic polynomials ... [elaborate, see p.218 of
  <inactive|<cite|solin>>]

  So we use cubic Hermite elements, which are <math|H<rsup|2>> conforming in
  <math|<with|math-font|Bbb|R>>.
</body>

<\initial>
  <\collection>
    <associate|item-vsep|<macro|0>>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?|../../../.TeXmacs/texts/scratch/no_name_140.tm>>
    <associate|auto-2|<tuple|1|?|../../../.TeXmacs/texts/scratch/no_name_140.tm>>
    <associate|auto-3|<tuple|2|?|../../../.TeXmacs/texts/scratch/no_name_140.tm>>
    <associate|auto-4|<tuple|1|?|../../../.TeXmacs/texts/scratch/no_name_140.tm>>
    <associate|auto-5|<tuple|2|?|../../../.TeXmacs/texts/scratch/no_name_140.tm>>
    <associate|auto-6|<tuple|3|?|../../../.TeXmacs/texts/scratch/no_name_140.tm>>
    <associate|auto-7|<tuple|4|?|../../../.TeXmacs/texts/scratch/no_name_140.tm>>
    <associate|eq:euler-bernoulli|<tuple|1|?|../../../.TeXmacs/texts/scratch/no_name_140.tm>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|2fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|font-size|<quote|1.19>|1<space|2spc>The
      Euler-Bernoulli beam model> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|1fn>

      <with|par-left|<quote|1tab>|1<space|2spc>Derivation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1tab>|2<space|2spc>Weak formulation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|4tab>|Essential boundary conditions:
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Natural boundary conditions:
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5><vspace|0.15fn>>

      <with|par-left|<quote|1tab>|3<space|2spc>Existence and uniqueness
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|<quote|1tab>|4<space|2spc>Discretization
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>
    </associate>
  </collection>
</auxiliary>