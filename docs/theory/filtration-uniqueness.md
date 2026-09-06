# When the polygons do not determine the filtration, and what does

*One of the derivations behind* [*artin*](https://claude.ai/README.md)*; the* [*index*](https://claude.ai/chat/README.md) *lists them all and names the code that implements each.*

**Setting.** As in [the filtration](https://claude.ai/chat/ramification-filtration.md): $D=\mathrm{Gal}(L/\mathbb{Q}*\ell)$, $I=D\_0$, $P=D\_1$ the wild subgroup, $i\_L$, $\nu=v*\ell(\mathrm{disc}f)$, $e=e(L)$, $f=f(L)$, $\delta(L)=\sum\_{\sigma\ne1}i\_L(\sigma)\le e\nu$. The *polygon data* are the coset-sum functions $J\_H(\sigma)=\sum\_{\tau\in\sigma H}i\_L(\tau)$ for $H$ in a family $\mathcal{H}$ of subgroups, the stabilizers $H\_\beta$ of the roots (factor polygons) and the intersections $H\_\beta\cap H\_{\beta'}$ (pairwise compositum polygons), known by [the filtration](https://claude.ai/chat/ramification-filtration.md) Lemma 3.1 from the ramification polygons of $L^H$. A *filtration function* is $j\:D\to\mathbb{Z}*{\ge0}\cup{\infty}$ with $j(1)=\infty$, $j=0$ off $I$, $j=1$ on $I\setminus P$, and $D\_i(j):={j\ge i+1}$ ($i\ge1$) a decreasing chain of normal subgroups of $D$ inside $P$ with elementary abelian successive quotients; $j$ is determined by the chain $(D\_i(j))*{i\ge1}$, and $i\_L$ is one such function. Write $\mathcal{C}(\mathcal{H})$ for the set of filtration functions $j$ with $J^{(j)}\_H=J\_H$ for all $H\in\mathcal{H}$, where $J^{(j)}*H(\sigma)=\sum*{\tau\in\sigma H}j(\tau)$.

## 1. What a polygon sees

**Lemma 1.1.** For a filtration function $j$ with chain $D\_i=D\_i(j)$ and any $H\le D$,

$$J^{(j)}*H(\sigma)=\sum*{i\ge0}|H\cap D\_i|\cdot[\sigma\in D\_iH]\qquad(\sigma\notin H),$$

and for $\sigma\in P\setminus H$, with $W=H\cap P$,
$$J^{(j)}*H(\sigma)=|H\cap I|+\sum*{i\ge1}|W\cap D\_i|\cdot[\sigma\in D\_iW].$$

*Proof.* $j(\tau)=#{i\ge0:\tau\in D\_i}$, so $J^{(j)}\_H(\sigma)=\sum\_i|\sigma H\cap D\_i|$. Since $D\_i\trianglelefteq D$, $\sigma H\cap D\_i$ is empty unless $\sigma\in D\_iH$, and then it is a coset of $H\cap D\_i$. For $\sigma\in P$ and $i\ge1$: $D\_i\le P$, so by the modular law $D\_iH\cap P=D\_i(H\cap P)=D\_iW$ and $H\cap D\_i=W\cap D\_i$; the $i=0$ term is $|H\cap I|$ because $\sigma\in P\subseteq IH$. $\square$

So $J\_H$ on $P$ is a step function along the descending chain $P\supseteq D\_1W\supseteq D\_2W\supseteq\cdots\supseteq W$, with value $|H\cap I|+\sum\_{k\le i}|W\cap D\_k|$ on $D\_iW\setminus D\_{i+1}W$. Two consequences:

**Proposition 1.2 (only the wild part matters).** With $\mathcal{W}={H\cap P\:H\in\mathcal{H}}$, $\mathcal{C}(\mathcal{H})=\mathcal{C}(\mathcal{W})$: a candidate is consistent with the polygon of $L^H$ iff it is consistent with the polygon of $L^{H\cap P}$. Consistency with $H$ means: the chains $(D\_i(j),W)*{i\ge1}$ and $(D\_i,W)*{i\ge1}$ coincide as chains of subgroups of $P$ *and* the level sums $\sum\_{k\le i}|W\cap D\_k(j)|$ and $\sum\_{k\le i}|W\cap D\_k|$ coincide at every level.

*Proof.* Lemma 1.1: off $P$ every filtration function is $0$ or $1$ as prescribed, and on $P$ the formula involves $H$ only through $W$ and $|H\cap I|$, the latter being the same for $j$ and $i\_L$. $\square$

**Theorem 1.3 (exact description of the undetermined cases).**

1. $\mathcal{C}(\mathcal{H})$ is finite: its elements are chains of normal subgroups of $D$ inside $P$ with elementary abelian quotients, with jumps at positions $\le b\_{\max}\le e\nu$, subject to the equalities of Proposition 1.2.
2. If some $H\in\mathcal{H}$ satisfies $H\cap P=1$, then $\mathcal{C}(\mathcal{H})={i\_L}$.
3. In general the elements of $\mathcal{C}(\mathcal{H})$ agree on every coset sum $\sum\_{\tau\in\sigma W}j(\tau)$, $W\in\mathcal{W}$, and need not agree pointwise anywhere on $P\setminus1$ (Example 1.4: two candidates differ on all of $B\cup B'$ minus $1$).
4. Any two elements $j\ne j'$ of $\mathcal{C}(\mathcal{H})$ differ at some $\sigma\in P$ and, for every subgroup $S$ of the explicit list

$$\mathcal{S}:={S\le D:\ S\cap P=1}$$

(the subgroups whose intersection with inertia has order prime to $\ell$), they differ on the fixed field $L^S$: $J^{(j)}\_S\ne J^{(j')}\_S$. Indeed for $S\in\mathcal{S}$ and $\sigma\in P\setminus{1}$,
$$J^{(j)}\_S(\sigma)=|S\cap I|+j(\sigma)-1 .$$

*Proof.* 1. Finitely many chains, and $b\_{\max}+1\le\delta(L)\le e\nu$ ([the filtration](https://claude.ai/chat/ramification-filtration.md) Lemma 2.3) bounds the positions. 2. With $W=1$ the second formula of Lemma 1.1 reads $J\_H(\sigma)=|H\cap I|+j(\sigma)-1$ on $P$, so $j$ is recovered. 3. The coset-sum agreement is the definition of $\mathcal{C}(\mathcal{H})$; Example 1.4 shows pointwise disagreement. 4. $S\cap P=1$ gives $|S\cap D\_i|=1$ for $i\ge1$ and $D\_iS\cap P=D\_i$ (modular law), so $[\sigma\in D\_iS]=[\sigma\in D\_i]$ for $\sigma\in P$; substitute in Lemma 1.1. Injectivity of $j\mapsto J\_S^{(j)}$ on $P$ follows. $\square$

**Example 1.4 (both layers of ambiguity).** $P=D\_1\cong C\_\ell^2$, $D\_2=B$ a line, $D\_{b+1}=1$ for the jump $b\ge2$, $\mathcal{W}={W\_1,W\_2}$ two lines different from $B$ (two factors whose wild stabilizers are these lines). For a candidate with $D\_2(j)=B'$ any line $\ne W\_1,W\_2$ that is normal in $D$ (all lines, when $D$ acts trivially on $P$ by conjugation) and the same jumps, $D\_i(j)W\_k=P$ for $2\le i\le b$ and $|W\_k\cap D\_i(j)|=1$ there, exactly as for $B$; so $\mathcal{C}(\mathcal{H})$ contains all $\ell-1$ lines $B'\ne W\_1,W\_2$. All candidates have the same orders $|D\_i|$, hence the same Herbrand function, upper breaks and different; they differ in the subgroups, hence in $f\_\ell(\chi)$ for characters $\chi$ of $D$ nontrivial on $P$ but trivial on some line. Any $S\in\mathcal{S}$ separates them by Theorem 1.3(4); an *unlabeled* polygon of $L^S$ does not, since the multisets ${J\_S(\sigma)}$ coincide (§2).

## 2. The discriminating resolvent

**Choice of $S$.** $\mathcal{S}$ is nonempty ($1\in\mathcal{S}$) and contains every $\ell'$-Hall subgroup of the solvable group $D$ (index $|D|*\ell$). The index $[D\:S]$ is at least $|P|$ for every $S\in\mathcal{S}$ (as $|SP|=|S||P|\le|D|$), with equality iff $S$ is a complement of $P$ in $D$, which exists whenever $\ell\nmid[D\:P]$ (Schur to Zassenhaus), i.e. whenever $\ell\nmid f\cdot|I/P|$, i.e. $\ell\nmid f$. Take $S\in\mathcal{S}$ of minimal index (a maximal member under inclusion); then $L^S$ has degree $[D\:S]\in[|P|,|D|*\ell]$, much smaller than $|D|$ when the tame and unramified parts are large.

**Construction 2.1.** With $L$ and $D$ explicit ([local Galois groups](https://claude.ai/chat/local-galois-groups.md)): compute a uniformizer $\varpi\_S$ and a Teichmüller generator $\omega\_S$ of $L^S={x\in L:\ sx=x\ \forall s\in S}$ (fixed field by linear algebra over $\mathbb{Q}*\ell$ in the tower basis; $\varpi\_S$ as an element of $L^S$ of minimal positive $v\_L$, $\omega\_S$ as the Teichmüller lift of a generator of the residue field of $L^S$), set $y\_S=\varpi\_S+\omega\_S$, so that $\mathcal{O}*{L^S}=\mathbb{Z}\_\ell[y\_S]$, and form

$$R\_S(x)=\prod\_{\tau\in D/S}\big(x-\tau y\_S\big)\in\mathbb{Z}\_\ell[x],\qquad \deg R\_S=[D\:S],$$

the minimal polynomial of $y\_S$ over $\mathbb{Q}\_\ell$, i.e. the resolvent of the pair $(D,S)$ at the $\ell$-adic roots with the invariant chosen so that its value generates the ring of integers of the fixed field.

**Theorem 2.2 (what $R\_S$ selects).**

1. (Unlabeled: Newton polygon.) The Newton polygon of $R\_S(Y+y\_S)/Y$ over $L^S$ (computable inside $L^S$, without splitting $R\_S$) has slopes $-J\_S(\sigma)/e(L/L^S)$ in $v\_{L^S}$-units, i.e. $-J\_S(\sigma)$ in $v\_L$-units, with multiplicities. The cosets meeting $P$ are in bijection with $P$ (as $S\cap P=1$), the cosets in $IS\setminus PS$ carry the value $|S\cap I|$ and those outside $IS$ the value $0$; by Theorem 1.3(4) the values exceeding $|S\cap I|$ form the multiset ${|S\cap I|+i\_L(\sigma)-1:\ \sigma\in P\setminus1}$, which determines the orders $|D\_i|$ for all $i\ge1$, hence the Herbrand function $\varphi\_L$, all upper breaks, $\delta(L)$, and the sequence of jump positions. It selects, among the candidates of $\mathcal{C}(\mathcal{H})$, exactly those with the true order sequence $(|D\_i|)\_i$.
2. (Labeled: root valuations.) The roots $\tau y\_S$ ($\tau\in D/S$) are elements of $L$ known in the tower, and $v\_L(\tau y\_S-y\_S)=J\_S(\tau)$ for each coset. By Theorem 1.3(4), for $\tau\in P\setminus1$, $i\_L(\tau)=J\_S(\tau)-|S\cap I|+1$, so the labeled valuations determine $i\_L$ on $P$, hence the chain $(D\_i)$ itself: they select the unique true element of $\mathcal{C}(\mathcal{H})$.
3. Both computations are exact at $\ell$-adic precision $\nu+1$: every $J\_S(\sigma)$ with $\sigma\notin S$ satisfies $J\_S(\sigma)\le\delta(L)\le e\nu$, so the valuations in $v\_L$-units are $\<e(\nu+1)$, and the Newton polygon of $R\_S(Y+y\_S)/Y$ involves only valuations $\<e(\nu+1)$ in $v\_L$-units, i.e. it is determined by the coefficients of $R\_S(Y+y\_S)$ modulo $\ell^{\nu+1}$. This is the precision $k\_\ell\ge\nu+1$ of the policy of [the precision policy](https://claude.ai/chat/precision-policy.md) §5 ($\ell^{k\_\ell}>\Delta\ge\ell^{\nu}$), so no new precision requirement arises; $y\_S$ generating $\mathcal{O}\_{L^S}$ is certified by [the ramified primes](https://claude.ai/chat/ramified-primes.md)'s $\ell$-maximality test at the same precision.

*Proof.* 1 to 2. [the filtration](https://claude.ai/chat/ramification-filtration.md) Lemma 3.1 with $H=S$ and $y=y\_S$ (a generator of $\mathcal{O}*{L^S}$), and Theorem 1.3(4). The Newton polygon of $\prod*{\tau S\ne S}(Y-(\tau y\_S-y\_S))$ has slopes the negatives of the valuations of the roots; over $L^S$ the roots are not individually available but the polygon is, as it is determined by the valuations of the coefficients, which lie in $L^S$. 3. [the filtration](https://claude.ai/chat/ramification-filtration.md) Lemma 2.3 and the bound on coset sums by the total $\delta(L)$; the sum of all slopes is $v\_L(\mathfrak{D}*{L^S/\mathbb{Q}*\ell})\le\delta(L)$. $\square$

**Remark 2.3.** When $L$ has been constructed, the labeled computation is equivalent to [the filtration](https://claude.ai/chat/ramification-filtration.md) §2 ($i\_L$ from $\varpi\_L$), and $R\_S$ is a certificate rather than a discovery. The polygon route is the cheap one when $L$ itself is avoided: constructing $L^S$ alone (a field of degree $[D\:S]$, obtainable as the compositum of the fixed fields of $S$ inside the factor fields when $S$ is chosen inside a stabilizer, or as the splitting field of a suitable resolvent of degree $[D\:S]$ over $\mathbb{Q}\_\ell$ whose roots need only be found in $L^S$'s own Galois closure) yields, through its polygon, the complete numerical filtration (orders and breaks) and, through its labeled roots in the Galois closure of $L^S$, which is $L^{\mathrm{core}\_D(S)}$, of degree $|D|/|\mathrm{core}\_D(S)|$, the full filtration on $D/\mathrm{core}\_D(S)$; since $\mathrm{core}\_D(S)\cap P=1$, this is the full wild filtration.

## 3. The two mandatory structural checks

Every filtration accepted, whether from [the filtration](https://claude.ai/chat/ramification-filtration.md) §2, from Theorem 2.2, or from the candidate enumeration of Theorem 1.3, must pass the following two checks, and failure of either rejects the local computation ([local Galois groups](https://claude.ai/chat/local-galois-groups.md)) as a whole.

**Check 1 (Hasse to Arf on abelian subquotients).** For every pair $a\<b$ with $D\_b\trianglelefteq D\_a$ and $D\_a/D\_b$ abelian, in particular for every consecutive pair $(D\_i,D\_{i+1})$, whose quotient is elementary abelian, and for $(D\_0,D\_1)$, whose quotient is cyclic, the upper jumps of the abelian extension $L^{D\_b}/L^{D\_a}$ are integers. Concretely: the lower filtration of $\mathrm{Gal}(L/L^{D\_a})=D\_a$ is $(D\_a)*u=D\_a\cap D\_u=D*{\max(a,u)}$ (compatibility of lower numbering with subgroups), its Herbrand function is $\varphi\_a(u)=\int\_0^u dt/[D\_a\cap D\_0\:D\_a\cap D\_t]$, its upper filtration is $(D\_a)^v=(D\_a)\_{\psi\_a(v)}$, and the upper jumps of $D\_a/D\_b$ are the $v$ at which $(D\_a)^vD\_b/D\_b$ changes; each such $v$ must lie in $\mathbb{Z}$. For $a=0$ this is Hasse to Arf for the abelian quotients of $I$ itself. The check is a finite computation on the chain and its orders, and by [the filtration](https://claude.ai/chat/ramification-filtration.md) Theorem 4.2(1) it is a function of the computed data only. Its power: in Example 1.4 with $\ell\nmid$ nothing further, it does not separate the $\ell-1$ candidates (they share the orders), but it rejects every candidate order sequence whose implied upper jumps are non-integral, which eliminates most wrong jump positions in Theorem 1.3(1).

**Check 2 (inertia orbits against the local invariants).** For each irreducible factor $g\_j$ with invariants $(e\_j,f\_j)$ from [the ramified primes](https://claude.ai/chat/ramified-primes.md) and wild part $\ell^{a\_j},|,e\_j$: the orbits of $(D\_\ell)*0=I$ on the roots of $g\_j$ are $f\_j$ orbits of size $e\_j$; the orbits of $(D*\ell)*1=P$ on the roots of $g\_j$ have size $\ell^{a\_j}$; and the orbits of $D*\ell$ are the root sets of the $g\_j$. These are read off the permutation groups $I,P\le S\_n$ in either numbering and compared with the integers $e\_j,f\_j$ computed independently by the $\ell$-maximal-order algorithm of [the ramified primes](https://claude.ai/chat/ramified-primes.md); equality is required for every $j$.

*Proof of the required equalities.* $K\_j\otimes\mathbb{Q}\_\ell^{\mathrm{ur}}$ is a product of $f\_j$ totally ramified extensions of degree $e\_j$, so $\mathrm{Gal}(L/L^{\mathrm{ur}})=I$ has $f\_j$ orbits of size $e\_j$ on the roots of $g\_j$; over the maximal tamely ramified extension $T$, each totally ramified piece of degree $e\_j=e'\_j\ell^{a\_j}$ ($\ell\nmid e'\_j$) splits into $e'\_j$ totally wildly ramified pieces of degree $\ell^{a\_j}$, and $\mathrm{Gal}(L/L\cap T)=P$ has orbits of that size. $\square$

Both checks are also the natural certificates for the choices of §2: Check 2 is insensitive to the choice of $S$, and Check 1 must hold for the filtration read off $R\_S$; a filtration produced by Theorem 2.2 that fails either check signals an error in $L$, in $D$, or in the precision, not in the filtration theory.

## 4. Summary of the decision procedure for the filtration at $\ell$

1. From [local Galois groups](https://claude.ai/chat/local-galois-groups.md): $D$, $L$ in the tower, $\nu$.
2. From [the filtration](https://claude.ai/chat/ramification-filtration.md) §1: $I$, $P$, $\Phi$.
3. Either [the filtration](https://claude.ai/chat/ramification-filtration.md) §2 ($i\_L$ from $\varpi\_L$), or, when $L$ is avoided: the polygons of the factor fields and pairwise composita ($\mathcal{H}$), the candidate set $\mathcal{C}(\mathcal{H})$ of Theorem 1.3 filtered by Check 1, and, if $|\mathcal{C}(\mathcal{H})|>1$, the resolvent $R\_S$ of §2 for a minimal-index $S\in\mathcal{S}$: its Newton polygon fixes the orders, its labeled roots (in $L^{\mathrm{core}\_D(S)}$) fix the subgroups.
4. Check 1 and Check 2 on the result; then the Herbrand conversion of [the filtration](https://claude.ai/chat/ramification-filtration.md) §4 and the transport to $G$ of [matching](https://claude.ai/chat/matching.md) to [well-definedness](https://claude.ai/chat/well-definedness.md).

All steps run at $\ell$-adic precision $\nu+1\le k\_\ell$.