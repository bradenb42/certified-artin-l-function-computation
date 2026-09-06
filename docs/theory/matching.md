# Matching the $\ell$-adic numbering with the global numbering

**Setting.** Global numbering $\alpha\_1,\dots,\alpha\_n$ (at $p$), Galois group $G\le S\_n$ in it, certificate chain $S\_n=G\_0>G\_1>\dots>G\_t=G$ with invariants $F\_i\in\mathbb{Z}[x\_1,\dots,x\_n]$ (Tschirnhaus transforms folded in), $\mathrm{Stab}*{G*{i-1}}(F\_i)=G\_i$, indices $m\_i=[G\_{i-1}\:G\_i]$, recorded integers $c\_i=F\_i(\alpha)$, and the recorded fact that the $m\_i$ values $F\_i(\tau\alpha)$, $\tau\in G\_{i-1}/G\_i$, are pairwise distinct. Local numbering $\beta\_1,\dots,\beta\_n\in L$ from [local Galois groups](https://claude.ai/chat/local-galois-groups.md), with local group $D\le S\_n$ in that numbering ($\phi(\beta\_i)=\beta\_{d(i)}$ for $d\in D$), inertia $I\trianglelefteq D$ and the higher ramification groups likewise in the $\beta$-numbering. Conventions: $(\rho P)(x)=P(x\_{\rho(1)},\dots,x\_{\rho(n)})$, a left action; for a tuple $\beta$ and $\rho\in S\_n$, $\beta^\rho:=(\beta\_{\rho(1)},\dots,\beta\_{\rho(n)})$, so $F(\beta^\rho)=(\rho F)(\beta)$; $g(P(\alpha))=(gP)(\alpha)$ for $g\in G$.

## 1. The matching problem

**Definition.** For an embedding $\iota\:N\to\bar{\mathbb{Q}}*\ell$ let $\tau*\iota\in S\_n$ be defined by $\beta\_i=\iota(\alpha\_{\tau\_\iota(i)})$. A permutation $\tau$ is a *matching* if $\tau=\tau\_\iota$ for some $\iota$.

**Lemma 1.1.** (a) The matchings form the left coset $G\tau\_0$ of any one matching $\tau\_0$. (b) For a matching $\tau$, $\tau D\tau^{-1}=D\_\iota\le G$, the decomposition group in the global numbering of the prime of $N$ determined by $\iota$; likewise $\tau I\tau^{-1}=I\_\iota$ and $\tau G\_{\mathfrak{L},u}\tau^{-1}$ are the inertia and higher ramification groups in the global numbering. (c) For every polynomial $F$ and every $\rho\in S\_n$, $F(\beta^\rho)=\iota\big((\tau\rho F)(\alpha)\big)$; in particular $\beta^\rho$ is the numbering induced by an embedding iff $\tau\rho\in G$.

*Proof.* (a) $\iota'=\iota\circ g$ gives $\beta\_i=\iota'(g^{-1}\alpha\_{\tau(i)})=\iota'(\alpha\_{g^{-1}\tau(i)})$, so $\tau\_{\iota'}=g^{-1}\tau$; all embeddings are of this form. (b) For $d\in D$ let $\phi$ be the corresponding automorphism of $L$ and $g\in G$ with $\phi\circ\iota=\iota\circ g$ (exists since $\iota(N)$ is stable under $\phi$). Then $\beta\_{d(i)}=\phi(\beta\_i)=\phi\iota(\alpha\_{\tau(i)})=\iota(\alpha\_{g\tau(i)})$, so $\tau d=g\tau$ and $g=\tau d\tau^{-1}$; the map $d\mapsto\tau d\tau^{-1}$ is the isomorphism $D\to D\_\iota$ that carries $I$ and the filtration along. (c) $F(\beta^\rho)=F(\beta\_{\rho(1)},\dots)=F(\iota\alpha\_{\tau\rho(1)},\dots)=\iota((\tau\rho F)(\alpha))$; and $\beta^\rho\_i=\iota(\alpha\_{\tau\rho(i)})$ is $\iota g(\alpha\_i)$ for $g=\tau\rho\in G$, while for $\tau\rho\notin G$ it is not induced by any embedding by (a). $\square$

So the problem is: find $\rho\in S\_n$ with $\tau\rho\in G$ (equivalently $\rho\in\tau^{-1}G$, a right coset), without knowing $\tau$. Any such $\rho$ gives $\rho^{-1}D\rho=(\tau\rho)^{-1}(\tau D\tau^{-1})(\tau\rho)\le G$, which is $D\_\ell$ up to conjugacy in $G$, and the outputs of [local Galois groups](https://claude.ai/chat/local-galois-groups.md) are transported by $\rho$. In the normalization "$\tau D\tau^{-1}\le G$" the matchings are the $\tau$ above; the output $\rho$ of the search satisfies $\rho^{-1}=g^{-1}\tau$ for some $g\in G$, so $\rho^{-1}$ is itself a matching.

**Proposition 1.2 (characterization).** $\tau\rho\in G$ iff $F\_i(\beta^\rho)=\iota(c\_i)$ for $i=1,\dots,t$. More precisely, if $\tau\rho\in G\_{i-1}$ then $\tau\rho\in G\_i$ iff $F\_i(\beta^\rho)=\iota(c\_i)$.

*Proof.* By Lemma 1.1(c), $F\_i(\beta^\rho)=\iota(F\_i((\tau\rho)\alpha))$, and [the precision policy](https://claude.ai/chat/precision-policy.md) Proposition 4.1 with $\rho\_{\text{there}}=\tau\rho$: if $\tau\rho\in G\_{i-1}$, then $F\_i((\tau\rho)\alpha)=c\_i=F\_i(\alpha)$ iff $\tau\rho G\_i=G\_i$, by the recorded distinctness over $G\_{i-1}/G\_i$. Induction from $G\_0=S\_n$. $\square$

## 2. The search along the chain

The candidate space is the set of cosets $\rho G$ (classes of $\rho$ under right multiplication by $G$: $\rho$ and $\rho g$ are both solutions or both not), i.e. $S\_n/G$. It is fibered by the chain:

$$S\_n/G=S\_n/G\_t\ \to\ S\_n/G\_{t-1}\ \to\ \dots\ \to\ S\_n/G\_1\ \to\ S\_n/G\_0={\ast},$$

each fiber of $S\_n/G\_i\to S\_n/G\_{i-1}$ having $m\_i$ elements, indexed by a transversal $T\_i$ of $G\_i$ in $G\_{i-1}$ ($\rho G\_{i-1}\mapsto{\rho\tau'G\_i:\tau'\in T\_i}$). The true class $\tau^{-1}G$ projects to a unique element at each level, and Proposition 1.2 recognizes it. The search descends one fiber at a time.

**Algorithm 2.1.**
Input: $\beta$ at $\ell$-adic precision $k\_\ell$ ([the precision policy](https://claude.ai/chat/precision-policy.md) §5), $D$ with generators $d\_1,\dots,d\_b$, the chain data.
$\pi\leftarrow\mathrm{id}$.
For $i=1,\dots,t$:

1. (structural filter) $C\_i\leftarrow{\tau'\in T\_i:\ (\pi\tau')^{-1}d\_k(\pi\tau')\in G\_i\text{ for }k=1,\dots,b}$.
2. (invariant test) $S\_i\leftarrow{\tau'\in C\_i:\ F\_i(\beta^{\pi\tau'})\equiv c\_i\pmod{\ell^{k\_\ell}}}$, the congruence meaning that the difference has valuation $\ge k\_\ell$ in $L$.
3. If $S\_i={\tau'}$, $\pi\leftarrow\pi\tau'$; if $|S\_i|>1$, branch (retain all); if $S\_i=\emptyset$, fail. Output $\rho=\pi$ (or the set of surviving $\pi$'s).

Step 1 uses only permutation arithmetic (sifting $(\pi\tau')^{-1}d\_k(\pi\tau')$ through a base and strong generating set of $G\_i$); step 2 is $|C\_i|$ evaluations of $F\_i$ in $\mathcal{O}*L/\ell^{k*\ell}$.

**Theorem 2.2 (no true matching is discarded).** Let $\pi$ satisfy $\tau\pi\in G\_{i-1}$ (true for $\pi=\mathrm{id}$ and $i=1$). Then the unique $\tau'*\ast\in T\_i$ with $\tau\pi\tau'*\ast\in G\_i$ satisfies $\tau'*\ast\in C\_i$ and $\tau'*\ast\in S\_i$, at every precision $k\ge1$, not only at $k\_\ell$. Hence, by induction, the branch of $\pi$'s containing the true class is never pruned: at each level the set of retained $\pi$'s contains one with $\tau\pi\in G\_i$, and after level $t$ it contains a $\rho$ with $\tau\rho\in G$.

*Proof.* Existence and uniqueness of $\tau'*\ast$: $\tau\pi\in G*{i-1}$ and $T\_i$ is a transversal of $G\_i$ in $G\_{i-1}$, so exactly one $\tau'$ has $\tau\pi\tau'\in G\_i$ (namely $\tau'G\_i=(\tau\pi)^{-1}G\_i$). Filter: $\tau D\tau^{-1}\le G\le G\_i$ (Lemma 1.1(b)), so $(\pi\tau'*\ast)^{-1}D(\pi\tau'*\ast)=(\tau\pi\tau'*\ast)^{-1}(\tau D\tau^{-1})(\tau\pi\tau'*\ast)\le(\tau\pi\tau'*\ast)^{-1}G\_i(\tau\pi\tau'*\ast)=G\_i$ since $\tau\pi\tau'*\ast\in G\_i$; so every generator passes. Invariant: by Proposition 1.2, $F\_i(\beta^{\pi\tau'*\ast})=\iota(c\_i)$ exactly in $L$; an exact equality is a congruence at every precision. $\square$

**Theorem 2.3 (completeness at the** [**the precision policy**](https://claude.ai/chat/precision-policy.md) **precision).** With $k\_\ell$ from [the precision policy](https://claude.ai/chat/precision-policy.md) §5 and $(F\_i,m\_i)\in\mathcal{F}$: for every $\tau'\in T\_i$ with $\tau'\ne\tau'*\ast$, $F\_i(\beta^{\pi\tau'})\not\equiv c\_i\pmod{\ell^{k*\ell}}$. Hence $S\_i={\tau'*\ast}$ at every level, no branching occurs, and the output $\rho=\tau'*{\ast,1}\cdots\tau'*{\ast,t}$ satisfies $\tau\rho\in G$; consequently $\rho^{-1}D\rho\le G$ is the decomposition group in the global numbering up to $G$-conjugacy, and $\rho^{-1}I\rho$, $\rho^{-1}G*{\mathfrak{L},u}\rho$ are the inertia and higher ramification groups.

*Proof.* Put $\gamma=F\_i((\tau\pi\tau')\alpha)-c\_i\in N$. Since $\tau\pi\tau'\in G\_{i-1}$ for every $\tau'\in T\_i$ (true or not), the conjugates $F\_i(g\tau\pi\tau'\alpha)-c\_i$, $g\in G$, depend only on the coset $g\tau\pi\tau'G\_i\subseteq G\_{i-1}$, so $\gamma$ has degree $\le m\_i$ over $\mathbb{Q}$ and all its conjugates are bounded by $2B\_{F\_i}$ ([the precision policy](https://claude.ai/chat/precision-policy.md) Lemma 1.1). It is nonzero for $\tau'\ne\tau'*\ast$ by Proposition 1.2.* [*the precision policy*](https://claude.ai/chat/precision-policy.md) *Lemma 2.1 gives $v(\gamma)\le m\_i\log*\ell(2B\_{F\_i})\<k\_\ell$, because $\ell^{k\_\ell}>M\_{F\_i}\ge(2B\_{F\_i})^{m\_i}$. The residue of $F\_i(\beta^{\pi\tau'})$ at precision $k\_\ell$ is exact (polynomial with integer coefficients in integral quantities known to that precision), so the congruence fails. The final statements are Lemma 1.1(b) transported by $\rho$: $\tau\rho=\:g\in G$, so $\rho^{-1}D\rho=g^{-1}(\tau D\tau^{-1})g$. $\square$

At an unramified prime with $\ell$ outside the excluded sets $E\_i$ of [the direct route](https://claude.ai/chat/direct-route.md), the same holds with $k\_\ell$ replaced by $1$ (residue field arithmetic in $\mathbb{F}\_{\ell^r}$), by [the direct route](https://claude.ai/chat/direct-route.md) Theorem 3.1; at ramified primes the [the precision policy](https://claude.ai/chat/precision-policy.md) precision is required in general because $\ell\mid\mathrm{disc}f$ and the resolvent roots $F\_i(\tau\alpha)$ may be congruent modulo $\mathfrak{L}$.

**Certificate of the output.** $\rho$ together with the checks $\rho^{-1}d\_k\rho\in G$ for all generators $d\_k$, and $F\_i(\beta^{\rho})\equiv c\_i$ for all $i$ at precision $k\_\ell$, is a certificate that $\rho$ is a matching; the first check is the target relation $\rho^{-1}D\rho\le G$ itself, the second is Theorem 2.3's characterization.

## 3. Bounds on the number of candidates

**Proposition 3.1 (without the structural filter).** At the [the precision policy](https://claude.ai/chat/precision-policy.md) precision the number of invariant evaluations is exactly $\sum\_{i=1}^{t}|T\_i|=\sum\_im\_i$ if every fiber is scanned, and at most that with early exit after the first success (expected $\sum\_i(m\_i+1)/2$ under a uniform position of $\tau'\_\ast$ in $T\_i$). By contrast the candidate space has $\prod\_im\_i=[S\_n\:G]$ classes, so the chain replaces a product by a sum.

**Proposition 3.2 (with the structural filter).** Let $D^{(i-1)}=\pi^{-1}D\pi\le G\_{i-1}$ be the local group in the numbering reached at level $i-1$. The number of invariant evaluations at level $i$ is

$$|C\_i|=\mathrm{fix}*{G*{i-1}/G\_i}\big(D^{(i-1)}\big)=#{\tau'G\_i\in G\_{i-1}/G\_i:\ D^{(i-1)}\le\tau'G\_i\tau'^{-1}}=\frac{#{\tau'\in G\_{i-1}:\tau'^{-1}D^{(i-1)}\tau'\le G\_i}}{|G\_i|}\ \le\ m\_i,$$

the mark of $D^{(i-1)}$ on $G\_{i-1}/G\_i$. When $D$ is cyclic, $D^{(i-1)}=\langle d\rangle$ (every unramified prime; at tamely ramified primes $D$ is metacyclic and the general formula applies with two generators), this is

$$|C\_i|=\frac{|C\_{G\_{i-1}}(d)|\cdot|d^{,G\_{i-1}}\cap G\_i|}{|G\_i|}=m\_i,\frac{|d^{,G\_{i-1}}\cap G\_i|}{|d^{,G\_{i-1}}|},$$

the fixed-point formula of [rational classes](https://claude.ai/chat/rational-classes.md) Lemma 1.1 in the group $G\_{i-1}$. The total is $\sum\_i\mathrm{fix}*{G*{i-1}/G\_i}(D^{(i-1)})$, which depends on $\ell$ only through the conjugacy class of $D\_\ell$ in $G$.

*Proof.* $(\pi\tau')^{-1}d\_k(\pi\tau')\in G\_i$ for all $k$ iff $\tau'^{-1}D^{(i-1)}\tau'\le G\_i$ iff $D^{(i-1)}$ fixes the coset $\tau'G\_i$; the count of fixed cosets is the number of $\tau'\in G\_{i-1}$ with the property divided by $|G\_i|$. The cyclic case is [rational classes](https://claude.ai/chat/rational-classes.md) Lemma 1.1. Independence of $\ell$: $D^{(i-1)}$ is conjugate in $G\_{i-1}$ to $\tau'*{\ast,i-1}\cdots$-transports of $D*\ell$, and marks are conjugation invariant. $\square$

The filter is decisive at the first level, where $m\_1=[S\_n\:G\_1]$ can be large (a primitive maximal $G\_1$, e.g. $[S\_{12}:\mathrm{PGL}*2(11)]=362880$): for $D=\langle d\rangle$ the count $|C\_1|=|C*{S\_n}(d)|,|d^{S\_n}\cap G\_1|/|G\_1|$ is the number of conjugates of $G\_1$ containing $d$, typically a small integer (for $G\_1=M\_{12}$ and $d$ an $11$-cycle, $|C\_1|=11\cdot17280/95040=2$; for $G\_1=A\_n$, $|C\_1|\in{0,2}$). Enumerating $C\_i$ costs $m\_i$ sifts of a permutation through the BSGS of $G\_i$, $O(m\_i,n^2)$ elementary operations, negligible against an invariant evaluation in $\mathcal{O}*L/\ell^{k*\ell}$; when $m\_i$ is very large, $C\_i$ is instead generated directly as the set of conjugates of $G\_i$ in $G\_{i-1}$ containing $D^{(i-1)}$ by a backtrack that places the generators $d\_k$ into $G\_i$, whose cost is bounded by $|C\_i|$ times the cost of one conjugacy-into-subgroup search.

**Proposition 3.3 (branching without separation).** If the test is run at a precision $k$ below the [the precision policy](https://claude.ai/chat/precision-policy.md) bound, Theorem 2.2 still guarantees that the true branch survives, and the number of $\pi$'s retained after level $i$ is at most $\prod\_{j\le i}s\_j$ where $s\_j$ is the number of $\tau'\in C\_j$ whose value $F\_j(\beta^{\pi\tau'})$ is congruent to $c\_j$ modulo $\ell^{k}$; every false survivor is eventually eliminated at a later level or by the final certificate $\rho^{-1}D\rho\le G$ failing, and raising the precision to $k\_\ell$ collapses every $s\_j$ to $1$. The cost of running at the [the precision policy](https://claude.ai/chat/precision-policy.md) precision from the start is therefore never worse than the branching search plus verification, and the branching variant is useful only when $k\_\ell$ digits are unavailable for the roots.

*Proof.* Retention count is the product of fiber survivors. A false final $\rho$ has $\tau\rho\notin G$, so its certificate fails: either some later exact test fails, or $\rho^{-1}D\rho\not\le G$; the latter is not guaranteed to fail (some $\rho$ with $\tau\rho\notin G$ may still conjugate $D$ into $G$), which is why the invariant checks at full precision are part of the certificate. $\square$

## 4. Cost and placement

Per prime, after [local Galois groups](https://claude.ai/chat/local-galois-groups.md): $\sum\_i|C\_i|$ evaluations of $F\_i$ at cost $c\_i$ operations each in $\mathcal{O}*L/\ell^{k*\ell}$, i.e. $\tilde O\big(\sum\_i|C\_i|,c\_i\cdot|D|,k\_\ell\log\ell\big)$ bit operations, plus $O(\sum\_im\_in^2)$ for the filters and $O(bn)$ for the final membership certificate. The evaluations are at the roots of [local Galois groups](https://claude.ai/chat/local-galois-groups.md) lifted to precision $k\_\ell$ by Hensel in the tower (cost $\tilde O(n|D|k\_\ell\log\ell)$).

The matching is computed once per ramified prime (giving $D\_\ell,I\_\ell$ and the filtration in both numberings) and, by [the direct route](https://claude.ai/chat/direct-route.md), once per unramified prime on which the direct route is taken (giving $\mathrm{Frob}*\ell$ exactly, with $\mathbb{F}*{\ell^r}$ in place of $\mathcal{O}*L/\ell^{k*\ell}$ and $D=\langle\phi\_\beta\rangle$). In both uses the same chain, transversals $T\_i$ and invariants $F\_i$ are reused, and the structural filter is applied with the local group available at that prime.