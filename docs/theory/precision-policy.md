# Precision management: one policy for every prime

**Setting.** $f\in\mathbb{Z}[x]$ monic separable of degree $n$, roots $\alpha\_1,\dots,\alpha\_n\in N\subset\bar{\mathbb{Q}}$, $\Delta=|\mathrm{disc}f|$, $G=\mathrm{Gal}(N/\mathbb{Q})\le S\_n$ with respect to the numbering fixed by the embedding $\iota\_p\:N\to\bar{\mathbb{Q}}*p$ at the numbering prime $p\nmid\Delta$. A root bound is any $R\ge1$ with $|\sigma(\alpha\_j)|\le R$ for all $j$ and all complex embeddings $\sigma$; Cauchy's $R=1+\max\_i|a\_i|$ for $f=x^n+\sum a\_ix^{n-i}$ is one. For $F\in\mathbb{Z}[x\_1,\dots,x\_n]$ write $|F|1$ for the sum of the absolute values of its coefficients and, for $\tau\in S\_n$, $F(\tau\alpha)=F(\alpha{\tau(1)},\dots,\alpha*{\tau(n)})$. If invariants are evaluated on Tschirnhaus-transformed roots $T(\alpha\_j)$, replace $F$ by $F\circ T$ throughout; $|F\circ T|\_1\le|F|\_1|T|\_1^{\deg F}$ and $\deg(F\circ T)=\deg F\deg T$.

Valuations at a prime $\ell$ are normalized by $v(\ell)=1$ and extended to $\bar{\mathbb{Q}}*\ell$; "precision $k$ at $\ell$" means every quantity is known modulo elements of valuation $\ge k$. In a field $L\supseteq\mathbb{Q}*\ell$ with ramification index $e$ this is $ek$ digits in a uniformizer.

## 1. Height bounds

**Lemma 1.1.** Let $F\in\mathbb{Z}[x\_1,\dots,x\_n]$, $d=\deg F$, $B\_F=|F|\_1R^{d}$. For every $\tau\in S\_n$ and every complex embedding $\sigma$, $|\sigma(F(\tau\alpha))|\le B\_F$, and $F(\tau\alpha)$ is an algebraic integer.

*Proof.* Each monomial of total degree $\le d$ evaluated at numbers of modulus $\le R\ge1$ has modulus $\le R^d$; sum over monomials with the coefficient moduli. Integrality: polynomial with integer coefficients in algebraic integers. $\square$

**Lemma 1.2 (resolvents).** Let $H\le G$, $m=[G\:H]$, and $F$ an $H$-invariant. The resolvent $\mathrm{Res}*{F,H}(x)=\prod*{\tau\in G/H}(x-F(\tau\alpha))$ lies in $\mathbb{Z}[x]$, its coefficient of $x^{m-j}$ has absolute value $\le\binom mjB\_F^{,j}$, so $|\mathrm{Res}*{F,H}|1\le(1+B\_F)^m$, and if it is squarefree, $1\le|\mathrm{disc},\mathrm{Res}{F,H}|\le(2B\_F)^{m(m-1)}$. The same holds for $\prod*{\tau\in G\_{i-1}/G\_i}(x-F\_i(\tau\alpha))$ for any group $G\_{i-1}\supseteq G$ containing the Galois group and any $G\_i$-invariant $F\_i$, with $m=[G\_{i-1}\:G\_i]$.

*Proof.* The multiset ${F(\tau\alpha)}$ is stable under $G$ because $F$ is $H$-invariant and $G$ permutes $G/H$; so the coefficients are rational and integral. The coefficient bound is the elementary symmetric function bound; the discriminant is a nonzero integer equal to $\prod\_{a\<b}(\beta\_a-\beta\_b)^2$ with $|\beta\_a-\beta\_b|\le2B\_F$. $\square$

The certificate records, for a chain $S\_n=G\_0>G\_1>\dots>G\_t=G$ with $G\_i$-invariants $F\_i$ and indices $m\_i=[G\_{i-1}\:G\_i]$, the integers $c\_i=F\_i(\alpha)$ and the fact that the $m\_i$ values $F\_i(\tau\alpha)$, $\tau\in G\_{i-1}/G\_i$, are pairwise distinct. By Lemma 1.1, $|c\_i|\le B\_{F\_i}$. The auxiliary resolvents used later are $\mathrm{Res}\_{F,H}$ for pairs $(F,H)$ with $H\le G$ chosen from $G$ alone; their data are $(B\_F,[G\:H])$.

## 2. The separation lemma and its transport to every prime

**Lemma 2.1 (separation).** Let $\gamma\ne0$ be an algebraic integer of degree $D$ over $\mathbb{Q}$ with $|\sigma(\gamma)|\le C$ for all complex embeddings $\sigma$. Then for every prime $\ell$ and every extension $v$ of $v\_\ell$ to $\bar{\mathbb{Q}}$,
$$v(\gamma)\le v\_\ell\big(N\_{\mathbb{Q}(\gamma)/\mathbb{Q}}(\gamma)\big)\le D\log\_\ell C.$$

*Proof.* $N(\gamma)=\prod\_\sigma\sigma(\gamma)$ over the $D$ embeddings of $\mathbb{Q}(\gamma)$; $v(N(\gamma))=\sum\_\sigma v(\sigma\gamma)$ where each $v\circ\sigma$ is a valuation extending $v\_\ell$, hence $\ge0$ on algebraic integers, and the identity embedding contributes $v(\gamma)$. $N(\gamma)$ is a nonzero integer with $|N(\gamma)|\le C^D$, so $v\_\ell(N(\gamma))\le\log\_\ell C^D$. $\square$

**Corollary 2.2 (transport).** Let $\gamma$ be as in Lemma 2.1 (or $\gamma=0$) and $k$ an integer with $\ell^k>C^D$. Then $\gamma\equiv0$ modulo valuation $\ge k$ iff $\gamma=0$. In particular a single integer $M\ge C^D$ fixes, at every prime $\ell$, the precision $k\_\ell(M)=\lfloor\log\_\ell M\rfloor+1$ at which the test is correct, and $k\_\ell(M)\log\_2\ell\le\log\_2M+\log\_2\ell$: the same bit budget at every prime, distributed into $\ell$-adic digits.

The tests that occur, and the pair $(C,D)$ for each:

- (T1) *Integer recovery.* $c\in\mathbb{Z}$ with $|c|\le B$ (an invariant value or a resolvent coefficient) is the centered representative of its residue once $\ell^k>2B$: two such integers differ by $\gamma$ with $|\gamma|\le2B$, $D=1$.
- (T2) *Root of a resolvent against an integer.* Whether $F(\tau\alpha)=c$ for $c\in\mathbb{Z}$, when $\tau\in G\_{i-1}$ and $F=F\_i$ (or $\tau\in G$ and $F$ an $H$-invariant): $\gamma=F(\tau\alpha)-c$, $C=2B\_F$, and the conjugates of $\gamma$ are $F(\sigma\tau\alpha)-c$, $\sigma\in G$, which take at most $m$ distinct values since they depend only on the coset $\sigma\tau G\_i$ (resp. $\sigma\tau H$); so $D\le m$.
- (T3) *Two resolvent roots.* Whether $F(\tau\alpha)=F(\tau'\alpha)$: by Lemma 1.2, $2v(\beta\_a-\beta\_b)\le v\_\ell(\mathrm{disc},\mathrm{Res})\le m(m-1)\log\_\ell(2B\_F)$, i.e. $(C,D)=(2B\_F,\binom m2)$.
- (T4) *Two roots of $f$.* $\gamma=\alpha\_i-\alpha\_j$: $2v(\gamma)\le v\_\ell(\Delta)$, i.e. $\ell^k>\ell^{v\_\ell(\Delta)/2}$ suffices, and $\ell^k>\Delta$ suffices without knowing $v\_\ell(\Delta)$.
- (T5) *Root counting for a monic $h\in\mathbb{Z}[x]$ in a local field* (§3): correct once $k>v\_\ell(\mathrm{disc},h)$; for $h=\mathrm{Res}\_{F,H}$ this is implied by $\ell^k>(2B\_F)^{m(m-1)}$, for $h=f$ by $\ell^k>\Delta$.

All five are covered by $\ell^k>M\_F:=(2B\_F+2)^{\mu(m)}$ with $\mu(1)=1$, $\mu(m)=m(m-1)$ for $m\ge2$, since $(2B\_F+2)^{\mu(m)}\ge\max{2(1+B\_F)^m,(2B\_F)^{m(m-1)}}$, together with $\ell^k>\Delta$.

## 3. Local precision: Hensel and Krasner

Let $f=\prod\_{j}g\_j$ over $\mathbb{Z}*\ell$ with $g\_j$ monic irreducible of degree $d\_j$, and for a root $\alpha$ of $g\_j$ let $\rho\_j=\max*{\alpha'\ne\alpha}v(\alpha-\alpha')$ over the conjugates $\alpha'$ of $\alpha$.

**Lemma 3.1 (root distances).** (a) $\rho\_j\le v(g\_j'(\alpha))=v\_\ell(\mathrm{disc},g\_j)/d\_j$. (b) For roots $\alpha$ of $g\_a$ and $\beta$ of $g\_b$, $a\ne b$: $v(\alpha-\beta)\le v\_\ell(\mathrm{Res}(g\_a,g\_b))$. (c) $\sum\_jv\_\ell(\mathrm{disc},g\_j)+2\sum\_{a\<b}v\_\ell(\mathrm{Res}(g\_a,g\_b))=v\_\ell(\Delta)$; in particular every quantity in (a),(b) is $\le v\_\ell(\Delta)$, and every root distance is $\le v\_\ell(\Delta)/2$.

*Proof.* (a) $g\_j'(\alpha)=\prod\_{\alpha'\ne\alpha}(\alpha-\alpha')$, all factors of valuation $\ge0$; $\mathrm{disc},g\_j=\pm\prod\_{\alpha}g\_j'(\alpha)$ and all conjugates of $g\_j'(\alpha)$ have equal valuation. (b) $\mathrm{Res}(g\_a,g\_b)=\prod\_{\alpha,\beta}(\alpha-\beta)$. (c) $\mathrm{disc}f=\prod\_j\mathrm{disc},g\_j\prod\_{a\<b}\mathrm{Res}(g\_a,g\_b)^2$. $\square$

**Lemma 3.2 (Hensel).** Let $g,h\in\mathbb{Z}*\ell[x]$ be monic with $f\equiv gh\pmod{\ell^k}$ and $r=v*\ell(\mathrm{Res}(g,h))\<k/2$. Then there are unique monic $g^*,h^*\in\mathbb{Z}*\ell[x]$ with $f=g^h^$ and $g^\equiv g$, $h^\equiv h\pmod{\ell^{k-r}}$; moreover $v*\ell(\mathrm{Res}(g^*,h^*))=r$, and $g^*,h^*$ can be computed to any precision from $g,h$ and the exact $f$.

This is the standard quadratic Hensel lemma; the lifting iteration is $\ell$-adic Newton on the map $(g,h)\mapsto gh$, whose derivative is invertible modulo $\ell^{,\cdot}$ up to the factor $\mathrm{Res}(g,h)$. Since $f$ is exact, the lift is available at any precision at the cost of the iteration, independently of the precision at which $(g,h)$ were first found.

**Lemma 3.3 (Krasner stability).** Let $g\in\mathbb{Z}*\ell[x]$ be monic irreducible of degree $d$, $\tilde g\in\mathbb{Z}*\ell[x]$ monic of degree $d$ with $\tilde g\equiv g\pmod{\ell^k}$, and $k>d\rho\_g$; by Lemma 3.1(a), $k>v\_\ell(\mathrm{disc},g)$ suffices. Then $\tilde g$ is irreducible, and there is a bijection $\tilde\alpha\mapsto\alpha'$ from the roots of $\tilde g$ to the roots of $g$, equivariant for $\mathrm{Gal}(\bar{\mathbb{Q}}*\ell/\mathbb{Q}*\ell)$, with $v(\tilde\alpha-\alpha')>\rho\_g$ and $\mathbb{Q}*\ell(\tilde\alpha)=\mathbb{Q}*\ell(\alpha')$.

*Proof.* $\tilde\alpha$ is integral (root of a monic integral polynomial) and $g(\tilde\alpha)=-(\tilde g-g)(\tilde\alpha)$ has valuation $\ge k$. Since $g(\tilde\alpha)=\prod\_{\alpha'}(\tilde\alpha-\alpha')$ over the $d$ roots, some factor has $v(\tilde\alpha-\alpha')\ge k/d>\rho\_g$. That $\alpha'$ is unique: two such would give $v(\alpha'-\alpha'')>\rho\_g$. Krasner's lemma gives $\mathbb{Q}*\ell(\alpha')\subseteq\mathbb{Q}*\ell(\tilde\alpha)$, and $[\mathbb{Q}*\ell(\tilde\alpha):\mathbb{Q}*\ell]\le d=[\mathbb{Q}*\ell(\alpha'):\mathbb{Q}*\ell]$ forces equality and the irreducibility of $\tilde g$. The map is Galois-equivariant by construction; its image is Galois-stable, and Galois is transitive on the roots of $g$, so the map is surjective between sets of size $d$, hence bijective. $\square$

**Proposition 3.4 (local policy).** Let $k\_\ell>v\_\ell(\Delta)$. Then, from $f$ modulo $\ell^{k\_\ell}$ and computations carried at precision $k\_\ell$:

1. Any factorization $f\equiv\prod\_jg\_j\pmod{\ell^{k\_\ell}}$ into monic factors with pairwise $v\_\ell(\mathrm{Res}(g\_a,g\_b))\<k\_\ell/2$ lifts uniquely (Lemma 3.2) to the true factorization into products of irreducibles, and the true factorization reduced modulo $\ell^{k\_\ell}$ satisfies this hypothesis (Lemma 3.1(c)).
2. After Hensel-lifting each candidate to precision $k\_\ell$, its irreducibility, and the fields $\mathbb{Q}*\ell(\alpha)$ of its roots (hence $e\_j$, $f\_j$, $d\_j=e\_jf\_j$), are determined by its residue class (Lemma 3.3, as $k*\ell>v\_\ell(\mathrm{disc},g\_j)$): every lift has the same root fields. So they may be read from the $\ell$-maximal order of $\mathbb{Z}\_\ell[x]/(\tilde g)$ for any lift $\tilde g$ ([the ramified primes](ramified-primes.md), §2), or from a Newton-polygon factorization, applied to the residue.
3. Distinct roots of $f$ in $\bar{\mathbb{Q}}*\ell$ are distinct at precision $k*\ell$ (Lemma 3.1(c)), so root-equality and root-counting tests on $f$ over any local field are correct.
4. The same holds for a monic $h\in\mathbb{Z}[x]$ in place of $f$ once $k\_\ell>v\_\ell(\mathrm{disc},h)$.

If $v\_\ell(\Delta)$ is not known ([the ramified primes](ramified-primes.md)'s factorization of $\Delta$ not yet done), $\ell^{k\_\ell}>\Delta$ implies $k\_\ell>v\_\ell(\Delta)$.

## 4. Tying the numbering at $p$ to the computations at $\ell$

An *$\ell$-numbering* is a list $\beta\_1,\dots,\beta\_n$ of the roots of $f$ in $\bar{\mathbb{Q}}*\ell$ (computed at precision $k*\ell$ inside the fields of Proposition 3.4). A permutation $\pi\in S\_n$ is *admissible* if there is an embedding $\iota\:N\to\bar{\mathbb{Q}}*\ell$ with $\iota(\alpha\_i)=\beta*{\pi(i)}$ for all $i$. The admissible permutations form one coset $\pi\_0G$: embeddings differ by automorphisms of $N$, i.e. by $G$. Transporting any object defined through the $p$-numbering (the elements of $G$, the subgroups $D\_\ell,I\_\ell$ to be located, Frobenius) to the $\ell$-side is conjugation by an admissible $\pi$, and any two choices differ by an element of $G$, which is why all outputs are stated up to conjugacy in $G$.

**Proposition 4.1 (alignment).** With the certificate chain $(G\_i,F\_i,c\_i)$ of §1: $\pi$ is admissible iff $F\_i(\beta\_\pi)=c\_i$ for $i=1,\dots,t$, where $\beta\_\pi=(\beta\_{\pi(1)},\dots,\beta\_{\pi(n)})$. Moreover, given $\pi$ with $F\_j(\beta\_\pi)=c\_j$ for $j\<i$, exactly one coset representative $\tau\in G\_{i-1}/G\_i$ satisfies $F\_i(\beta\_{\pi\tau})=c\_i$.

*Proof.* Fix an embedding $\iota\_0$ with permutation $\pi\_0$, so $\beta\_{\pi(i)}=\iota\_0(\alpha\_{\rho(i)})$ with $\rho=\pi\_0^{-1}\pi$, and $F(\beta\_\pi)=\iota\_0(F(\rho\alpha))$. If $\pi$ is admissible then $\rho\in G$ and $F\_i(\rho\alpha)=F\_i(\alpha)=c\_i$. Conversely, induct on $i$: if $\rho\in G\_{i-1}$ and $F\_i(\rho\alpha)=c\_i=F\_i(\alpha)$, the recorded distinctness of the values $F\_i(\tau\alpha)$ over $\tau\in G\_{i-1}/G\_i$ forces $\rho G\_i=G\_i$. So $\rho\in G\_t=G$, and $\iota=\iota\_0\circ\rho$ is an embedding with $\iota(\alpha\_i)=\beta\_{\pi(i)}$. For the second statement, $F\_i(\beta\_{\pi\tau})=\iota\_0(F\_i(\rho\tau\alpha))$ equals $c\_i$ iff $\rho\tau\in G\_i$ iff $\tau G\_i=\rho^{-1}G\_i$, which picks out exactly one representative. $\square$

So an admissible $\pi$ is found by $\sum\_im\_i$ evaluations of invariants at the $\beta$'s, one step of the Stauduhar descent per level, each step being a test of type (T2) with $(C,D)=(2B\_{F\_i},m\_i)$; the evaluations are polynomials with integer coefficients in integral quantities known at precision $k\_\ell$, so their residues are known at precision $k\_\ell$, and by Corollary 2.2 each test is decided correctly once $\ell^{k\_\ell}>(2B\_{F\_i})^{m\_i}$.

At an unramified $\ell$, the $\beta$'s lie in unramified extensions and Frobenius acts as $\beta\mapsto\beta^{\ell}$ modulo $\ell$, lifted uniquely; the resulting permutation $\phi$ of ${1,\dots,n}$ is exact (roots are simple modulo $\ell$), and the Frobenius class is that of $\pi^{-1}\phi\pi\in G$. Precision at $\ell$ is then needed only for the alignment tests. At a ramified $\ell$, the $\beta$'s lie in the fields of Proposition 3.4 and the tests locating $D\_\ell$ and $I\_\ell$ are root-equality tests (T4) and invariant tests (T2), (T3).

## 5. The policy and its correctness

**Data fixed before any local computation.** $R$; $\Delta$; the list $\mathcal{F}$ of pairs $(F,m)$: the chain invariants $(F\_i,m\_i)$ and every auxiliary $(F,[G\:H])$ that any later step will use. For each pair, $B\_F=|F|*1R^{\deg F}$ and $M\_F=(2B\_F+2)^{\mu(m)}$ with $\mu(1)=1$, $\mu(m)=m(m-1)$ for $m\ge2$. Set
$$M^\*=\max*{(F,m)\in\mathcal{F}}M\_F,\qquad\mathbf{M}=M^\*\cdot\Delta .$$

**Policy.** At every prime $\ell$, including $p$, compute at precision

$$k\_\ell=\lfloor\log\_\ell\mathbf{M}\rfloor+1,\qquad\text{i.e. the least }k\text{ with }\ell^{k}>\mathbf{M}.$$

In bits, $k\_\ell\log\_2\ell<\log\_2\mathbf{M}+\log\_2\ell$, with

$$\log\_2\mathbf{M}\le n\log\_2n+(2n-2)\log\_2|f|*2+\max*{(F,m)}\mu(m)\big(\log\_2|F|\_1+\deg F\log\_2R+2\big),$$

using $\Delta\le n^n|f|*2^{2n-2}$. Sharper, when $v*\ell(\Delta)$ is known: $k\_\ell=\max{\lfloor\log\_\ell M^*\rfloor+1,\ v\_\ell(\Delta)+1}$; at $p\nmid\Delta$ this is $\lfloor\log\_pM^*\rfloor+1$.

**Theorem 5.1.** Under the policy, every test of types (T1) to (T5) performed at any prime, the local factorization and root-field data of Proposition 3.4 at every prime, and the alignment of Proposition 4.1 at every prime, are decided correctly from residues at precision $k\_\ell$.

*Proof.* $\ell^{k\_\ell}>\mathbf{M}\ge M\_F$ for every pair in $\mathcal{F}$ and $\ell^{k\_\ell}>\Delta\ge\ell^{v\_\ell(\Delta)}$. (T1): $M\_F\ge2(1+B\_F)^m\ge2B\_F$. (T2): $M\_F\ge(2B\_F)^m$ and Corollary 2.2 with $D\le m$. (T3): $M\_F\ge(2B\_F)^{\binom m2}$ and Lemma 1.2. (T4) and Proposition 3.4: $k\_\ell>v\_\ell(\Delta)$. (T5): $M\_F\ge(2B\_F)^{m(m-1)}\ge\ell^{v\_\ell(\mathrm{disc},\mathrm{Res}*{F,H})}$ by Lemma 1.2, then Proposition 3.4(4). Alignment: each step is a (T2) test with $(2B*{F\_i},m\_i)$, and by Proposition 4.1 the outcome of the exact test at step $i$ is what the inductive hypothesis ($\rho\in G\_{i-1}$) requires, so the degree bound $D\le m\_i$ applies at each step. The residues fed to every test are exact at precision $k\_\ell$ because each is a polynomial with integer coefficients in integral elements known at that precision. $\square$

**What the policy does not cover.** Any invariant or subgroup not in $\mathcal{F}$ when $\mathbf{M}$ is fixed; the remedy is to add its $M\_F$ and recompute $\mathbf{M}$ before use, which only increases every $k\_\ell$ and invalidates nothing already decided. Values already computed at a smaller precision remain valid, since all tests are monotone in $k$.

**Certificate content for this step.** $R$, $\Delta$, the list $\mathcal{F}$ with $(\deg F,|F|*1,m)$ per entry, $M^\*$, $\mathbf{M}$, and for each prime used, $k*\ell$; the exact resolvents $\mathrm{Res}*{F\_i}$ with nonzero discriminants (the distinctness hypothesis of Proposition 4.1); and, per prime, the certified local factorization type. From these, a checker recomputes $k*\ell$ and confirms $\ell^{k\_\ell}>\mathbf{M}$ without repeating any local computation.