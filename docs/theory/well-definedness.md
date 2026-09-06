# Independence of $D\_\ell$ from the choices in [local Galois groups](https://claude.ai/chat/local-galois-groups.md) to [matching](https://claude.ai/chat/matching.md)

*One of the derivations behind* [*artin*](https://claude.ai/README.md)*; the* [*index*](https://claude.ai/chat/README.md) *lists them all and names the code that implements each.*

**Conventions (fixed for this note).** Global numbering $\alpha\_1,\dots,\alpha\_n$ with $G\le S\_n$; local numbering $\beta\_1,\dots,\beta\_n\in L$ from [local Galois groups](https://claude.ai/chat/local-galois-groups.md) with local group $D\le S\_n$ defined by $\phi(\beta\_i)=\beta\_{d(i)}$, inertia $I\trianglelefteq D$, higher ramification groups $D\_u\le I$ ($D\_0=I$), Frobenius coset $\Phi\in D/I$. A *matching* is $\tau\in S\_n$ with

$$\beta\_i=\iota(\alpha\_{\tau(i)})\quad(1\le i\le n)\qquad\text{for some embedding }\iota\:N\to\bar{\mathbb{Q}}\_\ell .$$

This is the $\tau$ of [matching](https://claude.ai/chat/matching.md) Lemma 1.1. [matching](https://claude.ai/chat/matching.md)'s search outputs $\rho$ with $\tau\rho\in G$, so $\rho^{-1}=g^{-1}\tau$ is itself a matching, and by [matching](https://claude.ai/chat/matching.md) Theorem 2.3 the outputs surviving at the [the precision policy](https://claude.ai/chat/precision-policy.md) precision are exactly the $\rho$ for which $\rho^{-1}$ is a matching. For $\kappa\in S\_n$ and a polynomial $F$, $(\kappa F)(x)=F(x\_{\kappa(1)},\dots,x\_{\kappa(n)})$, a left action, and $g(F(\alpha))=(gF)(\alpha)$ for $g\in G$. For a matching $\tau$ with embedding $\iota$ let $\mathfrak{L}*\iota$ be the prime of $N$ above $\ell$ determined by $\iota$, and $D*\iota,I\_\iota,(D\_\iota)*u,\mathrm{Fr}*\iota$ its decomposition group, inertia group, ramification groups and Frobenius coset in the global numbering. The *conjugation map* of $\tau$ is $c\_\tau\:D\to S\_n$, $d\mapsto\tau d\tau^{-1}$.

**Lemma 0.1 (**[**matching**](https://claude.ai/chat/matching.md) **Lemma 1.1, restated).** For a matching $\tau$ with embedding $\iota$: $c\_\tau(D)=\tau D\tau^{-1}=D\_\iota\le G$, $\tau I\tau^{-1}=I\_\iota$, $\tau D\_u\tau^{-1}=(D\_\iota)*u$, $\tau\Phi\tau^{-1}=\mathrm{Fr}*\iota I\_\iota$; and for every $\kappa\in S\_n$ and integral $F$, $F(\beta^\kappa)=\iota\big((\tau\kappa F)(\alpha)\big)$ where $\beta^\kappa=(\beta\_{\kappa(1)},\dots,\beta\_{\kappa(n)})$.

*Proof.* Let $d\in D$ with automorphism $\phi$ of $L$, and $g\in G$ with $\phi\circ\iota=\iota\circ g$ (it exists because $\iota(N)$ is $\phi$-stable, and $g$ is the element of $D\_\iota$ corresponding to $\phi|*{\iota(N)}$). Then $\beta*{d(i)}=\phi(\beta\_i)=\phi\iota(\alpha\_{\tau(i)})=\iota(\alpha\_{g\tau(i)})$ and $\beta\_{d(i)}=\iota(\alpha\_{\tau d(i)})$, so $\tau d=g\tau$, $g=\tau d\tau^{-1}$. The map $d\mapsto g$ is the isomorphism $D\cong D\_\iota$ given by $\phi\mapsto\phi|*{\iota(N)}$, which respects inertia, the ramification filtration and Frobenius. The last formula is $F(\beta*{\kappa(1)},\dots)=F(\iota\alpha\_{\tau\kappa(1)},\dots)=\iota(F(\alpha\_{\tau\kappa(1)},\dots))$. $\square$

## 1. Matchings form a $G$-coset and give conjugate subgroups

**Proposition 1.1.** Let $\tau$ be a matching with embedding $\iota$.

1. For $g\in G$, $g\tau$ is the matching of the embedding $\iota\circ g^{-1}$; every matching arises this way; the set of matchings is the left coset $G\tau$.
2. Left multiplication by $g\in G$ preserves every recorded invariant value: for $\kappa\in S\_n$ put $V\_i(\kappa)=(\kappa F\_i)(\alpha)$; then $V\_i(g\kappa)=g(V\_i(\kappa))$, so $V\_i(\kappa)=c\_i$ implies $V\_i(g\kappa)=c\_i$. In particular, the test of [matching](https://claude.ai/chat/matching.md) at a candidate numbering $\beta^{\pi}$ reads $F\_i(\beta^\pi)=\iota(V\_i(\tau\pi))$, and replacing $\tau$ by $g\tau$ (i.e. $\iota$ by $\iota\circ g^{-1}$) leaves every test outcome unchanged.
3. $c\_{g\tau}=\mathrm{inn}*g\circ c*\tau$. Hence $(g\tau)D(g\tau)^{-1}=g,(\tau D\tau^{-1}),g^{-1}$, and the same for $I$, each $D\_u$, and $\Phi$.

*Proof.* 1. With $\iota'=\iota\circ g^{-1}$: $\beta\_i=\iota(\alpha\_{\tau(i)})=\iota'(g\alpha\_{\tau(i)})=\iota'(\alpha\_{g\tau(i)})$, so $g\tau$ is the matching of $\iota'$. The embeddings of $N$ into $\bar{\mathbb{Q}}*\ell$ over $\mathbb{Q}$ form a torsor under $\mathrm{Aut}(N)=G$ by precomposition, and a matching determines its embedding ($\iota(\alpha\_j)=\beta*{\tau^{-1}(j)}$ on generators), so the matchings are exactly ${g\tau\:g\in G}$. 2. $V\_i(g\kappa)=(g\kappa F\_i)(\alpha)=g\big((\kappa F\_i)(\alpha)\big)=g(V\_i(\kappa))$, and $g$ fixes the rational number $c\_i$. The test formula is Lemma 0.1 with $F=F\_i$. 3. $(g\tau)d(g\tau)^{-1}=g(\tau d\tau^{-1})g^{-1}$. $\square$

**Theorem 1.2.** For any two matchings $\tau,\tau'$ there is $k\in G$ with

$$\tau'D\tau'^{-1}=k(\tau D\tau^{-1})k^{-1},\quad \tau'I\tau'^{-1}=k(\tau I\tau^{-1})k^{-1},\quad \tau'D\_u\tau'^{-1}=k(\tau D\_u\tau^{-1})k^{-1}\ (u\ge0),\quad \tau'\Phi\tau'^{-1}=k(\tau\Phi\tau^{-1})k^{-1},$$

and $c\_{\tau'}=\mathrm{inn}*k\circ c*\tau$. Consequently the filtered subgroup with distinguished coset

$$\big(D\_\ell\supseteq I\_\ell\supseteq D\_{\ell,1}\supseteq\cdots;\ \mathrm{Fr}*\ell I*\ell\big):=\big(\tau D\tau^{-1}\supseteq\tau I\tau^{-1}\supseteq\tau D\_1\tau^{-1}\supseteq\cdots;\ \tau\Phi\tau^{-1}\big)$$

is well defined up to simultaneous conjugation in $G$.

*Proof.* By Proposition 1.1(1), $\tau'=k\tau$ with $k\in G$; apply 1.1(3). Independently of the coset structure: $\tau D\tau^{-1}=D\_\iota$ and $\tau'D\tau'^{-1}=D\_{\iota'}$ are the decomposition groups of the primes $\mathfrak{L}*\iota,\mathfrak{L}*{\iota'}$ of $N$ above $\ell$ (Lemma 0.1), any two of which are $G$-conjugate, $\mathfrak{L}*{\iota'}=k\mathfrak{L}*\iota$, so that $D\_{\iota'}=kD\_\iota k^{-1}$ with inertia, ramification groups and Frobenius conjugated along; both arguments give the same $k$ (the one with $\iota'=\iota\circ k^{-1}$, for which $\mathfrak{L}*{\iota'}=k\mathfrak{L}*\iota$). $\square$

## 2. The full ambiguity of the conjugation map

**Proposition 2.1.** Let $\tau$ be a matching and $\tau'\in S\_n$ arbitrary. The following are equivalent:

1. $c\_{\tau'}=\mathrm{inn}*g\circ c*\tau$ on $D$ for some $g\in G$;
2. $\tau'\in G,\tau,C\_{S\_n}(D)$.

Hence: (a) every $\tau'=g\tau c$ with $g\in G$, $c\in C\_{S\_n}(D)$ satisfies $\tau'D\tau'^{-1}=g(\tau D\tau^{-1})g^{-1}$, and likewise for $I$, $D\_u$, $\Phi$, whether or not $\tau'$ is a matching; (b) any two matchings differ by an element of $G$ on the left composed with an element of $C\_{S\_n}(D)$ on the right, in fact with $c=1$ (Theorem 1.2), the factor $C\_{S\_n}(D)$ being exactly the relabelings that leave the map $c\_\tau$ unchanged; (c) $g\tau c$ is itself a matching iff $\tau c\tau^{-1}\in G$, i.e. iff $c\in\tau^{-1}G\tau\cap C\_{S\_n}(D)$.

*Proof.* (1)$\Leftrightarrow$(2): $\tau'd\tau'^{-1}=g\tau d\tau^{-1}g^{-1}$ for all $d\in D$ iff $(g\tau)^{-1}\tau'$ centralizes $D$ iff $\tau'\in g\tau C\_{S\_n}(D)$. (a) is (1) applied to $D$ and to its subgroups and cosets. (b): Theorem 1.2 gives $\tau'=k\tau$. (c): by Proposition 1.1(1) the matchings are $G\tau$, and $g\tau c\in G\tau$ iff $\tau c\tau^{-1}\in G$. $\square$

So the conjugation map $D\to G$ is determined up to $\mathrm{Inn}(G)$ by any element of $G\tau C\_{S\_n}(D)$, the matchings are the subset $G\tau$, and the invariant test of [matching](https://claude.ai/chat/matching.md) is what cuts $G\tau C\_{S\_n}(D)$ down to $G\tau$; the part of $C\_{S\_n}(D)$ it removes changes nothing in the map and nothing downstream.

## 3. Independence from the choices

The choices are: in [local Galois groups](https://claude.ai/chat/local-galois-groups.md), the composition order of the factors, the roots adjoined in the towers, the embeddings used to identify roots, and the numbering of the roots; in [matching](https://claude.ai/chat/matching.md), the transversals $T\_i$, the scan order, and (in the branching variant) the surviving $\rho$ retained; in both, the prime $\mathfrak{L}\mid\ell$, i.e. the embedding $\iota$.

**Proposition 3.1 (relabelings from** [**local Galois groups**](https://claude.ai/chat/local-galois-groups.md)**).** Two runs of [local Galois groups](https://claude.ai/chat/local-galois-groups.md) list the same set of roots of $f$ in $\bar{\mathbb{Q}}\_\ell$ in numberings $\beta$ and $\beta'=\beta^\sigma$ for some $\sigma\in S\_n$, with local data $D'=\sigma^{-1}D\sigma$, $I'=\sigma^{-1}I\sigma$, $D'\_u=\sigma^{-1}D\_u\sigma$, $\Phi'=\sigma^{-1}\Phi\sigma$. If $\tau$ is a matching for $\beta$ then $\tau\sigma$ is a matching for $\beta'$ (same embedding), and

$$(\tau\sigma)D'(\tau\sigma)^{-1}=\tau D\tau^{-1},$$

identically for $I$, $D\_u$, $\Phi$: the transported filtered datum is unchanged, not merely conjugate.

*Proof.* The root set is intrinsic ([local Galois groups](https://claude.ai/chat/local-galois-groups.md) Theorem 4.1; it is determined by $f\bmod\ell^{\nu+1}$ by [local Galois groups](https://claude.ai/chat/local-galois-groups.md) Proposition 1.1), so the two numberings differ by a permutation: $\beta'*i=\beta*{\sigma(i)}$. Local group in the new numbering: $\phi(\beta'*i)=\phi(\beta*{\sigma(i)})=\beta\_{d\sigma(i)}=\beta'*{\sigma^{-1}d\sigma(i)}$, so $d'=\sigma^{-1}d\sigma$; the subgroups and the coset transform the same way. Matching: $\beta'i=\beta{\sigma(i)}=\iota(\alpha*{\tau\sigma(i)})$, so $\tau'=\tau\sigma$. Then $\tau'D'\tau'^{-1}=\tau\sigma,\sigma^{-1}D\sigma,\sigma^{-1}\tau^{-1}=\tau D\tau^{-1}$. $\square$

**Proposition 3.2 (choices in** [**matching**](https://claude.ai/chat/matching.md)**).** At the [the precision policy](https://claude.ai/chat/precision-policy.md) precision the numberings surviving [matching](https://claude.ai/chat/matching.md) are exactly those induced by embeddings ([matching](https://claude.ai/chat/matching.md) Theorem 2.3), so any two outputs correspond to matchings of the same $\beta$ and differ by Theorem 1.2; with a fixed transversal system the output is unique. In the branching variant below the [the precision policy](https://claude.ai/chat/precision-policy.md) precision, every retained candidate that passes the final certificate ($\rho^{-1}D\rho\le G$ and all invariant congruences at full precision) is a matching, and the same conclusion holds for it.

*Proof.* [matching](https://claude.ai/chat/matching.md) Theorems 2.2, 2.3 and Proposition 3.3; Theorem 1.2. $\square$

**Theorem 3.3 (independence of downstream invariants).** Let $\mathcal{I}$ be any function of a filtered subgroup of $G$ with a distinguished coset that is invariant under simultaneous conjugation by $G$. Then $\mathcal{I}(D\_\ell\supseteq I\_\ell\supseteq\cdots;\mathrm{Fr}*\ell I*\ell)$ is independent of every choice in [local Galois groups](https://claude.ai/chat/local-galois-groups.md) and [matching](https://claude.ai/chat/matching.md) and of the prime $\mathfrak{L}\mid\ell$. This includes: the conjugacy classes of $D\_\ell$, $I\_\ell$ and each $D\_{\ell,u}$ in $G$; $e=|I\_\ell|$, $f=[D\_\ell\:I\_\ell]$, the orders $|D\_{\ell,u}|$ and the lower and upper ramification breaks; and, for every class function $\chi$ of $G$, the conductor exponent $f\_\ell(\chi)$, the Euler factor $P\_\ell(\chi;T)$ and the local root number $\varepsilon\_\ell(\chi)$ of [the character table](https://claude.ai/chat/character-table.md) §5.

*Proof.* By Proposition 3.1 the [local Galois groups](https://claude.ai/chat/local-galois-groups.md) choices do not change the transported datum at all; by Proposition 3.2 and Theorem 1.2 the [matching](https://claude.ai/chat/matching.md) choices and the choice of $\mathfrak{L}$ change it by simultaneous conjugation in $G$. The listed quantities are invariant under such conjugation: orders and indices trivially; ramification breaks because the filtration is transported as a whole; $f\_\ell(\chi)=\sum\_u[D\_{\ell,0}\:D\_{\ell,u}]^{-1}\big(\chi(1)-\langle\chi|*{D*{\ell,u}},1\rangle\big)$ and $P\_\ell(\chi;T)=\exp\big(-\sum\_{k\ge1}\frac{T^k}{k},|I\_\ell|^{-1}\sum\_{\sigma\in I\_\ell}\chi(\mathrm{Fr}*\ell^{,k}\sigma)\big)$ are sums of a class function over subgroups and cosets, unchanged when subgroups, coset and summation variable are conjugated together; $\varepsilon*\ell(\chi)$ depends on $\chi|*{D*\ell}$ as a class function of $D\_\ell$ and on the local field structure carried by $D\_\ell\cong D$, both transported by the conjugation ([the character table](https://claude.ai/chat/character-table.md) Theorem 5.1(4)). Independence of the lift of $\mathrm{Fr}\_\ell$ within its coset is [the character table](https://claude.ai/chat/character-table.md) Theorem 5.1. $\square$

**What is not invariant, and is not an output.** The subgroup $D\_\ell$ as a specific subset of $S\_n$, the map $c\_\tau$, and the particular matching are choice-dependent, by exactly $G$ on the left and $C\_{S\_n}(D)$ on the right (Proposition 2.1). The tables record the conjugacy class of $D\_\ell$ inside $G$ (a canonical representative under a fixed enumeration of the subgroup classes of $G$, or the class-function invariants of Theorem 3.3), never a representative permutation. Comparisons across primes, whether $D\_\ell$ and $D\_{\ell'}$ are conjugate, whether $I\_\ell$ is conjugate into $D\_{\ell'}$, whether two Frobenius classes coincide, are comparisons of conjugacy classes, which Theorem 3.3 licenses; comparisons of representatives are meaningless and are never made.