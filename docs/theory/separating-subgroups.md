# Choosing and realising the separating coset actions

**Setting.** As in [rational classes](rational-classes.md): $G\le S\_n$ with the numbering at $p$, classes $C\_1,\dots,C\_r$ with representatives $g\_k$, power maps and rational classes from [the character table](character-table.md) to [rational classes](rational-classes.md). $G\_1$ is the stabilizer of the point $1$, so $\pi\_{G\_1}$ is the degree-$n$ permutation character. $R\_0$ denotes a root bound of $f$ ([the precision policy](precision-policy.md)). $\mathcal{S}$ denotes the set of conjugacy classes of subgroups of $G$ (or any sub-collection containing all cyclic subgroups, when the full lattice is not enumerated).

## 1. Blocks and sub-blocks

**Cycle-type blocks.** For a class $C$ let $\mathrm{ct}(C)$ be the cycle type of its elements on ${1,\dots,n}$, a partition $\lambda\vdash n$. The *block* of $\lambda$ is $\mathcal{B}*\lambda={C:\mathrm{ct}(C)=\lambda}$. By* [*rational classes*](rational-classes.md)*, Prop. 1.2(3), the factorization of $f$ modulo $\ell$ locates $\mathrm{Frob}*\ell$ in a block and no further.

**Lemma 1.1.** Each block is a union of rational classes, and any two distinct rational classes in the same block have non-conjugate cyclic subgroups of the same order, neither conjugate to a subgroup of the other.

*Proof.* $\sigma$ and $\sigma^k$ with $\gcd(k,o(\sigma))=1$ generate the same cyclic group and so have the same cycle type; hence $\mathrm{ct}$ is constant on rational classes. The order of a permutation is the lcm of its cycle lengths, so two classes in $\mathcal{B}\_\lambda$ have the same order $\mathrm{lcm}(\lambda)$; their cyclic subgroups are non-conjugate (distinct rational classes, [rational classes](rational-classes.md) Lemma 2.1(a)) and of equal order, and a subgroup of a cyclic group of the same order is the whole group. $\square$

**Sub-blocks.** Write $\mathcal{B}*\lambda=\rho*{\lambda,1}\sqcup\dots\sqcup\rho\_{\lambda,r\_\lambda}$ as a union of rational classes. By [rational classes](rational-classes.md) Theorem 4.1, coset actions can separate the $\rho\_{\lambda,j}$ from one another and cannot split any $\rho\_{\lambda,j}$ into its conjugacy classes. So the sub-blocks that coset actions must separate are exactly the rational classes, and the separation task is:

$$U=\Big{{\rho,\rho'}:\ \rho\ne\rho'\text{ rational classes in the same block}\Big},\qquad u:=|U|=\sum\_\lambda\binom{r\_\lambda}{2}.$$

Blocks with $r\_\lambda=1$ contribute nothing; if all $r\_\lambda=1$ the degree-$n$ action already determines the rational class ([rational classes](rational-classes.md) Cor. 4.3). The partition is computed by hashing the class representatives on cycle type and then applying the union to find of [rational classes](rational-classes.md) §2 inside each block.

**What a subgroup separates.** For $H\in\mathcal{S}$ define the *$H$-signature* of a rational class $\rho$ with generator $\sigma$ as the cycle type of $\sigma$ on $G/H$, i.e. the function $d\mapsto\pi\_H(\sigma^d)$ for $d\mid o(\sigma)$, computed from [rational classes](rational-classes.md) Lemma 1.1 and the power maps ($|H|$ class identifications per $H$). $H$ *covers* ${\rho,\rho'}\in U$ if the signatures differ; let $S\_H\subseteq U$ be the set of covered pairs. Since the full factorization type of a resolvent over $\mathbb{Q}\_\ell$ is available at the same cost as the count of degree-one factors ([rational classes](rational-classes.md) Prop. 1.2), the signature is the right notion; using only $\pi\_H(\sigma)$ gives a sub-instance with smaller $S\_H$ and every statement below still holds.

**Lemma 1.2 (cyclic subgroups cover).** Let $\rho\in\mathcal{B}*\lambda$ have cyclic subgroup $C*\rho$. Then $\pi\_{C\_\rho}(\rho)=[N\_G(C\_\rho)\:C\_\rho]\ge1$ and $\pi\_{C\_\rho}(\rho')=0$ for every other $\rho'\in\mathcal{B}*\lambda$. Hence $C*\rho$ covers every pair ${\rho,\rho'}$ with $\rho'$ in the block of $\rho$, and the family ${C\_\rho:\rho\in\mathcal{B}*\lambda\setminus{\rho*\lambda^{\max}}}$, taken over all blocks with $r\_\lambda\ge2$ and omitting in each block one class $\rho\_\lambda^{\max}$ of largest $[G\:C\_\rho]$, is a feasible cover of $U$.

*Proof.* The value on $\rho$ is [rational classes](rational-classes.md) Theorem 3.1's diagonal entry. $\pi\_{C\_\rho}(\rho')\ne0$ would require an element of $\rho'$ inside $C\_\rho$, hence $C\_{\rho'}$ conjugate to a subgroup of $C\_\rho$, excluded by Lemma 1.1. Omitting one class per block is harmless: a class whose signature is $0$ for every chosen $C\_{\rho}$ of its block is the omitted one. $\square$

The same argument shows that any $H$ with $\sigma\in H$ and $H\cap\rho'=\emptyset$ for all other $\rho'$ in the block covers the same pairs as $C\_\rho$; taking $H$ maximal with this property lowers the weight. All $[G\:C\_\rho]$ in a block equal $|G|/\mathrm{lcm}(\lambda)$.

## 2. Weighted set cover

**Instance.** Universe $U$; for each $H\in\mathcal{S}$ the set $S\_H$ and weight $w\_H=[G\:H]$, the degree of the resolvent that realizes the test; the per-prime cost of the test (one factorization or gcd of a degree-$w\_H$ polynomial over $\mathbb{F}*\ell$) is linear in $w\_H$ up to logarithmic factors, and this cost is paid at every prime up to $X$, which is why the resolvent degree is the weight. Find $\mathcal{H}\subseteq\mathcal{S}$ with $\bigcup*{H\in\mathcal{H}}S\_H=U$ minimizing $w(\mathcal{H})=\sum\_{H\in\mathcal{H}}[G\:H]$. Feasibility is Lemma 1.2 (which is why $\mathcal{S}$ must contain the cyclic subgroups).

**Greedy.** $\mathcal{H}\leftarrow\emptyset$, $D\leftarrow\emptyset$ (covered pairs). While $D\ne U$: choose $H\in\mathcal{S}$ with $S\_H\setminus D\ne\emptyset$ minimizing $w\_H/|S\_H\setminus D|$; $\mathcal{H}\leftarrow\mathcal{H}\cup{H}$, $D\leftarrow D\cup S\_H$. Cost: $|\mathcal{S}|$ signature computations, then at most $u$ rounds of $|\mathcal{S}|$ ratio evaluations.

**Theorem 2.1 (greedy guarantee).** Let $\mathrm{OPT}$ be the minimum weight of a cover of $U$ by sets from $\mathcal{S}$ and $H\_u=\sum\_{j=1}^u1/j\le1+\ln u$. Then $w(\mathcal{H}\_{\mathrm{greedy}})\le H\_u\cdot\mathrm{OPT}$.

*Proof.* When greedy selects $H$ with newly covered set $A=S\_H\setminus D$, charge each pair $x\in A$ the price $p(x)=w\_H/|A|$; then $w(\mathcal{H}*{\mathrm{greedy}})=\sum*{x\in U}p(x)$. Let $x\_1,\dots,x\_u$ be the pairs in the order they are covered. Just before $x\_j$ is covered, at least $u-j+1$ pairs are uncovered; an optimal cover $\mathcal{H}^*$ covers them with total weight $\le\mathrm{OPT}$, so some $H^*\in\mathcal{H}^*$ has $w\_{H^*}/|S\_{H^*}\setminus D|\le\mathrm{OPT}/(u-j+1)$ (otherwise summing $|S\_{H^*}\setminus D|\cdot w\_{H^*}/|S\_{H^*}\setminus D|$ over $\mathcal{H}^\*$ would exceed $\mathrm{OPT}$). Greedy's chosen ratio is at most this, so $p(x\_j)\le\mathrm{OPT}/(u-j+1)$. Summing over $j$ gives $\mathrm{OPT}\cdot H\_u$. $\square$

**Theorem 2.2 (explicit bound).**

$$\sum\_{H\in\mathcal{H}*{\mathrm{greedy}}}[G\:H]\ \le\ (1+\ln u)\sum*{\lambda:,r\_\lambda\ge2}(r\_\lambda-1),\frac{|G|}{\mathrm{lcm}(\lambda)},\qquad u=\sum\_\lambda\binom{r\_\lambda}{2}\le\binom{r\_{\mathbb{Q}}}{2},$$

where $r\_\mathbb{Q}$ is the number of rational classes of $G$. In particular the total is at most $(1+\ln u)(r\_\mathbb{Q}-1)|G|/2$, and, whenever every rational class of $G$ lies in a block by itself, it is $0$.

*Proof.* Theorem 2.1 with $\mathrm{OPT}\le w$ of the cover of Lemma 1.2, whose weight is $\sum\_\lambda(r\_\lambda-1)|G|/\mathrm{lcm}(\lambda)$ because every $C\_\rho$ in $\mathcal{B}\_\lambda$ has order $\mathrm{lcm}(\lambda)$; nontrivial blocks have $\mathrm{lcm}(\lambda)\ge2$. $\square$

The bound is a worst case; greedy will usually prefer non-cyclic subgroups of small index (maximal subgroups, stabilizers of small sets) whose signatures separate many pairs at once.

## 3. Realizing a test as a resolvent

**Construction.** Fix once and for all a base $t=(i\_1,\dots,i\_k)$ of $G$: a tuple of points whose pointwise stabilizer $G\_t$ is trivial (from the BSGS of [the character table](character-table.md); an irredundant base has $k\le\log\_2|G|$). Put $m\_t=\prod\_{j=1}^kx\_{i\_j}^{,j}$ and, for $H\le G$,

$$F\_H=\sum\_{h\in H}h\cdot m\_t=\sum\_{(i'*1,\dots,i'k)\in H\cdot t}\prod{j=1}^kx*{i'\_j}^{,j},$$

the orbit sum of $m\_t$ under $H$, indexed by the $H$-orbit of the base tuple.

**Proposition 3.1.** (a) $\mathrm{Stab}\_G(F\_H)=H$. (b) $\deg F\_H=\binom{k+1}{2}$, $|F\_H|*1=|H|$, and $F\_H$ has $|H|$ distinct monomials. (c) $F\_H$ has integer coefficients, so $B*{F\_H}=|H|R\_0^{\binom{k+1}2}$ in [the precision policy](precision-policy.md)'s notation.

*Proof.* For $\sigma\in S\_n$, $\sigma\cdot m\_t=m\_{\sigma t}$, and distinct $k$-tuples of distinct points give distinct monomials (the exponent $j$ identifies which point sits in position $j$). So $\sigma$ fixes $F\_H$ iff $\sigma$ permutes the set of tuples $H\cdot t$. For $g\in G$: if $gHt=Ht$ then $gt=ht$ for some $h\in H$, so $h^{-1}g\in G\_t=1$ and $g\in H$; conversely $H$ preserves its own orbit. Since $H\_t\subseteq G\_t=1$, $|Ht|=|H|$, which gives (b); (c) is clear. $\square$

When $H$ is the stabilizer in $G$ of a subset $\Omega$, the cheaper $F=\prod\_{i\in\Omega}x\_i$ (degree $|\Omega|$, $|F|\_1=1$) has the same stabilizer and should be used instead; Proposition 3.1 is the construction that works for every $H$.

**Resolvent.** $R\_{G,H,F}(x)=\prod\_{\tau\in G/H}(x-F(\tau\alpha))$, of degree $m=[G\:H]$, with $\mathbb{Z}$-coefficients bounded by [the precision policy](precision-policy.md) Lemma 1.2; it is computed exactly from the $p$-adic roots at the precision of [the precision policy](precision-policy.md) §5 (test T1).

**Squarefreeness.** The $m$ polynomials $\tau F$, $\tau\in G/H$, are pairwise distinct by Proposition 3.1(a), but their values at $\alpha$ may coincide. Let $P\_F=\prod\_{{\tau H\ne\tau'H}}(\tau F-\tau'F)\in\mathbb{Z}[x\_1,\dots,x\_n]$, nonzero of degree $\le\binom m2\deg F$; $R\_{G,H,F}$ is squarefree iff $P\_F(\alpha)\ne0$, decided exactly by $\mathrm{disc}(R\_{G,H,F})\ne0$ ([the precision policy](precision-policy.md), T1/T3).

**Lemma 3.2 (Tschirnhaus repair).** For $T=\sum\_{j\<n}c\_jx^j$ set $\alpha^T=(T(\alpha\_1),\dots,T(\alpha\_n))$. The function $Q(c)=P\_F(\alpha^T)$ is a nonzero polynomial in $c=(c\_0,\dots,c\_{n-1})$ of degree $\le\binom m2\deg F$; hence for $c$ drawn uniformly from ${0,\dots,S-1}^n$, $\Pr[Q(c)=0]\le\binom m2\deg F/S$, and $\alpha^T$ has the same Galois action as $\alpha$, so $R\_{G,H,F}$ formed with $\alpha^T$ is a squarefree resolvent for the same $H$.

*Proof.* The map $c\mapsto\alpha^T$ is the linear map $c\mapsto Vc$ with $V$ the Vandermonde matrix of the distinct $\alpha\_i$, invertible; a nonzero polynomial composed with an invertible linear map is nonzero, of the same degree. Schwartz to Zippel gives the probability. $G$ acts on $\alpha^T$ through the same permutations. $\square$

After the repair, $F$ is replaced by $F\circ T$ in [the precision policy](precision-policy.md)'s list $\mathcal{F}$: $\deg(F\circ T)=\deg F\cdot\deg T$ and $|F\circ T|\_1\le|F|\_1|T|\_1^{\deg F}$. Trying $T$ of small degree first keeps these small; Lemma 3.2 is the guarantee that a search over degree $\<n$ terminates.

**Theorem 3.3 (what the factorization reads off).** Let $R=R\_{G,H,F}$ be squarefree, $\beta=F(\alpha)$.

1. $R$ is irreducible over $\mathbb{Q}$, $\mathbb{Q}(\beta)=N^H$, and the splitting field of $R$ is $N^{\mathrm{core}\_G(H)}$ with $\mathrm{core}\_G(H)=\bigcap\_gg Hg^{-1}$.
2. For every prime $\ell\nmid\mathrm{disc}(R)$: $\ell$ is unramified in $N^{\mathrm{core}*G(H)}$, so $I*\ell\le\mathrm{core}*G(H)$ acts trivially on $G/H$ and $\pi\_H$ is constant on the Frobenius coset $\phi I*\ell$; and the number of degree-one factors of $R$ modulo $\ell$ equals $\pi\_H(\phi)$, while the factorization type of $R$ modulo $\ell$ is the cycle type of $\phi$ on $G/H$. This holds whether or not $\ell\mid\mathrm{disc}f$.
3. For $\ell\mid\mathrm{disc}(R)$ with $\ell\nmid\mathrm{disc}f$, the same quantities are the number of $\mathbb{Q}*\ell$-roots and the $\mathbb{Q}*\ell$-factorization type of $R$, computed at precision $k\_\ell$ ([the precision policy](precision-policy.md), T5), since $(F,m)\in\mathcal{F}$.

*Proof.* 1. $G$ permutes the roots $F(\tau\alpha)$ transitively (transitively on $G/H$) and they are distinct, so $R$ is irreducible; $\mathrm{Stab}*G(\beta)={g\:F(g\alpha)=F(\alpha)}=H$ by distinctness, so $\mathbb{Q}(\beta)=N^H$; the splitting field is the Galois closure, fixed by the largest normal subgroup of $G$ inside $H$. 2. $\ell\nmid\mathrm{disc}(R)=\mathrm{disc}(\mathbb{Z}[\beta])$ implies $\mathbb{Z}[\beta]$ is $\ell$-maximal and $\ell\nmid\mathrm{disc},\mathcal{O}*{N^H}$, so $\ell$ is unramified in $N^H$ ([the ramified primes](ramified-primes.md) Lemma 1.2), hence in its Galois closure $N^{\mathrm{core}*G(H)}$ (*[*the ramified primes*](ramified-primes.md) *Theorem 1.1), i.e. $I*\ell\le\mathrm{core}*G(H)$. Kummer to Dedekind for the $\ell$-maximal order $\mathbb{Z}[\beta]$: the factorization of $R$ modulo $\ell$ mirrors the splitting of $\ell$ in $N^H$, whose residue degrees are the orbit lengths of the Frobenius of $N^{\mathrm{core}}$ on the roots, i.e. of $\phi$ (any element of $\phi I*\ell$) on $G/H$; the orbits of length one are the fixed cosets. 3. [the precision policy](precision-policy.md) Prop. 3.4(4) with $h=R$, as $\ell^{k\_\ell}>M\_F\ge\ell^{v\_\ell(\mathrm{disc}R)}$. $\square$

**Proposition 3.4 (excluded primes).** The set $E\_H={\ell:\ell\mid\mathrm{disc}(R\_{G,H,F})}$ satisfies

$$|E\_H|\le\log\_2|\mathrm{disc}(R\_{G,H,F})|\le m(m-1)\log\_2(2B\_F)=m(m-1)\Big(1+\log\_2|F|\_1+\deg F\cdot\log\_2R\_0\Big),$$

which for the construction of Proposition 3.1 (no Tschirnhaus repair) is $m(m-1)\big(1+\log\_2|H|+\binom{k+1}2\log\_2R\_0\big)$, and every $\ell\in E\_H$ is at most $|\mathrm{disc}(R)|\le(2B\_F)^{m(m-1)}$.

*Proof.* A distinct primes' product dividing a nonzero integer $D$ is at most $|D|$, so their number is at most $\log\_2|D|$; the bound on $|\mathrm{disc}R|$ is [the precision policy](precision-policy.md) Lemma 1.2 and $B\_F=|F|\_1R\_0^{\deg F}$. $\square$

At the primes of $E\_H$ the test is not lost, only more expensive: Theorem 3.3(3) supplies the same value at precision $k\_\ell$, and $E\_H$ is computed exactly from the factorization of $\mathrm{disc}(R)$ when wanted, or bounded as above without factoring. For the loop over $\ell\le X$, the primes of $\bigcup\_{H\in\mathcal{H}}E\_H\cap[1,X]$ are the only ones where a resolvent of the separating family is handled $\ell$-adically rather than modulo $\ell$.

**Summary of the family attached to $G$.** Before the unramified-prime loop: blocks and sub-blocks (§1), signatures for $H\in\mathcal{S}$, greedy cover $\mathcal{H}$ (§2) with weight bounded by Theorem 2.2, invariants $F\_H$ (§3) repaired by Lemma 3.2 until every $R\_{G,H,F}$ is squarefree, the pairs $(F\_H,[G\:H])$ appended to $\mathcal{F}$ so that [the precision policy](precision-policy.md)'s $\mathbf{M}$ covers them, the exact resolvents and their discriminants, and the excluded sets $E\_H$. At each $\ell\le X$: factor $f$ modulo $\ell$ (block); factor each $R\_{G,H,F}$, $H\in\mathcal{H}$, modulo $\ell$ or $\ell$-adically according to $\ell\notin E\_H$ or $\ell\in E\_H$ (rational class, by [rational classes](rational-classes.md) Theorem 4.1); then [the precision policy](precision-policy.md) §4 for the conjugacy class when the rational class is not a single class.