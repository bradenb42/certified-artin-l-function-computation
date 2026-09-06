# Character table, Schur indices, matrix models, and model-independence of the L-function data

**Conventions.** $G\le S\_n$ is given by generators acting on the fixed numbering of the $p$-adic roots. From the generators compute a base and strong generating set, $|G|$, the exponent $e=\exp(G)$, the classes $C\_1={1},\dots,C\_r$ with representatives $g\_k$, sizes $|C\_k|$, and the involution $k\mapsto k'$ with $g\_k^{-1}\in C\_{k'}$. Class identification (deciding which $C\_k$ contains a given $x\in G$) is by cycle-type fingerprint followed by conjugacy backtrack. The descent certificate is not consumed in this step beyond having proved that $G$ is the Galois group. $K=\mathbb{Q}(\chi)$ throughout; $\zeta\_m$ is a primitive $m$-th root of unity.

## 1. Class algebra and structure constants

Let $K\_i=\sum\_{x\in C\_i}x\in Z(\mathbb{C}G)$ and $K\_iK\_j=\sum\_k a\_{ijk}K\_k$ with

$$a\_{ijk}=#{(x,y)\in C\_i\times C\_j:\ xy=g\_k}=#{x\in C\_i:\ x^{-1}g\_k\in C\_j}\in\mathbb{Z}\_{\ge0}.$$

For fixed $i$ and $k$, one enumeration of $C\_i$ with one class identification per element yields $a\_{ijk}$ for all $j$; so the class matrix $M\_i=(a\_{ijk})\_{j,k}$ costs $r|C\_i|$ identifications, and only the $M\_i$ with small $|C\_i|$ are ever computed (§2.4).

For $\chi\in\mathrm{Irr}(G)$ the central character $\omega\_\chi(K\_k)=|C\_k|\chi(g\_k)/\chi(1)$ lies in $\mathbb{Z}[\zeta\_e]$, and $\omega\_\chi\:Z(\mathbb{C}G)\to\mathbb{C}$ is an algebra homomorphism; the $r$ homomorphisms $Z(\mathbb{C}G)\to\mathbb{C}$ are exactly the $\omega\_\chi$. Writing $\omega\_\chi$ as the column vector $(\omega\_\chi(K\_k))*k$, the identity $\omega*\chi(K\_i)\omega\_\chi(K\_j)=\sum\_k a\_{ijk}\omega\_\chi(K\_k)$ reads

$$M\_i,\omega\_\chi=\omega\_\chi(K\_i),\omega\_\chi\qquad(1\le i\le r).$$

The $\omega\_\chi$ are the common eigenvectors of the commuting family ${M\_i}$, and the table is a simultaneous-eigenvector problem. Dixon to Schneider solves it modulo a prime and lifts.

## 2. Dixon to Schneider modulo $p$

### 2.1 Choice of the prime

**Proposition 2.1.** Let $p$ be a prime with $p\equiv1\pmod e$ and $p>2\sqrt{|G|}$, and let $z\in\mathbb{F}\_p^\times$ have order $e$ (it exists since $e\mid p-1$). Let $\theta:\mathbb{Z}[\zeta\_e]\to\mathbb{F}\_p$, $\zeta\_e\mapsto z$, with kernel a degree-one prime $\mathfrak{p}\mid p$. Then:

1. $p\nmid|G|$.
2. The reduced vectors $\bar\omega\_\chi=\theta(\omega\_\chi)\in\mathbb{F}\_p^r$, $\chi\in\mathrm{Irr}(G)$, form a basis of $\mathbb{F}*p^r$ consisting of common eigenvectors of the $\bar M\_i$, with pairwise distinct joint eigenvalue tuples $(\bar\omega*\chi(K\_i))\_i$. Hence every joint eigenspace of ${\bar M\_i}$ is one-dimensional.
3. $\chi(1)$ is recovered from $\bar\omega\_\chi$ as the unique positive integer $d\le\sqrt{|G|}$ with

$$d^2\equiv|G|\Big(\sum\_k\bar\omega\_\chi(K\_k)\bar\omega\_\chi(K\_{k'})/|C\_k|\Big)^{-1}\pmod p.$$

4. Every character value lifts uniquely: for $g$ of order $o$ and $j\in\mathbb{Z}/o$, $m\_j:=\frac1o\sum\_{t=0}^{o-1}\chi(g^t)\zeta\_o^{-jt}$ is the multiplicity of the eigenvalue $\zeta\_o^j$ of $\rho(g)$, so $0\le m\_j\le\chi(1)\<p$; $\theta(m\_j)=o^{-1}\sum\_t\theta(\chi(g^t)),z\_o^{-jt}$ with $z\_o=z^{e/o}$; $m\_j$ is the unique representative of $\theta(m\_j)$ in $[0,\chi(1)]$; and $\chi(g)=\sum\_j m\_j\zeta\_o^j$.

*Proof.* (1) Every prime dividing $|G|$ divides $e\<p$. (2) Let $X=(\chi(g\_k))*{\chi,k}$. Column orthogonality gives $\bar X^{T}X=\mathrm{diag}(|C\_G(g\_k)|)$ and $\bar X$ is $X$ with rows permuted, so $(\det X)^2=\pm\prod\_k|C\_G(g\_k)|$, an integer prime to $p$; hence $\det X\notin\mathfrak{p}$. The matrix $(\omega*\chi(K\_k))*{\chi,k}$ is $X$ row-scaled by $\chi(1)^{-1}$ and column-scaled by $|C\_k|$, both prime to $p$, so its reduction is invertible: the $\bar\omega*\chi$ are $r$ linearly independent vectors, a basis. The eigenvalue of $\bar M\_i$ on $\bar\omega\_\chi$ is the $i$-th coordinate of $\bar\omega\_\chi$, so the joint eigenvalue tuple of $\bar\omega\_\chi$ is the vector itself; distinct vectors give distinct tuples. A basis of joint eigenvectors with pairwise distinct tuples forces each joint eigenspace to be the span of one of them. (3) Substituting $\chi(g\_k)=\chi(1)\omega\_\chi(K\_k)/|C\_k|$ and $|C\_{k'}|=|C\_k|$ into $\sum\_k|C\_k|\chi(g\_k)\chi(g\_k^{-1})=|G|$ gives $\chi(1)^2=|G|/\sum\_k\omega\_\chi(K\_k)\omega\_\chi(K\_{k'})/|C\_k|$; the denominator equals $|G|/\chi(1)^2$, nonzero mod $p$. If $d\ne d'$ in $[1,\sqrt{|G|}]$ had $d^2\equiv d'^2$, then $p\mid(d-d')(d+d')$ with $0<|d-d'|\<d+d'<2\sqrt{|G|}\<p$. (4) The multiplicity statement is the inverse discrete Fourier transform on $\langle g\rangle$; $\theta(\chi(g^t))=\chi(1)\bar\omega\_\chi(K\_{k(t)})/|C\_{k(t)}|$ where $C\_{k(t)}\ni g^t$ is given by the power map; the interval $[0,\chi(1)]$ has length $\<p$. $\square$

Any prime $p\equiv 1\pmod e$ beyond $2\sqrt{|G|}$ works. Replacing $z$ by $z^t$ composes $\theta$ with $\sigma\_t:\zeta\_e\mapsto\zeta\_e^t$ and permutes the rows of the lifted table by the Galois action, nothing else.

### 2.2 The computation

1. Power maps: for each class $k$ and each $0\le t\<o(g\_k)$, identify the class of $g\_k^t$. (It suffices to do $t\mid o(g\_k)$ and $t$ coprime to $o(g\_k)$ and compose.)
2. Maintain a partition of $\mathbb{F}\_p^r$ into joint eigenspaces of the matrices used so far, initially ${\mathbb{F}\_p^r}$. While some part $V$ has $\dim V>1$: pick a class $i$ with $|C\_i|$ small and $\bar M\_i$ not yet used on $V$; $V$ is $\bar M\_i$-stable because the $M\_i$ commute; compute $\bar M\_i|*V$ (this needs the entries $a*{ijk}$ for all $j$ and for $k$ in the support of a basis of $V$, obtained as in §1); its characteristic polynomial splits over $\mathbb{F}\_p$ by Prop. 2.1(2); replace $V$ by the eigenspaces of $\bar M\_i|\_V$.
3. A one-dimensional part is a joint eigenspace of the full family (each part is a sum of full joint eigenspaces), hence spanned by some $\bar\omega\_\chi$; normalize by $\bar\omega\_\chi(K\_1)=1$. Recover $\chi(1)$ by 2.1(3), the values $\theta(\chi(g\_k))$, and lift by 2.1(4).

### 2.3 Galois orbits

$\sigma\_t\in\mathrm{Gal}(\mathbb{Q}(\zeta\_e)/\mathbb{Q})$ acts by $\chi^{\sigma\_t}(g)=\chi(g^t)$, so

$$\omega\_{\chi^{\sigma\_t}}(K\_k)=\frac{|C\_k|}{|C\_{k^{(t)}}|},\omega\_\chi(K\_{k^{(t)}}),\qquad g\_k^t\in C\_{k^{(t)}}.$$

Modulo $\mathfrak{p}$ this is a permutation-and-rational-rescaling of the coordinates of $\bar\omega\_\chi$, computable from the power maps. Hence step 3 is done once per Galois orbit: from one lifted $\chi$ the whole orbit ${\chi^{\sigma\_t}}$ and its reduced central characters follow. A part $V$ of dimension $v$ that already contains $v$ known reduced central characters is spanned by them and needs no further splitting; this is checked by linear algebra over $\mathbb{F}\_p$ after each orbit is found. The table is closed under $\mathrm{Gal}(\mathbb{Q}(\zeta\_e)/\mathbb{Q})$ by construction.

### 2.4 Schneider's economies

Only class matrices for the classes actually used in step 2 are formed, and each only on the parts $V$ it is applied to. For a basis vector $b$ of $V$, $(M\_ib)*j=\sum\_k a*{ijk}b\_k$ needs $a\_{ijk}$ for all $j$ and $k\in\mathrm{supp}(b)$, i.e. $|C\_i|\cdot|\mathrm{supp}(b)|$ identifications. Classes with small $|C\_i|$ (large centralizers) are chosen first. $M\_1$ is the identity and never splits anything; for a central $g\_i$, $M\_i$ is a permutation matrix and splits parts according to the value of the central character on $g\_i$ at cost $|C\_i|=1$ per column.

### 2.5 Certificate

The lifted table is proved correct by an exact check that uses only the class matrices already formed.

**Proposition 2.2.** Let $S\subseteq{1,\dots,r}$ be nonempty and let $\omega^{(1)},\dots,\omega^{(r)}\in\mathbb{Z}[\zeta\_e]^r$ be such that (i) $M\_i\omega^{(\nu)}=\omega^{(\nu)}\_i,\omega^{(\nu)}$ for all $i\in S$ and all $\nu$, exactly in $\mathbb{Z}[\zeta\_e]$; (ii) the tuples $(\omega^{(\nu)}*i)*{i\in S}$ are pairwise distinct; (iii) $\omega^{(\nu)}*1=1$. Then ${\omega^{(\nu)}}={\omega*\chi:\chi\in\mathrm{Irr}(G)}$.

*Proof.* The family ${M\_i}*{i\in S}$ is commuting and diagonalizable over $\mathbb{C}$ (the $\omega*\chi$ diagonalize it), so its joint eigenvalue tuples form a multiset of size $r$ independent of the diagonalizing basis. Joint eigenvectors with distinct tuples are linearly independent, so by (i),(ii) the $\omega^{(\nu)}$ are a basis and their $r$ distinct tuples exhaust the multiset with multiplicity one: every joint eigenspace of ${M\_i}*{i\in S}$ is a line spanned by one $\omega^{(\nu)}$. Each $\omega*\chi$ lies in the joint eigenspace of its own tuple, hence is a scalar multiple of some $\omega^{(\nu)}$, and (iii) together with $\omega\_\chi(K\_1)=1$ makes the scalar $1$. $\square$

The set $S$ of classes used in §2.2 satisfies (ii) automatically (splitting terminated), so the certificate costs $|S|,r^3$ multiplications in $\mathbb{Z}[\zeta\_e]$, or $|S|,r\cdot\mathrm{nnz}(M\_i)$ using sparsity. Characters follow from central characters by $\chi(g\_k)=\chi(1)\omega\_\chi(K\_k)/|C\_k|$ with $\chi(1)$ the positive square root of $|G|/\sum\_k\omega\_\chi(K\_k)\omega\_\chi(K\_{k'})/|C\_k|$, an exact rational computation that must return a perfect square. The lifted values of §2.2 and the values from the certificate must agree; both are exact.

## 3. Fields of values and Schur indices

$K=\mathbb{Q}(\chi)\subseteq\mathbb{Q}(\zeta\_e)$ is the fixed field of the stabilizer of the row $\chi$ under the action of §2.3; its conductor $c$ (least $c$ with $K\subseteq\mathbb{Q}(\zeta\_c)$) is read off from that stabilizer as a subgroup of $(\mathbb{Z}/e)^\times$.

The central idempotent $e\_\chi=\frac{\chi(1)}{|G|}\sum\_g\chi(g^{-1})g$ lies in $KG$; $A\_\chi:=e\_\chi KG$ is central simple over $K$ of dimension $\chi(1)^2$ (since $A\_\chi\otimes\_K\bar{\mathbb{Q}}=e\_\chi\bar{\mathbb{Q}}G\cong M\_{\chi(1)}(\bar{\mathbb{Q}})$), so $A\_\chi\cong M\_k(D\_\chi)$ with $D\_\chi$ a division algebra of index $m\_\chi$, $km\_\chi=\chi(1)$. The unique irreducible $KG$-module $M\_\chi$ with $e\_\chi M\_\chi=M\_\chi$ has character $m\_\chi\chi$, and $m\_\chi=m\_K(\chi)$ is the least $m$ such that $m\chi$ is the character of a $KG$-module. For $F\supseteq K$, $m\_F(\chi)$ is the index of $D\_\chi\otimes\_KF$; it divides $m\_K(\chi)$.

Standard facts used below:

- (S1) $m\_\chi\mid\chi(1)$; $m\_{\mathbb{Q}(\zeta\_e)}(\chi)=1$ (Brauer); hence, by Lemma 3.1, $m\_\chi\mid[\mathbb{Q}(\zeta\_e)\:K]$.
- (S2) $\zeta\_{m\_\chi}\in K$ (Benard to Schacher). So $m\_\chi$ divides the number of roots of unity in $K$; in particular $m\_\chi\le2$ when $K\subseteq\mathbb{R}$ (Brauer to Speiser).
- (S3) $m\_\chi=\mathrm{lcm}*v,m\_v$ over the places $v$ of $K$, $m\_v$ the index of $D*\chi\otimes K\_v$; $m\_v=1$ for finite $v$ of residue characteristic $\ell\nmid|G|$; $m\_v$ is constant on Galois orbits of places (Benard to Schacher), so $m\_\chi=\mathrm{lcm}(m\_\infty,,m\_\ell:\ell\mid|G|)$ with $m\_\ell$ the common value over $\ell$.
- (S4) $m\_\infty$: complex conjugation acts on $K$ as $\sigma\_{-1}$, so $\chi\ne\bar\chi$ iff $K\not\subseteq\mathbb{R}$, in which case $m\_\infty=1$. If $\chi=\bar\chi$ then $m\_\infty=1$ or $2$ according as the Frobenius to Schur indicator $\nu\_2(\chi)=|G|^{-1}\sum\_g\chi(g^2)$ is $+1$ or $-1$.
- (S5) (Brauer to Witt) For $\ell\mid|G|$ there exist an $\ell$-quasi-elementary $H\le G$ and $\psi\in\mathrm{Irr}(H)$ with $\langle\chi\_H,\psi\rangle\not\equiv0\pmod\ell$ and $[K(\psi)\:K]$ prime to $\ell$ such that the $\ell$-parts of $m\_K(\chi)$ and $m\_{K(\psi)}(\psi)$ coincide. For any such pair the divisibility of the $\ell$-part of $m\_K(\chi)$ by that of $m\_{K(\psi)}(\psi)$ is Lemma 3.1; the theorem is the existence of a pair with equality. $H$ is an M-group, $A\_\psi\otimes K(\psi)$ is a cyclotomic algebra, and its Hasse invariants at the places over $\ell$ are given by explicit formulas (Yamada; Riese to Schmid). This is what Unger's Schur index algorithm evaluates.

**Lemma 3.1 (upper bound).** For $H\le G$, $\psi\in\mathrm{Irr}(H)$, $L=K(\psi)$:
$$m\_K(\chi)\ \mid\ \langle\chi\_H,\psi\rangle\cdot m\_L(\psi)\cdot[L\:K].$$

*Proof.* Let $U$ be an $LH$-module with character $m\_L(\psi)\psi$. Since $e\_\chi\in KG\subseteq LG$, the $\chi$-isotypic part $e\_\chi\mathrm{Ind}*H^GU$ is an $LG$-module with character $m\_L(\psi)\langle\chi\_H,\psi\rangle\chi$ (Frobenius reciprocity); every $LG$-module with character $a\chi$ is a multiple of $M*{\chi,L}$, so $m\_L(\chi)\mid m\_L(\psi)\langle\chi\_H,\psi\rangle$. An $LG$-module of character $m\_L(\chi)\chi$ regarded as a $K$-space has character $\mathrm{Tr}\_{L/K}\circ(m\_L(\chi)\chi)=m\_L(\chi)[L\:K]\chi$ because $\chi$ is $K$-valued, so $m\_K(\chi)\mid[L\:K],m\_L(\chi)$. $\square$

With $\psi=\lambda$ linear, $m\_L(\lambda)=1$ and the bound is $\langle\chi\_H,\lambda\rangle[K(\lambda)\:K]$; a pair with $\langle\chi\_H,\lambda\rangle=1$ and $\mathbb{Q}(\lambda)\subseteq K$ certifies $m\_\chi=1$.

**Procedure.** Compute $m\_\infty$ by (S4). For each $\ell\mid|G|$ compute $m\_\ell$ by (S5); this is exact. Independently, search pairs $(H,\lambda)$ as above for an upper bound, which is what §4 consumes. Two facts govern how $m\_\chi$ is used later: any certified multiple $m'$ of $m\_\chi$ may replace $m\_\chi$ in §4 and §5 with every statement remaining true (outputs are then divided or rooted by $m'$ instead of $m\_\chi$); and the exact value $m\_\chi$ is an output for the tables but is consumed nowhere else.

## 4. A matrix model of $m\_\chi\chi$

**Target.** $d=m\_\chi\chi(1)$; a homomorphism $\rho\:G\to\mathrm{GL}*d(\mathcal{O}F)$ with $\mathrm{tr}\rho=m\chi\chi$, where $F$ is a cyclotomic field containing $K$: $F=\mathbb{Q}(\zeta\_c)$ when the descent of 4.3 is performed, $F=\mathbb{Q}(\zeta*{\mathrm{lcm}(c,o)})$ otherwise (4.2). Once $\rho$ exists over $F$ with any denominators, integrality is a lattice choice: for $\Lambda\_0=\mathcal{O}*F^d$, $\Lambda=\sum*{g\in G}\rho(g)\Lambda\_0$ is a full $G$-stable $\mathcal{O}\_F$-lattice. If $\Lambda$ is free (always for the monomial models of 4.2, and always when $h(F)=1$) a basis of $\Lambda$ gives $\rho(G)\subseteq\mathrm{GL}*d(\mathcal{O}F)$. If not, use a pseudo-basis (entries then lie in fractional ideals determined by the Steinitz class) or restrict scalars to $\mathbb{Z}$ (a free $\mathbb{Z}$-lattice of rank $d[F:\mathbb{Q}]$ with character $m\chi\mathrm{Tr}*{F/\mathbb{Q}}\chi$).

**Lemma 4.1 (projector).** Let $W\_0$ be an $FG$-module with character $\theta$ and let $E\_\chi=\frac{\chi(1)}{|G|}\sum\_g\chi(g^{-1})\rho\_0(g)$ act on $W\_0$. Then $E\_\chi$ is an idempotent commuting with $G$, and $W=E\_\chi W\_0$ is an $FG$-module with character $\langle\theta,\chi\rangle\chi$. If $W\_0=\mathrm{Ind}*H^G\lambda$ with $\lambda$ linear, $\langle\theta,\chi\rangle=\langle\chi\_H,\lambda\rangle$. In all cases $W\cong M*{\chi,F}^{,a/m\_F(\chi)}$ with $a=\langle\theta,\chi\rangle$.

*Proof.* $E\_\chi$ is the image of $e\_\chi$ and projects onto the $\chi$-isotypic component, whose character is $\langle\theta,\chi\rangle\chi$; Frobenius reciprocity gives the multiplicity for induced modules; $e\_\chi W=W$ makes $W$ a module over the simple algebra $A\_\chi\otimes\_KF$, all of whose modules are multiples of the simple one. $\square$

### 4.2 Monomial route

Search $H\le G$ of index $d$ and linear $\lambda\in\mathrm{Hom}(H/H',\mu\_e)$ with $\langle\chi\_H,\lambda\rangle=m\_\chi$. Then $\mathrm{Ind}*H^G\lambda$ has degree $d$ and contains $\chi$ with multiplicity $m*\chi$; $m\_\chi$ copies of $\chi$ already have degree $d$, so $\mathrm{Ind}*H^G\lambda=m*\chi\chi$. In the basis of cosets $G/H$, $\rho(g)e\_{xH}=\lambda(h)e\_{yH}$ where $gx=yh$: $\rho(g)$ is monomial with nonzero entries in $\mu\_o$, $o=\mathrm{ord}\lambda$, the lattice $\mathbb{Z}[\zeta\_o][G/H]$ is free and stable, and no denominators appear. The subgroup search uses the subgroup structure of $G$ (for solvable $G$, the normal series); the conditions $[G\:H]=d$ and $\langle\chi\_H,\lambda\rangle=m\_\chi$ are checked from the table.

If $\mathbb{Q}(\zeta\_o)\subseteq\mathbb{Q}(\zeta\_c)$ the model is over $\mathbb{Z}[\zeta\_c]$ already. Otherwise it is over $\mathbb{Z}[\zeta\_{\mathrm{lcm}(c,o)}]$, which is accepted as the target ring; a model over $\mathbb{Q}(\zeta\_c)$ of $m\_{\mathbb{Q}(\zeta\_c)}(\chi)\chi$ exists (since $m\_{\mathbb{Q}(\zeta\_c)}(\chi)\mid m\_\chi$) and is obtained by the descent step of 4.3 when wanted.

### 4.3 Projector-and-splitting route (Unger)

Choose $(H,\lambda)$ with $a=\langle\chi\_H,\lambda\rangle$ as small as possible and form $W=E\_\chi\mathrm{Ind}*H^G\lambda$ over $F=\mathbb{Q}(\zeta*{\mathrm{lcm}(c,o)})$ (Lemma 4.1): a basis of the column space of the matrix $E\_\chi$, then $\rho(g)$ is the action of $g$ in that basis. $\dim\_FW=a\chi(1)$ and $\mathrm{tr}=a\chi$.

- If $a=m\_F(\chi)$, $W$ is the model. ($H=1$ gives $a=\chi(1)$. If $F$ is enlarged to contain $\zeta\_e$ then $m\_F=1$ and a pair with $a=1$ yields a model of $\chi$ itself by projection alone; this is Dixon's 1993 construction.)
- If $a>m\_F(\chi)$, compute $E=\mathrm{End}*{FG}(W)$ as the solution space of $X\rho(g\_j)=\rho(g\_j)X$ over the generators; $\dim\_FE=(a/m\_F)^2m\_F^2=a^2$ and $E\cong M*{a/m\_F}(D\_\chi\otimes\_KF)$. A primitive idempotent $\varepsilon\in E$ has $\varepsilon W$ of character $m\_F(\chi)\chi$, the model. Finding $\varepsilon$ is the explicit-isomorphism problem for the central simple $F$-algebra $E$ of degree $a$; Unger's splitting (and Fieker's descent algorithms) solve it by locating zero divisors through subfields of $E$ and norm equations in the corresponding cyclic algebras. Splitting along $X\in E$ with reducible minimal polynomial $f\_1f\_2$ gives the $G$-stable decomposition $W=\ker f\_1(X)\oplus\ker f\_2(X)$; iterate until a summand of trace $m\_F(\chi)\chi$ appears.

**Lemma 4.2 (maximal subfields give models).** Let $L\subseteq E$ be a maximal subfield, so $[L\:F]=a$, generated by some $X\in E$ with minimal polynomial of degree $a$. Then $W$, regarded as an $L$-space via $X$, is an $LG$-module of $L$-dimension $\chi(1)$ with character $\chi$.

*Proof.* Maximal subfields of a central simple algebra of degree $a$ over an infinite field have degree $a$ (the centralizer of a smaller subfield $L\_0$ is central simple over $L\_0$ of degree $>1$ and contains a proper extension of $L\_0$). $G$ commutes with $L$, so $W$ is an $LG$-module, of $L$-dimension $a\chi(1)/a$. Let $\chi\_L$ be its $L$-character. Since $\mathrm{tr}*F=\mathrm{Tr}*{L/F}\circ\mathrm{tr}*L$, $\sum*{\tau}\chi\_L^{\tau}=a\chi$, the sum over the $a$ embeddings $\tau$ of $L$ over $F$. Every irreducible constituent $\eta$ of $\chi\_L$ therefore satisfies $\eta^\tau=\chi$ for some $\tau$, so $\eta=\chi^{\tau^{-1}}=\chi$ as $\chi$ is $F$-valued; degree $\chi(1)$ forces $\chi\_L=\chi$. $\square$

This is the alternative to splitting: a model of $\chi$ over an extension $L/F$ of degree $a$, obtained by linear algebra only. It also exhibits the Schur index: for $a=m\_F(\chi)$, $E=D\_\chi\otimes F$ and $L$ is a maximal subfield of a division algebra.

**Descent to $\mathbb{Q}(\zeta\_c)$.** Let the model $\rho$ of $m\chi$ sit over $F'\supsetneq F=\mathbb{Q}(\zeta\_c)$ with $m=m\_F(\chi)$, and write $F'/F$ as a tower of cyclic steps of prime degree $q$ with generator $\sigma$. At each step: $m\chi$ is $F$-valued, so $\rho^\sigma\cong\rho$; solve the linear system $A\rho(g\_j)=\rho^\sigma(g\_j)A$ for an invertible $A$. Define the $\sigma$-semilinear map $\Phi=A^{-1}\circ\sigma$ on $F'^{,d}$ (apply $\sigma$ coordinatewise, then $A^{-1}$). Then $\Phi\rho(g)=\rho(g)\Phi$ for all $g$, and $\Phi^q$ is $F'$-linear and commutes with $G$, so $\Phi^q=B\in E'^\times$, $E'=\mathrm{End}*{F'G}(\rho)$. An $F$-form of $\rho$ is the same as a $G$-commuting $\sigma$-semilinear $\Phi'$ with $\Phi'^q=1$ (its fixed points are an $F$-form by Galois descent, and a basis $P$ of the fixed points gives $P^{-1}\rho P$ with entries in $F$). Every such $\Phi'$ is $\Phi u$ with $u\in E'^\times$, and $(\Phi u)^q=B\cdot u^{\Phi^{q-1}}\cdots u^{\Phi}u$ with $u^{\Phi}=\Phi u\Phi^{-1}$; so the step is a twisted norm equation $N*\Phi(u)=B^{-1}$ in $E'$, the same kernel as the splitting above. It is solvable at every step because an $F$-model of $m\chi$ exists.

## 5. The outputs depend only on $\chi$

Fix a prime $\ell$, a prime $\mathfrak{L}\mid\ell$ of the splitting field $N$, decomposition group $D=D\_{\mathfrak{L}}\le G$, inertia $I=G\_0\trianglelefteq D$, lower ramification groups $G\_i\trianglelefteq D$, and any $\mathrm{Fr}\in D$ mapping to the arithmetic Frobenius of $D/I$. For a class function $\chi$ on $G$ define

$$P\_\ell(\chi;T)=\exp\Big(-\sum\_{k\ge1}\frac{T^k}{k}\cdot\frac1{|I|}\sum\_{\sigma\in I}\chi(\mathrm{Fr}^k\sigma)\Big),\qquad f\_\ell(\chi)=\sum\_{i\ge0}\frac{|G\_i|}{|G\_0|}\Big(\chi(1)-\frac1{|G\_i|}\sum\_{\sigma\in G\_i}\chi(\sigma)\Big),$$

$$a\_\pm(\chi)=\tfrac12\big(\chi(1)\pm\chi(c)\big),\quad c=\text{complex conjugation (an involution class of $G$).}$$

**Theorem 5.1.** Let $\rho\:G\to\mathrm{GL}(V)$ be any representation over a field of characteristic $0$ (embedded in $\mathbb{C}$ or $\bar{\mathbb{Q}}\_\ell$ as needed) with $\mathrm{tr}\rho=\chi$, $\chi$ any character, not necessarily irreducible. Then:

1. $\det(1-\rho(\mathrm{Fr})T\mid V^I)=P\_\ell(\chi;T)$; so $L\_\ell(s,\chi)=P\_\ell(\chi;\ell^{-s})^{-1}$.
2. The Artin conductor exponent of $\rho$ at $\ell$ is $f\_\ell(\chi)$.
3. The archimedean factor is $\Gamma\_{\mathbb{R}}(s)^{a\_+(\chi)}\Gamma\_{\mathbb{R}}(s+1)^{a\_-(\chi)}$.
4. The local $\varepsilon$-factors $\varepsilon\_\ell(\rho|*D,\psi*\ell)$, hence the global root number $W(\chi)$, depend on $\rho$ only through $\chi$.
5. $P\_\ell(\chi;T)\in\mathcal{O}*K[T]$, $f*\ell(\chi)\in\mathbb{Z}$, and $P\_\ell(\chi^\tau;T)=\tau(P\_\ell(\chi;T))$, $f\_\ell(\chi^\tau)=f\_\ell(\chi)$ for $\tau\in\mathrm{Gal}(K/\mathbb{Q})$.

None of the quantities in (1) to (4) depends on the choice of $\mathfrak{L}$ or of the lift $\mathrm{Fr}$.

*Proof.* For $H\trianglelefteq D$ and $x\in D$, $e\_H=|H|^{-1}\sum\_{h\in H}\rho(h)$ is the projector onto $V^H$ and commutes with $\rho(x)$; so $\rho(x)e\_H$ acts as $\rho(x)|*{V^H}$ on $V^H$ and as $0$ on $(1-e\_H)V$, giving $\mathrm{tr}(\rho(x)|*{V^H})=|H|^{-1}\sum\_h\chi(xh)$. (1): apply this with $H=I$, $x=\mathrm{Fr}^k$, in $\log\det(1-AT)=-\sum\_k\mathrm{tr}(A^k)T^k/k$ for $A=\rho(\mathrm{Fr})|*{V^I}$. (2): Artin's definition $f*\ell=\sum\_i[G\_0\:G\_i]^{-1}\dim V/V^{G\_i}$ with $\dim V^{G\_i}=\mathrm{tr}(e\_{G\_i})$ (case $x=1$). (3): $a\_\pm$ are the dimensions of the $\pm1$ eigenspaces of $\rho(c)$. (4): Deligne's theorem defines $\varepsilon(\cdot,\psi\_\ell,dx)$ on the Grothendieck group of representations of $\mathrm{Gal}(\bar{\mathbb{Q}}*\ell/\mathbb{Q}*\ell)$, additively and inductively in degree $0$; a representation with finite image factoring through $D$ is determined in the Grothendieck group by its character; the archimedean $\varepsilon$ depends only on $a\_-(\chi)$; $W(\chi)$ is the normalized product. Independence of $\mathrm{Fr}$: if $\mathrm{Fr}'=\mathrm{Fr}\tau$ with $\tau\in I$, then $\mathrm{Fr}'^k\in\mathrm{Fr}^kI$ because $I\trianglelefteq D$, so the inner sums in $P\_\ell$ coincide. Independence of $\mathfrak{L}$: replacing $\mathfrak{L}$ by $g\mathfrak{L}$ conjugates $(D,I,G\_i,\mathrm{Fr})$ by $g$, and $\chi$ is a class function. (5): the coefficients of $\det(1-\rho(\mathrm{Fr})T\mid V^I)$ are elementary symmetric functions of roots of unity, hence algebraic integers; each $|G\_i|^{-1}\sum\_{G\_i}\chi=\langle\chi\_{G\_i},1\rangle$ is an integer; the exponential formula has coefficients in the field generated by the values of $\chi$, so $P\_\ell\in K[T]$, and applying $\tau$ to that formula gives the equivariance. $\square$

**Corollary 5.2 (the Schur index never enters).** Conductors, Euler factors at all primes, gamma factors, and root numbers are functions of the row $\chi$ of the character table together with the subgroups $D,I,G\_i$ and the elements $\mathrm{Fr}$, $c$ of $G$. If a model $\rho$ of $m\chi$ is used instead ($m=m\_\chi$ or any multiple), it returns $P\_\ell(\chi;T)^m$, $m f\_\ell(\chi)$, $m,a\_\pm(\chi)$ and $W(\chi)^m$: the first three recover the invariants of $\chi$ uniquely ($1+TK[[T]]$ is uniquely $m$-divisible in characteristic $0$), the last only up to $\mu\_m$.

Frobenius classes (the next step) are likewise read from the action of $G$ on the roots and never from a model.

**Root numbers.** $W(\chi)$ is a function of $\chi$ by 5.1(4), so the residual $\mu\_m$ ambiguity of a model-based computation is a defect of the model, not of the output, and is removed as follows.

- If $\chi=\bar\chi$ and $\nu\_2(\chi)=+1$ (orthogonal $\chi$; this covers every real $\chi$ with $m\_\infty=1$, whatever $m\_\chi$ is), then $W(\chi)=1$ (Fröhlich to Queyrut 1973; Deligne 1976 for the local constants).
- Otherwise the candidates are $w\_0\mu\_m$ where $w\_0^m$ is the model value, and the functional equation selects one:

**Lemma 5.3.** Assume $\Lambda(s,\chi)=\mathfrak{f}(\chi)^{s/2}\Gamma\_\mathbb{R}(s)^{a\_+}\Gamma\_\mathbb{R}(s+1)^{a\_-}L(s,\chi)$ continues meromorphically. Then at most one $w\in\mathbb{C}$ satisfies $\Lambda(s,\chi)=w,\Lambda(1-s,\bar\chi)$ identically.

*Proof.* Two values $w\ne w'$ force $\Lambda(1-s,\bar\chi)\equiv0$, hence $\Lambda(s,\bar\chi)\equiv0$, but $L(s,\bar\chi)\ne0$ for $\mathrm{Re},s>1$ (absolutely convergent Euler product). $\square$

The numerical verification is against this finite candidate set, at a real argument where the smoothed sum is provably nonzero, with truncation bounds from $|a\_n(\chi)|\le d\_{\chi(1)}(n)$: the coefficient of $T^k$ in $\prod\_{i\le\chi(1)}(1-\lambda\_iT)^{-1}$ with $|\lambda\_i|\le1$ is at most $\binom{k+\chi(1)-1}{\chi(1)-1}$, and multiplicativity gives $d\_{\chi(1)}(n)$, the number of ordered factorizations of $n$ into $\chi(1)$ positive factors. The computed discrepancy for the wrong candidates must exceed the error bound. The same check is the functional-equation verification demanded of the final output. The exact alternative, Deligne's inductive formula for $\varepsilon\_\ell$ on the solvable group $D$ via linear characters of its subgroups and Tate's local constants in the completions $N\_{\mathfrak{L}}^H$, belongs to the ramified-prime step and again uses only $\chi|\_D$.

Two consequences for the test families: for abelian $G$ the formulas above collapse to the Dirichlet-character formulas, since $D$ is cyclic and the $|I|$-averages pick out the characters trivial on inertia; and for the weight-one comparison, the coefficients $a\_n(\chi)$ lie in $\mathcal{O}\_K$ by 5.1(5) and are Galois-equivariant in $\chi$, which is what makes coefficient-by-coefficient comparison with a newform whose Hecke field is $K$ meaningful without any choice of model.