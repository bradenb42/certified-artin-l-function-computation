# The primes ramified in $N/\mathbb{Q}$

**Setting.** $f=\prod\_{i=1}^{s}f\_i$ with $f\_i\in\mathbb{Z}[x]$ monic irreducible, $n\_i=\deg f\_i$, $K\_i=\mathbb{Q}[x]/(f\_i)=\mathbb{Q}(\alpha\_i)$, $\mathcal{O}\_i$ the maximal order of $K\_i$, $A=\mathbb{Q}[x]/(f)\cong\prod\_iK\_i$ with maximal order $\mathcal{O}\_A=\prod\_i\mathcal{O}\_i$, $\mathrm{disc}(\mathcal{O}\_A)=\prod\_i\mathrm{disc}(\mathcal{O}*i)$. $N$ is the splitting field of $f$, the compositum inside $\bar{\mathbb{Q}}$ of all conjugates $K\_i^\sigma=\sigma(K\_i)$, $\sigma\in\mathrm{Hom}(K\_i,\bar{\mathbb{Q}})$. A prime $\ell$ is unramified in a number field $L$ if every prime of $L$ above $\ell$ has ramification index $1$; equivalently $\mathcal{O}L/\ell\mathcal{O}L$ is a product of fields; equivalently $L\otimes\mathbb{Q}\mathbb{Q}\ell$ is a product of unramified extensions of $\mathbb{Q}*\ell$. It is *unramified in the étale algebra $A$* if it is unramified in every $K\_i$.

## 1. The criterion

**Theorem 1.1.** For a prime $\ell$ the following are equivalent:

1. $\ell$ is unramified in $N$.
2. $\ell$ is unramified in $K\_i$ for every $i$ (i.e. in $A$).
3. $\ell\nmid\mathrm{disc}(\mathcal{O}\_i)$ for every $i$, i.e. $\ell\nmid\mathrm{disc}(\mathcal{O}\_A)$.

Hence the ramified primes of $N$ are exactly the primes dividing $\prod\_i\mathrm{disc}(\mathcal{O}\_i)$.

The proof uses three lemmas.

**Lemma 1.2 (Dedekind).** $\ell$ is unramified in $L$ iff $\ell\nmid\mathrm{disc}(\mathcal{O}\_L)$.

*Proof.* Let $\bar{\mathcal{O}}=\mathcal{O}\_L/\ell\mathcal{O}*L$, an $\mathbb{F}*\ell$-algebra of dimension $[L:\mathbb{Q}]$. The discriminant of the trace form of $\mathcal{O}*L$ on a $\mathbb{Z}$-basis reduces mod $\ell$ to the discriminant of the trace form of $\bar{\mathcal{O}}$ over $\mathbb{F}*\ell$ on the reduced basis, so $\ell\nmid\mathrm{disc}(\mathcal{O}\_L)$ iff the trace form of $\bar{\mathcal{O}}$ is nondegenerate. Write $\ell\mathcal{O}\_L=\prod\_j\mathfrak{p}\_j^{e\_j}$, so $\bar{\mathcal{O}}\cong\prod\_j\mathcal{O}\_L/\mathfrak{p}\_j^{e\_j}$ and the trace form is the orthogonal sum of the trace forms of the factors. If some $e\_j\ge2$, the factor $\mathcal{O}\_L/\mathfrak{p}*j^{e\_j}$ has a nonzero nilpotent $y$, and $\mathrm{Tr}(yx)=0$ for all $x$ since $yx$ is nilpotent, so the form is degenerate. If all $e\_j=1$, each factor is a finite field, separable over $\mathbb{F}*\ell$, whose trace form is nondegenerate. $\square$

**Lemma 1.3 (subfields).** If $\ell$ is unramified in $L$ and $L'\subseteq L$, then $\ell$ is unramified in $L'$.

*Proof.* For $\mathfrak{P}\mid\mathfrak{p}\mid\ell$ with $\mathfrak{P}$ a prime of $L$ and $\mathfrak{p}=\mathfrak{P}\cap\mathcal{O}\_{L'}$, $e(\mathfrak{P}/\ell)=e(\mathfrak{P}/\mathfrak{p}),e(\mathfrak{p}/\ell)$. $\square$

**Lemma 1.4 (conjugates and compositum).** (a) If $\ell$ is unramified in $L$ then it is unramified in every conjugate $\sigma(L)$. (b) If $\ell$ is unramified in $L\_1$ and in $L\_2$ (both inside $\bar{\mathbb{Q}}$) then it is unramified in $L\_1L\_2$.

*Proof.* (a) $\sigma\:L\to\sigma(L)$ is a $\mathbb{Q}$-isomorphism carrying $\mathcal{O}*L$ to $\mathcal{O}*{\sigma(L)}$; discriminants agree. (b) Multiplication $L\_1\otimes\_\mathbb{Q}L\_2\to L\_1L\_2$ is a surjective algebra map, so after $\otimes\_\mathbb{Q}\mathbb{Q}*\ell$, $L\_1L\_2\otimes\mathbb{Q}*\ell$ is a quotient of $(L\_1\otimes\mathbb{Q}*\ell)\otimes*{\mathbb{Q}*\ell}(L\_2\otimes\mathbb{Q}*\ell)=\prod\_{v,w}L\_{1,v}\otimes\_{\mathbb{Q}*\ell}L*{2,w}$. Each $L\_{1,v},L\_{2,w}$ is an unramified extension $\mathbb{Q}*{\ell^a},\mathbb{Q}*{\ell^b}$ of $\mathbb{Q}*\ell$, and $\mathbb{Q}*{\ell^a}\otimes\_{\mathbb{Q}*\ell}\mathbb{Q}*{\ell^b}\cong\mathbb{Q}*{\ell^{\mathrm{lcm}(a,b)}}^{\gcd(a,b)}$ (the unramified extension of degree $a$ is $\mathbb{Q}*\ell(\zeta\_{\ell^a-1})$, and its minimal polynomial factors over $\mathbb{Q}*{\ell^b}$ into $\gcd(a,b)$ factors of degree $a/\gcd(a,b)$). A quotient of a finite product of fields is a product of a subset of them. So $L\_1L\_2\otimes\mathbb{Q}*\ell$ is a product of unramified fields. $\square$

*Proof of Theorem 1.1.* (2)$\Leftrightarrow$(3) is Lemma 1.2 applied to each $K\_i$. (1)$\Rightarrow$(2): $K\_i\subseteq N$, Lemma 1.3. (2)$\Rightarrow$(1): $N$ is the compositum of the finitely many fields $K\_i^\sigma$; each is unramified at $\ell$ by Lemma 1.4(a); induction on the number of factors with Lemma 1.4(b). $\square$

**Corollary 1.5.** $\mathrm{disc}(f)=\prod\_i\mathrm{disc}(f\_i)\cdot\prod\_{i\<j}\mathrm{Res}(f\_i,f\_j)^2$ and $\mathrm{disc}(f\_i)=[\mathcal{O}*i:\mathbb{Z}[\alpha\_i]]^2,\mathrm{disc}(\mathcal{O}i)$. Hence every ramified prime divides $\prod\_i\mathrm{disc}(f\_i)$, and a prime dividing $\mathrm{disc}(f)$ only through the resultants $\mathrm{Res}(f\_i,f\_j)$ is unramified. Moreover $v\ell(\mathrm{disc},\mathcal{O}i)\equiv v\ell(\mathrm{disc},f\_i)\pmod 2$, so $\ell$ ramifies in $K\_i$, hence in $N$, whenever $v*\ell(\mathrm{disc},f\_i)$ is odd.

So the candidate set is $\mathcal{C}=\bigcup\_i{\ell:\ell\mid\mathrm{disc}(f\_i)}$, obtained from the factorizations of the $s$ integers $\mathrm{disc}(f\_i)=\pm\mathrm{Res}(f\_i,f\_i')$, not of $\mathrm{disc}(f)$; the resultants $\mathrm{Res}(f\_i,f\_j)$ are never factored. Candidates with odd valuation are ramified without further work. The remaining candidates, those with $v\_\ell(\mathrm{disc}f\_i)\ge2$ even for every $i$ with $\ell\mid\mathrm{disc}f\_i$, are decided by the procedure of §2. Complete factorization of each $\mathrm{disc}(f\_i)$ is required: for a prime in an unfactored cofactor the procedure cannot be run, and no other criterion here decides it.

## 2. Decision at a candidate prime by the $\ell$-maximal order

Fix $i$, write $n=n\_i$, $\alpha=\alpha\_i$, $\mathcal{O}=\mathbb{Z}[\alpha]$, $v=v\_\ell(\mathrm{disc}f\_i)$. An order $\mathcal{O}'\supseteq\mathbb{Z}[\alpha]$ of $K\_i$ is *$\ell$-maximal* if $\ell\nmid[\mathcal{O}*i:\mathcal{O}']$. Then $v*\ell(\mathrm{disc},\mathcal{O}')=v\_\ell(\mathrm{disc},\mathcal{O}\_i)$ and $\mathcal{O}'/\ell\mathcal{O}'\cong\mathcal{O}\_i/\ell\mathcal{O}\_i$, so by Lemma 1.2:

**Proposition 2.1.** If $\mathcal{O}'$ is $\ell$-maximal then $\ell$ is unramified in $K\_i$ iff $\mathcal{O}'/\ell\mathcal{O}'$ is reduced (its nilradical is $0$), iff $v\_\ell(\mathrm{disc},\mathcal{O}')=0$, where $v\_\ell(\mathrm{disc},\mathcal{O}')=v-2v\_\ell([\mathcal{O}':\mathbb{Z}[\alpha]])$.

### 2.1 Fast path: Dedekind's criterion

Factor $\bar f\_i=\prod\_j\bar g\_j^{e\_j}$ in $\mathbb{F}*\ell[x]$ with $\bar g\_j$ distinct monic irreducible, lift to monic $g\_j\in\mathbb{Z}[x]$, and set $h=(f\_i-\prod\_jg\_j^{e\_j})/\ell\in\mathbb{Z}[x]$. Then $\mathbb{Z}[\alpha]$ is $\ell$-maximal iff $\gcd(\bar h,\bar g\_j)=1$ in $\mathbb{F}*\ell[x]$ for every $j$ with $e\_j\ge2$ (Dedekind's criterion). If so, $\ell\mathcal{O}*i=\prod\_j(\ell,g\_j(\alpha))^{e\_j}$ and $\ell$ is ramified in $K\_i$ iff some $e\_j\ge2$. Cost: one factorization in $\mathbb{F}*\ell[x]$ of degree $n$ and one gcd computation per $j$, $O(n^2\log\ell+n^3)$ operations in $\mathbb{F}\_\ell$. This settles every candidate for which $\mathbb{Z}[\alpha\_i]$ is already $\ell$-maximal, in particular every $\ell$ with $v=1$ (redundantly with Corollary 1.5).

### 2.2 The $\ell$-maximal order by Round 2

For an order $\mathcal{O}$ let $I\_\ell(\mathcal{O})={x\in\mathcal{O}\:x^m\in\ell\mathcal{O}\text{ for some }m}$, the preimage of the nilradical of $\bar{\mathcal{O}}=\mathcal{O}/\ell\mathcal{O}$, and let $\mathcal{O}'=[I\_\ell\:I\_\ell]={x\in K\_i\:xI\_\ell\subseteq I\_\ell}$, the multiplier ring.

**Theorem 2.2 (Pohst to Zassenhaus).** $\mathcal{O}\subseteq\mathcal{O}'\subseteq\mathcal{O}\_i$, $\ell\mathcal{O}'\subseteq\mathcal{O}$, and $\mathcal{O}$ is $\ell$-maximal iff $\mathcal{O}'=\mathcal{O}$.

*Proof.* $I\_\ell$ is a full-rank $\mathcal{O}$-submodule of $\mathcal{O}$ (it contains $\ell\mathcal{O}$), so its multiplier ring is an order containing $\mathcal{O}$, hence inside $\mathcal{O}*i$. Since $\ell\in I*\ell$, $x\in\mathcal{O}'$ gives $\ell x\in I\_\ell\subseteq\mathcal{O}$; so $[\mathcal{O}':\mathcal{O}]$ divides $\ell^n$, and if $\mathcal{O}$ is $\ell$-maximal then $\mathcal{O}'=\mathcal{O}$.

Conversely suppose $\mathcal{O}$ is not $\ell$-maximal, and let $\mathcal{O}''$ be an order with $\mathcal{O}\subsetneq\mathcal{O}''\subseteq\mathcal{O}*i$ and $[\mathcal{O}'':\mathcal{O}]$ a power of $\ell$. Two preliminary facts. First, $I*\ell^{,n}\subseteq\ell\mathcal{O}$: the nilradical of the $n$-dimensional $\mathbb{F}*\ell$-algebra $\bar{\mathcal{O}}$ satisfies $\mathrm{nil}(\bar{\mathcal{O}})^n=0$. Second, choose $k\ge1$ minimal with $\ell^k\mathcal{O}''\subseteq\mathcal{O}$ and pick $x\in\ell^{k-1}\mathcal{O}''\setminus\mathcal{O}$; then $x\in\mathcal{O}''$ and $\ell x\in\mathcal{O}$. Let $j$ be minimal with $xI*\ell^{,j}\subseteq\mathcal{O}$; $j\le n$ because $xI\_\ell^{,n}\subseteq x\ell\mathcal{O}\subseteq\mathcal{O}$, and $j\ge1$ because $x\notin\mathcal{O}$. So there is $y\in I\_\ell^{,j-1}$ with $x':=xy\notin\mathcal{O}$, while $x'I\_\ell\subseteq xI\_\ell^{,j}\subseteq\mathcal{O}$. Now $x'\in\mathcal{O}''$, and for $z\in I\_\ell$, $x'z\in\mathcal{O}$ with $(x'z)^{n(k+1)}=x'^{,n(k+1)}(z^n)^{k+1}\in\mathcal{O}''\cdot\ell^{k+1}\mathcal{O}\subseteq\ell\cdot\ell^k\mathcal{O}''\subseteq\ell\mathcal{O}$; thus $x'z\in I\_\ell$. Hence $x'\in\mathcal{O}'\setminus\mathcal{O}$. $\square$

**Algorithm 2.3 (Round 2 at $\ell$).** Input $f\_i$, $\ell$. Set $\mathcal{O}\leftarrow\mathbb{Z}[\alpha]$. Repeat: compute $I\_\ell(\mathcal{O})$ and $\mathcal{O}'=[I\_\ell\:I\_\ell]$; if $\mathcal{O}'=\mathcal{O}$ stop, else $\mathcal{O}\leftarrow\mathcal{O}'$. Output $\mathcal{O}$, which is $\ell$-maximal by Theorem 2.2, together with $I\_\ell(\mathcal{O})$ from the last iteration. Then by Proposition 2.1:

$$\ell\text{ unramified in }K\_i\iff I\_\ell(\mathcal{O})=\ell\mathcal{O}\iff v-2v\_\ell([\mathcal{O}:\mathbb{Z}[\alpha]])=0 .$$

Both tests are read off the final iteration: the first from the computed radical, the second from the index, which is the $\ell$-adic valuation of the determinant of the basis of $\mathcal{O}$ in terms of the power basis.

**Number of iterations.** Each non-terminal iteration strictly enlarges $\mathcal{O}$ inside $\mathcal{O}*i$ by an $\ell$-power index, so the number of iterations is at most $1+v*\ell([\mathcal{O}*i:\mathbb{Z}[\alpha]])=1+\tfrac12\big(v-v*\ell(\mathrm{disc},\mathcal{O}\_i)\big)\le1+\lfloor v/2\rfloor$.

**Representation and precision.** Every intermediate order satisfies $\mathbb{Z}[\alpha]\subseteq\mathcal{O}\subseteq\ell^{-t}\mathbb{Z}[\alpha]$ with $t=v\_\ell([\mathcal{O}:\mathbb{Z}[\alpha]])\le\lfloor v/2\rfloor$. Represent $\mathcal{O}$ by an $\ell$-adic basis: $\mathcal{O}\otimes\mathbb{Z}*\ell$ is determined by the image of $\ell^t\mathcal{O}$ in $(\mathbb{Z}/\ell^{t}\mathbb{Z})^n$ (coordinates in the power basis), so all linear algebra is over $\mathbb{Z}/\ell^{m}$ with $m\le\lfloor v/2\rfloor+1$, on words of $O(v\log\ell)$ bits; the $\ell$-maximal order over $\mathbb{Z}$ (a global order) is recovered at the end as $\mathcal{O}*{i,\ell}=\mathbb{Z}[\alpha]+\ell^{-t}\Lambda$ where $\Lambda\subseteq\mathbb{Z}[\alpha]$ is any lift of the computed lattice, but it is not needed for the ramification decision.

**Cost of one iteration** (operations in $\mathbb{Z}/\ell^m$ or $\mathbb{F}\_\ell$; $\omega$ the matrix-multiplication exponent, $\omega\le3$):

1. *Multiplication table of $\mathcal{O}$.* Products of basis elements via the power basis reduced mod $f\_i$: $n^2$ products at $O(n^2)$ each, $O(n^4)$; or $O(n^3)$ using that $\mathcal{O}=\ell^{-t}\Lambda$ with $\Lambda$ a lattice in $\mathbb{Z}[\alpha]$ so that only $n$ multiplications-by-$\omega\_j$ matrices are needed, each $O(n^2)$ to form and $O(n^\omega)$ to change basis: $O(n^{\omega+1})$.
2. *Radical.* If $\ell>n$: $I\_\ell/\ell\mathcal{O}$ is the kernel of the trace form of $\bar{\mathcal{O}}$ (Lemma 1.2's argument shows the radical of the trace form contains the nilradical, and equals it when $\ell$ exceeds every local length, which is at most $n$): $O(n^3)$ to form the Gram matrix from the multiplication table, $O(n^\omega)$ for the kernel. If $\ell\le n$: the nilradical of $\bar{\mathcal{O}}$ is the kernel of $F^{r}$ where $F\:x\mapsto x^\ell$ is $\mathbb{F}*\ell$-linear and $r=\lceil\log*\ell n\rceil$; the matrix of $F$ costs $n$ exponentiations of $O(\log\ell)$ multiplications at $O(n^2)$ each, $O(n^3\log\ell)$, then $O(n^\omega\log n)$ for $F^{r}$ and its kernel.
3. *Multiplier ring.* Since $\ell\mathcal{O}'\subseteq\mathcal{O}$, $\mathcal{O}'=\ell^{-1}U$ with $U={x\in\mathcal{O}\:xI\_\ell\subseteq\ell I\_\ell}$, and $U/\ell\mathcal{O}$ is the kernel of the $\mathbb{F}*\ell$-linear map $\bar{\mathcal{O}}\to\mathrm{End}*{\mathbb{F}*\ell}(I*\ell/\ell I\_\ell)$, $x\mapsto(\text{multiplication by }x)$: an $n^2\times n$ system, $O(n^{\omega+1})$. Lifting to a basis of $\mathcal{O}'$ and putting it in Hermite form over $\mathbb{Z}/\ell^m$: $O(n^3)$.
4. *Termination test* $\mathcal{O}'=\mathcal{O}$: compare indices, $O(1)$ once the Hermite form is known.

So one iteration costs $O(n^{\omega+1}+n^3\log\ell)$ operations on $O(v\log\ell)$-bit words, and

**Proposition 2.4 (cost).** Deciding whether $\ell$ ramifies in $K\_i$, and computing $v\_\ell(\mathrm{disc},\mathcal{O}\_i)$ and the $\ell$-maximal order, costs

$$O\Big(\big(1+\lfloor v\_\ell(\mathrm{disc}f\_i)/2\rfloor\big)\cdot\big(n\_i^{\omega+1}+n\_i^{3}\log\ell\big)\Big)$$

operations on integers of $O(v\_\ell(\mathrm{disc}f\_i)\log\ell)$ bits. With $\omega=3$ this is $O\big(v\_\ell(\mathrm{disc}f\_i),n\_i^{4}+v\_\ell(\mathrm{disc}f\_i),n\_i^3\log\ell\big)$.

The factor $v\_\ell(\mathrm{disc}f\_i)/2$ is the Round 2 iteration count; Round 4 (Ford to Pauli to Roblot, via Montes' Newton-polygon and residual-polynomial refinement) computes the same $\ell$-maximal order with the number of refinement steps bounded by $O(\log v)$-type quantities in the worst case and independent of $v$ in the typical case, at the same per-step polynomial cost, and is the algorithm to use when $v$ is large; the correctness of the ramification decision (Proposition 2.1) is the same for either.

### 2.3 The full procedure

1. For each $i$, compute $\mathrm{disc}(f\_i)=\pm\mathrm{Res}(f\_i,f\_i')$ ($O(n\_i^2)$ big-integer operations, or a subresultant computation) and factor it completely.
2. $\mathcal{C}=\bigcup\_i{\ell:\ell\mid\mathrm{disc}(f\_i)}$. For $\ell\in\mathcal{C}$ and each $i$ with $\ell\mid\mathrm{disc}(f\_i)$:

- if $v\_\ell(\mathrm{disc}f\_i)$ is odd, $\ell$ ramifies in $K\_i$ (Corollary 1.5);
- else run Dedekind's criterion (2.1); if $\mathbb{Z}[\alpha\_i]$ is $\ell$-maximal, decide by the exponents $e\_j$;
- else run Algorithm 2.3 and decide by $I\_\ell(\mathcal{O})=\ell\mathcal{O}$.

3. $\ell$ is ramified in $N$ iff it is ramified in at least one $K\_i$ (Theorem 1.1). As soon as one $i$ shows ramification, the remaining $i$ need not be examined for that $\ell$.

Total cost, beyond the factorizations in step 1, is $\sum\_i\sum\_{\ell\mid\mathrm{disc}f\_i}$ of the bound in Proposition 2.4, and the sum of $v\_\ell(\mathrm{disc}f\_i)\log\ell$ over $\ell$ is $\log|\mathrm{disc}f\_i|=O(n\_i\log n\_i+n\_i\log|f\_i|\_\infty)$, so the whole decision is polynomial in $\sum\_in\_i$ and the bit size of $f$, given the factorizations.

**Output kept for later steps.** For each $i$ and each ramified $\ell$: the $\ell$-maximal order $\mathcal{O}*{i,\ell}$ (as an $\ell$-adic lattice), $v*\ell(\mathrm{disc},\mathcal{O}*i)$, and the decomposition $\mathcal{O}*{i,\ell}/\ell\mathcal{O}\_{i,\ell}\cong\prod\_j\mathcal{O}*i/\mathfrak{p}j^{e\_j}$, which gives the number of primes of $K\_i$ above $\ell$, the residue degrees $f\_j$ (dimensions of the simple factors of the semisimple quotient $\bar{\mathcal{O}}/\mathrm{nil}$) and the ramification indices $e\_j$ (local lengths); these are the inputs, together with the $\ell$-adic factorization of $f\_i$, for locating $D\ell$ and $I*\ell$ inside $G$.