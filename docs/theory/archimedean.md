# Complex conjugation, the Gamma factor and the archimedean root number

**Setting.** $N\subset\bar{\mathbb{Q}}$ the splitting field with $G\le S\_n$ in the $p$-adic numbering; for an embedding $\iota\_\infty\:N\to\mathbb{C}$, complex conjugation restricts to an automorphism $c\_{\iota\_\infty}\in G$; embeddings differ by elements of $G$, so the class $c^G$ is well defined. $D\_\infty=\langle c\rangle$ is the decomposition group at the infinite place; $c^2=1$.

## 1. The class of $c$

**Proposition 1.1.** Let $r$ be the number of real roots of $f$. Then (a) $n-r$ is even, and $c$ has cycle type $1^r2^{(n-r)/2}$ in the degree-$n$ action; (b) $c=1$ iff $r=n$ iff $N$ is totally real, and otherwise $N$ is totally complex; (c) the rational class of $c$ is its conjugacy class; (d) for every $H\le G$ and every squarefree resolvent $R\_{G,H,F}\in\mathbb{Z}[x]$, $\mathrm{fix}*{G/H}(c)$ equals the number of real roots of $R*{G,H,F}$; (e) the coset-action tests of [separating subgroups](separating-subgroups.md) for the block of cycle type $1^r2^{(n-r)/2}$, with $\mathrm{fix}\_{G/H}(c)$ supplied by (d), determine the class of $c$ exactly.

*Proof.* (a) Under $\iota\_\infty$ the roots are $\iota\_\infty(\alpha\_i)$, permuted by conjugation: real roots are fixed, non-real roots are exchanged in conjugate pairs. (b) $c=1$ iff $\iota\_\infty(N)\subseteq\mathbb{R}$ iff all roots real; $N/\mathbb{Q}$ being Galois, all embeddings have the same image, so $N$ is totally real or totally complex. (c) The powers of $c$ prime to its order are $c$ itself ([rational classes](rational-classes.md) §2). (d) $c$ acts on the roots $F(\tau\alpha)$ of $R$ by $F(\tau\alpha)\mapsto F(c\tau\alpha)=\overline{F(\tau\alpha)}$ under $\iota\_\infty$ (the invariant has rational coefficients), so a root is fixed iff it is real, and the roots are distinct. (e) [rational classes](rational-classes.md) Theorem 4.1: the coset actions determine the rational class, which by (c) is the class. $\square$

**Certified counts.** The number of real roots of a squarefree $g\in\mathbb{Z}[x]$ is computed exactly by Sturm's theorem (sign variations of the Sturm sequence at $\pm\infty$, rational arithmetic), or by certified real-root isolation (Descartes/Vincent bisection with exact interval endpoints); no floating-point estimate enters. Consistency: $\mathrm{fix}\_{G/H}(c)\equiv[G\:H]\pmod2$ since the non-real roots pair off.

**Procedure.** Count $r$ (Sturm on $f$). If $r=n$, $c=1$. Otherwise locate the block $\mathcal{B}*\lambda$, $\lambda=1^r2^{(n-r)/2}$, from the table (classes of involutions of that cycle type, by the power maps of* [*the character table*](character-table.md)*); if $\mathcal{B}*\lambda$ is a single class, done; else apply the separating family $\mathcal{H}*\lambda$ of* [*separating subgroups*](separating-subgroups.md) *with $\mathrm{fix}*{G/H}(c)$ = real-root count of $R\_{G,H,F}$ (Sturm on the exact resolvents already computed at $p$), and read the class by [rational classes](rational-classes.md) Theorem 4.1(1). No archimedean numerics of the roots themselves and no matching at the infinite place are needed. As a check, $c$ must satisfy $\chi(c)\in\mathbb{Z}$ with $\chi(c)\equiv\chi(1)\pmod 2$ for every $\chi$ (§2).

## 2. Eigenvalues of $\rho\_\chi(c)$

**Proposition 2.1.** For $\chi\in\mathrm{Irr}(G)$, $\rho\_\chi(c)$ is an involution with eigenvalues $\pm1$ of multiplicities

$$a\_\chi=\frac{\chi(1)+\chi(c)}{2},\qquad b\_\chi=\frac{\chi(1)-\chi(c)}{2},$$

nonnegative integers with $a\_\chi+b\_\chi=\chi(1)$, $a\_\chi-b\_\chi=\chi(c)\in\mathbb{Z}$; $\det\rho\_\chi(c)=(-1)^{b\_\chi}$; and $a\_\chi=\dim V\_\chi^{D\_\infty}=\langle\chi|*{\langle c\rangle},1\rangle$, $b*\chi=\langle\chi|\_{\langle c\rangle},\mathrm{sgn}\rangle$. These are the multiplicities $m\_0,m\_1$ of [the character table](character-table.md) Prop. 2.1(4) for $o=2$.

*Proof.* $\rho\_\chi(c)^2=1$ so $V\_\chi=V^+\oplus V^-$; the trace is $\dim V^+-\dim V^-$ and the dimension their sum; the determinant is $(-1)^{\dim V^-}$. $\square$

## 3. The Gamma factor

**Definition.** $\Gamma\_\mathbb{R}(s)=\pi^{-s/2}\Gamma(s/2)$, $\Gamma\_\mathbb{C}(s)=2(2\pi)^{-s}\Gamma(s)=\Gamma\_\mathbb{R}(s)\Gamma\_\mathbb{R}(s+1)$ (Legendre's duplication formula).

**Proposition 3.1.** The archimedean factor of $L(s,\chi)$ is

$$\gamma\_\chi(s)=\Gamma\_\mathbb{R}(s)^{a\_\chi},\Gamma\_\mathbb{R}(s+1)^{b\_\chi},$$

and $\Lambda(s,\chi)=\mathfrak{f}(\chi)^{s/2}\gamma\_\chi(s)L(s,\chi)$ is the completed $L$-function of [the character table](character-table.md) Lemma 5.3.

*Proof.* The local factor at the real place is defined for representations of $\mathrm{Gal}(\mathbb{C}/\mathbb{R})=\langle c\rangle$ by additivity from the two characters, with $L\_\infty(s,\mathbf{1})=\Gamma\_\mathbb{R}(s)$ and $L\_\infty(s,\mathrm{sgn})=\Gamma\_\mathbb{R}(s+1)$ (Tate's local factors for $\mathbb{R}^\times$, matching the Gamma factors of even and odd Dirichlet $L$-functions); $V\_\chi|*{\langle c\rangle}=a*\chi\mathbf{1}+b\_\chi,\mathrm{sgn}$. $\square$

For $c=1$ every $\gamma\_\chi(s)=\Gamma\_\mathbb{R}(s)^{\chi(1)}$; for $\chi(c)=-\chi(1)$, $\gamma\_\chi(s)=\Gamma\_\mathbb{R}(s+1)^{\chi(1)}$; and $\gamma\_\chi=\gamma\_{\bar\chi}=\gamma\_{\chi^\tau}$ for all Galois conjugates, as $\chi(c)\in\mathbb{Z}$.

## 4. The archimedean root number

**Proposition 4.1.** With the normalization in which the functional equation reads $\Lambda(s,\chi)=W(\chi)\Lambda(1-s,\bar\chi)$, $W(\chi)=\varepsilon\_\infty(\chi)\prod\_\ell\varepsilon\_\ell(\chi)$ with $|\varepsilon\_\ell|=1$ after the standard normalization, and the finite root number of a primitive Dirichlet character $\psi$ of conductor $q$ equal to $\tau(\psi)/\sqrt q$ with $\tau(\psi)=\sum\_{a\bmod q}\psi(a)e^{2\pi ia/q}$:
$$\varepsilon\_\infty(\chi)=i^{-b\_\chi}.$$

*Proof.* Root numbers are multiplicative in direct sums, so $\varepsilon\_\infty(\chi)=\varepsilon\_\infty(\mathbf{1})^{a\_\chi}\varepsilon\_\infty(\mathrm{sgn})^{b\_\chi}$; the values are fixed by the Dirichlet case: for primitive $\psi$ with $\psi(-1)=(-1)^\kappa$, $\Lambda(s,\psi)=q^{s/2}\Gamma\_\mathbb{R}(s+\kappa)L(s,\psi)$ satisfies $\Lambda(s,\psi)=\dfrac{\tau(\psi)}{i^{\kappa}\sqrt q},\Lambda(1-s,\bar\psi)$, so $\varepsilon\_\infty(\mathbf{1})=1$ and $\varepsilon\_\infty(\mathrm{sgn})=i^{-1}$. $\square$

Hence $\varepsilon\_\infty(\chi)\in{1,-i,-1,i}$ according to $b\_\chi\bmod4$, $\varepsilon\_\infty(\bar\chi)=\varepsilon\_\infty(\chi)$, and $\varepsilon\_\infty(\chi)\varepsilon\_\infty(\bar\chi)=(-1)^{b\_\chi}=\det\rho\_\chi(c)$, the sign of the nebentypus at $-1$ in the modular case.

## 5. The odd two-dimensional characters

**Proposition 5.1.** For $\chi\in\mathrm{Irr}(G)$ of degree $2$ the following are equivalent: $\chi$ is odd ($\det\rho\_\chi(c)=-1$); $b\_\chi=1$; $\chi(c)=0$. For an even $\chi$ of degree $2$, $\chi(c)=\pm2$ and $\rho\_\chi(c)=\pm1$. For odd $\chi$:

$$\gamma\_\chi(s)=\Gamma\_\mathbb{R}(s)\Gamma\_\mathbb{R}(s+1)=\Gamma\_\mathbb{C}(s)=2(2\pi)^{-s}\Gamma(s),\qquad\varepsilon\_\infty(\chi)=-i,$$

which are the Gamma factor and archimedean sign of a holomorphic newform of weight one, and $\lambda\_\chi=\det\rho\_\chi$ is an odd Dirichlet character (the nebentypus), $\lambda\_\chi(c)=(-1)^{b\_\chi}=-1$, in agreement with [the global conductor](global-conductor.md) Prop. 2.1 evaluated at $c$.

*Proof.* $a\_\chi+b\_\chi=2$, $a\_\chi-b\_\chi=\chi(c)$; $\det\rho\_\chi(c)=(-1)^{b\_\chi}$ is $-1$ iff $b\_\chi=1$ iff $\chi(c)=0$. The Gamma factor is Proposition 3.1 with $a=b=1$ and the duplication formula; $\varepsilon\_\infty=i^{-1}$. $\square$

**Selection of family 3.** From the table and the class of $c$: the odd two-dimensional characters are the rows with $\chi(1)=2$ and $\chi(c)=0$; among these, those with solvable $G$ (equivalently, by Langlands to Tunnell in the two-dimensional case, those whose projective image is dihedral, tetrahedral or octahedral) are the ones for which the coefficient-by-coefficient comparison with a weight-one newform of level $\mathfrak{f}(\chi)$ and nebentypus $\lambda\_\chi$ (Deligne to Serre, with existence by Langlands to Tunnell) is the correctness test of the whole pipeline; icosahedral cases lie outside family 3 as specified. The data consumed by that comparison are all now in hand: $\mathfrak{f}(\chi)$ ([conductor exponents](conductor-exponents.md) to [the global conductor](global-conductor.md)), $\lambda\_\chi$ ([the global conductor](global-conductor.md)), the Euler factors at every prime ([Euler factors](euler-factors.md) to [the Euler identities](euler-identities.md) and [rational classes](rational-classes.md) to [the direct route](direct-route.md)), and $\gamma\_\chi$, $\varepsilon\_\infty(\chi)$ (this note). If $c=1$ there are no odd characters and family 3 is empty for $f$.

## 6. Cost

One Sturm sequence for $f$ and, if the involution block has several classes, one per resolvent of the separating family, each a computation in exact rational arithmetic of degree $[G\:H]$; the multiplicities $a\_\chi,b\_\chi$ are read from the table; nothing else is archimedean. No numerical root of $f$ is ever computed to decide the class of $c$, which is why the archimedean data carry no precision requirement in [the precision policy](precision-policy.md)'s policy.