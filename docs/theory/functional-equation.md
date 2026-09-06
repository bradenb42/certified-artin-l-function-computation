# The smoothed functional-equation test

**Setting.** For $\chi\in\mathrm{Irr}(G)$: $d=\chi(1)$, $a=a\_\chi$, $b=b\_\chi$ ([the archimedean place](https://claude.ai/chat/archimedean.md)), $\mathfrak{f}=\mathfrak{f}(\chi)$ ([the global conductor](https://claude.ai/chat/global-conductor.md)), $\gamma\_\chi(s)=\Gamma\_\mathbb{R}(s)^a\Gamma\_\mathbb{R}(s+1)^b$, $L(s,\chi)=\sum\_{m\ge1}a\_m(\chi)m^{-s}=\prod\_\ell P\_\ell(\chi;\ell^{-s})^{-1}$ ([Euler factors](https://claude.ai/chat/euler-factors.md) to [the Euler identities](https://claude.ai/chat/euler-identities.md), [rational classes](https://claude.ai/chat/rational-classes.md) to [the direct route](https://claude.ai/chat/direct-route.md)), $W=W(\chi)$ ([the global root number](https://claude.ai/chat/global-root-numbers.md)), $\Lambda(s,\chi)=\mathfrak{f}^{s/2}\gamma\_\chi(s)L(s,\chi)$. The coefficients $a\_m(\chi)\in\mathcal{O}\_{\mathbb{Q}(\chi)}$ are exact; they are embedded into $\mathbb{C}$ through the fixed embedding $\mathbb{Q}(\zeta\_e)\hookrightarrow\mathbb{C}$ used for the table. Throughout, $\chi\ne1$ unless stated; $\chi=1$ is treated in §6.

**Hypothesis (H$\_\chi$).** $\Lambda(s,\chi)$ is entire of finite order. This is Artin's conjecture for $\chi$; it is a theorem when $\chi$ is monomial (Hecke), when $\chi$ is a character of a nilpotent or supersolvable group (M-groups), and for every $\chi$ of degree $2$ with solvable $G$ (Langlands to Tunnell), hence for all of family 3. Without (H$\_\chi$), $\Lambda$ is meromorphic (Brauer) with functional equation $\Lambda(s,\chi)=W\Lambda(1-s,\bar\chi)$, and its possible poles lie in $0<\mathrm{Re},s<1$.

## 1. The kernel $g\_\chi$

**Definition.** $g\_\chi(x)=\dfrac{1}{2\pi i}\displaystyle\int\_{(c)}\gamma\_\chi(s),x^{-s},ds$ for any $c>0$, $x>0$; equivalently $\gamma\_\chi(s)=\int\_0^\infty g\_\chi(x)x^{s},\dfrac{dx}{x}$ for $\mathrm{Re},s>0$.

**Lemma 1.1.** (a) $g\_\chi$ is the Mellin convolution of $a$ copies of $h\_0(x)=2e^{-\pi x^2}$ and $b$ copies of $h\_1(x)=2xe^{-\pi x^2}$; hence $g\_\chi>0$, and $g\_\chi$ is real-analytic on $(0,\infty)$. (b) Closed forms: $d=1$: $g=h\_0$ ($b=0$) or $h\_1$ ($b=1$); $d=2$: $g=4K\_0(2\pi x)$ for $(a,b)=(2,0)$, $g=2e^{-2\pi x}$ for $(1,1)$, $g=4xK\_0(2\pi x)$ for $(0,2)$. (c) For every $c>0$, $g\_\chi(x)\le M\_\chi(c),x^{-c}$ with $M\_\chi(c)=\frac{1}{2\pi}\int\_\mathbb{R}|\gamma\_\chi(c+iy)|,dy$; explicitly, using $|\Gamma\_\mathbb{R}(1+iy)|=\cosh(\pi y/2)^{-1/2}$ and $|\Gamma\_\mathbb{R}(2+iy)|=\pi^{-1}\big(\tfrac{\pi y/2}{\sinh(\pi y/2)}\big)^{1/2}$,

$$M\_\chi(1)\le\frac{2^{(a+b)/2}}{\pi^{1+b}}\Big[\frac{4}{\pi d}+\pi^{b/2},\Gamma!\Big(\frac b2+1\Big)\Big(\frac{4}{\pi d}\Big)^{b/2+1}\Big]=\:M^\*\_\chi .$$

(d) Exponential decay: for all $x>0$,
$$g\_\chi(x)\le K\_\chi,x^{-1}\exp!\Big(-\frac{\pi d}{2},x^{2/d}\Big),\qquad K\_\chi:=2^{(b+d)/2}M^\*\_\chi .$$

*Proof.* (a) $\Gamma\_\mathbb{R}(s)=\int\_0^\infty 2e^{-\pi x^2}x^{s}d^\times x$ and $\Gamma\_\mathbb{R}(s+1)=\int 2xe^{-\pi x^2}x^{s}d^\times x$; products of Mellin transforms are transforms of Mellin convolutions; positivity is inherited. (b) Standard transforms ($\Gamma(s)\leftrightarrow e^{-x}$ with $\Gamma\_\mathbb{C}(s)=2(2\pi)^{-s}\Gamma(s)$; $\Gamma\_\mathbb{R}(s)^2\leftrightarrow4K\_0(2\pi x)$; a shift $s\mapsto s+1$ multiplies the kernel by $x$). (c) Bound the inverse Mellin integral on $\mathrm{Re},s=c$; for $c=1$ use $\cosh u\ge e^{u}/2$ and $u/\sinh u=\frac{2ue^{-u}}{1-e^{-2u}}\le(2u+1)e^{-u}$ ($u\ge0$, from $\frac{u}{e^{2u}-1}\le\frac12$) to get $|\Gamma\_\mathbb{R}(1+iy)|\le\sqrt2e^{-\pi|y|/4}$ and $|\Gamma\_\mathbb{R}(2+iy)|\le\pi^{-1}\sqrt{\pi|y|+1},e^{-\pi|y|/4}$, then $(\pi y+1)^{b/2}\le2^{b/2}\big(1+(\pi y)^{b/2}\big)$ and integrate over $y\ge0$ (doubling for $y<0$). (d) In the convolution integral over $u\_1\cdots u\_d=x$, write $e^{-\pi\sum u\_i^2}=e^{-\pi\sum u\_i^2/2}e^{-\pi\sum u\_i^2/2}$ and use $\sum u\_i^2\ge d,x^{2/d}$ (AM to GM); the remaining integrand is the convolution of the $\tilde h\_i(u)=2u^{k\_i}e^{-\pi u^2/2}$, whose Mellin transform is $2^{(ds+b)/2}\gamma\_\chi(s)$, so it equals $2^{b/2}g\_\chi(2^{-d/2}x)\le2^{b/2}M^\*\_\chi2^{d/2}x^{-1}$ by (c). $\square$

The kernel is evaluated with certified error (interval arithmetic on the closed forms for $d\le2$, on the convolution integral or the Mellin integral otherwise); that evaluation error is added to the bounds below and is not written explicitly.

## 2. The theta function and the identity

**Definition.** $\Theta\_\chi(t)=\sum\_{m\ge1}a\_m(\chi),g\_\chi!\big(mt/\sqrt{\mathfrak f}\big)$, $t>0$; and for $X\ge1$ the finite sums $S\_\chi(t,X)=\sum\_{m\le X}a\_m(\chi)g\_\chi(mt/\sqrt{\mathfrak f})$, tails $T\_\chi(t,X)=\Theta\_\chi(t)-S\_\chi(t,X)$.

**Lemma 2.1 (coefficient bound).** $|a\_m(\chi)|\le d\_{,d}(m)$, the number of ordered factorizations of $m$ into $d=\chi(1)$ positive factors; in particular $|a\_m|\le\tau(m)$ for $d\le2$, and $d\_d(m)\le\tau(m)^{d-1}\le(2\sqrt m)^{d-1}$. The bound $|a\_m|\le\chi(1)\tau(m)$ holds for $d\le2$ but not in general: for $d\ge3$ and $m=p^k$, $a\_{p^k}$ can equal $\binom{k+d-1}{d-1}$, which exceeds $d(k+1)$ for large $k$.

*Proof.* At every prime $P\_\ell(\chi;T)^{-1}=\prod\_{i\le d'}(1-\lambda\_iT)^{-1}$ with $d'\le d$ and $|\lambda\_i|\le1$ (roots of unity), so $|a\_{\ell^k}|\le\binom{k+d'-1}{d'-1}\le\binom{k+d-1}{d-1}=d\_d(\ell^k)$; multiplicativity. $d\_d=\tau\*d\_{d-1}$ gives $d\_d(p^k)=\binom{k+d-1}{d-1}=\prod\_{i=1}^{d-1}\frac{k+i}{i}\le(k+1)^{d-1}$; $\tau(m)\le2\sqrt m$. Equality $a\_{p^k}=\binom{k+d-1}{d-1}$ occurs when $\rho(\mathrm{Fr}\_p)=1$. $\square$

**Theorem 2.2 (the identity).** Assume the data $a\_m(\chi)$, $\mathfrak f$, $\gamma\_\chi$, $W$ are correct and (H$\_\chi$) holds. Then for every $t>0$:

$$\Theta\_\chi(t)=W,\Theta\_{\bar\chi}(1/t),\qquad\text{hence}\qquad S\_\chi(t,X)=W,S\_{\bar\chi}(1/t,X)+\mathcal{R}*\chi(t,X),\quad\mathcal{R}*\chi(t,X)=T\_\chi(t,X)-W,T\_{\bar\chi}(1/t,X).$$

If instead $\Lambda(s,\chi)$ is meromorphic with finitely many poles $s\_0$, all in $0<\mathrm{Re},s<1$, and of polynomial growth in vertical strips, the same holds with the extra term $\mathcal{P}*\chi(t)=\sum*{s\_0}\mathrm{Res}\_{s=s\_0}\big(\Lambda(s,\chi)t^{-s}\big)$ on the right of the first identity.

*Proof.* For $\mathrm{Re},s>1$, $\int\_0^\infty\Theta\_\chi(t)t^{s}d^\times t=\sum\_ma\_m(\sqrt{\mathfrak f}/m)^s\gamma\_\chi(s)=\Lambda(s,\chi)$, the interchange being justified by absolute convergence ($|a\_m|\le d\_d(m)$, Lemma 1.1(d)). Mellin inversion: $\Theta\_\chi(t)=\frac{1}{2\pi i}\int\_{(2)}\Lambda(s,\chi)t^{-s}ds$. Under (H$*\chi$), $\Lambda$ is entire and, by Phragmén to Lindelöf between $\mathrm{Re},s=2$ and $\mathrm{Re},s=-1$ where it is bounded (on $\mathrm{Re},s=-1$ by the functional equation), of polynomial growth in the strip, so the contour moves to $\mathrm{Re},s=-1$ without residues; there $\Lambda(s,\chi)=W\Lambda(1-s,\bar\chi)$ and the substitution $s'=1-s$ gives $\frac{W}{2\pi i}\int*{(2)}\Lambda(s',\bar\chi)t^{-(1-s')}ds'=W\Theta\_{\bar\chi}(1/t)$. With poles, the residues appear; a pole of order $r$ contributes $t^{-s\_0}$ times a polynomial of degree $r-1$ in $\log t$. $\square$

**Remark 2.3 (converse).** If, for the given data, $\Theta\_\chi(t)=W\Theta\_{\bar\chi}(1/t)$ holds for all $t>0$, then $\Lambda(s,\chi)=\int\_1^\infty\Theta\_\chi(t)t^sd^\times t+W\int\_1^\infty\Theta\_{\bar\chi}(t)t^{1-s}d^\times t$ is entire and satisfies the functional equation (Hecke). So, with correct data, the identity is *equivalent* to (H$*\chi$); a numerical verification at finitely many $t$ cannot establish it, but its failure is a proof that either some datum is wrong or (H$*\chi$) fails. Everything below is the contrapositive.

## 3. The tail bound

**Proposition 3.1.** For $t>0$, $X\ge1$, put $\alpha=\alpha(t)=\frac{\pi d}{2}\big(t/\sqrt{\mathfrak f}\big)^{2/d}$ and $\kappa=\frac{d(d-1)}{4}$. If $\alpha X^{2/d}\ge\max{\tfrac{d(d-3)}{4},,2\kappa,,1}$, then

$$|T\_\chi(t,X)|\ \le\ \mathcal{E}*\chi(t,X):=2^{d-1}K*\chi,\frac{\sqrt{\mathfrak f}}{t}\cdot\frac d2,\alpha^{-\kappa},\Gamma\big(\kappa,\alpha X^{2/d}\big)\ \le\ 2^{d}K\_\chi,\frac{\sqrt{\mathfrak f}}{t}\cdot\frac d2,X^{(d-1)/2-2/d},\alpha^{-1},e^{-\alpha X^{2/d}},$$

and consequently $|\mathcal{R}*\chi(t,X)|\le\mathcal{E}*\chi(t,X)+\mathcal{E}\_{\bar\chi}(1/t,X)=:\mathcal{E}(X,\chi,t)$ ($\bar\chi$ has the same $d,a,b,\mathfrak f$).

*Proof.* By Lemmas 1.1(d) and 2.1, $|a\_mg\_\chi(mt/\sqrt{\mathfrak f})|\le(2\sqrt m)^{d-1}K\_\chi\frac{\sqrt{\mathfrak f}}{mt}e^{-\alpha m^{2/d}}=2^{d-1}K\_\chi\frac{\sqrt{\mathfrak f}}{t},m^{(d-3)/2}e^{-\alpha m^{2/d}}$. The function $u\mapsto u^{(d-3)/2}e^{-\alpha u^{2/d}}$ is decreasing for $\alpha u^{2/d}\ge d(d-3)/4$, so the sum over $m>X$ is at most $\int\_X^\infty u^{(d-3)/2}e^{-\alpha u^{2/d}}du=\frac d2\alpha^{-\kappa}\Gamma(\kappa,\alpha X^{2/d})$ (substitute $w=\alpha u^{2/d}$). The elementary bound $\Gamma(\kappa,y)\le2y^{\kappa-1}e^{-y}$ for $y\ge\max(2\kappa,1)$ gives the last expression. $\square$

For $d\le2$ the monotonicity condition is vacuous and $\kappa\in{0,\tfrac12}$ ($\Gamma(0,y)=E\_1(y)\le e^{-y}/y$, $\Gamma(\tfrac12,y)=\sqrt\pi,\mathrm{erfc}(\sqrt y)\le e^{-y}/\sqrt y$). To make $\mathcal{E}(X,\chi,t)\le\epsilon$ it suffices to take

$$X\ \ge\ \Big(\frac{\sqrt{\mathfrak f}}{\min(t,1/t)}\Big)\Big(\frac{2}{\pi d}\log\frac{C\_\chi}{\epsilon\min(t,1/t)}\Big)^{d/2}\quad(\text{with }C\_\chi\text{ explicit from the constants above}),$$

i.e. $X$ of size $\sqrt{\mathfrak f}\cdot(\log\frac1\epsilon)^{d/2}$, which is the standard cost of the test: the coefficients $a\_m(\chi)$ are needed for $m\le X$, so the Frobenius computations of [rational classes](https://claude.ai/chat/rational-classes.md) to [the direct route](https://claude.ai/chat/direct-route.md) must reach $\max(t,1/t)X\cdot$, in practice $t$ is kept in $[\tfrac12,2]$ and the bound $X$ of the pipeline is chosen accordingly.

## 4. The contrapositive, with thresholds

Fix $X$ and $t$; the computable discrepancy is $\Delta\_\chi(t,X;W'):=S\_\chi(t,X)-W'S\_{\bar\chi}(1/t,X)$ for the value $W'$ under test. **Decision rule:** the data are rejected if $|\Delta\_\chi(t,X;W')|>\mathcal{E}(X,\chi,t)$ (plus the certified evaluation error). By Theorem 2.2 and Proposition 3.1, correct data with $W'=W$ and (H$\_\chi$) never trigger the rule.

**Theorem 4.1 (wrong root number).** Suppose all data are correct and (H$*\chi$) holds, but $W'\ne W$ is used. Then $\Delta*\chi(t,X;W')=(W-W'),\Theta\_{\bar\chi}(1/t)-T\_\chi(t,X)+W'T\_{\bar\chi}(1/t,X)$, so the rule triggers at every $t$ with

$$|W-W'|\cdot|S\_{\bar\chi}(1/t,X)|\ >\ \mathcal{E}(X,\chi,t)+|W-W'|,\mathcal{E}\_{\bar\chi}(1/t,X),$$

a condition decided from the computed $S\_{\bar\chi}(1/t,X)$. The set of $t>0$ at which $\Theta\_{\bar\chi}(1/t)=0$ is discrete (a nonzero real-analytic function), so for $X$ large enough the rule triggers at every $t$ outside a discrete set; for self-dual $\chi$ with $W=-1$, $t=1$ belongs to that set ($\Theta\_\chi(1)=W\Theta\_\chi(1)$ forces $\Theta\_\chi(1)=0$), so for the candidates $W'=\pm1$ the test must be run at $t$ bounded away from $1$, and there $|W-W'|=2$. At most one $W'$ passes at a $t$ with $2|S\_{\bar\chi}(1/t,X)|>\mathcal{E}+2\mathcal{E}\_{\bar\chi}$, in agreement with [the character table](https://claude.ai/chat/character-table.md) Lemma 5.3.

*Proof.* Substitute the true identity into $\Delta$; the triggering condition is $|\Delta|>\mathcal{E}$ with $|(W-W')\Theta\_{\bar\chi}(1/t)|\ge|W-W'|(|S\_{\bar\chi}|-\mathcal{E}\_{\bar\chi})$ and the other terms bounded by $\mathcal{E}$. $\square$

**Theorem 4.2 (conductor wrong at one prime).** Suppose $\mathfrak f'=\mathfrak f,\ell^{k}$, $k\ne0$, is used with everything else correct. Then the computed sums are $S'*\chi(t,X)=S*\chi(t\ell^{-k/2},X)$, and

$$\Delta'*\chi(t,X;W)=W\Big(\Theta*{\bar\chi}\big(\ell^{k/2}/t\big)-\Theta\_{\bar\chi}\big(\ell^{-k/2}/t\big)\Big)+O(\mathcal{E}'),$$

where $\mathcal{E}'$ is the tail bound at the rescaled arguments. The bracket, as a function of $u=1/t$, is not identically zero: if $\Theta\_{\bar\chi}(u\lambda)=\Theta\_{\bar\chi}(u/\lambda)$ for all $u$ with $\lambda=\ell^{k/2}\ne1$, the Mellin transform gives $\lambda^{-s}\Lambda(s,\bar\chi)=\lambda^{s}\Lambda(s,\bar\chi)$, so $\Lambda\equiv0$. Hence for each $t$ outside a discrete set there is a finite $X\_0(\ell,k,t)$, the least $X$ with $\mathcal{E}'(X)<\tfrac12|\Theta\_{\bar\chi}(\ell^{k/2}/t)-\Theta\_{\bar\chi}(\ell^{-k/2}/t)|$, decided from the computed sums, beyond which the rule triggers; a priori, the discrepancy is a difference of two theta values at arguments in ratio $\ell^{k}$, so it is bounded below by the variation of $\Theta\_{\bar\chi}$ on $[\ell^{-|k|/2}/t,\ell^{|k|/2}/t]$, which for $t$ near $1$ is of the order of $\Theta\_{\bar\chi}$ itself once $\ell^{|k|/2}\ge2$.

*Proof.* $g\_\chi(mt/\sqrt{\mathfrak f'})=g\_\chi(m,t\ell^{-k/2}/\sqrt{\mathfrak f})$; then use the true identity at $t\ell^{-k/2}$. $\square$

**Theorem 4.3 (ramified Euler factor wrong at one prime).** Suppose at a ramified $\ell$ the factor $P'*\ell\ne P*\ell$ is used (with $P'*{\bar\chi,\ell}=\overline{P'*\ell}$), everything else correct. Write $P\_\ell(T)/P'*\ell(T)=\sum*{j\ge0}c\_jT^j$ ($c\_0=1$, $c\_j$ not all zero for $j\ge1$, $|c\_j|\le2^{\deg P\_\ell},d\_{\deg P'*\ell}(\ell^j)$). Then $a'm=\sum{\ell^j\mid m}c\_ja*{m/\ell^j}$, $\Theta'*\chi(t)=\sum*{j\ge0}c\_j\Theta\_\chi(\ell^jt)$, and

$$\Delta'*\chi(t,X;W)=W\sum*{j\ge1}\Big(c\_j,\Theta\_{\bar\chi}(\ell^{-j}/t)-\bar c\_j,\Theta\_{\bar\chi}(\ell^{j}/t)\Big)+O(\mathcal{E}'),$$

whose main term is not identically zero in $t$ (its Mellin transform is $\sum\_{j\ge1}(c\_j\ell^{js}-\bar c\_j\ell^{-js})\Lambda(s,\bar\chi)$, and a nonzero Laurent polynomial in $\ell^{s}$ does not vanish identically). Hence for $t$ outside a discrete set there is a finite $X\_0(\ell,\chi,t)$ beyond which the rule triggers, decided from the computed sums as in Theorem 4.2; the leading term is $c\_{j\_0}$ times a theta value at argument $\ell^{-j\_0}/t$ with $j\_0$ the least $j\ge1$ with $c\_j\ne0$, so the discrepancy is of the size of $\Theta\_{\bar\chi}$ at a point displaced from $1/t$ by the factor $\ell^{j\_0}$, and $X\_0$ is governed by $\mathcal{E}$ at that point.

*Proof.* Multiplicativity of the coefficients and the true identity at the arguments $\ell^jt$. $\square$

The same argument applies to a wrong coefficient at an unramified prime $\ell\le X$ (Frobenius class wrong), with $c\_j$ the coefficients of the ratio of the two Euler factors; there the displacement is again by powers of $\ell$, and a wrong class at a large prime $\ell$ is detected only when $g\_\chi(\ell t/\sqrt{\mathfrak f})$ is not below the tail bound, i.e. when $\ell t/\sqrt{\mathfrak f}$ is at most of the order of $(\log\tfrac1\epsilon)^{d/2}$: errors at primes near $X$ are the last to be seen and require a smaller $t$ or a larger $X$.

## 5. What the test cannot detect

1. Anything at $m>X$: a wrong Frobenius class at a prime $\ell>X$, or at a prime $\ell\le X$ whose kernel value $g\_\chi(\ell t/\sqrt{\mathfrak f})$ is below $\mathcal{E}$ at every $t$ used.
2. A consistent Galois relabeling: if every datum of $\chi$ is replaced by that of $\chi^\tau$ (Frobenius classes, Euler factors, $W$) the identity holds. In particular $\chi\leftrightarrow\bar\chi$ with $W\leftrightarrow\overline W$. The test verifies an $L$-function; it does not verify that the $L$-function is the one attached to the row $\chi$ of the table rather than to a conjugate row. That identification is the matching of [the precision policy](https://claude.ai/chat/precision-policy.md) §4/[matching](https://claude.ai/chat/matching.md) and is certified there, and, for family 3, by the comparison with a newform given together with an embedding of its coefficient field.
3. Errors below the resolution: any change of the data whose effect on $\Delta$ is smaller than $\mathcal{E}(X,\chi,t)$ at all $t$ used; the resolution is improved only by increasing $X$ (and the precision of $g$).
4. The distinction between a wrong datum and a failure of (H$*\chi$): a pole contributes $\mathcal{P}*\chi(t)$ (Theorem 2.2), a wrong datum contributes the expressions of §4; both are nonzero, and only the shape in $t$ distinguishes them (a pole gives $t^{-s\_0}$-type terms with $0<\mathrm{Re},s\_0<1$, a conductor or Euler-factor error gives combinations of $\Theta$ at rescaled arguments). When (H$\_\chi$) is a theorem (family 3, monomial $\chi$), this ambiguity is absent.
5. Errors in $\gamma\_\chi$ (i.e. in the class of $c$, [the archimedean place](https://claude.ai/chat/archimedean.md)) are detected in principle (the identity with a wrong $(a,b)$ cannot hold for all $t$: by Remark 2.3 it would give a second functional equation with a different Gamma factor, and $\gamma'(s)/\gamma(s)$ is not invariant under $s\mapsto1-s$ up to a constant), but no a priori threshold is given here; [the archimedean place](https://claude.ai/chat/archimedean.md)'s checks ($\chi(c)\equiv\chi(1)\bmod2$, the class of $c$ from exact real-root counts) are the primary safeguard.
6. Compensating simultaneous errors that happen to produce the data of some other genuine $L$-function satisfying a functional equation of the same shape; by strong multiplicity one for Artin $L$-functions (Chebotarev), two Artin $L$-functions with the same coefficients at all but finitely many primes are equal, so this can only be the relabeling of item 2 combined with item 1.

The test is therefore a rejection test: passing it at several $t$ with small $\mathcal{E}$ is strong evidence, and it is treated as the acceptance criterion for the *numerical* consistency of the assembled data, but the algebraic certificates of [matching](https://claude.ai/chat/matching.md) to [the global root number](https://claude.ai/chat/global-root-numbers.md) are what prove the individual components.

## 6. The trivial character

For $\chi=1$: $L(s,1)=\zeta(s)$, $\mathfrak f=1$, $a=1$, $b=0$, $g\_1(x)=2e^{-\pi x^2}$, $W=1$, and $\Lambda(s)=\Gamma\_\mathbb{R}(s)\zeta(s)$ has simple poles at $s=1$ (residue $1$) and $s=0$ (residue $\Gamma\_\mathbb{R}$'s pole times $\zeta(0)=-\tfrac12$, i.e. $-1$). Theorem 2.2's polar version gives

$$\Theta\_1(t)=t^{-1}\Theta\_1(1/t)+t^{-1}-1,\qquad\Theta\_1(t)=\sum\_{m\ge1}2e^{-\pi m^2t^2}=\vartheta(t)-1,$$

which is Jacobi's $\vartheta(t)=t^{-1}\vartheta(1/t)$ for $\vartheta(t)=\sum\_{m\in\mathbb{Z}}e^{-\pi m^2t^2}$; the finite version is $S\_1(t,X)=t^{-1}S\_1(1/t,X)+t^{-1}-1+\mathcal{R}\_1(t,X)$ with $|\mathcal{R}*1|\le\mathcal{E}(X,1,t)$, $\mathcal{E}$ from Proposition 3.1 with $d=1$, $\mathfrak f=1$, $K\_1=2\sqrt2M^\*1$. The data of $\chi=1$ contain nothing computed by the pipeline; the identity is included because (i) it exercises the kernel evaluation, tail bounds and decision rule on a case with known answer, (ii) the polar terms $t^{-1}-1$ are exactly the terms that must be added in the product identity for $\zeta\_N=\prod\chi L(s,\chi)^{\chi(1)}$ (*[*Euler factors*](https://claude.ai/chat/euler-factors.md) *check 4 at the level of theta functions: $\prod$ over $\chi$ of the identities, with the single pole of $\zeta\_N$ at $s=1$ coming from $\chi=1$ alone when all other $L(s,\chi)$ are entire), and (iii) a nontrivial $\chi$ for which $L(s,\chi)$ had a pole at $s=1$ would show, in the test of §2, precisely a term $\mathrm{Res}*{s=1}\Lambda(s,\chi)\cdot t^{-1}$, which is excluded for irreducible $\chi\ne1$ by the theorem that $L(s,\chi)$ is holomorphic and nonvanishing on $\mathrm{Re},s\ge1$; so a $t^{-1}$-shaped discrepancy for $\chi\ne1$ is always a data error, most likely a wrong multiplicity of $\mathbf{1}$ in the assembled character.