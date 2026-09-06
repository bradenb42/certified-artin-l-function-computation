# The conductor $\mathfrak{f}(\chi)$ and its cross-checks against $\det\rho\_\chi$ and $\chi\bar\chi$

**Setting.** $f\_\ell(\chi)$ from [conductor exponents](conductor-exponents.md) for every class function $\chi$ of $G$ and every ramified $\ell$ (the list of [the ramified primes](ramified-primes.md)); $f\_\ell(\chi)=0$ at unramified $\ell$ since $D\_0=1$ there. For $\chi\in\mathrm{Irr}(G)$ with representation $\rho\_\chi$ on $V$, $\lambda\_\chi:=\det\rho\_\chi$ is a linear character of $G$, hence a character of $G^{ab}$; $\widehat{G^{ab}}$ and the Dirichlet characters $\chi\_\psi$ of conductor $c\_\psi$ attached to $\psi\in\widehat{G^{ab}}$ are those of [cyclotomic refinement](cyclotomic-refinement.md) §1. For a linear character $\lambda$ of order $o=o(\lambda)$ let $N\_\lambda=N^{\ker\lambda}$, cyclic of degree $o$ over $\mathbb{Q}$, with discriminant $d\_{N\_\lambda}$.

## 1. The global conductor

**Definition and Proposition 1.1.** $\mathfrak{f}(\chi):=\prod\_{\ell}\ell^{f\_\ell(\chi)}$, a finite product over the ramified primes. For characters, $\mathfrak{f}(\chi)\in\mathbb{Z}*{\ge1}$ (*[*conductor exponents*](conductor-exponents.md) *Theorem 2.3); $\mathfrak{f}$ is multiplicative in direct sums, $\mathfrak{f}(\chi^\tau)=\mathfrak{f}(\chi)$ for $\tau\in\mathrm{Gal}(\mathbb{Q}(\chi)/\mathbb{Q})$, $\mathfrak{f}(\bar\chi)=\mathfrak{f}(\chi)$, and $\mathfrak{f}(\chi)=1$ iff $\chi$ is trivial on every inertia group iff the representation factors through the maximal unramified quotient, which for $N/\mathbb{Q}$ means $\chi=\chi(1)\cdot1$ (Minkowski). $|d\_N|=\prod*{\chi\in\mathrm{Irr}(G)}\mathfrak{f}(\chi)^{\chi(1)}$ ([conductor exponents](conductor-exponents.md) Theorem 4.1 at each prime).

**Inflation.** If $\chi$ is a character of $G/M$ inflated to $G$ ($M\trianglelefteq G$), then $f\_\ell(\chi)$ computed in $G$ equals $f\_\ell(\chi)$ computed in $G/M$ with the filtration of $N^M$. *Proof.* [the filtration](ramification-filtration.md) Prop. 4.1(b): $f\_\ell(\chi)=\int\_{-1}^{\infty}\mathrm{codim},V^{D^v},dv$, and $V^{D^v}=V^{D^vM/M}$ for inflated $V$, while $D^vM/M$ is the upper filtration of the decomposition group of $N^M$ (Herbrand, [the filtration](ramification-filtration.md) Prop. 4.1(c)). $\square$ So the conductors of the linear characters below may be computed either in $G$ ([conductor exponents](conductor-exponents.md)) or in the small quotient $\mathrm{Gal}(N\_\lambda/\mathbb{Q})$, and the two must agree.

## 2. The determinant from the character table

**Proposition 2.1.** For $\sigma\in G$ of order $o$, with $m\_j(\sigma)=\frac1o\sum\_{t=0}^{o-1}\chi(\sigma^t)\zeta\_o^{-jt}$ the eigenvalue multiplicities of [the character table](character-table.md) Prop. 2.1(4),

$$\lambda\_\chi(\sigma)=\det\rho\_\chi(\sigma)=\prod\_{j\in\mathbb{Z}/o}\zeta\_o^{,j,m\_j(\sigma)}=\zeta\_o^{\sum\_jj,m\_j(\sigma)} .$$

$\lambda\_\chi$ is a class function computable from the row $\chi$ and the power maps alone, with values in $\mathbb{Q}(\zeta\_{e})$; it is a homomorphism $G\to\mu\_{e}$, so it is identified with the element $\psi\in\widehat{G^{ab}}$ agreeing with it on the generators of $G$ ([cyclotomic refinement](cyclotomic-refinement.md) §1), and $o(\lambda\_\chi)$ and $\ker\lambda\_\chi$ follow.

*Proof.* The determinant is the product of the eigenvalues with multiplicity. $\square$

## 3. Cross-check I: $\mathfrak{f}(\det\rho\_\chi)$ divides $\mathfrak{f}(\chi)$

**Theorem 3.1.** For every character $\chi$ and every $\ell$: $f\_\ell(\lambda\_\chi)\le f\_\ell(\chi)$, and more precisely the tame parts satisfy $[\lambda\_\chi|*{D\_0}\ne1]\le\mathrm{codim},V^{D\_0}$ and the Swan parts $\mathrm{sw}*\ell(\lambda\_\chi)\le\mathrm{sw}*\ell(\chi)$. Hence $\mathfrak{f}(\det\rho*\chi)\mid\mathfrak{f}(\chi)$.

*Proof.* If $\lambda\_\chi|*{D\_i}\ne1$ then $D\_i$ acts nontrivially on $V$, so $\mathrm{codim},V^{D\_i}\ge1$. With $c=\max{i:\lambda*\chi|*{D\_i}\ne1}$ (or no term if $\lambda*\chi|*{D\_0}=1$),* [*conductor exponents*](conductor-exponents.md) *Theorem 2.1 gives $f*\ell(\lambda\_\chi)=\sum\_{i=0}^{c}\frac{|D\_i|}{|D\_0|}\le\sum\_{i=0}^{c}\frac{|D\_i|}{|D\_0|}\mathrm{codim},V^{D\_i}\le f\_\ell(\chi)$; the same inequality restricted to $i=0$ and to $i\ge1$ gives the tame and Swan statements. $\square$

## 4. Cross-check II: $\mathfrak{f}(\det\rho\_\chi)$ against the abelian data

**Theorem 4.1.** For every linear character $\lambda$ of $G$, $\mathfrak{f}(\lambda)=c\_{\chi\_\lambda}$, the conductor of the Dirichlet character $\chi\_\lambda$ with $\chi\_\lambda(\ell)=\lambda(\mathrm{Frob}*\ell)$ (*[*cyclotomic refinement*](cyclotomic-refinement.md) *§1). Locally: $f*\ell(\lambda)$ is the least $n\ge0$ such that the character $\lambda\circ\mathrm{rec}*\ell$ of $\mathbb{Q}*\ell^\times$ is trivial on $U^{(n)}=1+\ell^n\mathbb{Z}*\ell$ ($U^{(0)}=\mathbb{Z}*\ell^\times$).

*Proof.* By §1 the conductor may be computed in $\mathrm{Gal}(N\_\lambda/\mathbb{Q})$. For the cyclic extension $N\_\lambda$ with decomposition group $H$ at $\ell$ and $\lambda$ faithful on $H$: $f\_\ell(\lambda)=1+\varphi(c)$ with $c$ the largest lower index where $\lambda|*{H\_c}\ne1$, i.e. $\varphi(c)$ is the largest upper jump of $H$ (*[*conductor exponents*](conductor-exponents.md) *Theorem 2.1). Local class field theory identifies the upper filtration with the unit filtration: $\mathrm{rec}*\ell(U^{(v)})=H^v$ for all $v\ge0$ (Serre XV §2 Thm 1'), so $\lambda\circ\mathrm{rec}\_\ell$ is trivial on $U^{(n)}$ iff $H^n=1$ iff $n>\varphi(c)$, i.e. iff $n\ge\varphi(c)+1$. The conductor of a Dirichlet character is the product of these local exponents. $\square$

The check: for every $\chi\in\mathrm{Irr}(G)$, compute $\lambda\_\chi$ by Proposition 2.1, identify it in $\widehat{G^{ab}}$, and require $\prod\_\ell\ell^{f\_\ell(\lambda\_\chi)}=c\_{\lambda\_\chi}$ with $f\_\ell(\lambda\_\chi)$ from [conductor exponents](conductor-exponents.md) and $c\_{\lambda\_\chi}$ from [cyclotomic refinement](cyclotomic-refinement.md). Since [cyclotomic refinement](cyclotomic-refinement.md) obtained the conductors from the abelian subfields and the kernels from rational Frobenius data, and [conductor exponents](conductor-exponents.md) obtained the exponents from the ramification filtration, the two sides are independent. For an odd two-dimensional $\chi$ (family 3 of the project), $\mathfrak{f}(\chi)$ is the level and $\lambda\_\chi$ the nebentypus of the weight-one form, so this check is the first half of the comparison with the newform.

## 5. Cross-check III: discriminants of the cyclic fields $N\_\lambda$

**Theorem 5.1.** For a linear character $\lambda$ of order $o$,

$$|d\_{N\_\lambda}|=\prod\_{k=0}^{o-1}\mathfrak{f}(\lambda^k),\qquad v\_\ell(d\_{N\_\lambda})=\sum\_{k=0}^{o-1}f\_\ell(\lambda^k)\ \text{ for every }\ell,$$

where the $f\_\ell(\lambda^k)$ are computed in $G$ by [conductor exponents](conductor-exponents.md), and $d\_{N\_\lambda}$ is computed independently by [the ramified primes](ramified-primes.md) applied to a resolvent of degree $o$ (an invariant with stabilizer $\ker\lambda$, [separating subgroups](separating-subgroups.md) §3; $N\_\lambda=\mathbb{Q}(\beta\_\lambda)$ for its value $\beta\_\lambda$, and $\mathrm{disc},\mathcal{O}*{N*\lambda}$ from the $\ell$-maximal orders at the primes dividing $\mathrm{disc}$ of the resolvent). Equality is required for every $\lambda\in\Lambda\_\chi:={\lambda\_\chi}\cup{\lambda\text{ linear}:\langle\chi\bar\chi,\lambda\rangle\ne0}$, and for every $\chi$.

*Proof.* [conductor exponents](conductor-exponents.md) Theorem 4.1 for the field $N\_\lambda$ with group $C=\mathrm{Gal}(N\_\lambda/\mathbb{Q})\cong\langle\lambda\rangle^{\vee}$, whose irreducible characters are the $\lambda^k$ (all of degree $1$), together with §1 (inflation). $\square$

## 6. Cross-check IV: the linear constituents of $\chi\bar\chi$

**Lemma 6.1.** For $\chi\in\mathrm{Irr}(G)$ and $\lambda$ linear, $\langle\chi\bar\chi,\lambda\rangle=\langle\chi,\lambda\chi\rangle\in{0,1}$, and it is $1$ iff $\lambda\chi=\chi$. The set $\mathrm{Stab}(\chi)={\lambda:\lambda\chi=\chi}$ is a subgroup of $\widehat{G^{ab}}$; every $\lambda\in\mathrm{Stab}(\chi)$ satisfies $\ker\lambda\supseteq{\sigma:\chi(\sigma)\ne0}$, and $|\mathrm{Stab}(\chi)|$ divides $\chi(1)$.

*Proof.* $\langle\chi\bar\chi,\lambda\rangle=\frac{1}{|G|}\sum\chi\bar\chi\bar\lambda=\langle\chi,\lambda\chi\rangle$, and $\lambda\chi$ is irreducible. If $\lambda\chi=\chi$ then $\chi(\sigma)(\lambda(\sigma)-1)=0$ for all $\sigma$. For the last claim let $K\_0=\bigcap\_{\lambda\in\mathrm{Stab}(\chi)}\ker\lambda$, so $G/K\_0$ is abelian and $\mathrm{Stab}(\chi)\subseteq\mathrm{Irr}(G/K\_0)$. Isaacs, *Character Theory* 6.17: if $\chi|*{K\_0}=e\sum*{i=1}^t\theta\_i$ with distinct irreducible $\theta\_i$, then the number of $\mu\in\mathrm{Irr}(G/K\_0)$ with $\mu\chi=\chi$ is $et$. Hence $|\mathrm{Stab}(\chi)|=et$, which divides $\chi(1)=et,\theta\_1(1)$. For a single $\lambda$ of order $o$ with $K=\ker\lambda$, the same theorem with $K$ in place of $K\_0$ gives $o=et$, so $\chi(1)=o,\theta(1)$ for a constituent $\theta$ of $\chi|\_K$; as $\mathrm{Ind}\_K^G\theta$ has degree $o,\theta(1)$ and contains $\chi$ (Frobenius reciprocity), $\chi=\mathrm{Ind}\_K^G\theta$. $\square$

**Theorem 6.2 (global induction formula).** For $K\le G$, $F=N^K$, and a character $\theta$ of $K$:

$$\mathfrak{f}(\mathrm{Ind}*K^G\theta)=|d\_F|^{\theta(1)}\cdot N*{F/\mathbb{Q}}\big(\mathfrak{f}\_F(\theta)\big),$$

where $\mathfrak{f}*F(\theta)$ is the Artin conductor of $\theta$ over $F$ (an ideal of $F$) and $N$ the ideal norm. Per prime: $f*\ell(\mathrm{Ind}*K^G\theta)=\theta(1),v*\ell(d\_F)+\sum\_{\mathfrak{p}\mid\ell}f(\mathfrak{p}/\ell),f\_{\mathfrak{p}}(\theta)$.

*Proof.* Mackey: $\mathrm{Res}*D\mathrm{Ind}K^G\theta=\sum{DgK}\mathrm{Ind}*{D\cap gKg^{-1}}^D(\theta^g)$, the double cosets corresponding to the primes $\mathfrak{p}$ of $F$ above $\ell$, with $D\cap gKg^{-1}$ the Galois group of $L/F\_{\mathfrak{p}}$. [conductor exponents](conductor-exponents.md) Proposition 2.2 for each term gives $f(\mathfrak{p}/\ell)f\_{\mathfrak{p}}(\theta^g)+v\_\ell(\mathfrak{d}*{F*{\mathfrak{p}}/\mathbb{Q}*\ell})\theta(1)$, and $\sum*{\mathfrak{p}\mid\ell}v\_\ell(\mathfrak{d}*{F*{\mathfrak{p}}/\mathbb{Q}*\ell})=v*\ell(d\_F)$. $\square$

**Theorem 6.3.** Let $\chi\in\mathrm{Irr}(G)$ and $\lambda\in\mathrm{Stab}(\chi)$ of order $o$, $K=\ker\lambda$, so $\chi=\mathrm{Ind}\_K^G\theta$ with $\theta(1)=\chi(1)/o$. Then for every $\ell$:

1. $f\_\ell(\chi)\ \ge\ \dfrac{\chi(1)}{o},v\_\ell(d\_{N\_\lambda})$, i.e. $|d\_{N\_\lambda}|^{\chi(1)/o}$ divides $\mathfrak{f}(\chi)$;
2. $f\_\ell(\chi)-\dfrac{\chi(1)}{o},v\_\ell(d\_{N\_\lambda})$ is a nonnegative multiple of the residue degree $f\_\ell(N\_\lambda)$ of $\ell$ in $N\_\lambda$ (equal to $f(\mathfrak{p}/\ell)$ for every $\mathfrak{p}\mid\ell$, the field being Galois), and vanishes iff $\theta$ is unramified at every prime of $N\_\lambda$ above $\ell$;
3. combining with Theorem 5.1, $f\_\ell(\chi)\ge\dfrac{\chi(1)}{o}\sum\_{k=0}^{o-1}f\_\ell(\lambda^k)$, an inequality between quantities all computed by [conductor exponents](conductor-exponents.md).

*Proof.* Theorem 6.2 with $F=N\_\lambda$, $\theta(1)=\chi(1)/o$; the primes above $\ell$ in a Galois extension have a common residue degree, so the sum $\sum\_{\mathfrak{p}}f(\mathfrak{p}/\ell)f\_{\mathfrak{p}}(\theta)$ is a nonnegative multiple of $f\_\ell(N\_\lambda)$, zero iff all $f\_{\mathfrak{p}}(\theta)=0$. Item 3 substitutes Theorem 5.1. $\square$

**The check.** For every $\chi$ and every $\lambda\in\mathrm{Stab}(\chi)$ (found from the table by testing $\lambda\chi=\chi$ over $\widehat{G^{ab}}$; nontrivial only when $\chi$ vanishes on a coset of a proper normal subgroup of index $o(\lambda)$), require items 1 to 2 of Theorem 6.3 at every ramified $\ell$, with $v\_\ell(d\_{N\_\lambda})$ from [the ramified primes](ramified-primes.md) (independent of the filtration) and $f\_\ell(N\_\lambda)=[D\ker\lambda\:D\_0\ker\lambda]$, the residue degree of the image of $D$ in $G/\ker\lambda$ (inertia maps onto inertia in quotients), equivalently the common degree of the irreducible factors of the degree-$o$ resolvent over $\mathbb{Q}\_\ell$.

## 7. A bound on the conductor of $\chi\bar\chi$

**Proposition 7.1.** $f\_\ell(\chi\bar\chi)\le2\chi(1)f\_\ell(\chi)$, and $f\_\ell(\psi)\le f\_\ell(\chi\bar\chi)$ for every irreducible constituent $\psi$ of $\chi\bar\chi$; in particular $f\_\ell(\lambda)\le2\chi(1)f\_\ell(\chi)$ for $\lambda\in\mathrm{Stab}(\chi)$, and $\det(\rho\otimes\rho^\*)=1$ so $\mathfrak{f}(\det(\chi\bar\chi))=1$.

*Proof.* With $a\_i=\dim V^{D\_i}$, $d=\chi(1)$: $(V\otimes V^*)^{D\_i}\supseteq V^{D\_i}\otimes(V^*)^{D\_i}$ has dimension $\ge a\_i^2$, so $\mathrm{codim}(V\otimes V^\*)^{D\_i}\le d^2-a\_i^2\le2d(d-a\_i)$; sum with the weights $|D\_i|/|D\_0|$. Conductors of constituents are bounded by the conductor of the sum because every term of [conductor exponents](conductor-exponents.md) Proposition 1.1 is nonnegative. $\square$

This is weaker than §6 and is recorded only as a sanity bound; the substantive constraints from $\chi\bar\chi$ are those of Theorem 6.3 through its linear constituents.

## 8. Summary of the cross-checks on the exponents of [conductor exponents](conductor-exponents.md)

For each $\chi\in\mathrm{Irr}(G)$ and each ramified $\ell$:

| check statement independent input                                    |                                                                                                                                         |                                                                                                                          |
| -------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| I                                                                    | $f\_\ell(\det\rho\_\chi)\le f\_\ell(\chi)$, tame and Swan parts separately                                                              | table only                                                                                                               |
| II                                                                   | $\prod\_\ell\ell^{f\_\ell(\det\rho\_\chi)}=c\_{\lambda\_\chi}$                                                                          | Dirichlet conductors of [cyclotomic refinement](cyclotomic-refinement.md)                         |
| III                                                                  | $v\_\ell(d\_{N\_\lambda})=\sum\_kf\_\ell(\lambda^k)$ for $\lambda\in\Lambda\_\chi$                                                      | discriminants of degree-$o(\lambda)$ resolvent fields ([the ramified primes](ramified-primes.md)) |
| IV                                                                   | $f\_\ell(\chi)-\frac{\chi(1)}{o}v\_\ell(d\_{N\_\lambda})\in f\_\ell(N\_\lambda),\mathbb{Z}\_{\ge0}$ for $\lambda\in\mathrm{Stab}(\chi)$ | same discriminants, residue degrees                                                                                      |
| [conductor exponents](conductor-exponents.md) | $\sum\_\chi\chi(1)f\_\ell(\chi)=v\_\ell(d\_N)$; $f\_\ell(\pi\_{H\_j})=v\_\ell(\mathrm{disc},\mathcal{O}\_{K\_j})$                       | local different; [the ramified primes](ramified-primes.md)                                        |

Checks I, III and IV are consequences of the induction formula and of Hasse to Arf, and fail only if the filtration or the class identification is wrong; Check II ties the filtration to class field theory and fails if the pinning of [cyclotomic refinement](cyclotomic-refinement.md) or the ramification data disagree. All are decided by integer arithmetic on the outputs of [the ramified primes](ramified-primes.md), [cyclotomic refinement](cyclotomic-refinement.md) and [conductor exponents](conductor-exponents.md), at no additional $\ell$-adic cost beyond the degree-$o(\lambda)$ resolvents of Check III, whose fields are small ($o(\lambda)\le\exp G^{ab}$).