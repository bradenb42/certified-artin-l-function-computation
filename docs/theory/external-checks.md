# External correctness checks

The artifact produced by [the character table](character-table.md) to [the verifier](verifier.md) for a monic separable $f\in\mathbb{Z}[x]$ with proven Galois group $G$ must satisfy the following three checks, each against an independent description of the same outputs. A check applies to $f$ only when $f$ lies in the corresponding family; the artifact records which checks applied and whether they passed.

## Check 1. Abelian $G$, class field theory

Applies when $G$ is abelian. Then $N\subseteq\mathbb{Q}(\zeta\_c)$ for the conductor $c$ of $N$, every $\chi\in\mathrm{Irr}(G)$ is one-dimensional and corresponds to a primitive Dirichlet character $\psi\_\chi$ of conductor $\mathfrak{f}(\chi)$, and every output has a closed form:

- $\mathfrak{f}(\chi)=$ the conductor of $\psi\_\chi$;
- for every unramified $\ell\le X$, $\mathrm{Frob}\_\ell$ is the element of $G$ corresponding to $\ell\bmod c$ under the fixed isomorphism $G\cong(\mathbb{Z}/c)^\times/U$;
- for every $\ell$, $P\_\ell(\chi;T)=1-\psi\_\chi(\ell)T$ ($\psi\_\chi(\ell)=0$ when $\ell\mid\mathfrak{f}(\chi)$);
- $W(\chi)=\tau(\psi\_\chi)/\big(i^{\kappa}\sqrt{\mathfrak{f}(\chi)}\big)$ with $\psi\_\chi(-1)=(-1)^\kappa$;
- $\gamma\_\chi(s)=\Gamma\_\mathbb{R}(s+\kappa)$.

The artifact's conductors, Frobenius classes, Euler factors and root numbers must agree with these, character by character, for every $\chi$ and every prime up to $X$.

## Check 2. $\deg f\le5$, explicit splitting field

Applies when $n=\deg f\le5$. Then $N$ is constructed explicitly (a defining polynomial of degree $[N:\mathbb{Q}]=|G|\le120$ and its maximal order), and for every ramified prime $\ell$ the decomposition group $D\_\ell$, its inertia group, the full lower and upper ramification filtration, and the Artin conductor exponents $f\_\ell(\chi)$ for every $\chi$ are read off the factorization of $\ell$ in $\mathcal{O}\_N$ and the ramification groups of the primes above $\ell$ computed directly in $N$.

The artifact's decomposition groups (up to $G$-conjugacy), filtrations (with their Herbrand conversions) and conductor exponents must agree with these, prime by prime and character by character.

## Check 3. Odd irreducible two-dimensional $\chi$ with $G$ solvable, weight-one newforms

Applies to every $\chi\in\mathrm{Irr}(G)$ with $\chi(1)=2$, $\chi(c)=0$ and $G$ solvable. Then there is a holomorphic newform $g\_\chi\in S\_1(\Gamma\_0(\mathfrak{f}(\chi)),\det\rho\_\chi)$ with $L(s,g\_\chi)=L(s,\chi)$ (Deligne to Serre; existence by Langlands to Tunnell). The newform is computed independently (as a weight-one eigenform of the stated level and nebentypus, with its coefficient field embedded compatibly with $\mathbb{Q}(\chi)\hookrightarrow\mathbb{C}$).

The artifact's coefficients $a\_m(\chi)$ for all $m\le X$ and its root number $W(\chi)$ must equal the Fourier coefficients $a\_m(g\_\chi)$ and the Atkin to Lehner root number of $g\_\chi$, for every such $\chi$.

## Outside the three families

For $f$ and $\chi$ not covered by Checks 1 to 3, the artifact returns the computed data (decomposition groups and filtrations at the ramified primes, Frobenius classes up to $X$, conductors, Euler factors, Gamma factors, root numbers, and the functional-equation defects against their bounds) together with the certificate of [the certificate](certificate.md) and the verdict of the verifier of [the verifier](verifier.md). It certifies nothing beyond what that verdict states: acceptance proves the algebraic outputs conditional on the certificate for $G$ ([the verifier](verifier.md) Theorem 2.1, Corollary 2.2, Theorem 3.1), and the functional-equation test is a rejection test ([the functional-equation test](functional-equation.md)).