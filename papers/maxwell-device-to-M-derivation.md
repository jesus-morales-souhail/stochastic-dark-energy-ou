# Maxwell + device geometry → transfer operator $M$

Jesús Morales Souhail · [github.com/jesus-morales-souhail](https://github.com/jesus-morales-souhail)  
July 2026 · hygiene / pedagogy · **not a DESI claim** · not peer reviewed  

Code: `scripts/derive_M_maxwell_device.py` (checks pass)  
Related: `papers/optics-ou-analogies-and-limits.md` (Q2), `papers/no-go-superoscillation-tesseract.md`, `papers/EXPLORATORY_BOUNDARY.md`

---

## What this note is

I close Q2 of `optics-ou-analogies-and-limits.md` — “How does a tesseract projection become an ABCD matrix on an SLM?” — by deriving the transfer operator of a real optical device from Maxwell, and then asking whether a 4-cube / Coxeter group $B_{4}$ / “8 cubic cells” structure can appear in that derivation.

It cannot, except as a *constraint* that reduces freedom. The derivation is independent of any polytope.  
This is lab-optics hygiene, **not** a BAO residual result and **not** a DESI claim.

---

## 1. From Maxwell to the paraxial wave operator

Source-free Maxwell in a non-magnetic dielectric $\varepsilon(\mathbf{r})$, monochromatic field
$\mathfrak{E}(\mathbf{r},t)=\mathrm{Re}[E(\mathbf{r})e^{-i\omega t}]$, gives the Helmholtz equation

$$
(\nabla^{2}+k^{2}n^{2}(\mathbf{r}))E=0,\qquad k=\omega/c.
$$

Write $E(\mathbf{r}_{\perp},z)=e^{ikz}u(\mathbf{r}_{\perp},z)$ and assume paraxial propagation
($|\partial_{z}^{2}u|\ll k|\partial_{z}u|$). The slow envelope obeys the paraxial wave equation
(Schrödinger form, $z$ as “time”, mass $k$):

$$
\partial_{z}u=\frac{i}{2k}\nabla_{\perp}^{2}u+(\text{index terms}).
$$

Free space ($n=1$) is the first term alone.

---

## 2. Free propagation = the Fresnel operator $P_{d}$

The solution is a linear operator on $L^{2}(\mathbb{R}^{2})$:

$$
(P_{d}u)(\mathbf{r}'_{\perp})
=\frac{k}{2\pi i d}
\int\mathrm{d}^{2}r_{\perp}\,
\exp\!\Big(i k|\mathbf{r}'_{\perp}-\mathbf{r}_{\perp}|^{2}/(2d)\Big)\,
u(\mathbf{r}_{\perp}).
$$

$P_{d}$ is unitary. It is the **metaplectic** quantisation of the ray matrix
$\bigl(\begin{smallmatrix}1&d\\0&1\end{smallmatrix}\bigr)$. This is free-space transfer;
nothing beyond Maxwell + the paraxial approximation was used.

---

## 3. Thin phase element (SLM / lens) = a diagonal operator $L_{\Phi}$

A thin element of optical thickness adds a phase $\Phi(\mathbf{r}_{\perp})$:

$$
(L_{\Phi}u)(\mathbf{r}_{\perp})=e^{i\Phi(\mathbf{r}_{\perp})}u(\mathbf{r}_{\perp}).
$$

A thin lens of focal length $f$ is the quadratic specialisation
$\Phi_{f}=-k|\mathbf{r}_{\perp}|^{2}/(2f)$. A general phase-only SLM is an *arbitrary*
$\Phi(\mathbf{r}_{\perp})$ — **not**, in general, a $2\times 2$ matrix.

---

## 4. The device transfer operator $M$

A device is a sequence of phase screens and free propagations:

$$
M_{\mathrm{device}}
=P_{d_{n}}L_{\Phi_{n-1}}\cdots P_{d_{2}}L_{\Phi_{1}}P_{d_{1}},
$$

a composition of Fresnel kernels with phase insertions — i.e. an **integral operator** on the field, not a matrix, unless every $\Phi$ is quadratic.

### Quadratic phase ⇒ ABCD (metaplectic / Collins integral)

If every $\Phi$ is quadratic, the operator acts linearly on the ray phase space
$(\mathbf{r}_{\perp},\boldsymbol{\theta})$ and reduces to a symplectic matrix
$M=\bigl(\begin{smallmatrix}A&B\\C&D\end{smallmatrix}\bigr)\in\mathrm{Sp}(2,\mathbb{R})$,
with the closed-form kernel (Collins integral)

$$
u_{\mathrm{out}}(\mathbf{r}'_{\perp})
=\frac{1}{i\lambda B}
\int\mathrm{d}^{2}r_{\perp}\,
\exp\!\Big[
\frac{i}{2B}\big(
A|\mathbf{r}_{\perp}|^{2}
-2\mathbf{r}_{\perp}\!\cdot\!\mathbf{r}'_{\perp}
+D|\mathbf{r}'_{\perp}|^{2}
\big)
\Big]
u_{\mathrm{in}}.
$$

Free propagation $d$ and thin lens $f$ give the textbook factors
$P_{d}=\bigl(\begin{smallmatrix}1&d\\0&1\end{smallmatrix}\bigr)$,
$L_{f}=\bigl(\begin{smallmatrix}1&0\\-1/f&1\end{smallmatrix}\bigr)$;
the imaging condition $B_{\mathrm{tot}}=0$ is Descartes $1/s_{o}+1/s_{i}=1/f$.
For two transverse planes $(x,\theta_{x},y,\theta_{y})$ the matrix is $4\times 4$ symplectic:

$$
M\in\mathrm{Sp}(4,\mathbb{R}).
$$

That is the only legitimate “4D” object I accept from this derivation.

**Verification** (`scripts/derive_M_maxwell_device.py`, check 1). A Gaussian waist
propagated through $P_{d_{1}}L_{f}P_{d_{2}}$ by the integral operators of §§2–3 is
compared to the analytic ABCD law $q_{\mathrm{out}}=(A q_{\mathrm{in}}+B)/(C q_{\mathrm{in}}+D)$.
Result:

```
q_ABCD = 0.11554 + 0.08677 i
q_num  = 0.11554 + 0.08677 i -> relative error 2.2e-13
```

The device operator *is* the ABCD action on the mode, to machine precision.

---

## 5. Controllable degrees of freedom = hardware count, not group dimension

The space of achievable transfer operators is the space of achievable $\Phi$ functions.
A phase-only SLM with $N_{\mathrm{pix}}$ independently addressable pixels gives
$\Phi\in[0,2\pi)^{N_{\mathrm{pix}}}$: **$N_{\mathrm{pix}}$ real parameters**.
Adding spectral pulse shaping with $N_{\omega}$ bins multiplies this by $N_{\omega}$.
This is the actual controllable dimension. It is set by the device, not by the dimension of an abstract group.

**Verification (check 2).** A $1920\times 1080$ SLM: $N_{\mathrm{pix}}=2{,}073{,}600$; with 128 spectral bins, $\sim 2.65\times 10^{8}$. No “8”, “16”, “24”, or “32” appears from polytope combinatorics.

---

## 6. Where a 4-cube / $B_{4}$ could enter — and what it actually does

The tesseract is the 4-cube $\{|x_{i}|\le 1\}\subset\mathbb{R}^{4}$; its boundary has 8 cubic cells.
Its Coxeter group $B_{4}$ (hyperoctahedral, signed permutations of four coordinates) has order

$$
\lvert $B_{4}$\rvert = 2^{4}\cdot 4! = 384.
$$

The question is whether $B_{4}$ enters §§1–4. It can, in exactly one way: **as a discrete symmetry of the optical phase space**
$(x,\theta_{x},y,\theta_{y})\cong\mathbb{R}^{4}$, i.e. as a subgroup of $\mathrm{Sp}(4,\mathbb{R})$.
A device invariant under such a symmetry is *constrained* by it — its transfer matrix must commute with the group —
so the symmetry reduces the available operator to the group-invariant subalgebra. It cannot create channels.

**Verification (check 3).** Enumerate all 384 signed-permutation matrices; keep only those preserving the optical symplectic form
$J=\mathrm{blockdiag}(J_{2},J_{2})$:

$$
\lvert $B_{4}$\cap\mathrm{Sp}(4,\mathbb{R})\rvert = 32
\quad\text{(only 32 of 384 preserve the symplectic structure).}
$$

```
|B4| = 384
|B4 ∩ Sp(4,R)| = 32
dim sp(4,R) = 10
dim(commutant) = 1
```

A $B_{4}$-type symmetry collapses the 10-dimensional algebra $\mathfrak{sp}(4,\mathbb{R})$ to a
**1-dimensional** invariant subspace. The 4-cube geometry *removes* nine parameters of freedom; it adds none.
The combinatorial counts of the 4-cube (8 cells / 16 vertices / 24 faces / 32 edges) appear nowhere as a rank or degree of freedom.

---

## 7. Conclusion (closes Q2)

Starting from Maxwell and a physical device (phase screen + propagation, optional spectral shaping), the derived transfer operator is

$$
M_{\mathrm{device}}
=P_{d_{n}}L_{\Phi_{n-1}}\cdots L_{\Phi_{1}}P_{d_{1}}
\in\mathcal{U}\big(L^{2}(\mathbb{R}^{2})\big),
$$

with the quadratic (first-order) specialisation

$$
M\in\mathrm{Sp}(4,\mathbb{R})
$$

and the controllable dimension $N_{\mathrm{pix}}\,[\times N_{\omega}]$.

**A “tesseract / $B_{4}$ / 8-cube” structure is not produced by this derivation.**  
It enters only if inserted by hand as an ansatz; inserted as a *symmetry* it reduces freedom ($10\to 1$);
inserted as a *parametrisation* of free parameters using 4-cube combinatorics it is decoration
(the same “undeclared physical power” pattern as $\ln 4\equiv\omega_{R}$, flagged in `no-go-superoscillation-tesseract.md`).
There is no derived map “8 cubes → 8 optical channels.”

I keep the honest ABCD/engineering programme (§4). The tesseract programme does not stand until someone supplies a Maxwell-and-geometry derivation that inserts the 4-cube — which, by the symmetry-count result above, can only *limit* the device, not extend it.

This note is **not** a DESI claim. Companion code: `scripts/derive_M_maxwell_device.py`.
