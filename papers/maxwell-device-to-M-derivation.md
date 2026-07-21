# Maxwell + device geometry → M : a parameter-free derivation

**Author:** derivation note (companion to `optics-ou-analogies-and-limits.md`, Q2)
**Date:** July 2026
**Status:** Formal derivation + verified numerics. No tesseract ad hoc free-parameter tuning.
**Code:** `scripts/derive_M_maxwell_device.py` (all three checks pass)

This note closes Q2 of `optics-ou-analogies-and-limits.md` ("How does a tesseract
projection become an ABCD matrix on an SLM?") by *deriving* the transfer operator
of a real optical device from Maxwell, and then asking — once the derivation is on
the table — whether a 4-cube / Coxeter B₄ / "8 cubic cells" structure can appear
in it. It cannot, except as a *constraint* that reduces freedom. The derivation is
independent of any polytope.

---

## 1. From Maxwell to the paraxial wave operator

Source-free Maxwell in a non-magnetic dielectric ε(r), monochomorphic field
𝔈(r,t)=Re[E(r)e^{−iωt}], gives the Helmholtz equation

\[
(\nabla^2 + k^2 n^2(\mathbf r)) E = 0,\qquad k=\omega/c .
\]

Write \(E(\mathbf r_\perp,z)=e^{ikz}u(\mathbf r_\perp,z)\) and assume paraxial
propagation (\(|\partial_z^2 u|\ll k|\partial_z u|{ }). The slow envelope obeys the
paraxial wave equation (Schrödinger form, \(z\) as "time", mass \(k\)):

\[
\partial_z u = \frac{i}{2k}\nabla_\perp^2 u + (\text{index terms}).
\]

Free space (\(n=1\)) is the first term alone.

## 2. Free propagation = the Fresnel operator \(P_d\)

The solution is a linear operator on \(L^2(\mathbb R^2)\):

\[
(P_d u)(\mathbf r'_\perp)=\frac{k}{2\pi i d}\int d^2r_\perp 
e^{i k|\mathbf r'_\perp-\mathbf r_\perp|^2/(2d)} u(\mathbf r_\perp).
\]

\(P_d\) is unitary. It is the **metaplectic** quantisation of the ray matrix
\(\bigl(\begin{smallmatrix}1&d\\0&1\end{smallmatrix}\bigr)\). This is the free-space
transfer; nothing beyond Maxwell + the paraxial approximation was used.

## 3. Thin phase element (SLM / lens) = a diagonal operator \(L_\Phi\)

A thin element of optical thickness adds a phase \(\Phi(\mathbf r_\perp)\):

\[
(L_\Phi u)(\mathbf r_\perp)=e^{i\Phi(\mathbf r_\perp)} u(\mathbf r_\perp).
\]

A thin lens of focal length \(f\) is the quadratic specialisation
\(\Phi_f=-k|\mathbf r_\perp|^2/(2f)\). A general phase-only SLM is an *arbitrary*
\(\Phi(\mathbf r_\perp)\) — **not**, in general, a 2×2 matrix.

## 4. The device transfer operator \(M\)

A device is a sequence of phase screens and free propagations:

\[
\boxed{ M_{\rm device}=P_{d_n} L_{\Phi_{n-1}}\cdots P_{d_2} L_{\Phi_1} P_{d_1} }
\]

a composition of Fresnel kernels with phase insertions — i.e. an **integral
operator** on the field, not a matrix, unless every \(\Phi\) is quadratic.

### Quadratic phase ⇒ ABCD (metaplectic / Collins integral)

If every \(\Phi\) is quadratic the operator acts linearly on the ray phase space
\((\mathbf r_\perp,\boldsymbol\theta)\) and reduces to a symplectic matrix
\(M=\bigl(\begin{smallmatrix}A&B\\C&D\end{smallmatrix}\bigr)\in Sp(2,\mathbb R)\),
with the closed-form kernel (Collins integral)

\[
u_{\rm out}(\mathbf r'_\perp)=\frac{1}{i\lambda B}\int d^2r_\perp 
\exp\!\Big[\frac{i}{2B}\big(A|\mathbf r_\perp|^2-2\mathbf r_\perp\!\cdot\!\mathbf r'_\perp+D|\mathbf r'_\perp|^2\big)\Big]u_{\rm in}.
\]

Free propagation \(d\) and thin lens \(f\) give the textbook factors
\(P_d=\bigl(\begin{smallmatrix}1&d\\0&1\end{smallmatrix}\bigr)\),
\(L_f=\bigl(\begin{smallmatrix}1&0\\-1/f&1\end{smallmatrix}\bigr)\); the imaging
condition \(B_{\rm tot}=0\) is Descartes \(1/s_o+1/s_i=1/f\). For two transverse
planes \((x,\theta_x,y,\theta_y)\) the matrix is \(4\times4\) symplectic,
\(M\in Sp(4,\mathbb R)\). This is the only legitimate "4D" object.

**Verification (`scripts/derive_M_maxwell_device.py`, check 1).** A Gaussian waist
propagated through \(P_{d_1}L_fP_{d_2}\) by the *integral operators* of §2–§3 is
compared to the analytic ABCD law \(q_{\rm out}=(Aq_{\rm in}+B)/(Cq_{\rm in}+D)\)
applied to the Gaussian \(q\)-parameter. Result:

```
q_ABCD = 0.11554 + 0.08677 i
q_num = 0.11554 + 0.08677 i -> relative error 2.2e-13
```

The device operator *is* the ABCD action on the mode, to machine precision.

## 5. Controllable degrees of freedom = hardware count, not group dimension

The space of achievable transfer operators is the space of achievable \(\Phi\)
functions. A phase-only SLM with \(N_{\rm pix}\) independently addressable pixels
gives \(\Phi\in[0,2\pi)^{N_{\rm pix}}\): **\(N_{\rm pix}\) real parameters**.
Adding spectral pulse shaping with \(N_\omega\) bins multiplies this by
\(N_\omega\). This is the actual controllable dimension. It is set by the device,
not by the dimensionality of an abstract group.

**Verification (check 2).** A 1920×1080 SLM: \(N_{\rm pix}=2{,}073{,}600\); with 128
spectral bins, \(2.65\times10^8\). No "8", "16", "24", or "32" appears.

## 6. Where a 4-cube / B₄ could enter — and what it actually does

The tesseract is the 4-cube \(\{|x_i|\le1\}\subset\mathbb R^4\); its boundary has 8
cubic cells. Its Coxeter group \(B_4\) (hyperoctahedral, signed permutations of
four coordinates, \(|B_4|=2^4\cdot4!=384\)) is real. The question is whether it
enters §1–§4. It can, in exactly one way: **as a discrete symmetry of the optical
phase space** \((x,\theta_x,y,\theta_y)\cong\mathbb R^4\), i.e. as a subgroup of
\(Sp(4,\mathbb R)\). A device invariant under such a symmetry is *constrained* by
it — its transfer matrix must commute with the group — so the symmetry reduces the
available operator to the group-invariant subalgebra. It cannot create channels.

**Verification (check 3).** Enumerate all 384 signed-permutation matrices; keep only
those preserving the optical symplectic form \(J={\rm blockdiag}(J_2,J_2)\):

```
|B4| = 384
|B4 ∩ Sp(4,R)| = 32 (only 32 of 384 preserve the symplectic structure)
dim sp(4,R) = 10
dim(commutant) = 1 (B4-symmetric devices live in a 1-dim subspace)
```

A \(B_4\)-type symmetry collapses the 10-dimensional \(sp(4,\mathbb R)\) to a
**1-dimensional** invariant subspace. The 4-cube geometry *removes* nine
parameters of freedom; it adds none. The combinatorial counts of the 4-cube
(8 cells / 16 vertices / 24 faces / 32 edges) appear nowhere as a rank or DOF.

## 7. Conclusion (closes Q2)

Starting from Maxwell and a physical device (phase screen + propagation, optional
spectral shaping), the derived transfer operator is

\[
M_{\rm device}=P_{d_n}L_{\Phi_{n-1}}\cdots L_{\Phi_1}P_{d_1}\in{\cal U}(L^2(\mathbb R^2)),
\]

with the quadratic (first-order) specialisation \(M\in Sp(4,\mathbb R)\) and the
controllable dimension \(N_{\rm pix} [\times N_\omega]\). **A "tesseract / B₄ /
8-cube" structure is not produced by this derivation.** It enters only if inserted
by hand as an ansatz; inserted as a *symmetry* it reduces freedom
(\(10\to1\)); inserted as a *parametrisation* of \(M\)'s free parameters using
4-cube combinatorics it is decoration (the same "undeclared physical power"
pattern as \(\ln4\equiv\omega_R\), flagged in `no-go-superoscillation-tesseract.md`).
There is no derived map "8 cubes → 8 optical channels."

The honest ABCD/engineering programme (§4) stands; the tesseract programme does not,
until someone supplies a Maxwell-and-geometry derivation that inserts the 4-cube —
which, by the symmetry-count result above, can only *limit* the device, not extend
it.

---

*End of derivation note. Companion code: `scripts/derive_M_maxwell_device.py`.*
