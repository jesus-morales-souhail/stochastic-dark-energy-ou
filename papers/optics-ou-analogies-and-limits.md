# Three questions: quantum lens, tesseract→ABCD, and “focal length” in the OU model

**Author:** Jesús Morales Souhail · [github.com/jesus-morales-souhail](https://github.com/jesus-morales-souhail)  
**Date:** July 2026  
**Status:** Formal mapping with explicit **yes / analogy only / no** (hygiene / pedagogy). **Not a DESI claim.** Not peer reviewed.  
**Related:** `no-go-superoscillation-tesseract.md`, `car-drone-pupil-newton-einstein.md`, `stochastic-dark-energy-desi-dr2.md`, `papers/EXPLORATORY_BOUNDARY.md`  
**Repo $\theta$, $\sigma$:** OU kernel $C(\Delta x)=\sigma_X^2 e^{-\theta\lvert \Delta x \rvert}$; MLE drives $\sigma_X\to 0$, $\theta$ poorly constrained near boundary

---

## What this note is

I answer three questions that sit at the boundary of optics and the stochastic DE repo. The goal is **analogy discipline**: what transfers, what does not.  
This file is method hygiene. It is **not** part of the peer-review abstract of the main manuscript. Please do **not** cite it as a DESI result.

---

## Q1. Can we derive a “quantum mirror” equation that replaces Descartes for a single photon?

### What is legitimate

In **classical paraxial wave optics**, a thin lens (or quadratic-phase “mirror”) multiplies the field by

$$
M(\mathbf{r}_\perp)=A(\mathbf{r}_\perp) \exp\!\left(-i\frac{k}{2f} \lvert \mathbf{r}_\perp \rvert^2\right).
$$

Propagating with the Fresnel kernel from object plane $z=-s_o$ to image plane $z=s_i$, the quadratic phase at the pupil collects the factor

$$
\Delta R^{-1}\equiv \frac{1}{s_o}+\frac{1}{s_i}-\frac{1}{f}.
$$

When $\Delta R^{-1}=0$, the pupil phase is flat (to quadratic order) and one recovers the usual **imaging condition of Descartes / thin-lens formula**. The field at the image plane is then essentially the Fourier transform of the aperture (Airy for a circular stop), not a Dirac delta.

For a **single photon in a pure spatial mode**, the detection probability is

$$
P(\mathbf{r},s_i)=\lvert \psi(\mathbf{r},s_i) \rvert^2,
$$

with $\psi$ the same mode function (Born). So I may write, schematically,

$$
P(\mathbf{r},s_i)=\Big|\mathcal{F}\!\left\{\psi_{\mathrm{in}} A \exp\!\Big(i\frac{k}{2} |\mathbf{r}_\perp\lvert ^2 \Delta R^{-1}\Big)\right\}\Big \rvert^2
$$

(up to known Fresnel prefactors and normalizations). This is **not** a new law of physics: it is Fresnel + quadratic lens + Born.

### What does *not* replace Descartes

| Claim | Verdict |
|-------|---------|
| Descartes fails for one photon | **No.** The imaging condition $\Delta R^{-1}=0$ still marks the plane of best focus of the **mode**. |
| Probability concentrates to a point when $\Delta R^{-1}=0$ | **No.** Focused mode still has finite Airy / diffraction size set by $A$ and $\lambda$. |
| “Perfect phase” forces detection at the geometric image point | **No.** Born uses the full $\lvert \psi \rvert^2$; side lobes remain. |

**Answer to Q1:** Yes to a **careful** single-photon *mode* statement:  
**Descartes remains the condition $\Delta R^{-1}=0$ for quadratic phase cancellation; probability is $\lvert \psi \rvert^2$ after Fresnel propagation.**  
No to replacing Descartes with a formula that gives a deterministic sub-diffraction hit for one photon.

---

## Q2. How does a tesseract projection become an ABCD matrix on an SLM?

### What is legitimate

- Ray optics in the paraxial regime: **ABCD** (or 4×4 for $x,\theta_x,y,\theta_y$) is standard and symplectic.
- An ideal thin lens is $A=D=1$, $B=0$, $C=-1/f$.
- A **phase-only SLM** implements approximately a pure phase $\Phi(x,y)$, i.e. a position-dependent optical path — locally like a thin element of type “$C$-only” (plus possible pixelation and $2\pi$ wraps).
- Any smooth quadratic $\Phi$ can be read as a combination of focus / astigmatism / tilt; higher-order holograms can approximate more general transfers **within the limits of a 2D phase screen**.

### What is **not** derived

The proposed map

$$
M_{\mathrm{tesseract}}=\begin{pmatrix}A&B\\C&D\end{pmatrix} \quad\text{with}\quad A=\cos\phi_w I, B\propto d_w\sin\phi_w R(\theta), \ldots
$$

and

$$
\Phi_{\mathrm{SLM}}\propto (x^2\cos\phi_w+y^2\sin\phi_w)+\gamma_4 \mathrm{atan2}(y,x) \ln(r^2/r_0^2)
$$

**asserts** that:

1. SO(4) / Coxeter $B_4$ / “8 cubes” define $\phi_w,d_w,f_w,\gamma_4$;  
2. those parameters equal a physical 4D geometry of light;  
3. this exceeds ordinary SLM + pulse shaping.

None of (1)–(3) is derived from Maxwell + a physical tesseract cavity in the text. A 4D hypercube has 8 cubic cells as a **mathematical boundary**; that does **not** automatically supply eight independent optical channels beyond a 2D phase mask.

A **real** SLM programme would be:

1. Choose a target transfer (e.g. focus $f$, astigmatism, OAM charge $\ell$).  
2. Compute $\Phi(x,y)$ from standard Fourier / Fresnel / GS algorithms.  
3. Measure PSF and $E_{\mathrm{core}}/E_{\mathrm{total}}$.  

No $B_4$ required. OAM $\propto \ell\phi$ is standard and **not** a tesseract theorem.

Full Maxwell → device derivation: see `papers/maxwell-device-to-M-derivation.md`.

**Answer to Q2:**  
**I implement quadratic (and higher) phases on an SLM as usual.**  
**I do *not* get a derived “tesseract ABCD” until someone builds Maxwell + device geometry → $M$.** Until then, I treat tesseract language as **undeclared physical power** (same pattern as $\ln 4 \equiv \omega_R$).

---

## Q3. What is the analogue of focal length $f$ in the stochastic DE repo?

### What $\theta$ actually is in the repository

Ornstein–Uhlenbeck (Axiom A3 / scripts):

$$
\mathrm{d}X = -\theta X \mathrm{d}x + \sigma \mathrm{d}W_x, \qquad x=\ln a,
$$

with stationary variance $\mathrm{Var}(X)=\sigma^2/(2\theta)$ (when $\theta>0$) and correlation

$$
\mathrm{Cov}\big(X(x_i),X(x_j)\big) =\frac{\sigma^2}{2\theta} e^{-\theta\lvert x_i-x_j \rvert}.
$$

Physical reading in the paper: $\theta$ is an **effective mean-reversion / damping rate** (related to Hubble friction $3H(1+w)$ in the continuity story), **not** a geometric focal length.

Empirically (DESI DR2 BAO residual MLE): $\sigma_X\to 0$; $\theta$ is not a well-measured “focus of the cosmos.”  
(That empirical fact belongs to the main pipeline; this *analogy* note is still not a DESI claim.)

### Honest analogy table (formal, not identity)

| Property | Thin lens $f$ | OU $\tau=1/\theta$ |
|----------|-----------------|----------------------|
| Role | Sets **spatial** quadratic phase curvature of a wavefront | Sets **correlation length in $x=\ln a$** (memory / return rate) |
| Limit $\to\infty$ | Flat phase (no lens) | No mean reversion (Wiener / undamped) |
| Limit $\to 0$ | Infinite optical power (formal) | Instant snap to mean (if $\sigma$ fixed, variance formula singular) |
| Observable in repo | — | Enters residual **covariance between BAO bins**, not a PSF on the sky |

I may say, **as a metaphor only**:

$$
\tau \equiv \frac{1}{\theta} \quad\text{“plays a similar structural role to $f$”}
$$

in the sense that both are **single parameters controlling how strongly a quadratic (lens) or linear (OU drift) restoring mechanism curves a profile** (phase vs correlation).

I **must not** say:

$$
f_{\mathrm{cosmo}}\equiv\frac{1}{\theta} \quad\text{is the cosmological focal length of spacetime}
$$

without a derived map $\theta \leftrightarrow$ metric curvature / wavefront radius of the universe. That map is **not** in the repo.

### Desqueezing side note

In open-system notes, $t_{1/2}=\ln 2/\gamma$ with $\gamma\leftrightarrow\theta H_0$ is a **relaxation time**, again analogous to a rate, not to an optical $f$ in metres.

**Answer to Q3:**  
The closest **repo parameter** to a “restoring strength” is **$\theta$** (or $\tau=1/\theta$).  
I call it a **formal analogy** to $1/f$, not a replacement of Descartes by cosmology.

---

## Unified verdict (“¿vas por aquí?”)

| Path | Go / stop |
|------|-----------|
| Fresnel + quadratic phase + Born as “quantum thin lens” | **Go** (standard, careful wording) |
| SLM phase masks for focus / astigmatism / OAM | **Go** (engineering) |
| Tesseract / $B_4$ as superior ABCD engine | **Stop** (not derived) |
| $\tau=1/\theta$ as poetic “focal length of OU” | **Go only as analogy**, never as optical identity |
| Using that analogy to claim DE “focuses” probability on the sky | **Stop** (wrong operator / undeclared power) |

**Direction for this project:** keep **scale + operator** discipline (BAO residuals, model kills, Option 0 slip with amplitude honesty). Do **not** merge tesseract optics into the DESI/SDiff papers.

---

## Minimal equations to keep (if writing a pedagogy note)

**Optics (OK):**

$$
\frac{1}{s_o}+\frac{1}{s_i}-\frac{1}{f}=0,\qquad P=\lvert \psi \rvert^2.
$$

**OU (OK, repo):**

$$
\mathrm{d}X=-\theta X \mathrm{d}x+\sigma \mathrm{d}W,\qquad C_{ij}\propto e^{-\theta\lvert x_i-x_j \rvert}.
$$

**Bridge (analogy only):**

$$
f \longleftrightarrow \frac{1}{\theta} \quad\text{as “restoring-parameter” roles, not units-equivalent physics.}
$$

---

*End of hygiene note. Not a DESI claim.*
