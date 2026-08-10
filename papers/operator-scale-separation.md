# Operator and scale separation (claim hygiene)

Jesús Morales Souhail · [github.com/jesus-morales-souhail](https://github.com/jesus-morales-souhail) · August 2026

**Purpose.** Keep three research threads that are often conflated in conversation **strictly separate**. This note is hygiene and non-claims, not a new detection and not a derivation of \(\ell_\ast\).

**Sister theory map:** [measurable-stochastic-vacuum `NEXT_QUESTION.md`](https://github.com/jesus-morales-souhail/measurable-stochastic-vacuum/blob/main/papers/r1_kernel/NEXT_QUESTION.md) · three-gate lock `results/r1_three_gate/`.

---

## 0. Three threads (do not fuse)

| Thread | What it is | Home |
|:-------|:-----------|:-----|
| **A. Residual DE (this repo)** | Phenomenological \(X\equiv\delta\Omega_\Lambda\) on a fixed smooth background; OU/QNM and rank-1 tachyonic residual kernels on public DESI DR2 BAO \(\alpha\) | `stochastic-dark-energy-ou` |
| **B. Counting cell / R1 (theory sister)** | \(N_{\mathrm{eff}}=(L/\ell_\ast)^d\), \(\sigma_{0,\mathrm{eff}}=(\ell_\ast/L)^{d/2}\); three-gate compatibility if \(\ell_\ast=R_{\mathrm{nl}}\) | `measurable-stochastic-vacuum` |
| **C. Strong-field light (not this programme)** | VLBI / EHT: photon sphere, Schwarzschild shadow, mm-wave interferometry | Outside the claim set |

Shared attitude: measure before inventing amplifiers; do not mix operators or scales.  
Shared physics slogan (“GR affects light”) is **not** a map between kernels.

---

## 1. Thread A — residual DE (DESI BAO)

### 1.1 What the pipeline actually implements

OU in logarithmic scale factor \(x=\ln a\):

$$
\mathrm{d}X=-\theta X\,\mathrm{d}x+\sigma\,\mathrm{d}W_x,\qquad \mathrm{Var}(X)=\frac{\sigma^{2}}{2\theta}\quad(\theta>0).
$$

Likelihood kernel (code: `scripts/desi_dr2_data.py` → `add_ou_kernel`):

$$
C_{ij}=C^{\mathrm{meas}}_{ij}+S_{i}S_{j}\,\sigma_X^{2}\,e^{-\theta\lvert\Delta x_{ij}\rvert}.
$$

Here \(\sigma_X\) is the **residual amplitude** that multiplies the sensitivity kernel \(S\). It is the quantity profiled against DESI. It is **not** the baryon density \(\omega_b\), and it is **not** the Sorkin seed \(\sigma_0\sim 10^{-61}\) without a map.

### 1.2 Published residual ceiling (claim C2)

$$
\sigma_X < 2.5\times 10^{-2}\quad\text{(95\% CL, free \(\theta\), 7 \(\alpha\) bins, full meas. cov)}.
$$

Best fit: \(\sigma_X\to 0\). Former working target \(1.5\times 10^{-4}\) is **not** this profile.

### 1.3 Rank-1 tachyonic growth (separate operator)

Coherent growing mode uses \(C=C^{\mathrm{meas}}+\mathbf{u}\mathbf{u}^{\mathsf{T}}\), not the stationary OU form. Active growth excluded; preferred \(t_c\to\infty\). See `quantum-fluid-instabilities-desi-dr2.md`.

### 1.4 Symbol hygiene (do not mix)

| Symbol | Meaning in this corpus | Not equal to |
|:-------|:-----------------------|:-------------|
| \(\sigma_X\) | Residual amplitude in \(C_{\mathrm{OU}}\) | \(\omega_b=\Omega_b h^2\approx 0.022\) |
| \(2.5\times 10^{-2}\) | DESI free-\(\theta\) residual ceiling | Baryon budget / BBN number |
| \(\sigma\) (SDE diffusion) | Langevin noise strength | \(\sigma_X\) unless \(\sigma_X=\sigma/\sqrt{2\theta}\) is stated |
| \(\theta\) | Mean reversion of \(X\) in \(\ln a\) | \(k\)-dependent BAO wiggle damping, Jeans scale, or \(t_c\) |
| \(\sigma_0\sim 10^{-61}\) | Motivational holographic seed | Observed residual amplitude |

**Illegal move:** identify \(\sigma_{\mathrm{free}}\) from counting with the SDE diffusion coefficient and then divide again by \(\sqrt{2\theta}\). Counting \(\sigma_{\mathrm{free}}\) is already compared to \(\sigma_X\) as an amplitude of the same kind (three-gate G1).

---

## 2. Thread B — counting cell and three-gate (sister repo)

Under the counting hypothesis with \(d=3\), \(L=L_H\), and the **proposal** \(\ell_\ast=R_{\mathrm{nl}}\):

$$
\sigma_{\mathrm{free}}=\Bigl(\frac{R_{\mathrm{nl}}}{L_H}\Bigr)^{3/2}\approx 8.5\times 10^{-5}
$$

(\(R_{\mathrm{nl}}\approx 8.61\,\mathrm{Mpc}\), \(L_H\approx 4448\,\mathrm{Mpc}\); `r1_three_gate_lock.py`).

| Gate | Statement | Status |
|:-----|:----------|:-------|
| G1 | \(\sigma_{\mathrm{free}} < 2.5\times 10^{-2}\) (headroom \(\approx 294\times\)) | PASS |
| G2 | CF4 \(r_e\) in sandwich band around \(R_{\mathrm{nl}}\) (primary \(L=20\)) | PASS |
| G3 | Block \(\eta=\sigma_v/(H_0 L)=\mathcal{O}(1)\) near \(R_{\mathrm{nl}}\) | PASS |

**What G1 is:** compatibility — if \(\ell_\ast=R_{\mathrm{nl}}\), the free residual amplitude sits well below the DESI residual ceiling.  
**What G1 is not:** a derivation that \(\ell_\ast\) must equal \(R_{\mathrm{nl}}\); a soft amplifier; a new OU prediction from \(\sqrt{\sigma^2/2\theta}\) with \(\theta\sim 1\).

**Still open (R1):** a principle that *forces* a mesoscopic counting cell for the DE residual sector. See sister `r1-open-kernel.md` and `r1-principle-nonlinear-matter.md`.

Ordered next steps (theory side, already listed in `NEXT_QUESTION.md`): path-RMS / slip floors with \(\ell_\ast=R_{\mathrm{nl}}\) **fixed** (no retune); then principle-hunting only after that.

---

## 3. Thread C — VLBI, Schwarzschild, photon sphere (outside claim set)

| Object | Definition | Role here |
|:-------|:-----------|:----------|
| \(R_s=2GM/c^2\) | Event horizon (Schwarzschild) | Strong-field GR |
| \(1.5\,R_s=3GM/c^2\) | Unstable **photon** circular orbit (photon sphere) | Light, not massive particles |
| \(3\,R_s=6GM/c^2\) | Innermost stable circular orbit (**massive** particles) | Not the photon sphere |
| VLBI / EHT | Earth-scale radio baselines; mm image of shadow / photon ring | Different instrument, different scale |

Pedagogical parallels only: instrument resolution vs predicted feature size; unstable orbits vs growth rates; “do not mix emission model with geometry.”

**Non-claims:**

- No map \(\ell_\ast \leftrightarrow R_s\) or residual OU \(\leftrightarrow\) photon ring.
- No paper claim that baryonic feedback “is” the photon sphere.
- No use of EHT data in the DESI residual likelihood.

If a strong-field note is ever written, it belongs in exploratory notes or a separate GR project — not under residual BAO claims.

---

## 4. Baryons (related science, different equation)

| Quantity | Approx. value | Role |
|:---------|:--------------|:-----|
| \(\omega_b=\Omega_b h^2\) | \(\approx 0.022\) | Total ordinary matter (BBN/CMB) |
| DESI residual ceiling \(\sigma_X\) | \(2.5\times 10^{-2}\) | Stochastic residual amplitude on BAO \(\alpha\) |
| BAO itself | scale \(\sim 150\,\mathrm{Mpc}\) | Early plasma sound horizon; standard DESI product |

Baryonic feedback, WHIM, cluster gas, and kSZ can **bias** large-scale structure and distances. They are real systematics and real astrophysics. They are **not** the same as \(\sigma_X\), and the numerical coincidence “0.025 vs 0.022” is **not** a physical identification.

---

## 5. What is publishable here vs what is not

| Content | Status |
|:--------|:-------|
| DESI residual null + ceiling C2; tachyonic rank-1 exclusion; amplification walls | Claims (this repo) |
| Three-gate compatibility under \(\ell_\ast=R_{\mathrm{nl}}\) | Sister result; not a derivation of \(\ell_\ast\) |
| This separation / symbol hygiene note | Publishable as programme hygiene |
| “RG fixed point forces \(\ell_\ast=R_{\mathrm{nl}}\)” without microphysics | **Not** publishable as a result of this corpus |
| “\(\theta\sim 1\) closes a new gap beyond G1” | Restates G1; do not advertise as a new theorem |
| Unification residual DE \(\leftrightarrow\) EHT photon sphere | **Not** a claim of this programme |

---

## 6. One-line rules

1. \(\sigma_X\) residual \(\neq\) \(\omega_b\) baryons \(\neq\) \(\sigma_0\) Sorkin seed.  
2. OU \(\theta\) (in \(\ln a\)) \(\neq\) tachyonic \(t_c\) \(\neq\) BAO \(k\)-damping.  
3. G1 = sub-ceiling counting under a **proposal**; R1 = missing **principle**.  
4. VLBI / \(1.5\,R_s\) = strong-field light; keep it out of the residual claim set.

---

## Reproducibility (threads A–B only)

```bash
# Thread A
python scripts/profile_sigma_x_desi.py
python scripts/tachyonic_rank1_mle.py

# Thread B (sister repo)
python scripts/r1/r1_three_gate_lock.py
python scripts/r1/r1_cf4_re_jackknife.py
```
