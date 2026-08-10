# Scope and mix-ups I keep having to unpick

Jesús Morales Souhail · August 2026  
[github.com/jesus-morales-souhail](https://github.com/jesus-morales-souhail)

This is a short note for myself (and anyone reading the public repos) about what belongs where. It is not a new measurement and not a derivation of \(\ell_\ast\).

---

## Three different things

**1. This repository.** Late-time residual on top of a smooth dark-energy background. I model \(X \equiv \delta\Omega_\Lambda\) as an Ornstein–Uhlenbeck process in \(x = \ln a\), fit public DESI DR2 isotropic \(\alpha\), and exclude a coherent rank-1 growing mode. That is the claim set.

The process and the kernel used in the code are


$$
\mathrm{d}X=-\theta X\,\mathrm{d}x+\sigma\,\mathrm{d}W_x,\qquad \mathrm{Var}(X)=\frac{\sigma^{2}}{2\theta}\quad(\theta>0).
$$


$$
C_{ij}=C^{\mathrm{meas}}_{ij}+S_{i}S_{j}\,\sigma_X^{2}\,e^{-\theta\lvert\Delta x_{ij}\rvert}.
$$


The number I report is the residual amplitude \(\sigma_X\) in that second equation (see `scripts/desi_dr2_data.py`, `add_ou_kernel`). Free-\(\theta\) profile, 7 bins, full measurement covariance:


$$
\sigma_X < 2.5\times 10^{-2}\quad\text{(95\% CL)}.
$$


Best fit still drives \(\sigma_X\to 0\). The old \(1.5\times 10^{-4}\) was a working target, not this profile.

The tachyonic rank-1 test is a *different* covariance (\(C^{\mathrm{meas}}+\mathbf{u}\mathbf{u}^{\mathsf{T}}\)). Preferred limit is no growth. Details in `quantum-fluid-instabilities-desi-dr2.md`.

**2. Sister theory repository.** Counting of a residual cell: \(N_{\mathrm{eff}}=(L/\ell_\ast)^d\) and \(\sigma_{0,\mathrm{eff}}=(\ell_\ast/L)^{d/2}\). If one *proposes* \(\ell_\ast = R_{\mathrm{nl}}\) with \(d=3\), the free amplitude is about \(8.5\times 10^{-5}\), well under the DESI ceiling above (three-gate G1, headroom ~294). That is compatibility, not a proof that the cell must be \(R_{\mathrm{nl}}\). What still missing is a principle that fixes \(\ell_\ast\). Map: [NEXT_QUESTION.md](https://github.com/jesus-morales-souhail/measurable-stochastic-vacuum/blob/main/papers/r1_kernel/NEXT_QUESTION.md).

**3. Black holes and VLBI.** Photon sphere at \(1.5\,R_s\), ISCO for massive particles at \(3\,R_s\), EHT images of the shadow — strong-field light, millimetre interferometry. Same general relativity in a different regime. I do not map \(\ell_\ast\) or \(\sigma_X\) onto \(R_s\), and I do not put EHT data in the BAO residual likelihood. If I ever write that up, it will not sit under residual-DE claims.

---

## Numbers that look alike but are not the same

| Number | What it is |
|:-------|:-----------|
| \(\sigma_X < 2.5\times 10^{-2}\) | Residual amplitude ceiling on DESI BAO \(\alpha\) |
| \(\omega_b \approx 0.022\) | Physical baryon density from BBN/CMB |
| \(\sigma_{\mathrm{free}} \approx 8.5\times 10^{-5}\) | Counting amplitude if \(\ell_\ast=R_{\mathrm{nl}}\), \(d=3\) |
| \(\sigma_0 \sim 10^{-61}\) | Motivational Sorkin / holographic seed |

\(\sigma_X\) in the likelihood is the residual amplitude. It is not \(\omega_b\). It is not the Sorkin seed unless someone builds a map. In the OU SDE, the diffusion coefficient \(\sigma\) and the reported \(\sigma_X\) only match after \(\sigma_X = \sigma/\sqrt{2\theta}\); the three-gate number \(\sigma_{\mathrm{free}}\) is already compared to \(\sigma_X\) as an *amplitude*, so I do not divide it by \(\sqrt{2\theta}\) again.

\(\theta\) in the OU is mean reversion of \(X\) in \(\ln a\). It is not the BAO wiggle damping scale, not a Jeans wavenumber, and not the tachyonic \(t_c\).

Baryonic feedback, WHIM, clusters, kSZ matter for large-scale structure and systematics. They are real work. They are not the residual \(\sigma_X\) object.

---

## What I will not claim from a chat

- That G1 “becomes a theorem” by assuming \(\theta\sim 1\). G1 is already the sub-ceiling comparison.
- That an RG story without microphysics forces \(\ell_\ast = R_{\mathrm{nl}}\).
- That residual dark energy and the EHT photon ring are the same project.

Next real calculation on the theory side remains path-RMS / slip with \(\ell_\ast=R_{\mathrm{nl}}\) held fixed — no retune to DESI.

```bash
python scripts/profile_sigma_x_desi.py
python scripts/tachyonic_rank1_mle.py
# sister: python scripts/r1/r1_three_gate_lock.py
```
