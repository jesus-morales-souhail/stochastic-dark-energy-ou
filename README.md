# stochastic-dark-energy-ou

Jesús Morales Souhail
[ORCID 0009-0000-7637-1818](https://orcid.org/0009-0000-7637-1818) · jmskjym@gmail.com
[github.com/jesus-morales-souhail](https://github.com/jesus-morales-souhail)

July 2026 · independent research

I use this repo for the **data** side: public DESI DR2 BAO, residual models (OU / QNM), what is excluded, and how large the soft-amplification gap is if you start from a pure Sorkin seed.

---

## Start here

| Document | What it is |
|:---------|:-----------|
| [`START_HERE.md`](START_HERE.md) | short map |
| [`manuscript/PREPRINT.md`](manuscript/PREPRINT.md) | main write-up |
| [`manuscript/CLAIMS.md`](manuscript/CLAIMS.md) | claims vs non-claims |
| [`papers/HONEST_HEADLINES.md`](papers/HONEST_HEADLINES.md) | allowed vs forbidden wording |
| [`papers/EXPLORATORY_BOUNDARY.md`](papers/EXPLORATORY_BOUNDARY.md) | what is not part of the paper |

Theory: [measurable-stochastic-vacuum](https://github.com/jesus-morales-souhail/measurable-stochastic-vacuum)

Exploratory notes (not cosmology claims): [stochastic-de-exploratory-notes](https://github.com/jesus-morales-souhail/stochastic-de-exploratory-notes)

---

## What the project is asking

1. **Data.** Do DESI DR2 BAO residual statistics require a stationary stochastic (OU/QNM) piece on top of a smooth background?
2. **Geometry.** Can SDiff / unimodular structure explain why late-time vacuum residuals stay small?
3. **Method.** Right scale and right operator for each claim — no undeclared free parameters.

### Data (BAO residuals)

| Result | Detail |
|:-------|:-------|
| CPL background (BAO-only) | $w_0\approx -0.99$, $w_a\approx -0.02$ (near $\Lambda$CDM) |
| OU / QNM residuals | amplitude driven to the floor when the background is free |
| Working upper limit | $\sigma_X < 1.5\times 10^{-4}$ (95% CL) under the OU kernel I use |
| Tachyonic quantum fluid | coherent growing mode excluded ($\Delta\ln\mathcal{L}\approx -11.35$ vs $\Lambda$CDM) |

Notes: `papers/stochastic-dark-energy-desi-dr2.md`, `papers/quantum-fluid-instabilities-desi-dr2.md`, `papers/resume.txt`.

### Smoothness (SDiff)

Unimodular / SDiff projects out isotropic vacuum stress $T_{\mu\nu}\propto g_{\mu\nu}$. That is a candidate reason residuals are small — not a particle detection.

### Amplification gap

A pure Sorkin seed $\sigma_0\sim 10^{-61}$ needs roughly $G\sim 10^{56}$ to hit Euclid-scale residuals $\sim 10^{-5}$, and $\sim 10^{57}$ to the DESI ceiling $1.5\times 10^{-4}$. Soft desqueezing gives $\mathcal{O}(10)$ at $r\sim 1.5$, not $10^{56}$. See `papers/amplification-gap.md` and `papers/HONEST_HEADLINES.md`.

### Slip (Option 0)

SDiff does not cancel shear. Gravitational slip is the right *kind* of operator, but at $\sigma_X\sim 10^{-4}$ the signal is $\sim 10^{-4}$, far below Maus-class errors $\mathcal{O}(0.1)$. No shortcut.

### Wrong paths

Lab optics, tesseract slogans, and undeclared free powers live only in the exploratory repo — not under `papers/` here.

---

## Layout

```
manuscript/ preprint and claims
papers/ technical notes
scripts/ analysis code
data/ local DESI products (public provenance)
results/ numerical outputs
notes/ working notes
```

See [`STRUCTURE.md`](STRUCTURE.md) if present for more detail.
