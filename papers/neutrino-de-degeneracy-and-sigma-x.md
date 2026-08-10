# Neutrino mass, dark energy, and our $\sigma_X$ null

Jesús Morales Souhail
[github.com/jesus-morales-souhail](https://github.com/jesus-morales-souhail)
August 2026

Planning note. No new fit. Not a DESI result of mine.

See also: `manuscript/PREPRINT.md`, `manuscript/CLAIMS.md`, `scripts/joint_w0wa_sigma_desi.py`.

---

DESI’s own neutrino analysis (DR2 BAO + CMB) reports roughly:

| Setup | $\sum m_\nu$ (95% CL) |
|-------|------------------------|
| $\Lambda$CDM | $< 0.064\,\mathrm{eV}$ |
| $w_0 w_a$ CDM (dynamical DE) | $< 0.16\,\mathrm{eV}$ |

Lab oscillations put a floor near $\sum m_\nu \gtrsim 0.059\,\mathrm{eV}$ for normal ordering. In $\Lambda$CDM the cosmological ceiling sits almost on that floor; inverted ordering is under pressure. Open the dark-energy equation of state and the ceiling rises by about a factor of 2.5. That is a published degeneracy between $ (w_0,w_a)$ and $\sum m_\nu$, not something I invented.

---

My joint on public DESI DR2 BAO (`joint_w0wa_sigma_desi.py`) varies

$$
\{w_0,\,w_a,\,\theta,\,\sigma_X\}
$$

with an OU residual kernel on $\alpha_{\mathrm{obs}}-\alpha_{\mathrm{pred}}(w_0,w_a)$. It does **not** vary $\sum m_\nu$. The MLE I already have is $ (w_0,w_a)\approx(-0.99,-0.016)$ and $\sigma_X\to 0$.

So:

| Parameter | Meaning | Units | In the joint? |
|-----------|---------|-------|----------------|
| $\sigma_X$ | residual amplitude in the BAO covariance (OU) | dimensionless | yes — driven to zero |
| $w_0,w_a$ | CPL dark energy | dimensionless | yes |
| $\sum m_\nu$ | total neutrino mass | eV | **no** |

I have looked at how $\sigma_X$ trades with a free CPL background. I have **not** re-derived DESI’s $\sum m_\nu$–$ (w_0,w_a)$ result. BAO alone is a weak handle on neutrino mass anyway: distances, not the growth suppression that massive neutrinos leave. Dumping $\sum m_\nu$ into this likelihood without CMB or growth would mostly recover the prior. That is easy to dress up as a “bound.” I will not.

---

**What I will say:**

Under a free CPL background, DESI DR2 BAO (this pipeline) does not need a stationary OU residual: $\sigma_X\to 0$. DESI+CMB under $\Lambda$CDM put a tight ceiling on $\sum m_\nu$; that ceiling loosens when dark energy is dynamical. Neutrino mass is fixed out of my likelihood, so the $\sigma_X$ null is for the BAO-only, fixed-neutrino setup I actually ran.

**What I will not say:**

- that my $\sigma_X$ null constrains $\sum m_\nu$, or the reverse
- that $\sigma_X$ and $\sum m_\nu$ are the same operator (they are not; both just compete with the same smooth DE background)
- that this note is a likelihood result (nothing new was fitted)

---

If I extend the fit later: add $\sum m_\nu$ only with a public CMB or growth prior; keep $\sum m_\nu$ and $\sigma_X$ as separate parameters. Do not fold one into the other’s kernel.

