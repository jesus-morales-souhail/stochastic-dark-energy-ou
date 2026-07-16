# Magnetar extreme conditions and open-system vacuum recovery

**Author:** Jesús Morales Souhail  
**Date:** July 2026  
**ORCID:** [0009-0000-7637-1818](https://orcid.org/0009-0000-7637-1818)  
**Code:** `scripts/magnetar_vacuum_response.py`  
**Outputs:** `results/magnetar_vacuum_response/`  
**Status:** Research note — not peer reviewed  

---

## Abstract

Magnetars host magnetic fields of order \(10^{14}\)–\(10^{15}\,\mathrm{G}\) and surface gravity comparable to other neutron stars. Their giant flares are among the most intense electromagnetic transients in the Galaxy. Independently, this repository models late-time vacuum fluctuations with an Ornstein–Uhlenbeck kernel and an open-system desqueezing law \(t_{1/2}=\ln 2/\gamma_{\mathrm{eff}}\). Here we place both frameworks side by side: we compute gravitational time dilation at the neutron-star surface, a schematic QED vacuum birefringence \(\Delta n\sim(\alpha/4\pi)(B/B_{\mathrm{QED}})^2\), and desqueezing recovery curves with half-lives matched to magnetar timescales **as an analogy only**. The result is a quantitative scale separation: magnetar recovery rates, if written as \(\gamma\sim\ln 2/t_{\mathrm{flare}}\), are \(\sim 10^{-3}\)–\(10\,\mathrm{s}^{-1}\), while cosmological mean-reversion \(\gamma\sim\theta H_0\) is \(\sim 10^{-21}\)–\(10^{-18}\,\mathrm{s}^{-1}\). The shared idea is a dynamical vacuum; the physics is not the same, and no DESI–magnetar parameter matching is performed.

---

## 1. Motivation

Both magnetar magnetospheres and cosmological vacuum models raise the same qualitative question: how does the vacuum respond when conditions are extreme? In magnetars the extreme agent is a supercritical magnetic field and dense matter under strong gravity. In this repository the extreme agent is a controlled squeezed state (laboratory sector) or residual dark-energy fluctuations (BAO sector).

The comparison is useful only if it stays honest:

- **Do** compare dimensionless response measures and timescales.  
- **Do not** identify flare frequencies with BAO \(\omega_R\), or flare energy with \(\sigma_X\).

---

## 2. Magnetar anchors

| Quantity | Scale used | Role |
|----------|------------|------|
| \(B\) | \(10^{14}\)–\(10^{15}\,\mathrm{G}\) | Field strength |
| \(B_{\mathrm{QED}}\) | \(\simeq 4.41\times 10^{13}\,\mathrm{G}\) | Schwinger critical field (SI formula) |
| \(R\), \(M\) | \(12\,\mathrm{km}\), \(1.4\,M_\odot\) | Compactness / redshift |
| \(P\) | \(\sim 7.56\,\mathrm{s}\) | Spin (e.g. SGR 1806−20 order) |
| Flare spike / tail | \(\sim 0.2\,\mathrm{s}\) / \(\sim 400\,\mathrm{s}\) | Electromagnetic timescales |

Gravity is treated with a Schwarzschild exterior: at the surface,
\(\alpha=\sqrt{1-R_s/R}\) so that a proper interval \(\mathrm{d}\tau\) maps to distant time \(\mathrm{d}t_\infty=\mathrm{d}\tau/\alpha\). Magnetars remain outside their Schwarzschild radius; they have a surface, not an event horizon.

---

## 3. Vacuum response in strong \(B\)

A schematic low-energy Euler–Heisenberg measure of vacuum modification is

$$
\Delta n \sim \frac{\alpha}{4\pi}\left(\frac{B}{B_{\mathrm{QED}}}\right)^2.
$$

For \(B\sim 10^{14}\)–\(10^{15}\,\mathrm{G}\), \(B/B_{\mathrm{QED}}\sim\mathrm{few}\)–\(20\), so \(\Delta n\) is no longer negligible compared with laboratory QED tests. This is **local QED**, not a cosmological gauge argument. It does show that the vacuum is not a passive medium under extreme conditions.

Propagation near the star is further shaped by:

- magnetic channeling of radiation;  
- possible photon splitting and pair processes in supercritical fields;  
- gravitational light bending and redshift for a distant observer.

A full radiative-transfer simulation is outside the scope of this note; the script reports order-of-magnitude \(\Delta n(B)\) and time dilation only.

---

## 4. Open-system desqueezing (repository law)

From the Lindblad analysis in this project (QuTiP sector; analytic limit used here):

$$
\lvert\langle a^2\rangle\rvert(t)
=
\lvert\langle a^2\rangle\rvert_0\,
e^{-\gamma(1+2n_{\mathrm{th}})t},
\qquad
t_{1/2}=\frac{\ln 2}{\gamma(1+2n_{\mathrm{th}})}.
$$

If one asks, purely as a **dictionary exercise**, what \(\gamma\) would make \(t_{1/2}\) equal a magnetar timescale \(T\), then

$$
\gamma=\frac{\ln 2}{T}.
$$

That maps spikes, spin periods, and tails onto open-system rates. It does **not** mean that a giant flare is a squeezed cosmological vacuum.

---

## 5. Cosmological dictionary in this repository

The OU continuity structure used for DESI BAO is

$$
\Gamma_{\mathrm{phys}}(z)=\theta\,H(z),
\qquad
\gamma \leftrightarrow \theta H_0
$$

at the present epoch. With \(H_0\simeq 67.4\,\mathrm{km\,s^{-1}\,Mpc^{-1}}\):

| \(\theta\) | \(\gamma\sim\theta H_0\) | \(t_{1/2}=\ln 2/\gamma\) |
|------------|--------------------------|---------------------------|
| \(1\) | \(\sim 2\times 10^{-18}\,\mathrm{s}^{-1}\) | \(\sim 10\,\mathrm{Gyr}\) |
| \(10^{-3}\) (numerical floor scale) | \(\sim 2\times 10^{-21}\,\mathrm{s}^{-1}\) | \(\sim 10^{4}\,\mathrm{Gyr}\) |

These numbers come from the repository’s cosmological mapping, not from magnetars.

---

## 6. Results of the comparison

Running `scripts/magnetar_vacuum_response.py` produces:

1. **Gravity:** for \(M=1.4\,M_\odot\), \(R=12\,\mathrm{km}\), \(\alpha\simeq 0.76\), so a \(0.2\,\mathrm{s}\) local spike maps to \(\sim 0.26\,\mathrm{s}\) at infinity (order-unity dilation, not horizon-like).  
2. **QED:** \(\Delta n\) rises steeply through the magnetar \(B\) range as \((B/B_{\mathrm{QED}})^2\).  
3. **Desqueezing analogy:** \(\gamma\sim 3\,\mathrm{s}^{-1}\) for a \(0.2\,\mathrm{s}\) half-life; \(\gamma\sim 0.09\,\mathrm{s}^{-1}\) for \(7.56\,\mathrm{s}\); \(\gamma\sim 2\times 10^{-3}\,\mathrm{s}^{-1}\) for \(400\,\mathrm{s}\).  
4. **Scale gap:** magnetar-associated \(\gamma\) exceeds cosmological \(\theta H_0\) by roughly \(10^{18}\)–\(10^{22}\).

Figure: `figures/magnetar_vacuum_desqueezing.png`.

---

## 7. What this contributes to the repository

| Claim | Status |
|-------|--------|
| Vacuum can respond non-trivially under extreme conditions | Supported (QED + open systems) |
| Desqueezing half-life law is portable mathematics | Supported (same formula as QuTiP sector) |
| Magnetar flares measure DESI \(\sigma_X\) or SDiff | **Not claimed** |
| A single \(\gamma\) unifies magnetars and cosmology | **Ruled out by scale separation** |

The scientific use of magnetars for this project is therefore **contrast and calibration of ideas**, not a multi-messenger joint fit. Future work that remains non-numerological would need a first-principles open-system derivation of an effective \(\gamma(B)\) from QED (or magnetospheric kinetics), then ask whether any residual couples to infrared vacuum modes—without equating microphysical rates to \(\theta H_0\).

---

## 8. How to reproduce

```bash
cd stochastic-dark-energy-ou
python3 scripts/magnetar_vacuum_response.py
```

Requires only `numpy` and `matplotlib` (no QuTiP). Analytic desqueezing matches the half-life law already documented in `notes/desqueezing-relaxation-vacuum-fluctuations-note.md`.

---

## References (selection)

1. DESI Collaboration, arXiv:2503.14738 (2025).  
2. Magnetar and giant-flare literature (e.g. SGR 1806−20, spin \(\sim 7.5\,\mathrm{s}\), \(B\sim 10^{15}\,\mathrm{G}\)).  
3. Euler–Heisenberg effective QED; Schwinger critical field.  
4. This repository: desqueezing note; DESI OU/QNM analysis; Euclid protocol.
