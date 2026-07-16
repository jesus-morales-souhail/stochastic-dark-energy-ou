# Literature anchors: magnetars and strong-field QED

**Purpose.** Curated, citable numbers for any future analysis that needs magnetar / QED scales.  
**Rule.** These are **published astrophysical and QED values only**. They are **not** free parameters of the DESI BAO analysis and are **not** to be matched by hand to \(\sigma_X\), \(\omega_R\), or \(\ln 4\).

**Primary catalog:** [McGill Online Magnetar Catalog](https://www.physics.mcgill.ca/~pulsar/magnetar/main.html) (Olausen & Kaspi 2014 and updates).

---

## 1. Class properties (order of magnitude)

| Quantity | Value | Notes | Source |
|----------|-------|-------|--------|
| Surface \(B\) (magnetar range) | \(\sim 10^{13}\)–\(10^{15}\,\mathrm{G}\) (\(\sim 10^{9}\)–\(10^{11}\,\mathrm{T}\)) | Dipole estimate from spin-down; local multipoles can be larger | Wikipedia Magnetar; standard reviews |
| Ordinary NS dipole | \(\sim 10^{12}\,\mathrm{G}\) | Contrast | literature |
| Earth surface \(B\) | \(\sim 0.25\)–\(0.65\,\mathrm{G}\) | Contrast | geomagnetism |
| Radius | \(\sim 10\)–\(12\,\mathrm{km}\) (diameter \(\lesssim 20\,\mathrm{km}\)) | NS equation of state | standard |
| Mass | \(\sim 1.4\,M_\odot\) | Typical NS | standard |
| Spin period (magnetars) | \(\sim 2\)–\(12\,\mathrm{s}\) | Slow rotators vs radio pulsars | McGill catalog |
| Active lifetime | \(\sim 10^{4}\,\mathrm{yr}\) (order) | Field decay | reviews |
| Galactic sample | \(\sim 30\) confirmed / candidates | Growing | McGill catalog; OA-Roma notes |

---

## 2. Schwinger / QED critical field

| Quantity | Value | Formula / note | Source |
|----------|-------|----------------|--------|
| \(B_{\mathrm{QED}}\) | \(4.41\times 10^{9}\,\mathrm{T} = 4.41\times 10^{13}\,\mathrm{G}\) | \(B_{\mathrm{c}}=m_e^{2}c^{2}/(e\hbar)\) (SI) | Schwinger limit (Wikipedia); QED texts |
| \(E_{\mathrm{QED}}\) (electric) | \(\sim 1.3\times 10^{18}\,\mathrm{V\,m^{-1}}\) | Pair-production scale | Schwinger 1951 / standard QED |
| Magnetar vs \(B_{\mathrm{QED}}\) | \(B/B_{\mathrm{QED}}\sim\mathrm{few}\)–\(20\) for \(B\sim 10^{14}\)–\(10^{15}\,\mathrm{G}\) | Supercritical for QED vacuum effects | literature |

Vacuum birefringence, photon splitting, and related processes become relevant when \(B\) approaches or exceeds \(B_{\mathrm{QED}}\). That is **local QED**, not a cosmological DE observable.

---

## 3. Reference objects (giant-flare sources)

### SGR 1806−20

| Quantity | Value | Source |
|----------|-------|--------|
| Spin period \(P\) | \(7.54773(2)\,\mathrm{s}\) (catalog); \(\simeq 7.56\,\mathrm{s}\) in flare tail | McGill catalog; Hurley et al. / RHESSI analyses |
| Dipole \(B\) (spin-down order) | \(\gtrsim 10^{15}\,\mathrm{G}\) (often quoted \(\sim 2\times 10^{15}\,\mathrm{G}\)) | McGill; Wikipedia |
| Distance | \(\sim 8.7\)–\(15\,\mathrm{kpc}\) (uncertain; older popular accounts \(\sim 15\,\mathrm{kpc}\) / \(\sim 50{,}000\,\mathrm{ly}\)) | McGill; literature scatter |
| Giant flare date | 2004 Dec 27 | Palmer et al. 2005, *Nature* |
| Isotropic flare energy | \(\sim 2\times 10^{46}\,\mathrm{erg}\) (total isotropic; \(\sim 100\times\) previous giant flares) | Palmer et al. 2005 |
| Spike duration | \(\sim 0.2\,\mathrm{s}\) | Hurley et al.; Rea 2006 review notes |
| Pulsating tail | \(\sim 400\)–\(600\,\mathrm{s}\), pulsed at \(P\simeq 7.56\,\mathrm{s}\) | RHESSI / multi-satellite |
| Tail energy | \(\sim 5\times 10^{43}\,d_{15}^{2}\,\mathrm{erg}\) (order; comparable across giant flares) | Rea 2006 and refs |

### SGR 1900+14

| Quantity | Value | Source |
|----------|-------|--------|
| \(P\) | \(5.19987(7)\,\mathrm{s}\) | McGill catalog |
| Giant flare | 1998 Aug 27 | Hurley et al. 1999 |
| Isotropic energy | \(\sim 10^{44}\)–\(10^{45}\,\mathrm{erg}\) (order; smaller than 1806−20) | giant-flare literature |

### SGR 0526−66

| Quantity | Value | Source |
|----------|-------|--------|
| \(P\) | \(8.0544(2)\,\mathrm{s}\) | McGill catalog |
| Giant flare | 1979 Mar 5 | Mazets et al. |
| Location | LMC | catalog |

---

## 4. Energies and magnetic reservoir (orders of magnitude)

| Quantity | Scale | Comment |
|----------|-------|---------|
| Giant-flare isotropic energy | \(10^{44}\)–\(10^{46}\,\mathrm{erg}\) | Spike-dominated; beaming may reduce true energy |
| Magnetic energy in crust / multipoles | \(\sim 10^{46}\)–\(10^{47}\,\mathrm{erg}\) (model-dependent) | Enough to power a giant flare; see e.g. crust-field simulations |
| Persistent \(L_X\) | \(\sim 10^{31}\)–\(10^{35}\,\mathrm{erg\,s^{-1}}\) | Quiescent / outburst states |
| Short bursts | \(\sim 10^{40}\,\mathrm{erg}\) scale | Much weaker than giant flares |

**Dipole field energy (order-of-magnitude estimate only):**
\[
E_B \sim \frac{B^{2}}{8\pi}\,R^{3}
\sim 10^{47}\left(\frac{B}{10^{15}\,\mathrm{G}}\right)^{2}
\left(\frac{R}{12\,\mathrm{km}}\right)^{3}\,\mathrm{erg}
\]
(coefficient \(\mathcal{O}(1)\); not a precise interior integral).

---

## 5. Gravity vs magnetism (physical statement)

- Escape speed from a \(1.4\,M_\odot\), \(12\,\mathrm{km}\) NS is a large fraction of \(c\); the star is held by gravity, not by \(B\).
- Magnetic stresses dominate the **crustal mechanics** and power bursts; they do not replace the gravitational binding of the star.
- Magnetars have a surface; they are not event horizons. Photons and \(B\)-field structure remain exterior observables.

---

## 6. Use inside this repository

| Allowed | Not allowed |
|---------|-------------|
| Use these numbers as **external physical context** | Fit DESI BAO to magnetar \(P\) or \(B\) |
| Compare **timescales** to open-system \(t_{1/2}=\ln 2/\gamma\) with explicit analogy label | Identify flare energy with \(\sigma_X\) |
| Cite McGill / Palmer / QED for methods sections | Invent “total EM field of the universe” from one magnetar |

Script that consumes a subset of these anchors: `scripts/magnetar_vacuum_response.py`.  
Machine-readable copy: `literature/magnetar_literature_anchors.json`.

---

## 7. Key references

1. Olausen, S. A. & Kaspi, V. M., “The McGill Magnetar Catalog,” ApJS 212, 6 (2014); online update: https://www.physics.mcgill.ca/~pulsar/magnetar/main.html  
2. Palmer, D. M. et al., “A giant γ-ray flare from the magnetar SGR 1806−20,” *Nature* 434, 1107 (2005).  
3. Hurley, K. et al., multi-satellite giant-flare papers (1999, 2005).  
4. Schwinger limit / \(B_{\mathrm{QED}}\): standard QED; Wikipedia “Schwinger limit” and textbooks.  
5. Kaspi, V. M. & Beloborodov, A. M., magnetar reviews (ARA&A).
