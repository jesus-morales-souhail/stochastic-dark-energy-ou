# Car + drone @ 120 km/h: pupil diffraction under Newton and Einstein

**Author:** Jesús Morales Souhail · [github.com/jesus-morales-souhail](https://github.com/jesus-morales-souhail)
**Date:** July 2026
**Status:** Quantitative thought experiment (method notes). **Not a DESI claim.**
**Script:** `scripts/car_drone_pupil_newton_einstein.py`
**Related:** `papers/self-shielding-vs-untestability.md`, `papers/no-go-superoscillation-tesseract.md`, `papers/EXPLORATORY_BOUNDARY.md`

---

## What this note is

I use this note to pin down a simple setup: highway speed, a millimetre aperture, Newton vs Einstein. I am **not** reporting a cosmological detection, and I am **not** claiming that a car or a pupil measures dark energy. The point is scale and operator discipline — the same hygiene I apply when I refuse wrong scale claims in the BAO programme.

**Not included:** hypercubes, Coxeter $B_{4}$, or undeclared 4D optical engines.

---

## Setup

- Car on a highway at **constant** $v = 120 \mathrm{km/h}$.
- Drone flies **parallel at the same velocity** (relative velocity car–drone = 0).
- Inside the car: diffraction of light ($\lambda = 550\,\mathrm{nm}$) through an aperture **$D = 1\,\mathrm{mm}$** (“pupila 1 mm”), also compared to a 4 mm human pupil.
- Question I ask myself: does the Airy pattern change under Newton vs Einstein, and does cosmic expansion show up?

---

## Premises that hold

1. **EM rigidity.** Electromagnetic binding in atoms and solids is enormously stronger than Hubble tidal effects on human / car scales. Local bound systems do **not** expand with the Hubble flow.
2. **Inertial frames.** Constant velocity (no acceleration) ⇒ Galilean / special-relativistic inertial physics; a closed car is approximately a lab at rest for the passenger.

---

## Results (from the script)

### A) Expansion vs local binding


H_0 \sim 2\times 10^{-18} \mathrm{s}^{-1} \qquad\Rightarrow\qquad v_H = H_0 L


| $L$ | $v_H = H_0 L$ |
|-------|------------------|
| $10^{-10} \mathrm{m}$ (atom) | $\sim 10^{-28} \mathrm{m/s}$ |
| $1 \mathrm{m}$ | $\sim 10^{-18} \mathrm{m/s}$ |
| $10^6 \mathrm{m}$ | $\sim 10^{-12} \mathrm{m/s}$ |

In 60 s, stretch of 1 m: $\Delta L \sim 10^{-16} \mathrm{m}$ — irrelevant next to atomic sizes fixed by EM.

### B) Diffraction (car rest frame)

Airy angle to first minimum: $\theta \approx 1.22 \lambda/D$.

| $D$ | $\theta$ | Airy radius on $f\approx 17 \mathrm{mm}$ |
|-------|------------|-----------------------------------------------|
| 1 mm | $\sim 6.7\times 10^{-4} \mathrm{rad}$ | $\sim 11 \mu\mathrm{m}$ |
| 4 mm | $\sim 1.7\times 10^{-4} \mathrm{rad}$ | $\sim 2.9 \mu\mathrm{m}$ |

### C) Newton / Galileo

Co-moving car and drone: **relative velocity = 0**.
Passenger and drone see the **same** local diffraction geometry.
I do not need an absolute-space “diagonal light + $v$” story between them; they share a rest frame.

### D) Einstein (special relativity)


\beta = \frac{v}{c} \approx 1.1\times 10^{-7},\qquad \gamma-1 \approx 6\times 10^{-15}.


- **Correct setup** (matched velocities): relative $\beta=0$ ⇒ **no** contraction/Doppler between car and drone; identical Airy pattern.
- **Wrong setup** (drone fixed to road): aberration $\sim\beta$ rad $\approx 0.023''$, while Airy angle $\sim 138''$ for 1 mm ⇒ ratio $\beta/\theta_{\mathrm{Airy}}\sim 10^{-4}$. **Undetectable** as a “tesseract-like” distortion of the pupil pattern.

### E) When SR would matter optically

Aberration comparable to Airy (1 mm): $v \sim \theta c \approx 200 \mathrm{km/s}$ ($\sim 6000\times$ 120 km/h).
$\gamma-1\sim 1\%$: $v\sim 0.14 c$.

### F) Cosmological expansion during the drive

Not a 4D optical effect in the pupil. I measure expansion with **redshift of free photons over Gpc**, not with a 1 mm aperture at highway speed.

---

## Link to “the question the mechanism blinds”

| Local car/pupil | Cosmology (SDiff / $\sigma_X$) |
|-----------------|----------------------------------|
| EM binds matter against Hubble flow | SDiff projects isotropic vacuum stress $\propto g_{\mu\nu}$ |
| Wrong scale for expansion | Wrong operator for generic “does the vacuum tremble?” |
| Measure free light over cosmic baselines | Measure residuals / slip / specific models with amplitudes |

The “shield” is **local physics at the right scale**, not magic immunity forever.
The blinded question is **mis-scaled**, not “physics fails inside the car.”

---

## Conclusion

1. **Newton (co-moving):** same diffraction inside and for the drone.
2. **Einstein (co-moving):** same; relative rest.
3. **Einstein (drone on road):** SR corrections $\sim 10^{-15}$–$10^{-7}$, negligible vs Airy.
4. **No tesseract** appears in any consistent calculation at 120 km/h.
5. Same methodological rule as Option 0 / no-go optics: **derive the map, quote the amplitude, or do not claim the effect.**

Run: `python scripts/car_drone_pupil_newton_einstein.py`

---

## Status: CLOSED

**Result locked:** this setup tests **local EM diffraction + inertial SR at** $\beta\sim 10^{-7}$.
It does **not** test cosmic expansion or any tesseract optical engine.
It is **not** a DESI result.

| Diagnosis | Content |
|-----------|---------|
| Wrong **scale** | Hubble / cosmic geometry on millimetre–metre / seconds |
| Wrong **operator** | Pupil Airy pattern (and non-derived 4D “projection”) |

**Next experiments:** only those with **correct scale + correct operator** — see `papers/scale-operator-experiment-map.md`.

*End of hygiene note.*
