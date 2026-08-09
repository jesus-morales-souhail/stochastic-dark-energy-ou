# Sensitivity Kernel S(z) for DESI BAO Tracers

Jesús Morales Souhail · github.com/jesus-morales-souhail · July 2026

---

I define the sensitivity kernel as


$$
S(z) = \frac{\partial \ln D_V(z)}{\partial \Omega_\Lambda}.
$$


I compute it numerically for a flat ΛCDM fiducial with $\Omega_m = 0.315$ and $H_0 = 67.4$ km/s/Mpc.

| z_eff | Tracer | S(z) | |S(z)| |
|-------|--------------|--------|--------|
| 0.295 | BGS | −0.284 | 0.284 |
| 0.510 | LRG1 | −0.462 | 0.462 |
| 0.706 | LRG2 | −0.595 | 0.595 |
| 0.934 | LRG3+ELG1 | −0.719 | 0.719 |
| 1.321 | ELG2 | −0.870 | 0.870 |
| 1.484 | QSO | −0.917 | 0.917 |
| 2.330 | Lyα | −1.070 | 1.070 |

**Key bin.** At $z = 0.934$(LRG3+ELG1) the predicted OU floor exceeds the DESI DR2 observational error by ~2.8σ. With the present data that is the cleanest diagnostic bin I have.

## Rayleigh cosmological criterion

Survey range $\Delta x = \ln(1+z_{\max}) - \ln(1+z_{\min})$:

| Survey | z range | Δx | ω_R_min | Max testable ω_R |
|--------------|--------------|------|---------|-----------------|
| DESI DR1 | 0.295–2.33 | 0.94 | 6.68 | BAO floor test |
| DESI DR2 | 0.295–2.33 | 0.94 | 6.68 | BAO floor test |
| Euclid DR1 | ~0.9–1.8 | ~0.5 | ~12.6 | BAO floor test |

For $\omega_R < \omega_{R,\min}$, the oscillatory kernel H1 (QNM) cannot be distinguished from the pure OU kernel H0, no matter how good the photometry is. That is a geometric limit set by survey range, not a noise limit.

