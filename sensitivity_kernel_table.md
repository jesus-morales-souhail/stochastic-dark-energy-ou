# Sensitivity Kernel S(z) for DESI BAO Tracers

S(z) = ∂ ln D_V(z) / ∂ Ω_Λ

Computed numerically using flat ΛCDM fiducial (Ω_m = 0.315, H0 = 67.4 km/s/Mpc).

| z_eff | Tracer       | S(z)   | |S(z)| |

|-------|--------------|--------|--------|

| 0.295 | BGS          | −0.284 | 0.284  |

| 0.510 | LRG1         | −0.462 | 0.462  |

| 0.706 | LRG2         | −0.595 | 0.595  |

| 0.934 | LRG3+ELG1    | −0.719 | 0.719  |

| 1.321 | ELG2         | −0.870 | 0.870  |

| 1.484 | QSO          | −0.917 | 0.917  |

| 2.330 | Lyα          | −1.070 | 1.070  |

**Key bin:** z = 0.934 (LRG3+ELG1) is the diagnostic bin where the predicted 
OU floor exceeds the DESI DR2 observational error by ~2.8σ. This is the 
cleanest signal available with current data.

## Rayleigh cosmological criterion

Survey range Δx = ln(1+z_max) − ln(1+z_min):

| Survey       | z range      | Δx   | ω_R_min | Max testable ω_R |
|--------------|--------------|------|---------|-----------------|
| DESI DR1     | 0.295–2.33   | 0.94 | 6.68    | BAO floor test  |
| DESI DR2     | 0.295–2.33   | 0.94 | 6.68    | BAO floor test  |
| Euclid DR1   | ~0.9–1.8     | ~0.5 | ~12.6   | BAO floor test  |

For ω_R < ω_R_min, the oscillatory kernel H1 (QNM) cannot be distinguished
from the pure OU kernel H0 regardless of photometric precision.
This is a geometric limit of the survey range, not a noise limit.
