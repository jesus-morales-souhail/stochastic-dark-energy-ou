# First-Principles Mapping Tables (repo anchors only)

**No illustrative \(A_0=0.01\).** Seed: \(\sigma_0\sim 10^{-61}\) (repo Axiom A2 / GPE paper).

Fiducial: \(H_0=67.4\), \(\Omega_m=0.315\), \(t_0=13.796\,\mathrm{Gyr}\).

## Table 1 — \(\theta\) from repo → \(\gamma\), \(t_{1/2}\), Sorkin residual

| Case | \(\theta\) | \(\gamma\) [1/Gyr] | \(t_{1/2}\) [Gyr] | \(\sigma_{\rm res}\) (DESI path) |
|------|----------:|-------------------:|------------------:|----------------------------------:|
| MLE floor (resume/paper) | 1.0000e-03 | 6.8930e-05 | 1.0056e+04 | 9.991e-62 |
| Old base THETA_BASE=1.2 (superseded) | 1.2000e+00 | 8.2716e-02 | 8.3798e+00 | 3.237e-62 |
| Half-decay over DESI path Dx=0.94 | 7.3739e-01 | 5.0829e-02 | 1.3637e+01 | 5.000e-62 |
| Half-decay since recombination | 9.9094e-02 | 6.8306e-03 | 1.0148e+02 | 9.111e-62 |
| Hubble-scale: t_1/2 = t_H0 (Map A) | 6.9315e-01 | 4.7779e-02 | 1.4507e+01 | 5.212e-62 |
| Age-scale: t_1/2 = t0 (Map A) | 7.2887e-01 | 5.0241e-02 | 1.3796e+01 | 5.040e-62 |

## Table 2 — \(A_0\) required to hit DESI/Euclid (vs Sorkin \(10^{-61}\))

| \(\theta\) case | path | target | \(A_0\) required | \(A_0/\sigma_0\) |
|----------------|------|--------|-----------------:|-----------------:|
| theta_MLE=0.001 | DESI_0.94 | DESI_limit_1.5e-4 | 1.501e-04 | 1.501e+57 |
| theta_MLE=0.001 | DESI_0.94 | Euclid_1e-5 | 1.001e-05 | 1.001e+56 |
| theta_MLE=0.001 | recomb | DESI_limit_1.5e-4 | 1.511e-04 | 1.511e+57 |
| theta_MLE=0.001 | recomb | Euclid_1e-5 | 1.007e-05 | 1.007e+56 |
| theta=1.2 (old) | DESI_0.94 | DESI_limit_1.5e-4 | 4.634e-04 | 4.634e+57 |
| theta=1.2 (old) | DESI_0.94 | Euclid_1e-5 | 3.089e-05 | 3.089e+56 |
| theta=1.2 (old) | recomb | DESI_limit_1.5e-4 | 6.630e-01 | 6.630e+60 |
| theta=1.2 (old) | recomb | Euclid_1e-5 | 4.420e-02 | 4.420e+59 |
| theta=ln2~0.693 (Hubble half-life) | DESI_0.94 | DESI_limit_1.5e-4 | 2.878e-04 | 2.878e+57 |
| theta=ln2~0.693 (Hubble half-life) | DESI_0.94 | Euclid_1e-5 | 1.919e-05 | 1.919e+56 |
| theta=ln2~0.693 (Hubble half-life) | recomb | DESI_limit_1.5e-4 | 1.913e-02 | 1.913e+59 |
| theta=ln2~0.693 (Hubble half-life) | recomb | Euclid_1e-5 | 1.275e-03 | 1.275e+58 |
| theta=0.099 (half over recomb path) | DESI_0.94 | DESI_limit_1.5e-4 | 1.646e-04 | 1.646e+57 |
| theta=0.099 (half over recomb path) | DESI_0.94 | Euclid_1e-5 | 1.098e-05 | 1.098e+56 |
| theta=0.099 (half over recomb path) | recomb | DESI_limit_1.5e-4 | 3.000e-04 | 3.000e+57 |
| theta=0.099 (half over recomb path) | recomb | Euclid_1e-5 | 2.000e-05 | 2.000e+56 |

## Table 3 — Page/scrambling with \(S=10^{122}\)

| Scale | T [Gyr] | \(\gamma\) | \(\theta_A\) |
|-------|--------:|-----------:|-------------:|
| t_H0 | 1.450742e+01 | 4.777881e-02 | 6.931472e-01 |
| t0_LCDM | 1.379634e+01 | 5.024137e-02 | 7.288725e-01 |
| t_since_recomb | 1.379587e+01 | 5.024311e-02 | 7.288978e-01 |
| t_scr = t_H ln S | 4.075357e+03 | 1.700826e-04 | 2.467459e-03 |
| t_Page ~ t_H sqrt(S) | 1.450742e+62 | 4.777881e-63 | 6.931472e-62 |
| t_Page ~ t_H S | 1.450742e+123 | 4.777881e-124 | 6.931472e-123 |
