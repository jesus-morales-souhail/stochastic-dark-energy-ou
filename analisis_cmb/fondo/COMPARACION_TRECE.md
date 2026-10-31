# Texto: cotas con los 13 números de DESI, antes / después

Números de `analisis_cmb/fondo/perfil_13.json`.

**Diferencia con lo pedido:** donde el encargo decía "by 11–15 %", "by 14–25 %" y "by 11–25 %", el texto dice "by a factor 1.11–1.15", "1.14–1.25" y "1.11–1.25". Esos números son cocientes de cotas (7 α / 13). Leídos como reducción porcentual de la cota, la reducción real es del 10–13 % (fijo o perfilado) y del 13–20 % (marginalizado).

## 1. PREPRINT §3.2

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
with the same θ → 0 mode degenerate with $\Omega_m$ ([`analisis_cmb/sn/SN.md`](../analisis_cmb/sn/SN.md)).
```

**Después:**

```text
with the same θ → 0 mode degenerate with $\Omega_m$ ([`analisis_cmb/sn/SN.md`](../analisis_cmb/sn/SN.md)). Using the 13 published numbers ($D_M$ and $D_H$ separately, $D_V$ for BGS) instead of the 7 isotropic α tightens every limit by a factor 1.11–1.15 with the background fixed or profiled ($2.3\times 10^{-2}$ at θ → 0, $1.3\times 10^{-2}$ at θ = 1) and by a factor 1.14–1.25 marginalised ($3.8\times 10^{-2}$ at θ = 1), without changing its meaning; the baseline pipeline keeps the 7 α ([`analisis_cmb/fondo/TRECE.md`](../analisis_cmb/fondo/TRECE.md)).
```

## 2. CLAIMS C2

Fichero: `manuscript/CLAIMS.md`

**Antes:**

```text
meaning: `analisis_cmb/fondo/FONDO.md`, `analisis_cmb/sn/SN.md`;
```

**Después:**

```text
meaning: `analisis_cmb/fondo/FONDO.md`, `analisis_cmb/fondo/TRECE.md`, `analisis_cmb/sn/SN.md`;
```

## 3. analisis_cmb/README.md

Fichero: `analisis_cmb/README.md`

**Antes:**

```text
| `fondo/` | The DESI σ_X profile with the background free (`FONDO.md`, `perfil_fondo.py`): the α reference is the DESI DR2 best fit, the θ → 0 limit is the Ω_m direction (ln det term), θ ≥ 1 limits by background treatment; σ_X needed to resolve θ with DESI and a Euclid-like survey (`RESOLVER_THETA.md`, `resolver_theta.py`) |
```

**Después:**

```text
| `fondo/` | The DESI σ_X profile with the background free (`FONDO.md`, `perfil_fondo.py`): the α reference is the DESI DR2 best fit, the θ → 0 limit is the Ω_m direction (ln det term), θ ≥ 1 limits by background treatment; σ_X needed to resolve θ with DESI and a Euclid-like survey (`RESOLVER_THETA.md`, `resolver_theta.py`) |
| `fondo/TRECE.md` | The same profile with the 13 DESI numbers (D_M, D_H separately, D_V for BGS) instead of the 7 isotropic α (`perfil_13.py`): limits tighter by a factor 1.11–1.25, same meaning |
```

## 4. Limitación 1

Fichero: `manuscript/PREPRINT.md`

**Antes:**

```text
1. **Summary statistics:** baseline residual analysis uses published BAO $\alpha$ with the official Gaussian BAO covariance projected to $\alpha$ (block-diagonal across bins). Cross-bin covariances, if released, may change the numerical factor in $\sigma_X$ while preserving the null preference if residuals remain consistent with noise.
```

**Después:**

```text
1. **Summary statistics:** baseline residual analysis uses published BAO $\alpha$ with the official Gaussian BAO covariance projected to $\alpha$ (block-diagonal across bins). Cross-bin covariances, if released, may change the numerical factor in $\sigma_X$ while preserving the null preference if residuals remain consistent with noise. Using the full 13-number vector tightens the limits by a factor 1.11–1.25 (§3.2).
```
