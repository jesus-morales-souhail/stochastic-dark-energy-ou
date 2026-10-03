# analisis_cmb — OU a varianza máxima DESI frente a Planck 2018

Pregunta (versión vigente): con σ_X = 2.506e-2 (límite 95 % del perfil DESI DR2 BAO,
`../results/profile_sigma_x/profile_sigma_x.json`), ¿alguna observable del CMB (TT a l<30,
C_L^φφ, θ_*) separa la banda 68 % de 200 realizaciones OU de la curva w0-wa más de lo que
permite Planck 2018? ¿Qué cota sobre σ_X impondría Planck? Método vigente: H0 de cada curva
resuelto para fijar θ_* al de ΛCDM; observables TT l<30, C_L^φφ y distribución de H0.

No toca los scripts de ajuste BAO del repo. Todas las decisiones, valores de procedencia y
problemas abiertos: `NOTAS_ANALISIS_CMB.md`.

```
venv/                 CAMB 2.0.4 (ignorado por git)
datos_planck/         ficheros del Planck Legacy Archive
ou_cmb.py             X(ln a) → ρ_DE → w(a) → CAMB (PPF)
planck_datos.py       carga de errores Planck (TT, bandpowers lensing, θ_*)
realizacion_unica.py  una realización θ=1 + prueba de convergencia 500/2000 puntos
diagnosticos/         pruebas de convergencia malla-camino / tabla, banda por Δ; obsoleto/ = versiones viejas
realizaciones.py      200 realizaciones por (θ, Δ, σ_X), H0 resuelto → resultados/<halofit>/sigma_<σ>/
estadisticos.py       Δχ², percentiles, H0 → resumen.json, tabla.md
resultados/           salidas (resultados/sigma_2.506e-02/ y tabla.py: método viejo, H0 fijo — obsoletos)
```

```bash
./venv/bin/python realizacion_unica.py
./venv/bin/python realizaciones.py --halofit takahashi            # σ_X = 2.506e-2
./venv/bin/python estadisticos.py resultados/takahashi/sigma_2.506e-02
```

Estado (2026-10-03): análisis terminado; resultado en `RESULTADO.md`.

## Volver a descargar los datos (`datos_planck/` está en .gitignore)

```bash
mkdir -p datos_planck
curl -L -o datos_planck/COM_PowerSpect_CMB-TT-full_R3.01.txt \
  "http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=COM_PowerSpect_CMB-TT-full_R3.01.txt"
# sha256: ccf3113604020536f6f13ccf51680a7316ad0f32da558eee7f625e613bdd5522
./venv/bin/pip install camb cobaya numpy scipy matplotlib
./venv/bin/cobaya-install planck_2018_lensing.native -p "$PWD/datos_planck/cobaya_packages"
```
