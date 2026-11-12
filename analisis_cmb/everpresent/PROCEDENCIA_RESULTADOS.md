# Procedencia de los resultados que respaldan "0/20000 y 0/500000" (PREPRINT §4.2, §8.3 y el abstract)

Comprobado el 2026-11-12 (sha256 sobre los ficheros en disco; nada recalculado). El PREPRINT (l. 239) cita como fuente
`EVERPRESENT_vs_DESI.md` §3a, §3c y §3e; los números de esa nota salen de los ficheros de `resultados/` listados aquí.

## Modelo 2 (Zwane et al. 2018, α = 0.8824, μ = 0.9804): 0/20000

| fichero (`resultados/`) | en git | sha256 | contenido |
|---|---|---|---|
| `bao_m2.json` | sí | `7ab4a000a7247c545afefe41fb1d641ff93b7df56a962cbf803cd181d5063ab6` | 20000 realizaciones con r_d de Planck: 1 pasa A y B (`passing_realisations`, 1 entrada) y falla el proxy CMB; malla de α |
| `bao_m2.txt` | sí | `06d5d70ab67851b694e5dc6460996702a32a62ed4580102c9e79d24b7216fdc4` | resumen legible |
| `bao_m2_bestfit.npz` | sí | `16e230f358ba8ba709dc82a2fe5bbe09c87ada54c6af26a21913978b48352ca8` | χ², H0 y Ω_de(rec) de las 20000 |
| `bao_m2_alphagrid.npz` | sí | `2d046dacd4fa97536cb6d63e056c81182c0e76b5273b975ffdc9cadb9e71eff0` | χ² en la malla de α |
| `rd_real.json` → `model2.c_rd_real_window_wm_2sigma` | sí | `0d71d3ef6559a9fee4b8fbf72f7efd336ae1e80a22f15477099a6687f14e4c99` | criterio final (c): n_total 20000, nA = 0, nB = 0 (χ²_min 230.7); sin ventana: 6/20000 |
| `rd_real_m2.npz` | sí | `0163cc3b8d8905955d881c068d7d970a48486c1e7f434c98ba30a8ad48fd1321` | r_d por realización y proxy en la ventana |

Scripts y comandos (desde `analisis_cmb/everpresent/`, con `../venv/bin/python` o `python3`):

    python3 bao_sim.py m2          # escribe bao_m2.{json,txt}, bao_m2_bestfit.npz, bao_m2_alphagrid.npz (semilla 20261016)
    python3 rd_real.py m2          # escribe rd_real_m2.npz y la parte model2 de rd_real.json (normales nuevas con semilla 20261017)

## Modelo 1 (Das et al. 2023, n = 10⁵): 0/500000 = 0/300000 (α = 0.01) + 0/200000 (α = 0.005)

| fichero (`resultados/`) | en git | sha256 | contenido |
|---|---|---|---|
| `bao_m1_alpha0.01.json` | sí | `18c9689d8bca43b51b87ee43c266741a4d811c3ba6ce1156e0d1331ea3dc7860` | n_seeds 300000, 82795 completan; `unconditional_over_seeds` nA = 0, nB = 0 |
| `bao_m1_alpha0.01.txt` | sí | `262c63a9eeaebf8d06093dad1e7d9bb74ddd48784de649b623064c35a0f9dea4` | resumen |
| `bao_m1_alpha0.01.npz` (40.4 MB) | sí | `15323a55e9df20ee7c138724152816261bea30ec60d51fda24a80cf5da03e311` | distancias, χ², H0, r_d parcial por semilla |
| `bao_m1_alpha0.005.json` | sí | `76cbab8fc82f204036184660a9335d2f077aa92023d1fffe58627654d039f331` | n_seeds 200000, 198858 completan; nA = 0, nB = 0 |
| `bao_m1_alpha0.005.txt` | sí | `347ce32601f3d895e3d563e4a9046aa7d785fd39543ede8eb980705a8ade550b` | resumen |
| `bao_m1_alpha0.005.npz` (100.9 MB) | **no** (`.gitignore` de esta carpeta; > 50 MB) | `2d4a41153a17cbaf6d5d9dc13704ace009b18bc2a2298562f3a6783de6b561c2` | ídem para α = 0.005; el blob `8c2fa0a7` sigue en el historial (commit a65c4f8) |
| `rd_real.json` → `model1_alpha0.01.c_rd_real_window_wm_2sigma` | sí | (arriba) | (c): n_total 300000, nA = 0, nB = 0 (χ²_min 11889); sin ventana 1/300000 |
| `rd_real.json` → `model1_alpha0.005.c_rd_real_window_wm_2sigma` | sí | (arriba) | (c): n_total 200000, nA = 0, nB = 0 (χ²_min 89.3); sin ventana 0/200000 |
| `rd_real_m1_alpha0.01.npz` (15.5 MB) | sí | `180b21605a4d8089eca9a2bcdc0c38c6ab040d447084b80e852b7032dcc74ed5` | r_d real por semilla |
| `rd_real_m1_alpha0.005.npz` (40.9 MB) | sí | `cab5c9a38e6e99890fef6b57e171150a03d6534d2a570c21d9b01e23c37c19bc` | ídem |
| `rd_planck.json` | sí | `dbd681b534ad848f1c9ff241cd856af188a79ba489376d23db7a142c10c345eb` | r_d de Planck y ω_m usados por `bao_sim.py` |

Scripts y comandos (`bao_sim.py m1 ALPHA NSEEDS`; bloques de 20000 semillas con semilla 20261016 + índice de bloque; tiempos del json):

    python3 bao_sim.py m1 0.01 300000     # 1534 s; escribe bao_m1_alpha0.01.{json,txt,npz}
    python3 bao_sim.py m1 0.005 200000    # 853 s;  escribe bao_m1_alpha0.005.{json,txt,npz}  (el npz es el que no está en git)
    python3 rd_real.py m1                 # lee los dos npz; escribe rd_real_m1_alpha*.npz y la parte model1_alpha* de rd_real.json

Prueba de identidad ya hecha (EVERPRESENT_vs_DESI.md, "Reproducción", 2026-10-17): regenerando `bao_m1_alpha0.005.npz` con el mismo
comando en otro directorio, json idéntico salvo `seconds` y npz idéntico bit a bit (mismo sha256 `2d4a41153a17cbaf…`).

Lo único que respalda los números y no está en git es `bao_m1_alpha0.005.npz`; sus resúmenes (`json`, `txt`) y el `rd_real.json` que
recoge el criterio final sí están. Ambos modelos tienen fichero de resultados: no falta ninguno.

Dependencias de datos: `desi_dr2_data.py` (copia local) lee el vector gaussiano BAO de 13 números y su covarianza y los 7 α_DV; su
procedencia y sha256 están en `data/PROCEDENCIA_DESI_DR2.md`.
