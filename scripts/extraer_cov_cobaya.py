import numpy as np
import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
data_dir = DATA_DIR / "dr2_data/dr2-bao-zenodo/cosmology_chains/bao_data/desi_bao_dr2"

# Leer la media (segunda columna) y también los comentarios para saber el orden
# Usamos usecols=(1,2) para leer valor y etiqueta
# Pero el archivo tiene comentarios con '#', así que usamos comentarios para saltar la primera línea
mean_data = np.loadtxt(os.path.join(data_dir, 'desi_gaussian_bao_ALL_GCcomb_mean.txt'),
                       comments='#', usecols=(1,2), dtype=str)
# mean_data es un array de strings (valor y etiqueta)
# Convertir valores a float
values = mean_data[:,0].astype(float)
labels = mean_data[:,1]  # strings como 'DV_over_rs', 'DM_over_rs', 'DH_over_rs'

# Leer la covarianza (13x13)
cov_dist = np.loadtxt(os.path.join(data_dir, 'desi_gaussian_bao_ALL_GCcomb_cov.txt'))
assert cov_dist.shape == (13,13), "Covarianza no es 13x13"

# Mapear índices según las etiquetas
# Definir el orden de los alphas: BGS, LRG1, LRG2, LRG3+ELG1, ELG2, QSO, Lya
# Para cada uno, necesitamos el índice de DM y DH (o DV para BGS)

# Buscar índices de cada parámetro
idx = {}
for i, label in enumerate(labels):
    idx[label] = i

# Ahora construimos el Jacobiano
alphas = np.array([1.0030, 0.9947, 1.0016, 0.9960, 1.0020, 0.9963, 1.0008])

J = np.zeros((7, 13))

# BGS: alpha = DV / DV_fid, donde DV es el valor en idx['DV_over_rs']
DV_BGS = values[idx['DV_over_rs']]
J[0, idx['DV_over_rs']] = alphas[0] / DV_BGS

# LRG1: alpha = (DM^2 * DH)^(1/3) / fid
# DM_LRG1 y DH_LRG1
DM_LRG1 = values[idx['DM_over_rs']]  # ¿cuál es el índice? Pero hay varios DM. Necesitamos distinguir por redshift.
# Los labels no incluyen redshift, solo 'DM_over_rs'. Eso es un problema: hay múltiples DM y DH.
# Necesitamos una forma de identificar cuál es cuál. El orden en el archivo es el mismo que en la media.
# El archivo de media tiene el orden: DV_BGS, DM_LRG1, DH_LRG1, DM_LRG2, DH_LRG2, DM_LRG3, DH_LRG3, DM_ELG2, DH_ELG2, DM_QSO, DH_QSO, DH_Lya, DM_Lya.
# Entonces podemos usar índices fijos, ya que el orden es conocido.

# Usar índices fijos basados en el orden del archivo (0 a 12)
# Índices fijos:
# 0: DV_BGS
# 1: DM_LRG1
# 2: DH_LRG1
# 3: DM_LRG2
# 4: DH_LRG2
# 5: DM_LRG3
# 6: DH_LRG3
# 7: DM_ELG2
# 8: DH_ELG2
# 9: DM_QSO
# 10: DH_QSO
# 11: DH_Lya
# 12: DM_Lya

# Asignamos los valores medios a variables para claridad
DM_LRG1_vals = values[1]
DH_LRG1_vals = values[2]
# etc.

# Construir J con índices fijos
J[0, 0] = alphas[0] / values[0]  # DV_BGS

J[1, 1] = (2/3) * alphas[1] / values[1]  # DM_LRG1
J[1, 2] = (1/3) * alphas[1] / values[2]  # DH_LRG1

J[2, 3] = (2/3) * alphas[2] / values[3]
J[2, 4] = (1/3) * alphas[2] / values[4]

J[3, 5] = (2/3) * alphas[3] / values[5]
J[3, 6] = (1/3) * alphas[3] / values[6]

J[4, 7] = (2/3) * alphas[4] / values[7]
J[4, 8] = (1/3) * alphas[4] / values[8]

J[5, 9] = (2/3) * alphas[5] / values[9]
J[5, 10] = (1/3) * alphas[5] / values[10]

J[6, 11] = (1/3) * alphas[6] / values[11]  # DH_Lya
J[6, 12] = (2/3) * alphas[6] / values[12]  # DM_Lya

# Calcular cov_alpha
cov_alpha = J @ cov_dist @ J.T

print("Matriz de covarianza de alpha_iso (7x7):")
print(cov_alpha)
print("\nDiagonal (varianzas):", np.diag(cov_alpha))
print("\nMatriz de correlaciones (rho_ij):")
std = np.sqrt(np.diag(cov_alpha))
corr = cov_alpha / np.outer(std, std)
print(corr)

# Guardar
out_path = DATA_DIR / 'desi_cov_alpha_iso.npy'
np.save(out_path, cov_alpha)
print(f"\n✅ Guardado como '{out_path}'")