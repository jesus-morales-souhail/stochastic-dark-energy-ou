#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Descarga y extrae la matriz de covarianza completa de DESI DR2
desde el repositorio Zenodo (DOI: 10.5281/zenodo.16644576).
"""

import os
from pathlib import Path
import requests

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(exist_ok=True)
import zipfile
import tarfile
import numpy as np
from astropy.io import fits
import json

# ------------------------------------------------------------------
# 1. Descargar el archivo complementario desde Zenodo
# ------------------------------------------------------------------
DOI = "10.5281/zenodo.16644576"
# Usamos la API de Zenodo para obtener el enlace de descarga
api_url = f"https://zenodo.org/api/records/{DOI.split('.')[-1]}"
response = requests.get(api_url)
if response.status_code != 200:
    raise RuntimeError("No se pudo acceder a Zenodo API")

data = response.json()
files = data.get('files', [])
if not files:
    raise RuntimeError("No se encontraron archivos en el registro de Zenodo")

# Tomamos el primer archivo (normalmente el más grande)
file_info = files[0]
download_url = file_info['links']['self']
filename = file_info['key']

print(f"Descargando {filename} desde Zenodo...")
r = requests.get(download_url, stream=True)
with open(filename, 'wb') as f:
    for chunk in r.iter_content(chunk_size=8192):
        f.write(chunk)
print("Descarga completada.")

# ------------------------------------------------------------------
# 2. Descomprimir según extensión
# ------------------------------------------------------------------
ext = os.path.splitext(filename)[1].lower()
if ext == '.zip':
    with zipfile.ZipFile(filename, 'r') as zf:
        zf.extractall('.')
        extracted_files = zf.namelist()
elif ext in ['.gz', '.tar'] or filename.endswith('.tar.gz'):
    with tarfile.open(filename, 'r:gz') as tf:
        tf.extractall('.')
        extracted_files = tf.getnames()
else:
    # Si es un solo archivo, lo usamos directamente
    extracted_files = [filename]

print("Archivos extraídos:", extracted_files)

# ------------------------------------------------------------------
# 3. Buscar la matriz de covarianza (7x7)
# ------------------------------------------------------------------
cov_matrix = None
for fname in extracted_files:
    # Buscamos archivos con nombres típicos
    if 'cov' in fname.lower() or 'matrix' in fname.lower() or 'correlation' in fname.lower():
        print(f"Posible archivo de covarianza: {fname}")
        # Intentar leer según extensión
        if fname.endswith('.fits'):
            try:
                hdul = fits.open(fname)
                # Inspeccionar extensiones
                for i, hdu in enumerate(hdul):
                    if hdu.data is not None and len(hdu.data.shape) == 2:
                        # Buscar matriz 7x7
                        if hdu.data.shape == (7, 7):
                            cov_matrix = hdu.data
                            print(f"Matriz 7x7 encontrada en extensión {i} de {fname}")
                            break
                        # A veces la matriz está en una tabla con columnas
                        if hasattr(hdu, 'columns'):
                            # Intentar reconstruir desde columnas
                            colnames = hdu.columns.names
                            if 'COV' in colnames or 'COVARIANCE' in colnames:
                                # Puede ser una lista de 49 elementos
                                arr = hdu.data['COV'][0]  # asumimos primera fila
                                if len(arr) == 49:
                                    cov_matrix = arr.reshape(7, 7)
                                    print(f"Matriz 7x7 reconstruida desde columna COV en {fname}")
                                    break
                hdul.close()
            except Exception as e:
                print(f"No se pudo leer FITS {fname}: {e}")
        elif fname.endswith('.csv'):
            try:
                cov_matrix = np.loadtxt(fname, delimiter=',')
                if cov_matrix.shape == (7, 7):
                    print(f"Matriz 7x7 cargada desde CSV {fname}")
            except:
                pass
        elif fname.endswith('.txt') or fname.endswith('.dat'):
            try:
                # Puede ser una lista de 49 números en una columna
                raw = np.loadtxt(fname)
                if raw.size == 49:
                    cov_matrix = raw.reshape(7, 7)
                    print(f"Matriz 7x7 reconstruida desde TXT {fname}")
                elif raw.shape == (7, 7):
                    cov_matrix = raw
                    print(f"Matriz 7x7 cargada desde TXT {fname}")
            except:
                pass
        if cov_matrix is not None:
            break

if cov_matrix is None:
    raise RuntimeError("No se pudo encontrar la matriz de covarianza 7x7 en los archivos extraídos.")

# ------------------------------------------------------------------
# 4. Guardar como .npy para reutilizar
# ------------------------------------------------------------------
out_path = DATA_DIR / 'desi_cov_full.npy'
np.save(out_path, cov_matrix)
print(f"Matriz de covarianza guardada como '{out_path}'")
print("Dimensiones:", cov_matrix.shape)
print("Diagonal (varianzas):", np.diag(cov_matrix))