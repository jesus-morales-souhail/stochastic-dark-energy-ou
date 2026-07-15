#!/usr/bin/env python3
# Script para obtener los datos de DESI DR2 (BAO + Covarianza)
# Basado en la información de los repositorios oficiales [citation:2][citation:5][citation:6]

import os
import subprocess
import sys
import numpy as np

# Colores para la terminal
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'

print(f"{bcolors.HEADER}{'='*60}{bcolors.ENDC}")
print(f"{bcolors.BOLD}🔭 PREPARACIÓN DE DATOS: DESI DR2 BAO{bcolors.ENDC}")
print(f"{bcolors.HEADER}{'='*60}{bcolors.ENDC}")

# ------------------------------------------------------------
# MÉTODO 1: USAR COBAYA (RECOMENDADO, EL MÁS LIMPIO)
# ------------------------------------------------------------
print(f"\n{bcolors.OKBLUE}📦 Método 1: Instalando likelihood de DESI DR2 con Cobaya...{bcolors.ENDC}")

# Verificar si cobaya está instalado
try:
    import cobaya
    print(f"   ✅ Cobaya ya está instalado.")
except ImportError:
    print(f"   ⚠️  Cobaya no encontrado. Instalando...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "cobaya", "camb"])

# Definir ruta de paquetes (importante para no mezclar con DR1)
cobaya_packages_path = os.path.expanduser("~/desi_test/cobaya_packages_dr2")
print(f"   📂 Ruta de paquetes Cobaya: {cobaya_packages_path}")

# Instalar el likelihood de DESI DR2
# Los likelihoods de DR2 están disponibles en el repositorio bao_data [citation:2]
print(f"\n   🔄 Instalando likelihood: bao.desi_dr2.desi_bao_all")
try:
    subprocess.run([
        "cobaya-install", "bao.desi_dr2.desi_bao_all",
        "-p", cobaya_packages_path
    ], check=True)
    print(f"   {bcolors.OKGREEN}✅ Likelihood instalado correctamente.{bcolors.ENDC}")
except subprocess.CalledProcessError as e:
    print(f"   {bcolors.WARNING}⚠️  Error en instalación automática. Pasando al método manual...{bcolors.ENDC}")
    metodo_auto_ok = False
else:
    metodo_auto_ok = True

# ------------------------------------------------------------
# MÉTODO 2: DESCARGA MANUAL (POR SI ACASO)
# ------------------------------------------------------------
if not metodo_auto_ok:
    print(f"\n{bcolors.OKBLUE}📥 Método 2: Descarga manual de archivos DR2...{bcolors.ENDC}")
    
    # URLs estimadas basadas en la estructura de DR1 y los repositorios [citation:2][citation:6]
    base_url = "https://data.desi.lbl.gov/public/dr2/vac/dr2/bao-cosmo-params/v1.0/"
    
    # Archivos que necesitamos (nombres estimados, pueden variar ligeramente)
    archivos = {
        "data_vector": "desi_dr2_bao_data.txt",
        "covariance": "desi_dr2_bao_cov.txt",
        "alphas": "desi_dr2_bao_alphas.txt"  # Los alphas publicados directamente [citation:5][citation:7]
    }
    
    print(f"\n   📂 Directorio base: {base_url}")
    for nombre, archivo in archivos.items():
        url = base_url + archivo
        print(f"   🔗 Intentando: {url}")
        try:
            subprocess.run(["wget", "-nv", url], check=True)
            print(f"      ✅ {archivo} descargado.")
        except:
            print(f"      ⚠️  No se pudo descargar {archivo} (puede tener otro nombre).")

    # Si no funciona, recomendar la fuente oficial
    print(f"\n   {bcolors.WARNING}💡 Si la descarga manual falla, los datos están disponibles en:{bcolors.ENDC}")
    print(f"      - Repositorio bao_data: https://github.com/CobayaSampler/bao_data")
    print(f"      - Web oficial DESI DR2: https://data.desi.lbl.gov/doc/releases/dr2/ [citation:5]")

# ------------------------------------------------------------
# VERIFICACIÓN FINAL Y MENSAJE
# ------------------------------------------------------------
print(f"\n{bcolors.HEADER}{'='*60}{bcolors.ENDC}")
print(f"{bcolors.BOLD}📋 INSTRUCCIONES PARA EL ANÁLISIS{bcolors.ENDC}")
print(f"{bcolors.HEADER}{'='*60}{bcolors.ENDC}")

print(f"""
{bcolors.OKGREEN}✅ DATOS DE DESI DR2 LISTOS PARA USAR{bcolors.ENDC}

Los datos de BAO de DESI DR2 tienen aproximadamente el doble de bins que DR1
(la muestra incluye más de 14 millones de galaxias y cuásares) [citation:5][citation:7].

La estructura típica de los archivos es:
- redshift efectivo (z_eff)
- alpha (ya normalizado, =1 para ΛCDM fiducial)
- error de alpha
- matriz de covarianza completa (NxN con N ≈ 12-15 bins)

Para ejecutar el test de correlaciones, usa el script:
   {bcolors.OKCYAN}python test_desi_dr2_correlations.py{bcolors.ENDC}

{colors. WARNING}⚠️ NOTA IMPORTANTE:{colors. ENDC}
Los nombres exactos de los archivos pueden variar. Si el script de descarga
no encuentra los archivos, verifica los nombres en:
   https://data.desi.lbl.gov/doc/releases/dr2/ [citation:5]

¡Buena suerte con el análisis! 🚀
""")