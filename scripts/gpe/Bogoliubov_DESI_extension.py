import numpy as np

def obtener_potencial(x, tipo="armonico", **kwargs):
    """
    Permite alternar dinámicamente entre diferentes geometrías de confinamiento.
    """
    if tipo == "armonico":
        omega = kwargs.get("omega", 0.01)
        return 0.5 * (omega * x**2)
        
    elif tipo == "doble_pozo":
        a = kwargs.get("a", 0.005)
        b = kwargs.get("b", 0.5)
        return a * x**4 - b * x**2
        
    elif tipo == "red_optica":
        V0 = kwargs.get("V0", 1.0)  # Profundidad de la red
        k_L = kwargs.get("k_L", 1.0) # Vector de onda del láser
        return V0 * np.sin(k_L * x)**2
        
    else:
        return np.zeros_like(x)

def obtener_operador_cinetico(k, dt, relacion_dispersion="masa_negativa_pura", **kwargs):
    """
    Soporta relaciones de dispersión no parabólicas extraídas de acoplamientos espín-órbita.
    """
    if relacion_dispersion == "masa_negativa_pura":
        m_tilde = kwargs.get("m_tilde", -1.0)
        return np.exp(-1j * (k**2 / (2 * m_tilde)) * dt)
        
    elif tipo == "espin_orbita_efectivo":
        # Relación de dispersión real con estructura de mínimos locales invertidos
        k_soc = kwargs.get("k_soc", 1.0)
        Omega = kwargs.get("Omega", 1.0) # Acoplamiento de Raman
        E_k = (k**2 / 2) - np.sqrt((k * k_soc)**2 + (Omega / 2)**2)
        return np.exp(-1j * E_k * dt)