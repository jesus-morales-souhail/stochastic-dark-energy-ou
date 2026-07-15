import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ============================================================
# Parámetros físicos y de simulación
# ============================================================
L = 200.0          # Longitud del dominio
N = 2048           # Número de puntos de la malla
dx = L / N
x = np.linspace(-L/2, L/2, N, endpoint=False)

# Masa efectiva negativa (adimensional)
m_tilde = -1.0     # m* < 0

# Parámetros de interacción
g = 1.0            # Interacción repulsiva (g > 0)

# Potencial externo: doble pozo o trampa armónica
V_ext = 0.5 * (0.01 * x**2)   # Trampa armónica suave

# Estado inicial: paquete gaussiano centrado en x0 con momento k0
x0 = -20.0
k0 = 0.0
sigma0 = 2.0
psi = (1.0 / (np.pi * sigma0**2)**0.25) * np.exp(-(x - x0)**2 / (2 * sigma0**2)) * np.exp(1j * k0 * x)

# Normalización
psi = psi / np.sqrt(np.sum(np.abs(psi)**2) * dx)

# ============================================================
# Parámetros de evolución temporal
# ============================================================
dt = 0.005
T_total = 50.0
N_steps = int(T_total / dt)

# ============================================================
# Espacio de momentos (frecuencias)
# ============================================================
dk = 2 * np.pi / L
k = np.fft.fftfreq(N, d=dx) * 2 * np.pi   # Orden: [0, 1, ..., N/2-1, -N/2, ..., -1]

# Operador cinético en el espacio de Fourier (con masa negativa)
T_k = np.exp(-1j * (k**2 / (2 * m_tilde)) * dt)

# ============================================================
# Función de evolución con split-step de segundo orden
# ============================================================
def evolve_step(psi, V_ext, g, dt, T_k):
    """
    Un paso de evolución usando Strang splitting.
    """
    # 1. Mitad del paso potencial (posición)
    V_half = V_ext + g * np.abs(psi)**2
    psi = np.exp(-1j * V_half * dt / 2) * psi
    
    # 2. Paso cinético completo (momento)
    psi_k = np.fft.fft(psi)
    psi_k = T_k * psi_k
    psi = np.fft.ifft(psi_k)
    
    # 3. Mitad del paso potencial (posición)
    V_half = V_ext + g * np.abs(psi)**2
    psi = np.exp(-1j * V_half * dt / 2) * psi
    
    return psi

# ============================================================
# Bucle de evolución y almacenamiento
# ============================================================
print(f"Evolucionando durante {N_steps} pasos con dt={dt}...")

# Almacenar cada `save_every` pasos para visualización
save_every = 50
N_saves = N_steps // save_every + 1
psi_history = np.zeros((N_saves, N), dtype=complex)
times = np.zeros(N_saves)

psi_history[0] = psi.copy()
times[0] = 0.0

for step in range(1, N_steps + 1):
    psi = evolve_step(psi, V_ext, g, dt, T_k)
    
    if step % save_every == 0:
        idx = step // save_every
        psi_history[idx] = psi.copy()
        times[idx] = step * dt
        
        if idx % 10 == 0:
            print(f"Paso {step}/{N_steps} completado")

# ============================================================
# Visualización de la densidad y la fase
# ============================================================
fig, axes = plt.subplots(2, 1, figsize=(12, 8))

# Densidad |ψ|²
ax1 = axes[0]
ax1.set_title("Evolución de la densidad |ψ(x,t)|²")
ax1.set_xlabel("x")
ax1.set_ylabel("Densidad")
ax1.set_xlim(-L/2, L/2)

# Fase de la función de onda
ax2 = axes[1]
ax2.set_title("Fase de ψ(x,t)")
ax2.set_xlabel("x")
ax2.set_ylabel("Fase (rad)")
ax2.set_xlim(-L/2, L/2)

# Elegir algunos tiempos para graficar
plot_indices = [0, N_saves//4, N_saves//2, 3*N_saves//4, N_saves-1]
colors = ['blue', 'green', 'orange', 'red', 'purple']

for idx, color in zip(plot_indices, colors):
    t = times[idx]
    dens = np.abs(psi_history[idx])**2
    fase = np.angle(psi_history[idx])
    
    ax1.plot(x, dens, color=color, label=f"t = {t:.1f}")
    ax2.plot(x, fase, color=color, label=f"t = {t:.1f}")

ax1.legend()
ax2.legend()
plt.tight_layout()
plt.savefig("gpe_negative_mass_evolution.png", dpi=150)
plt.show()

# ============================================================
# Análisis de la inestabilidad: crecimiento de la densidad pico
# ============================================================
peak_density = np.array([np.max(np.abs(psi_history[i])**2) for i in range(N_saves)])

fig2, ax = plt.subplots(figsize=(10, 5))
ax.plot(times, peak_density, 'b-', linewidth=2)
ax.set_xlabel("Tiempo t")
ax.set_ylabel("Densidad pico |ψ|²_max")
ax.set_title("Crecimiento de la densidad pico (inestabilidad dinámica)")
ax.grid(True, alpha=0.3)
plt.savefig("gpe_negative_mass_instability.png", dpi=150)
plt.show()

print("Simulación completada.")