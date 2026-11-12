"""Comprobación adicional de A por diferencias finitas en nodos especiales: nodo 0, los 7 nodos BAO, el último nodo."""
import json, numpy as np
from modelo import Modelo
m = Modelo(); p0 = np.load("salidas/p0.npy"); A = np.load("salidas/A.npy"); N = m.N; eps = 1e-4
nodos = [0] + sorted(set(m.k_bao.tolist())) + [N - 1]
res = []
for j in nodos:
    X = np.zeros(N)
    if j == 0: X[1:] = -eps
    else: X[j] = eps
    fd = (m.vector_modelo(p0, X) - m.vector_modelo(p0, -X)) / (2 * eps)
    an = A[:, j]; esc = np.max(np.abs(an)); mask = np.abs(an) > 1e-6 * esc
    rel = np.abs(fd[mask] - an[mask]) / np.abs(an[mask])
    i = np.where(mask)[0][np.argmax(rel)]
    res.append(dict(nodo=int(j), z=float(m.z[j]), n_filas=int(mask.sum()), max_rel=float(rel.max()), fila_peor=int(i),
                    A_an=float(an[i]), A_fd=float(fd[i]), max_abs=float(np.max(np.abs(fd - an))), ok_1pct=bool(rel.max() < 0.01)))
    print(res[-1])
json.dump(res, open("salidas/paso2b.json", "w"), indent=1)
