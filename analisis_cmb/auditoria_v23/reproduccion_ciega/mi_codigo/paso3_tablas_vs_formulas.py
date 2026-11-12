"""Paso 3 (tarea aparte, sección 13): tablas r_d(omega_m) y z*(omega_m) frente a fórmulas de ajuste escritas por mí.
   r_d: Eisenstein & Hu (1998), ec. 26:  s = 44.5 ln(9.83/omega_m) / sqrt(1 + 10 omega_b^{3/4})  Mpc
   z*:  Hu & Sugiyama (1996):  z* = 1048 [1 + 0.00124 omega_b^{-0.738}] [1 + g1 omega_m^{g2}],
        g1 = 0.0783 omega_b^{-0.238} / (1 + 39.5 omega_b^{0.763}),  g2 = 0.560 / (1 + 21.1 omega_b^{1.81})
   Diagnóstico adicional: r_d por la integral de EH98 (ecs. 2-6: z_eq, k_eq, z_d con b1,b2, R(z), s = 2/(3 k_eq) sqrt(6/R_eq) ln[...])."""
import json, numpy as np
omega_b = 0.02237
t = np.loadtxt("../tablas/rd_rs_zstar_vs_omega_m.csv", delimiter=",", comments="#")
om, rd_tab, rs_tab, zs_tab = t.T

def rd_EH98_ec26(om, ob=omega_b):
    return 44.5 * np.log(9.83 / om) / np.sqrt(1.0 + 10.0 * ob ** 0.75)

def zstar_HS96(om, ob=omega_b):
    g1 = 0.0783 * ob ** (-0.238) / (1.0 + 39.5 * ob ** 0.763)
    g2 = 0.560 / (1.0 + 21.1 * ob ** 1.81)
    return 1048.0 * (1.0 + 0.00124 * ob ** (-0.738)) * (1.0 + g1 * om ** g2)

def rd_EH98_integral(om, ob=omega_b, Theta27=2.7255 / 2.7):
    # EH98 ecs. 2-6 (horizonte de sonido en z_d con la fórmula cerrada)
    zeq = 2.50e4 * om * Theta27 ** (-4)
    keq = 7.46e-2 * om * Theta27 ** (-2)              # Mpc^-1
    b1 = 0.313 * om ** (-0.419) * (1.0 + 0.607 * om ** 0.674)
    b2 = 0.238 * om ** 0.223
    zd = 1291.0 * om ** 0.251 / (1.0 + 0.659 * om ** 0.828) * (1.0 + b1 * ob ** b2)
    R = lambda z: 31.5 * ob * Theta27 ** (-4) * (1e3 / z)
    Rd, Req = R(zd), R(zeq)
    return 2.0 / (3.0 * keq) * np.sqrt(6.0 / Req) * np.log((np.sqrt(1 + Rd) + np.sqrt(Rd + Req)) / (1 + np.sqrt(Req))), zd

res = {}
for nombre, tab, f in [("r_d_vs_EH98_ec26", rd_tab, rd_EH98_ec26), ("z_star_vs_HS96", zs_tab, zstar_HS96)]:
    v = f(om)
    dif = (tab - v) / v
    i = int(np.argmax(np.abs(dif)))
    j = int(np.argmin(np.abs(dif)))
    res[nombre] = dict(max_dif_rel=float(dif[i]), omega_m_max=float(om[i]), tabla=float(tab[i]), formula=float(v[i]),
                       min_abs_dif_rel=float(dif[j]), omega_m_min=float(om[j]),
                       dif_rel_en_om_0p08=float(dif[0]), dif_rel_en_om_0p14=float(dif[np.argmin(np.abs(om - 0.14))]),
                       dif_rel_en_om_0p26=float(dif[-1]), media_dif_rel=float(dif.mean()))
    print(nombre, res[nombre])
v_int, zd = rd_EH98_integral(om)
dif = (rd_tab - v_int) / v_int
i = int(np.argmax(np.abs(dif)))
res["r_d_vs_EH98_integral_cerrada_diag"] = dict(max_dif_rel=float(dif[i]), omega_m_max=float(om[i]), tabla=float(rd_tab[i]),
                                                formula=float(v_int[i]), dif_rel_en_om_0p14=float(dif[np.argmin(np.abs(om - 0.14))]),
                                                z_d_en_om_0p14=float(zd[np.argmin(np.abs(om - 0.14))]))
print("diag integral EH98:", res["r_d_vs_EH98_integral_cerrada_diag"])
# también r_s(z*) de la tabla frente a r_d: cociente
res["cociente_rs_rd_tabla"] = dict(min=float((rs_tab / rd_tab).min()), max=float((rs_tab / rd_tab).max()))
res["omega_b_usado"] = omega_b
with open("salidas/paso3.json", "w") as f:
    json.dump(res, f, indent=1)
