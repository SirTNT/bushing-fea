"""Analytical checks: beam theory for the static and modal cases, 1D radial heat conduction for the
thermal case, free thermal elongation for the thermal-structural case.

All beam checks treat the part above the flange as a cantilever clamped at the flange bottom face
(z = z_f0, the fixed support plane) and loaded at the conductor top (z = z_c1).
"""
import numpy as np
from scipy.linalg import eigh

from common import derived_geometry, load_config, read_results, update_results


def sections(cfg):
    d = derived_geometry(cfg)
    m = cfg["materials"]
    A_c = np.pi * (d["r_co"] ** 2 - d["r_ci"] ** 2)
    A_k = np.pi * (d["r_k"] ** 2 - d["r_co"] ** 2)
    I_c = np.pi / 4 * (d["r_co"] ** 4 - d["r_ci"] ** 4)
    I_k = np.pi / 4 * (d["r_k"] ** 4 - d["r_co"] ** 4)
    EI_c = m["copper"]["E"] * I_c
    EI_k = m["epoxy"]["E"] * I_k
    mu_c = m["copper"]["rho"] * A_c
    mu_k = m["epoxy"]["rho"] * A_k
    return dict(A_c=A_c, A_k=A_k, I_c=I_c, I_k=I_k, EI_c=EI_c, EI_k=EI_k, EI_eff=EI_c + EI_k,
                m_c=mu_c, m_k=mu_k, m_eff=mu_c + mu_k)


def cantilever(cfg):
    d = derived_geometry(cfg)
    s = sections(cfg)
    F = cfg["loads"]["cantilever_force"]
    L = d["z_c1"] - d["z_f0"]            # flange (fixed plane) to conductor top
    a = d["L_core"] - d["z_f0"]          # segment 1: core plus conductor
    b = d["z_c1"] - d["L_core"]          # segment 2: bare conductor
    delta_uniform = F * L ** 3 / (3 * s["EI_eff"])
    EI1, EI2 = s["EI_eff"], s["EI_c"]
    # segment 1 carries shear F and end moment F*b; segment 2 is a cantilever on the rotated end of segment 1
    d1 = F * a ** 3 / (3 * EI1) + F * b * a ** 2 / (2 * EI1)
    th1 = F * a ** 2 / (2 * EI1) + F * b * a / EI1
    d2 = F * b ** 3 / (3 * EI2)
    delta_two = d1 + th1 * b + d2
    # bending stress in the bare conductor at the core top (moment F*b), for reference
    sigma_cond = F * b * d["r_co"] / s["I_c"]
    return {"F_N": F, "L_m": L, "a_core_segment_m": a, "b_bare_conductor_m": b,
            "EI_conductor_Nm2": s["EI_c"], "EI_core_Nm2": s["EI_k"], "EI_eff_Nm2": s["EI_eff"],
            "delta_uniform_m": delta_uniform, "delta_two_segment_m": delta_two,
            "two_segment_parts_m": {"segment1_deflection": d1, "segment1_rotation_times_b": th1 * b,
                                    "segment2_cantilever": d2},
            "bending_stress_conductor_at_core_top_Pa": sigma_cond}


def cantilever_profile(cfg, s_pos):
    """Two-segment beam deflection w(s) under tip load F, s measured from the fixed plane."""
    c = cantilever(cfg)
    F, a, b = c["F_N"], c["a_core_segment_m"], c["b_bare_conductor_m"]
    EI1, EI2 = c["EI_eff_Nm2"], c["EI_conductor_Nm2"]
    L = a + b
    s_pos = np.asarray(s_pos, dtype=float)
    # segment 1: M(s) = F (L - s); w'' = M / EI1
    w1 = F / EI1 * (L * s_pos ** 2 / 2 - s_pos ** 3 / 6)
    wa = F / EI1 * (L * a ** 2 / 2 - a ** 3 / 6)
    ta = F / EI1 * (L * a - a ** 2 / 2)
    t = s_pos - a
    w2 = wa + ta * t + F / EI2 * (b * t ** 2 / 2 - t ** 3 / 6)
    return np.where(s_pos <= a, w1, w2)


N_BEAM_ELEMENTS = 200


def stepped_beam_frequencies(cfg, n_el=N_BEAM_ELEMENTS, n=3):
    """Euler-Bernoulli finite element model (Hermite elements) of the two-segment cantilever."""
    d = derived_geometry(cfg)
    s = sections(cfg)
    a = d["L_core"] - d["z_f0"]
    L = d["z_c1"] - d["z_f0"]
    z = np.linspace(0, L, n_el + 1)
    ndof = 2 * (n_el + 1)
    K = np.zeros((ndof, ndof))
    M = np.zeros((ndof, ndof))
    for e in range(n_el):
        le = z[e + 1] - z[e]
        zm = 0.5 * (z[e] + z[e + 1])
        EI, mu = (s["EI_eff"], s["m_eff"]) if zm < a else (s["EI_c"], s["m_c"])
        ke = EI / le ** 3 * np.array([[12, 6 * le, -12, 6 * le], [6 * le, 4 * le ** 2, -6 * le, 2 * le ** 2],
                                      [-12, -6 * le, 12, -6 * le], [6 * le, 2 * le ** 2, -6 * le, 4 * le ** 2]])
        me = mu * le / 420 * np.array([[156, 22 * le, 54, -13 * le], [22 * le, 4 * le ** 2, 13 * le, -3 * le ** 2],
                                       [54, 13 * le, 156, -22 * le], [-13 * le, -3 * le ** 2, -22 * le, 4 * le ** 2]])
        idx = slice(2 * e, 2 * e + 4)
        K[idx, idx] += ke
        M[idx, idx] += me
    free = slice(2, ndof)  # clamped at z = 0
    w2 = eigh(K[free, free], M[free, free], eigvals_only=True, subset_by_index=[0, n - 1])
    return np.sqrt(w2) / (2 * np.pi)


def modal(cfg):
    d = derived_geometry(cfg)
    s = sections(cfg)
    L = d["z_c1"] - d["z_f0"]
    f1_uniform = 1.875 ** 2 / (2 * np.pi) * np.sqrt(s["EI_eff"] / (s["m_eff"] * L ** 4))
    f_step = stepped_beam_frequencies(cfg)
    return {"L_m": L, "EI_eff_Nm2": s["EI_eff"], "m_eff_kg_per_m": s["m_eff"],
            "f1_uniform_Hz": float(f1_uniform), "f1_two_segment_beam_Hz": float(f_step[0]),
            "f2_bending_two_segment_beam_Hz": float(f_step[1]), "n_beam_elements": N_BEAM_ELEMENTS}


def joule_heating(cfg):
    s = sections(cfg)
    I = cfg["loads"]["current"]
    rho_e = cfg["materials"]["copper"]["resistivity"]
    q = I ** 2 * rho_e / s["A_c"] ** 2
    return {"q_W_per_m3": q, "Q_per_length_W_per_m": q * s["A_c"],
            "Q_total_W": q * s["A_c"] * cfg["geometry"]["conductor"]["length"]}


def radial_temperature(cfg, r, q_scale=1.0):
    """Long cylinder: uniform generation q in the copper tube (insulated bore), conduction through
    the epoxy, convection h to T_amb at the core surface. Returns T(r) in C."""
    d = derived_geometry(cfg)
    m = cfg["materials"]
    ld = cfg["loads"]
    q = joule_heating(cfg)["q_W_per_m3"] * q_scale
    kc, ke, h = m["copper"]["k"], m["epoxy"]["k"], ld["h_conv"]
    rci, rco, rk = d["r_ci"], d["r_co"], d["r_k"]
    Qp = q * np.pi * (rco ** 2 - rci ** 2)
    T_s = ld["T_ambient"] + Qp / (2 * np.pi * rk * h)
    T_co = T_s + Qp / (2 * np.pi * ke) * np.log(rk / rco)
    r = np.asarray(r, dtype=float)
    T = np.where(r >= rco, T_s + Qp / (2 * np.pi * ke) * np.log(rk / np.maximum(r, rco)),
                 T_co + q / (4 * kc) * (rco ** 2 - r ** 2) - q * rci ** 2 / (2 * kc) * np.log(rco / np.maximum(r, rci)))
    return T


def thermal(cfg):
    d = derived_geometry(cfg)
    T = radial_temperature(cfg, [d["r_ci"], d["r_co"], d["r_k"]])
    jh = joule_heating(cfg)
    return {**jh, "T_bore_C": float(T[0]), "T_conductor_surface_C": float(T[1]),
            "T_core_surface_C": float(T[2]),
            "dT_epoxy_K": float(T[1] - T[2]), "dT_convection_K": float(T[2] - cfg["loads"]["T_ambient"])}


def run(cfg):
    out = {"sections": sections(cfg), "cantilever": cantilever(cfg), "modal": modal(cfg), "thermal": thermal(cfg)}
    update_results("hand_calcs", out)
    return out


def pct(fea, ref):
    return 100.0 * (fea - ref) / ref


def verify():
    """Compare hand calculations with FEA results stored in results.json; writes a 'verification' table."""
    r = read_results()
    hc = r["hand_calcs"]
    rows = []
    fea = r["static"]["S1_cantilever"]["top_deflection_x_m"]
    rows.append({"check": "S1 top deflection, uniform composite cantilever", "quantity": "deflection",
                 "unit": "mm", "analytical": hc["cantilever"]["delta_uniform_m"] * 1e3, "fea": fea * 1e3})
    rows.append({"check": "S1 top deflection, two-segment cantilever", "quantity": "deflection",
                 "unit": "mm", "analytical": hc["cantilever"]["delta_two_segment_m"] * 1e3, "fea": fea * 1e3})
    f1 = r["modal"]["modes"][0]["frequency_Hz"]
    rows.append({"check": "f1, uniform composite cantilever", "quantity": "frequency", "unit": "Hz",
                 "analytical": hc["modal"]["f1_uniform_Hz"], "fea": f1})
    rows.append({"check": "f1, two-segment beam (1D beam elements)", "quantity": "frequency", "unit": "Hz",
                 "analytical": hc["modal"]["f1_two_segment_beam_Hz"], "fea": f1})
    th = r["thermal"]
    Ta = r["inputs"]["loads"]["T_ambient"] if "inputs" in r else load_config()["loads"]["T_ambient"]
    rows.append({"check": "Conductor surface temperature rise above ambient at core mid-height, 1D radial",
                 "quantity": "temperature rise", "unit": "K",
                 "analytical": hc["thermal"]["T_conductor_surface_C"] - Ta, "fea": th["mid_height"]["T_conductor_surface_C"] - Ta})
    rows.append({"check": "Core surface temperature rise above ambient at core mid-height, 1D radial",
                 "quantity": "temperature rise", "unit": "K",
                 "analytical": hc["thermal"]["T_core_surface_C"] - Ta, "fea": th["mid_height"]["T_core_surface_C"] - Ta})
    rows.append({"check": "Temperature rise across epoxy at mid-height, 1D radial", "quantity": "temperature difference",
                 "unit": "K", "analytical": hc["thermal"]["dT_epoxy_K"], "fea": th["mid_height"]["dT_epoxy_K"]})
    rows.append({"check": "Temperature rise across epoxy at mid-height, 1D radial with FEA local heat flow",
                 "quantity": "temperature difference", "unit": "K",
                 "analytical": th["mid_height"]["dT_epoxy_1D_with_local_heat_flow_K"], "fea": th["mid_height"]["dT_epoxy_K"]})
    tm = r["thermo_mech"]["elongation_check"]
    rows.append({"check": "Conductor top axial displacement, free copper expansion", "quantity": "displacement",
                 "unit": "mm", "analytical": tm["free_copper_elongation_m"] * 1e3, "fea": tm["fea_top_uz_m"] * 1e3})
    rows.append({"check": "Conductor top axial displacement, bonded composite expansion (EA-weighted) plus free copper overhang",
                 "quantity": "displacement", "unit": "mm",
                 "analytical": tm["free_core_governed_elongation_m"] * 1e3, "fea": tm["fea_top_uz_m"] * 1e3})
    for row in rows:
        row["difference_pct"] = pct(row["fea"], row["analytical"])
    update_results("verification", rows)
    return rows


if __name__ == "__main__":
    cfg = load_config()
    out = run(cfg)
    c, m, t = out["cantilever"], out["modal"], out["thermal"]
    print(f"EI_eff = {c['EI_eff_Nm2']:.4e} N m2, EI_cond = {c['EI_conductor_Nm2']:.4e}, EI_core = {c['EI_core_Nm2']:.4e}")
    print(f"delta uniform = {c['delta_uniform_m'] * 1e3:.4f} mm, two-segment = {c['delta_two_segment_m'] * 1e3:.4f} mm")
    print(f"f1 uniform = {m['f1_uniform_Hz']:.2f} Hz, two-segment beam = {m['f1_two_segment_beam_Hz']:.2f} Hz")
    print(f"q = {t['q_W_per_m3']:.4e} W/m3, Q' = {t['Q_per_length_W_per_m']:.2f} W/m, total {t['Q_total_W']:.2f} W")
    print(f"T conductor surface {t['T_conductor_surface_C']:.2f} C, core surface {t['T_core_surface_C']:.2f} C")
