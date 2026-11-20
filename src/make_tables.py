"""Generate every number used in the report, README, explainer and CV entry from results/results.json.

- report/generated/values.tex : one macro per value, used in report.tex as \\V{key}
- report/generated/*.tex      : tables
- README.md, docs/EXPLAINER.md, docs/CV_ENTRY.md : filled from docs/templates/*.md ({{key}} placeholders)

A placeholder without a value stops the script, so no hand-typed number can slip in.
"""
import datetime as dt
import re

from common import ROOT, read_results

GEN = ROOT / "report" / "generated"
TEMPLATES = ROOT / "docs" / "templates"


def f(x, nd=2):
    return f"{x:,.{nd}f}"


def values(r):
    V = {}
    inp, g = r["inputs"], r["geometry"]
    ld, mt = inp["loads"], inp["materials"]
    # inputs and geometry
    V["F-N"] = f(ld["cantilever_force"], 0)
    V["seis-g"] = f(ld["seismic_horizontal_g"], 1)
    V["I-A"] = f(ld["current"], 0)
    V["h-conv"] = f(ld["h_conv"], 0)
    V["T-amb"] = f(ld["T_ambient"], 0)
    V["T-ref"] = f(ld["T_reference"], 0)
    V["g"] = f(ld["g"], 2)
    for k, nd in [("r_ci", 0), ("r_co", 0), ("r_k", 0), ("r_f", 0), ("t_f", 0)]:
        V[f"{k.replace('_', '')}-mm"] = f(g[k] * 1e3, nd)
    for k in ["L_core", "L_cond", "z_f0", "z_f1", "z_c0", "z_c1"]:
        V[f"{k.replace('_', '')}-m"] = f(g[k], 2)
    V["Dco-mm"] = f(2 * g["r_co"] * 1e3, 0)
    V["Dci-mm"] = f(2 * g["r_ci"] * 1e3, 0)
    V["Dk-mm"] = f(2 * g["r_k"] * 1e3, 0)
    V["Df-mm"] = f(2 * g["r_f"] * 1e3, 0)
    V["excl-mm"] = f(inp["postprocessing"]["corner_exclusion"] * 1e3, 0)
    V["n-modes"] = str(inp["modal"]["n_modes"])
    # software
    for k, v in r["software"].items():
        V[f"ver-{k}"] = v
    # mesh
    m = r["mesh"]
    V["mesh-h-mm"] = f(m["h_global_m"] * 1e3, 0)
    V["mesh-hc-mm"] = f(m["h_conductor_m"] * 1e3, 0)
    V["mesh-hf-mm"] = f(m["h_flange_m"] * 1e3, 0)
    V["mesh-nel"] = f(m["n_elements"], 0)
    V["mesh-nnodes"] = f(m["n_nodes"], 0)
    V["mesh-qmin"] = f(m["quality_minSICN_min"], 2)
    V["mesh-qmean"] = f(m["quality_minSICN_mean"], 2)
    V["mesh-gmin"] = f(m["quality_gamma_min"], 2)
    V["ndofs"] = f(r["static"]["n_dofs"], 0)
    V["mass-kg"] = f(r["static"]["total_mass_kg"], 1)
    # hand calcs
    hc = r["hand_calcs"]
    c = hc["cantilever"]
    V["EIc"] = f"{c['EI_conductor_Nm2']:.3g}"
    V["EIk"] = f"{c['EI_core_Nm2']:.3g}"
    V["EIeff"] = f"{c['EI_eff_Nm2']:.3g}"
    V["L-beam-m"] = f(c["L_m"], 2)
    V["a-beam-m"] = f(c["a_core_segment_m"], 2)
    V["b-beam-m"] = f(c["b_bare_conductor_m"], 2)
    V["meff"] = f(hc["modal"]["m_eff_kg_per_m"], 1)
    V["mc"] = f(hc["sections"]["m_c"], 2)
    V["sig-cond-beam-MPa"] = f(c["bending_stress_conductor_at_core_top_Pa"] / 1e6, 1)
    p = c["two_segment_parts_m"]
    V["seg1-mm"] = f(p["segment1_deflection"] * 1e3, 3)
    V["seg1rot-mm"] = f(p["segment1_rotation_times_b"] * 1e3, 3)
    V["seg2-mm"] = f(p["segment2_cantilever"] * 1e3, 3)
    V["share-overhang-pct"] = f(100 * (p["segment1_rotation_times_b"] + p["segment2_cantilever"]) /
                                c["delta_two_segment_m"], 0)
    th_hc = hc["thermal"]
    V["q-Wm3"] = f"{th_hc['q_W_per_m3']:.3g}"
    V["Qp-Wm"] = f(th_hc["Q_per_length_W_per_m"], 1)
    V["Qtot-W"] = f(th_hc["Q_total_W"], 1)
    # static
    st = r["static"]
    for case, key in [("S1_cantilever", "S1"), ("S2_self_weight", "S2"), ("S3_seismic_static_equivalent", "S3")]:
        s = st[case]
        V[f"{key}-top-mm"] = f(s["top_deflection_x_m"] * 1e3, 3)
        V[f"{key}-umax-mm"] = f(s["max_displacement_m"] * 1e3, 3)
        V[f"{key}-umax-z"] = f(s["max_displacement_location"]["z"], 2)
        for b, v in s["von_mises"].items():
            V[f"{key}-{b}-max"] = f(v["max_Pa"] / 1e6, 2)
            V[f"{key}-{b}-away"] = f(v["max_away_from_corners_Pa"] / 1e6, 2)
            V[f"{key}-{b}-maxz"] = f(v["max_location"]["z"], 2)
            V[f"{key}-{b}-awayz"] = f(v["away_location"]["z"], 2)
    V["S1-uxcoretop-mm"] = f(st["S1_cantilever"]["ux_at_core_top_m"] * 1e3, 3)
    V["S2-umax-um"] = f(st["S2_self_weight"]["max_displacement_m"] * 1e6, 2)
    V["S3-top-um"] = f(st["S3_seismic_static_equivalent"]["top_deflection_x_m"] * 1e6, 1)
    # modal
    md = r["modal"]
    for i, mode in enumerate(md["modes"], 1):
        V[f"f{i}"] = f(mode["frequency_Hz"], 2)
        V[f"f{i}-effx"] = f(100 * mode["eff_mass_fraction"]["x"], 1)
        V[f"f{i}-effy"] = f(100 * mode["eff_mass_fraction"]["y"], 1)
    V["f12-diff-pct"] = f(md["f1_f2_relative_difference_pct"], 3)
    V["f1-uniform"] = f(hc["modal"]["f1_uniform_Hz"], 1)
    V["n-beam-el"] = str(hc["modal"]["n_beam_elements"])
    V["f1-twoseg"] = f(hc["modal"]["f1_two_segment_beam_Hz"], 1)
    # thermal
    th = r["thermal"]
    V["Tmax"] = f(th["T_max_C"], 1)
    V["Tmax-z"] = f(th["T_max_location"]["z"], 2)
    V["Tmax-rise"] = f(th["T_max_C"] - ld["T_ambient"], 1)
    for b, v in th["T_max_per_body_C"].items():
        V[f"Tmax-{b}"] = f(v, 1)
    V["Q-gen"] = f(th["heat_generated_W"], 2)
    V["Q-out"] = f(th["heat_convected_W"], 2)
    V["energy-err"] = f"{abs(th['energy_balance_error_pct']):.1e}"
    mh = th["mid_height"]
    V["mid-z"] = f(mh["z_m"], 2)
    V["mid-Tco-fea"] = f(mh["T_conductor_surface_C"], 1)
    V["mid-Tk-fea"] = f(mh["T_core_surface_C"], 1)
    V["mid-Tco-an"] = f(mh["T_conductor_surface_analytical_C"], 1)
    V["mid-Tk-an"] = f(mh["T_core_surface_analytical_C"], 1)
    V["mid-dT-fea"] = f(mh["dT_epoxy_K"], 2)
    V["mid-dT-local"] = f(mh["dT_epoxy_1D_with_local_heat_flow_K"], 2)
    V["mid-Qlocal"] = f(mh["Q_radial_local_W_per_m"], 1)
    V["mid-maxdiff"] = f(mh["max_abs_difference_K"], 1)
    V["mid-Qlocal-excess-pct"] = f(100 * (mh["Q_radial_local_W_per_m"] / mh["Q_generated_per_length_W_per_m"] - 1), 0)
    V["dT-epoxy-an"] = f(th_hc["dT_epoxy_K"], 2)
    # thermal-structural
    tm = r["thermo_mech"]
    for case, key in [("thermal_only", "TH"), ("thermal_plus_S1", "THS1"), ("thermal_only_sliding_support", "THSL")]:
        for b, v in tm[case]["von_mises"].items():
            V[f"{key}-{b}-max"] = f(v["max_Pa"] / 1e6, 1)
            V[f"{key}-{b}-away"] = f(v["max_away_from_corners_Pa"] / 1e6, 1)
            V[f"{key}-{b}-maxz"] = f(v["max_location"]["z"], 2)
            V[f"{key}-{b}-maxr"] = f(v["max_location"]["r"] * 1e3, 0)
        for iface, sides in tm[case]["interfaces"].items():
            for b, v in sides.items():
                V[f"{key}-{iface}-{b}-max"] = f(v["max_Pa"] / 1e6, 1)
                V[f"{key}-{iface}-{b}-away"] = f(v["max_away_from_corners_Pa"] / 1e6, 1)
    V["TH-uz-mm"] = f(tm["thermal_only"]["top_displacement_z_m"] * 1e3, 3)
    V["THS1-top-mm"] = f(tm["thermal_plus_S1"]["top_deflection_x_m"] * 1e3, 3)
    el = tm["elongation_check"]
    V["el-dTcond"] = f(el["mean_dT_conductor_K"], 1)
    V["el-dTcondseg"] = f(el["mean_dT_conductor_core_segment_K"], 1)
    V["el-dTcore"] = f(el["mean_dT_core_segment_K"], 1)
    V["el-dTover"] = f(el["mean_dT_conductor_overhang_K"], 1)
    V["el-free-mm"] = f(el["free_copper_elongation_m"] * 1e3, 3)
    V["el-comp-mm"] = f(el["free_core_governed_elongation_m"] * 1e3, 3)
    V["el-fea-mm"] = f(el["fea_top_uz_m"] * 1e3, 3)
    # verification rows
    for i, row in enumerate(r["verification"], 1):
        nd = 3 if row["unit"] == "mm" else 2
        V[f"ver{i}-an"] = f(row["analytical"], nd)
        V[f"ver{i}-fea"] = f(row["fea"], nd)
        V[f"ver{i}-diff"] = f"{row['difference_pct']:+.1f}"
        V[f"ver{i}-absdiff"] = f"{abs(row['difference_pct']):.1f}"
    # convergence
    if "convergence" in r:
        cv = r["convergence"]
        for i, row in enumerate(cv["rows"]):
            k = ["coarse", "medium", "fine"][i]
            V[f"cv-{k}-h"] = f(row["h_m"] * 1e3, 1)
            V[f"cv-{k}-nel"] = f(row["n_elements"], 0)
            V[f"cv-{k}-dofs"] = f(row["n_dofs"], 0)
            V[f"cv-{k}-top"] = f(row["top_deflection_m"] * 1e3, 4)
            V[f"cv-{k}-vm"] = f(row["max_von_mises_Pa"] / 1e6, 1)
            V[f"cv-{k}-vmaway"] = f(row["max_von_mises_away_from_corners_Pa"] / 1e6, 1)
            V[f"cv-{k}-f1"] = f(row["f1_Hz"], 2)
            for b, d in row["von_mises_per_body"].items():
                V[f"cv-{k}-{b}-max"] = f(d["max_Pa"] / 1e6, 1)
                V[f"cv-{k}-{b}-away"] = f(d["away_Pa"] / 1e6, 1)
        for k, v in cv["change_medium_to_fine"].items():
            V[f"cv-chg-{k.replace('_', '')}"] = f"{v:+.2f}"
        for b, d in cv["change_medium_to_fine_per_body"].items():
            V[f"cv-chg-{b}-max"] = f"{d['max_pct']:+.1f}"
            V[f"cv-chg-{b}-away"] = f"{d['away_pct']:+.1f}"
    V["alpha-cu"] = f(mt["copper"]["alpha"] * 1e6, 0)
    V["alpha-ep"] = f(mt["epoxy"]["alpha"] * 1e6, 0)
    sec = hc["sections"]
    V["EA-ratio"] = f(mt["epoxy"]["E"] * sec["A_k"] / (mt["copper"]["E"] * sec["A_c"]), 1)
    V["run-time-min"] = f(r.get("run_time_s", 0) / 60, 0)
    V["date"] = dt.date.today().strftime("%B %Y")
    return V


def tex_escape(s):
    return s.replace("%", r"\%").replace("_", r"\_").replace("&", r"\&")


def write_tex(r, V):
    GEN.mkdir(parents=True, exist_ok=True)
    lines = ["% Generated by src/make_tables.py from results/results.json. Do not edit.",
             r"\makeatletter",
             r"\newcommand{\V}[1]{\ifcsname val@#1\endcsname\csname val@#1\endcsname\else\errmessage{Missing value #1}\fi}"]
    for k, v in V.items():
        lines.append(rf"\expandafter\def\csname val@{k}\endcsname{{{tex_escape(str(v))}}}")
    lines.append(r"\makeatother")
    (GEN / "values.tex").write_text("\n".join(lines) + "\n")

    # Verification table
    rows = []
    for i, row in enumerate(r["verification"], 1):
        rows.append(f"{tex_escape(row['check'])} & {row['unit'].replace('C', r'\si{\celsius}')} & "
                    f"{V[f'ver{i}-an']} & {V[f'ver{i}-fea']} & {V[f'ver{i}-diff']} \\\\")
    (GEN / "verification_table.tex").write_text("\n".join([
        r"\begin{tabular}{p{10.2cm}lrrr}", r"\toprule",
        r"Check & Unit & Analytical & FEA & Diff. (\%) \\", r"\midrule", *rows, r"\bottomrule", r"\end{tabular}"]) + "\n")

    # Materials table from the inputs
    mt = r["inputs"]["materials"]
    names = {"copper": "Copper (conductor)", "epoxy": "Epoxy (core)", "aluminium": "Aluminium (flange)"}
    mrows = []
    for k, n in names.items():
        p = mt[k]
        res = f"{p['resistivity']:.3g}" if "resistivity" in p else "n/a"
        mrows.append(f"{n} & {p['E'] / 1e9:.0f} & {p['nu']:.2f} & {p['rho']:.0f} & {p['k']:g} & "
                     f"{p['alpha'] * 1e6:.0f} & {res} \\\\")
    (GEN / "materials_table.tex").write_text("\n".join([
        r"\begin{tabular}{lrrrrrr}", r"\toprule",
        r"Material & $E$ (GPa) & $\nu$ & $\rho$ (kg/m$^3$) & $k$ (W/(m K)) & $\alpha$ ($10^{-6}$/K) & $\rho_e$ ($\Omega$ m) \\",
        r"\midrule", *mrows, r"\bottomrule", r"\end{tabular}"]) + "\n")

    # Static results
    srows = []
    for key, name in [("S1", "S1 cantilever"), ("S2", "S2 self-weight"), ("S3", "S3 static-equivalent seismic")]:
        srows.append(f"{name} & {V[f'{key}-top-mm']} & {V[f'{key}-umax-mm']} & "
                     + " & ".join(f"{V[f'{key}-{b}-max']} ({V[f'{key}-{b}-away']})" for b in ["conductor", "core", "flange"])
                     + r" \\")
    (GEN / "static_table.tex").write_text("\n".join([
        r"\begin{tabular}{lrrrrr}", r"\toprule",
        r"Load case & Top defl. (mm) & Max $|u|$ (mm) & Conductor (MPa) & Core (MPa) & Flange (MPa) \\",
        r"\midrule", *srows, r"\bottomrule", r"\end{tabular}"]) + "\n")

    # Modal
    mrows = []
    for i, mode in enumerate(r["modal"]["modes"], 1):
        mrows.append(f"{i} & {V[f'f{i}']} & {tex_escape(mode['type'])} & {tex_escape(mode['dominant_region'])} & "
                     f"{V[f'f{i}-effx']} & {V[f'f{i}-effy']} \\\\")
    (GEN / "modal_table.tex").write_text("\n".join([
        r"\begin{tabular}{rrllrr}", r"\toprule",
        r"Mode & $f$ (Hz) & Type & Region & $m_{\mathrm{eff},x}$ (\%) & $m_{\mathrm{eff},y}$ (\%) \\",
        r"\midrule", *mrows, r"\bottomrule", r"\end{tabular}"]) + "\n")

    # Thermal-structural
    trows = []
    for key, name in [("TH", "Thermal, fixed flange face"), ("THSL", "Thermal, sliding flange face"),
                      ("THS1", "Thermal + S1, fixed flange face")]:
        trows.append(f"{name} & " + " & ".join(f"{V[f'{key}-{b}-max']} ({V[f'{key}-{b}-away']})"
                                               for b in ["conductor", "core", "flange"]) + r" \\")
    (GEN / "thermo_mech_table.tex").write_text("\n".join([
        r"\begin{tabular}{lrrr}", r"\toprule",
        r"Case & Conductor (MPa) & Core (MPa) & Flange (MPa) \\", r"\midrule", *trows,
        r"\bottomrule", r"\end{tabular}"]) + "\n")
    irows = []
    for iface, lab in [("flange_to_core", "Flange to core"), ("conductor_to_core", "Conductor to core")]:
        for b in r["thermo_mech"]["thermal_only"]["interfaces"][iface]:
            irows.append(f"{lab} & {b} side & {V[f'TH-{iface}-{b}-max']} ({V[f'TH-{iface}-{b}-away']}) & "
                         f"{V[f'THSL-{iface}-{b}-max']} ({V[f'THSL-{iface}-{b}-away']}) \\\\")
    (GEN / "interface_table.tex").write_text("\n".join([
        r"\begin{tabular}{llrr}", r"\toprule",
        r"Interface & Side & Fixed face (MPa) & Sliding face (MPa) \\", r"\midrule", *irows,
        r"\bottomrule", r"\end{tabular}"]) + "\n")

    # Convergence
    if "convergence" in r:
        crows = []
        for k in ["coarse", "medium", "fine"]:
            crows.append(f"{k.capitalize()} & {V[f'cv-{k}-h']} & {V[f'cv-{k}-nel']} & {V[f'cv-{k}-dofs']} & "
                         f"{V[f'cv-{k}-top']} & {V[f'cv-{k}-vm']} & {V[f'cv-{k}-vmaway']} & {V[f'cv-{k}-f1']} \\\\")
        crows.append(r"\midrule")
        crows.append(f"Change medium to fine (\\%) & & & & {V['cv-chg-topdeflectionpct']} & {V['cv-chg-maxvonmisespct']} & "
                     f"{V['cv-chg-maxvonmisesawayfromcornerspct']} & {V['cv-chg-f1pct']} \\\\")
        (GEN / "convergence_table.tex").write_text("\n".join([
            r"\begin{tabular}{lrrrrrrr}", r"\toprule",
            r"Mesh & $h$ (mm) & Elements & DOFs & Top defl. (mm) & Max VM (MPa) & VM away (MPa) & $f_1$ (Hz) \\",
            r"\midrule", *crows, r"\bottomrule", r"\end{tabular}"]) + "\n")


def md_verification_table(r, V):
    out = ["| Check | Analytical | FEA | Difference |", "|---|---:|---:|---:|"]
    for i, row in enumerate(r["verification"], 1):
        u = row["unit"]
        out.append(f"| {row['check']} | {V[f'ver{i}-an']} {u} | {V[f'ver{i}-fea']} {u} | {V[f'ver{i}-diff']} % |")
    return "\n".join(out)


def fill_templates(r, V):
    V = dict(V)
    V["verification-table-md"] = md_verification_table(r, V)
    targets = {"README.md": ROOT / "README.md", "EXPLAINER.md": ROOT / "docs" / "EXPLAINER.md",
               "CV_ENTRY.md": ROOT / "docs" / "CV_ENTRY.md"}
    for name, dest in targets.items():
        src = TEMPLATES / name
        if not src.exists():
            continue
        text = src.read_text()
        missing = sorted(set(k for k in re.findall(r"\{\{([^}]+)\}\}", text) if k not in V))
        if missing:
            raise KeyError(f"{name}: no value for {missing}")
        text = re.sub(r"\{\{([^}]+)\}\}", lambda m: str(V[m.group(1)]), text)
        header = f"<!-- Generated by src/make_tables.py from docs/templates/{name} and results/results.json. -->\n"
        dest.write_text(header + text)
        print("written", dest.relative_to(ROOT))


def main():
    r = read_results()
    V = values(r)
    write_tex(r, V)
    fill_templates(r, V)
    print(f"{len(V)} values written to report/generated/values.tex")


if __name__ == "__main__":
    main()
