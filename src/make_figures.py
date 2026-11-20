"""Stage 8: figures (PNG, 300 dpi). Fields are cut on a meridian plane with PyVista and drawn with
Matplotlib (no OpenGL needed)."""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import pyvista as pv
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch

from common import FIGURES_DIR, RESULTS_DIR, TAGS, derived_geometry, load_config, read_results
from hand_calcs import cantilever_profile

# Palette: categorical slots 1-3 and a single-hue sequential ramp (light to dark blue)
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#d9d8d4"
SEQ = LinearSegmentedColormap.from_list(
    "seq_blue", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"])
BODY_COLOURS = {"conductor": C2, "core": C1, "flange": C3}
DPI = 300

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "axes.edgecolor": INK2,
    "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
    "savefig.dpi": DPI, "savefig.bbox": "tight", "legend.frameon": False,
})


def cut(grid, normal="y"):
    """Slice on the plane through the axis; returns a triangulated PolyData."""
    s = grid.slice(normal=normal, origin=(0, 0, 0)).triangulate()
    return s


def tri_of(s, normal="y", disp=None, scale=0.0):
    P = np.array(s.points)
    if disp is not None:
        P = P + scale * np.asarray(s.point_data[disp])
    h = P[:, 0] if normal == "y" else P[:, 1]
    faces = s.faces.reshape(-1, 4)[:, 1:]
    return mtri.Triangulation(h, P[:, 2], faces)


def field_plot(ax, s, values, normal="y", disp=None, scale=0.0, cmap=SEQ, vmin=None, vmax=None, edges=False):
    tri = tri_of(s, normal, disp, scale)
    tpc = ax.tripcolor(tri, values, shading="gouraud", cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True)
    if edges:
        ax.triplot(tri, color="white", lw=0.15, alpha=0.6)
    ax.set_aspect("equal")
    return tpc


def outline(ax, cfg, normal_sign=1, colour=INK2, lw=0.6):
    """Undeformed outline of the three bodies in the r-z half planes (both sides of the axis)."""
    d = derived_geometry(cfg)
    for sgn in (1, -1):
        rects = [(d["r_ci"], d["z_c0"], d["r_co"] - d["r_ci"], d["L_cond"]),
                 (d["r_co"], 0.0, d["r_k"] - d["r_co"], d["L_core"]),
                 (d["r_k"], d["z_f0"], d["r_f"] - d["r_k"], d["t_f"])]
        for (r0, z0, w, hh) in rects:
            xs = sgn * np.array([r0, r0 + w, r0 + w, r0, r0])
            ax.plot(xs, [z0, z0, z0 + hh, z0 + hh, z0], color=colour, lw=lw, ls="--")


def colourbar(fig, mappable, ax, label):
    cb = fig.colorbar(mappable, ax=ax, shrink=0.85, pad=0.02)
    cb.set_label(label)
    cb.outline.set_visible(False)
    return cb


def body_values_on_cut(s):
    return np.asarray(s.cell_data["body"]) if "body" in s.cell_data else None


def fig_geometry(cfg, res):
    g = pv.read(RESULTS_DIR / "static_vm.vtu")
    s = cut(g)
    tri = tri_of(s)
    body = np.asarray(s.cell_data["body"])
    colours = {TAGS[b]: c for b, c in BODY_COLOURS.items()}
    fc = np.array([TAGS[b] for b in BODY_COLOURS])
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 6.2), gridspec_kw={"width_ratios": [1, 1.35]})
    for ax, (xl, zl) in zip(axes, [((-0.17, 0.17), (-0.25, 1.45)), ((0.0, 0.17), (0.33, 0.50))]):
        for tag in fc:
            m = body == tag
            ax.tripcolor(mtri.Triangulation(tri.x, tri.y, tri.triangles[m]),
                         facecolors=np.ones(m.sum()), cmap=matplotlib.colors.ListedColormap([colours[tag]]),
                         edgecolors="none", rasterized=True)
        ax.triplot(tri, color="white", lw=0.12 if ax is axes[0] else 0.3)
        ax.set_xlim(*xl)
        ax.set_ylim(*zl)
        ax.set_aspect("equal")
        ax.set_xlabel("x (m)")
        ax.set_ylabel("z (m)")
    d = derived_geometry(cfg)
    axes[0].annotate("fixed support\n(flange\nbottom face)", xy=(0.12, d["z_f0"]), xytext=(0.19, 0.15),
                     fontsize=7, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6), annotation_clip=False)
    axes[0].annotate("terminal\nload F", xy=(0.02, d["z_c1"]), xytext=(0.19, 1.30), fontsize=7,
                     arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6), annotation_clip=False)
    axes[0].set_xticks([-0.1, 0.1])
    axes[0].set_title("Axial cut (y = 0)")
    axes[1].set_title("Mesh cut at the flange")
    handles = [Patch(color=c, label=f"{b} ({m})") for (b, c), m in
               zip(BODY_COLOURS.items(), ["copper", "epoxy", "aluminium"])]
    fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.02))
    m = res["mesh"]
    fig.suptitle(f"Simplified bushing model: {m['n_elements']:,} quadratic tetrahedra, h = {m['h_global_m'] * 1e3:.0f} mm",
                 fontsize=10)
    fig.savefig(FIGURES_DIR / "fig1_geometry_mesh.png")
    plt.close(fig)


def fig_s1_displacement(cfg, res):
    g = pv.read(RESULTS_DIR / "static_u.vtu")
    s = cut(g)
    u = np.asarray(s.point_data["S1_cantilever"])
    mag = np.linalg.norm(u, axis=1) * 1e3
    scale = 50.0
    fig, ax = plt.subplots(figsize=(4.2, 6.2))
    tpc = field_plot(ax, s, mag, disp="S1_cantilever", scale=scale)
    outline(ax, cfg)
    colourbar(fig, tpc, ax, "displacement magnitude (mm)")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")
    r = res["static"]["S1_cantilever"]
    ax.set_title(f"S1 cantilever, F = {cfg['loads']['cantilever_force']:.0f} N\n"
                 f"deformation x{scale:.0f}, top deflection {r['top_deflection_x_m'] * 1e3:.3f} mm", fontsize=9)
    fig.savefig(FIGURES_DIR / "fig2_s1_displacement.png")
    plt.close(fig)


def fig_s1_stress(cfg, res):
    g = pv.read(RESULTS_DIR / "static_vm.vtu")
    s = cut(g)
    vm = np.asarray(s.point_data["S1_cantilever"]) / 1e6
    d = derived_geometry(cfg)
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 5.2), gridspec_kw={"width_ratios": [0.55, 1.45]})
    P = np.array(s.points)
    for ax, (xl, zl), lab in zip(axes, [((-0.17, 0.17), (-0.25, 1.45)), ((0.0, 0.16), (0.30, 0.56))],
                                 ["von Mises (MPa)", "von Mises (MPa), zoom scale"]):
        inside = (P[:, 0] >= xl[0]) & (P[:, 0] <= xl[1]) & (P[:, 2] >= zl[0]) & (P[:, 2] <= zl[1])
        tpc = field_plot(ax, s, vm, vmin=0, vmax=float(np.percentile(vm[inside], 99.9)))
        ax.set_xlim(*xl)
        ax.set_ylim(*zl)
        ax.set_xlabel("x (m)")
        ax.set_ylabel("z (m)")
        colourbar(fig, tpc, ax, lab)
    axes[0].set_xticks([-0.1, 0.1])
    outline(axes[1], cfg, lw=0.5)
    for (rr, zz) in [(d["r_k"], d["z_f0"]), (d["r_k"], d["z_f1"])]:
        axes[1].add_patch(plt.Circle((rr, zz), cfg["postprocessing"]["corner_exclusion"], fill=False,
                                     ec=C2, lw=0.8, ls=":"))
    axes[0].set_title("Whole bushing")
    axes[1].set_title("Zoom at the flange (tension side, x > 0)")
    axes[1].text(0.005, 0.305, "dotted circles: corner exclusion radius", color=INK2, fontsize=7)
    fig.suptitle("S1 cantilever: von Mises stress on the axial cut (y = 0)", fontsize=10)
    fig.savefig(FIGURES_DIR / "fig3_s1_von_mises_flange.png")
    plt.close(fig)


def fig_modes(cfg, res):
    g = pv.read(RESULTS_DIR / "modal_u.vtu")
    fig, axes = plt.subplots(1, 2, figsize=(6.4, 6.0), sharey=True)
    modes = res["modal"]["modes"]
    for k, (ax, normal) in enumerate(zip(axes, ["y", "x"])):
        s = cut(g, normal)
        key = f"mode{k + 1}"
        u = np.asarray(s.point_data[key])
        mag = np.linalg.norm(u, axis=1)
        tpc = field_plot(ax, s, mag, normal=normal, disp=key, scale=0.12, vmin=0, vmax=1)
        outline(ax, cfg)
        ax.set_xlabel("x (m)" if normal == "y" else "y (m)")
        f = modes[k]["frequency_Hz"]
        ax.set_title(f"Mode {k + 1}: {f:.2f} Hz\n{modes[k]['type']}", fontsize=9)
    axes[0].set_ylabel("z (m)")
    colourbar(fig, tpc, axes[1], "normalised displacement (-)")
    fig.suptitle(f"First two mode shapes (difference {res['modal']['f1_f2_relative_difference_pct']:.3f} %)",
                 fontsize=10)
    fig.savefig(FIGURES_DIR / "fig4_mode_shapes.png")
    plt.close(fig)


def fig_thermal(cfg, res):
    g = pv.read(RESULTS_DIR / "thermal_T.vtu")
    s = cut(g)
    T = np.asarray(s.point_data["T"])
    th = res["thermal"]
    fig, axes = plt.subplots(1, 2, figsize=(7.8, 5.4), gridspec_kw={"width_ratios": [0.9, 1.3], "wspace": 0.55})
    tpc = field_plot(axes[0], s, T)
    axes[0].set_xlim(-0.17, 0.17)
    axes[0].set_xlabel("x (m)")
    axes[0].set_ylabel("z (m)")
    axes[0].set_title(f"Temperature, I = {cfg['loads']['current']:.0f} A\nhot spot {th['T_max_C']:.1f} C")
    colourbar(fig, tpc, axes[0], "temperature (C)")
    mh = th["mid_height"]
    r = np.array(mh["r_m"]) * 1e3
    ax = axes[1]
    ax.plot(r, mh["T_fea_C"], color=C1, lw=2, label="FEA")
    ax.plot(r, mh["T_analytical_C"], color=C2, lw=2, ls="--", label="1D, local generation")
    ax.plot(r, mh["T_analytical_local_heat_flow_C"], color=C3, lw=2, ls=":", label="1D, FEA local heat flow")
    d = derived_geometry(cfg)
    for rr in (d["r_co"], d["r_k"]):
        ax.axvline(rr * 1e3, color=GRID, lw=0.8, zorder=0)
    ax.text(d["r_co"] * 1e3 + 1, ax.get_ylim()[0] + 1, "copper | epoxy", fontsize=7, color=INK2)
    ax.set_xlabel("radius r (mm)")
    ax.set_ylabel("temperature (C)")
    ax.set_title(f"Radial profile at core mid-height, z = {mh['z_m']:.2f} m")
    ax.legend(loc="upper right")
    ax.grid(axis="y", color=GRID, lw=0.5)
    fig.savefig(FIGURES_DIR / "fig5_temperature.png")
    plt.close(fig)


def fig_thermal_stress(cfg, res):
    g = pv.read(RESULTS_DIR / "thermo_mech_vm.vtu")
    s = cut(g)
    a = np.asarray(s.point_data["vm_thermal"]) / 1e6
    b = np.asarray(s.point_data["vm_thermal_sliding"]) / 1e6
    vmax = np.percentile(np.concatenate([a, b]), 99.9)
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 4.4), sharey=True)
    for ax, v, t in zip(axes, [a, b], ["fixed flange face (baseline)", "sliding flange face (sensitivity)"]):
        tpc = field_plot(ax, s, v, vmin=0, vmax=vmax)
        ax.set_xlim(0.0, 0.16)
        ax.set_ylim(0.30, 0.56)
        outline(ax, cfg, lw=0.5)
        ax.set_xlabel("x (m)")
        ax.set_title(t)
    axes[0].set_ylabel("z (m)")
    colourbar(fig, tpc, axes[1], "von Mises stress (MPa)")
    fig.suptitle(f"Thermal stress near the flange (T_ref = {cfg['loads']['T_reference']:.0f} C, no mechanical load)",
                 fontsize=10)
    fig.savefig(FIGURES_DIR / "fig6_thermal_stress_flange.png")
    plt.close(fig)


def fig_convergence(cfg, res):
    rows = res["convergence"]["rows"]
    n = np.array([r["n_elements"] for r in rows])
    fig, axes = plt.subplots(1, 3, figsize=(7.6, 2.8))
    panels = [
        ([np.array([r["top_deflection_m"] for r in rows]) * 1e3], ["top deflection"], "S1 top deflection (mm)"),
        ([np.array([r["max_von_mises_Pa"] for r in rows]) / 1e6,
          np.array([r["max_von_mises_away_from_corners_Pa"] for r in rows]) / 1e6],
         ["maximum (at a corner)", "maximum away from corners"], "S1 von Mises stress (MPa)"),
        ([np.array([r["f1_Hz"] for r in rows])], ["f1"], "first natural frequency (Hz)"),
    ]
    for ax, (ys, labels, yl) in zip(axes, panels):
        for y, lab, c in zip(ys, labels, [C1, C2]):
            ax.plot(n, y, color=c, lw=2, marker="o", ms=5, label=lab)
        ax.set_xscale("log")
        ax.set_xlabel("number of elements")
        ax.set_ylabel(yl)
        ax.grid(axis="y", color=GRID, lw=0.5)
        if len(ys) > 1:
            ax.legend(fontsize=7, loc="best")
        ax.set_xticks(n)
        ax.set_xticklabels([f"{v / 1e3:.0f}k" for v in n])
        ax.minorticks_off()
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "fig7_mesh_convergence.png")
    plt.close(fig)


def fig_deflection_profile(cfg, res):
    p = res["static"]["S1_cantilever"]["profile"]
    d = derived_geometry(cfg)
    z = np.array(p["z_m"])
    s_pos = z - d["z_f0"]
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    ax.plot(z, np.array(p["ux_m"]) * 1e3, color=C1, lw=2, label="FEA (conductor mid-wall, neutral plane)")
    ax.plot(z, cantilever_profile(cfg, s_pos) * 1e3, color=C2, lw=2, ls="--", label="two-segment beam")
    ax.axvline(d["L_core"], color=GRID, lw=0.8, zorder=0)
    ax.text(d["L_core"] - 0.01, 0.05, "core top", ha="right", fontsize=7, color=INK2)
    ax.set_xlabel("height z (m)")
    ax.set_ylabel("deflection u_x (mm)")
    ax.set_title("S1 deflection along the bushing above the flange")
    ax.legend(loc="upper left", fontsize=7)
    ax.grid(axis="y", color=GRID, lw=0.5)
    fig.savefig(FIGURES_DIR / "fig8_s1_deflection_profile.png")
    plt.close(fig)


def main():
    FIGURES_DIR.mkdir(exist_ok=True)
    cfg = load_config()
    res = read_results()
    for f in [fig_geometry, fig_s1_displacement, fig_s1_stress, fig_modes, fig_thermal, fig_thermal_stress,
              fig_deflection_profile]:
        f(cfg, res)
        print("written", f.__name__)
    if "convergence" in res:
        fig_convergence(cfg, res)
        print("written fig_convergence")


if __name__ == "__main__":
    main()
