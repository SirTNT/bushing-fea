"""Stage 6: thermal stress from the Stage 5 temperature field.

Thermal strain alpha (T - T_ref), T_ref = 20 C (stress-free), same fixed support, no mechanical load.
The combination with S1 uses superposition of displacement fields (the model is linear).
"""
import numpy as np
import ufl
from dolfinx import fem
from mpi4py import MPI

from common import (BODIES, TAGS, body_max, derived_geometry, load_config, location, singular_rings,
                    surface_average, update_results)


def interface_max(vm, ct, body, r0, z_lo, z_hi, rings=None, excl=0.0, tol=1e-6):
    """Maximum of a DG1 field on one side (body) of a cylindrical interface r = r0, z_lo <= z <= z_hi."""
    W = vm.function_space
    dofs = np.unique(W.dofmap.list[ct.find(TAGS[body])].flatten())
    X = W.tabulate_dof_coordinates()[dofs]
    r = np.hypot(X[:, 0], X[:, 1])
    sel = (np.abs(r - r0) < 1e-5) & (X[:, 2] > z_lo - tol) & (X[:, 2] < z_hi + tol)
    if rings:
        for (rr, zz) in rings:
            sel &= np.hypot(r - rr, X[:, 2] - zz) > excl
    vals = vm.x.array[dofs][sel]
    i = int(np.argmax(vals))
    return float(vals[i]), location(X[sel][i])


def mean_over(model, f, body, z_lo, z_hi):
    x = ufl.SpatialCoordinate(model.msh)
    ind = ufl.conditional(ufl.And(ufl.ge(x[2], z_lo), ufl.le(x[2], z_hi)), 1.0, 0.0)
    num = fem.assemble_scalar(fem.form(f * ind * model.dx(TAGS[body])))
    den = fem.assemble_scalar(fem.form(ind * model.dx(TAGS[body])))
    return float(model.msh.comm.allreduce(num, op=MPI.SUM) / model.msh.comm.allreduce(den, op=MPI.SUM))


def stress_summary(model, u, dT, cfg):
    d = derived_geometry(cfg)
    rings = singular_rings(cfg)
    excl = cfg["postprocessing"]["corner_exclusion"]
    vm = model.von_mises(u, dT)
    out = {"von_mises": {}, "interfaces": {}}
    for b in BODIES:
        s, loc = body_max(vm, model.ct, b)
        s2, loc2 = body_max(vm, model.ct, b, rings, excl)
        out["von_mises"][b] = {"max_Pa": s, "max_location": loc, "max_away_from_corners_Pa": s2, "away_location": loc2}
    ifs = {"flange_to_core": (d["r_k"], d["z_f0"], d["z_f1"], ["core", "flange"]),
           "conductor_to_core": (d["r_co"], 0.0, d["L_core"], ["conductor", "core"])}
    for name, (r0, z0, z1, bodies) in ifs.items():
        out["interfaces"][name] = {}
        for b in bodies:
            s, loc = interface_max(vm, model.ct, b, r0, z0, z1)
            s2, loc2 = interface_max(vm, model.ct, b, r0, z0, z1, rings, excl)
            out["interfaces"][name][b] = {"max_Pa": s, "max_location": loc,
                                          "max_away_from_corners_Pa": s2, "away_location": loc2}
    return out, vm


def run(mesh_data, cfg, T, model=None, u_s1=None, save=True):
    from static import ElasticModel, save_fields
    model = model or ElasticModel(mesh_data, cfg)
    d = derived_geometry(cfg)
    ld = cfg["loads"]
    m = cfg["materials"]
    dT = T - ld["T_reference"]
    u_th = model.solve(dT=dT)
    out = {"T_reference_C": ld["T_reference"]}
    th, vm_th = stress_summary(model, u_th, dT, cfg)
    th["top_displacement_z_m"] = surface_average(u_th[2], model.ds(TAGS["top_end"]), model.msh)
    th["top_displacement_x_m"] = surface_average(u_th[0], model.ds(TAGS["top_end"]), model.msh)
    out["thermal_only"] = th

    # Combination with S1 by superposition of displacements
    if u_s1 is None:
        A_end = np.pi * (d["r_co"] ** 2 - d["r_ci"] ** 2)
        u_s1 = model.solve(traction_top=(ld["cantilever_force"] / A_end, 0, 0))
    u_c = fem.Function(model.V)
    u_c.x.array[:] = u_th.x.array + u_s1.x.array
    comb, vm_c = stress_summary(model, u_c, dT, cfg)
    comb["top_deflection_x_m"] = surface_average(u_c[0], model.ds(TAGS["top_end"]), model.msh)
    out["thermal_plus_S1"] = comb

    # Sanity check: axial elongation from the fixed plane to the conductor top
    a = d["L_core"] - d["z_f0"]
    b = d["z_c1"] - d["L_core"]
    L = a + b
    dT_cond = mean_over(model, dT, "conductor", d["z_f0"], d["z_c1"])
    dT_cond_seg = mean_over(model, dT, "conductor", d["z_f0"], d["L_core"])
    dT_core_seg = mean_over(model, dT, "core", d["z_f0"], d["L_core"])
    dT_over = mean_over(model, dT, "conductor", d["L_core"], d["z_c1"])
    free_cu = m["copper"]["alpha"] * L * dT_cond
    EA_c = m["copper"]["E"] * np.pi * (d["r_co"] ** 2 - d["r_ci"] ** 2)
    EA_k = m["epoxy"]["E"] * np.pi * (d["r_k"] ** 2 - d["r_co"] ** 2)
    strain_seg = (EA_c * m["copper"]["alpha"] * dT_cond_seg + EA_k * m["epoxy"]["alpha"] * dT_core_seg) / (EA_c + EA_k)
    composite = strain_seg * a + m["copper"]["alpha"] * dT_over * b
    out["elongation_check"] = {
        "L_m": L, "mean_dT_conductor_K": dT_cond, "mean_dT_conductor_core_segment_K": dT_cond_seg,
        "mean_dT_core_segment_K": dT_core_seg, "mean_dT_conductor_overhang_K": dT_over,
        "free_copper_elongation_m": free_cu,
        "free_core_governed_elongation_m": composite,
        "fea_top_uz_m": th["top_displacement_z_m"],
        "note": "composite estimate: EA-weighted free strain of the bonded core segment plus free copper overhang",
    }
    # Sensitivity: same thermal load with a sliding support (axial restraint only at the flange face)
    sl = ElasticModel(mesh_data, cfg, support="sliding")
    u_sl = sl.solve(dT=dT)
    sl_sum, vm_sl = stress_summary(sl, u_sl, dT, cfg)
    out["thermal_only_sliding_support"] = sl_sum
    print("sliding support: " + ", ".join(f"{k} {v['max_Pa'] / 1e6:.2f} MPa (away {v['max_away_from_corners_Pa'] / 1e6:.2f}) at z {v['away_location']['z']}"
                                          for k, v in sl_sum["von_mises"].items()))
    print(f"thermal only: " + ", ".join(f"{k} {v['max_Pa'] / 1e6:.2f} MPa (away {v['max_away_from_corners_Pa'] / 1e6:.2f})"
                                       for k, v in th["von_mises"].items()))
    print(f"thermal + S1: " + ", ".join(f"{k} {v['max_Pa'] / 1e6:.2f} MPa (away {v['max_away_from_corners_Pa'] / 1e6:.2f})"
                                       for k, v in comb["von_mises"].items()))
    print(f"top uz FEA {th['top_displacement_z_m'] * 1e3:.4f} mm, free copper {free_cu * 1e3:.4f} mm, composite {composite * 1e3:.4f} mm")
    if save:
        save_fields(model, "thermo_mech", {"u_thermal": u_th, "u_thermal_plus_S1": u_c},
                    {"vm_thermal": vm_th, "vm_thermal_plus_S1": vm_c, "vm_thermal_sliding": vm_sl})
        update_results("thermo_mech", out)
    return out


if __name__ == "__main__":
    from geometry import create_mesh
    import thermal
    cfg = load_config()
    data, _ = create_mesh(cfg)
    T, _ = thermal.run(data, cfg, save=False)
    run(data, cfg, T)
