"""Stage 5: steady-state temperature from Joule heating in the conductor.

Heat source: q = I^2 rho_e / A^2 in the copper, resistivity at 20 C (no temperature dependence).
Boundary conditions:
  - convection h to T_ambient on the outer surfaces of core and flange (including the flange
    bottom face; the tank cover is not modelled),
  - conductor bore insulated,
  - exposed conductor surfaces outside the core (terminal ends) insulated, a conservative choice:
    their heat must leave through the core and flange.
Perfect thermal contact between bodies (conforming mesh, continuous temperature).
"""
import numpy as np
import pyvista as pv
import ufl
from dolfinx import fem, plot
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI

from common import (BODIES, RESULTS_DIR, TAGS, derived_geometry, eval_points, load_config, location,
                    material_field, update_results)
from hand_calcs import joule_heating, radial_temperature


def solve(mesh_data, cfg):
    msh, ct, ft = mesh_data.mesh, mesh_data.cell_tags, mesh_data.facet_tags
    ld = cfg["loads"]
    V = fem.functionspace(msh, ("Lagrange", cfg["mesh"]["element_order"]))
    k = material_field(msh, ct, cfg, "k")
    dx = ufl.Measure("dx", msh, subdomain_data=ct)
    ds = ufl.Measure("ds", msh, subdomain_data=ft)
    q = fem.Constant(msh, joule_heating(cfg)["q_W_per_m3"])
    h = fem.Constant(msh, ld["h_conv"])
    Tinf = fem.Constant(msh, ld["T_ambient"])
    T, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    ds_conv = ds(TAGS["convection"]) + ds(TAGS["fixed"])
    a = k * ufl.inner(ufl.grad(T), ufl.grad(v)) * dx + h * T * v * ds_conv
    L = q * v * dx(TAGS["conductor"]) + h * Tinf * v * ds_conv
    problem = LinearProblem(a, L, petsc_options={"ksp_type": "preonly", "pc_type": "cholesky",
                                                 "pc_factor_mat_solver_type": "mumps"},
                            petsc_options_prefix="thermal_")
    Th = problem.solve()
    if isinstance(Th, tuple):
        Th = Th[0]
    Th.name = "T"
    # Energy balance: generated heat vs heat leaving by convection
    Q_gen = fem.assemble_scalar(fem.form(q * dx(TAGS["conductor"])))
    Q_out = fem.assemble_scalar(fem.form(h * (Th - Tinf) * ds_conv))
    return Th, float(msh.comm.allreduce(Q_gen, op=MPI.SUM)), float(msh.comm.allreduce(Q_out, op=MPI.SUM))


def run(mesh_data, cfg, save=True):
    d = derived_geometry(cfg)
    Th, Q_gen, Q_out = solve(mesh_data, cfg)
    V = Th.function_space
    X = V.tabulate_dof_coordinates()
    i = int(np.argmax(Th.x.array))
    out = {"n_dofs": int(V.dofmap.index_map.size_global),
           "heat_generated_W": Q_gen, "heat_convected_W": Q_out,
           "energy_balance_error_pct": 100 * (Q_out - Q_gen) / Q_gen,
           "T_max_C": float(Th.x.array[i]), "T_max_location": location(X[i]),
           "T_min_C": float(Th.x.array.min())}
    # Max temperature per body
    ct = mesh_data.cell_tags
    out["T_max_per_body_C"] = {}
    for b in BODIES:
        dofs = np.unique(V.dofmap.list[ct.find(TAGS[b])].flatten())
        out["T_max_per_body_C"][b] = float(Th.x.array[dofs].max())

    # Axial line: the bore is empty, so the "axis" line is taken at the conductor mid-wall radius
    pp = cfg["postprocessing"]
    r_mid = 0.5 * (d["r_ci"] + d["r_co"])
    z = np.linspace(d["z_c0"], d["z_c1"], pp["axial_line_points"])
    Tz = eval_points(Th, np.column_stack([np.full_like(z, r_mid), np.zeros_like(z), z]))[:, 0]
    z_s = np.linspace(0.0, d["L_core"], pp["axial_line_points"])
    Ts = eval_points(Th, np.column_stack([np.full_like(z_s, d["r_k"] * (1 - 1e-6)), np.zeros_like(z_s), z_s]))[:, 0]
    out["axial"] = {"line": "conductor mid-wall radius, y = 0", "z_m": z.tolist(), "T_C": Tz.tolist(),
                    "core_surface_z_m": z_s.tolist(), "core_surface_T_C": Ts.tolist()}

    # Radial line at core mid-height, compared with the 1D analytical solution
    zm = d["L_core"] / 2
    r = np.linspace(d["r_ci"] * (1 + 1e-6), d["r_k"] * (1 - 1e-6), pp["radial_line_points"])
    Tr = eval_points(Th, np.column_stack([r, np.zeros_like(r), np.full_like(r, zm)]))[:, 0]
    Ta = radial_temperature(cfg, r)
    T_co = float(eval_points(Th, np.array([[d["r_co"], 0, zm]]))[0, 0])
    T_k = float(Tr[-1])
    T_co_a, T_k_a = radial_temperature(cfg, [d["r_co"], d["r_k"]])
    out["mid_height"] = {"z_m": zm, "r_m": r.tolist(), "T_fea_C": Tr.tolist(), "T_analytical_C": Ta.tolist(),
                         "T_conductor_surface_C": T_co, "T_core_surface_C": T_k, "dT_epoxy_K": T_co - T_k,
                         "max_abs_difference_K": float(np.max(np.abs(Tr - Ta))),
                         "T_conductor_surface_analytical_C": float(T_co_a),
                         "T_core_surface_analytical_C": float(T_k_a)}
    # Local radial heat flow at mid-height (from the convective flux at the core surface) compared with
    # the generation per unit length; the 1D conduction law is then checked with this local heat flow.
    m = cfg["materials"]
    Qp_local = 2 * np.pi * d["r_k"] * cfg["loads"]["h_conv"] * (T_k - cfg["loads"]["T_ambient"])
    out["mid_height"]["Q_radial_local_W_per_m"] = float(Qp_local)
    out["mid_height"]["Q_generated_per_length_W_per_m"] = joule_heating(cfg)["Q_per_length_W_per_m"]
    out["mid_height"]["T_analytical_local_heat_flow_C"] = radial_temperature(
        cfg, r, q_scale=Qp_local / joule_heating(cfg)["Q_per_length_W_per_m"]).tolist()
    out["mid_height"]["dT_epoxy_1D_with_local_heat_flow_K"] = float(
        Qp_local / (2 * np.pi * m["epoxy"]["k"]) * np.log(d["r_k"] / d["r_co"]))
    print(f"Q generated {Q_gen:.2f} W, convected {Q_out:.2f} W; T max {out['T_max_C']:.2f} C at {out['T_max_location']}")
    print(f"mid-height: conductor surface {T_co:.2f} C (1D {T_co_a:.2f}), core surface {T_k:.2f} C (1D {T_k_a:.2f})")
    if save:
        grid = pv.UnstructuredGrid(*plot.vtk_mesh(V))
        grid.point_data["T"] = Th.x.array
        grid.save(str(RESULTS_DIR / "thermal_T.vtu"))
        update_results("thermal", out)
    return Th, out


if __name__ == "__main__":
    from geometry import create_mesh
    cfg = load_config()
    data, _ = create_mesh(cfg)
    run(data, cfg)
