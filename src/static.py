"""Stage 3: linear static analysis.

S1 cantilever: horizontal force F (x direction) as a uniform traction on the top conductor end face.
S2 self-weight: gravity (-z).
S3 static-equivalent seismic: horizontal body acceleration a_h (x direction) plus gravity.
    Simplified static-equivalent check, not a seismic qualification.
Support: flange bottom face fully fixed. All bodies perfectly bonded (conforming mesh).
"""
import numpy as np
import pyvista as pv
import ufl
from dolfinx import fem, plot
from dolfinx.fem.petsc import apply_lifting, assemble_matrix, assemble_vector, set_bc
from petsc4py import PETSc

from common import (BODIES, RESULTS_DIR, TAGS, body_max, derived_geometry, eps, load_config,
                    material_field, nodal_max_displacement, project_dg1, sigma, singular_rings,
                    surface_average, update_results, von_mises)


class ElasticModel:
    """Assembles the stiffness matrix once (LU factorisation with MUMPS) and solves several load cases."""

    def __init__(self, mesh_data, cfg, support="fixed"):
        self.cfg = cfg
        self.msh = mesh_data.mesh
        self.ct = mesh_data.cell_tags
        self.ft = mesh_data.facet_tags
        msh = self.msh
        self.V = fem.functionspace(msh, ("Lagrange", cfg["mesh"]["element_order"], (3,)))
        self.E = material_field(msh, self.ct, cfg, "E")
        self.nu = material_field(msh, self.ct, cfg, "nu")
        self.rho = material_field(msh, self.ct, cfg, "rho")
        self.alpha = material_field(msh, self.ct, cfg, "alpha")
        self.dx = ufl.Measure("dx", msh, subdomain_data=self.ct)
        self.ds = ufl.Measure("ds", msh, subdomain_data=self.ft)
        u, v = ufl.TrialFunction(self.V), ufl.TestFunction(self.V)
        self.v = v
        fixed = self.ft.find(TAGS["fixed"])
        self.support = support
        if support == "fixed":
            zero = np.zeros(3, dtype=PETSc.ScalarType)
            self.bc = fem.dirichletbc(zero, fem.locate_dofs_topological(self.V, 2, fixed), self.V)
            self.bcs = [self.bc]
        elif support == "sliding":
            # Sensitivity case: flange bottom face restrained axially only (free radial expansion);
            # four points on the outer edge remove the rigid-body motions in x, y and rotation about z.
            self.bcs = self._sliding_bcs(fixed)
        else:
            raise ValueError(support)
        self.a = fem.form(ufl.inner(sigma(u, self.E, self.nu), eps(v)) * self.dx)
        self.A = assemble_matrix(self.a, bcs=self.bcs)
        self.A.assemble()
        self.ksp = PETSc.KSP().create(msh.comm)
        self.ksp.setOperators(self.A)
        self.ksp.setType("preonly")
        pc = self.ksp.getPC()
        pc.setType("cholesky")
        pc.setFactorSolverType("mumps")
        self.ksp.setFromOptions()
        self.ndofs = self.V.dofmap.index_map.size_global * 3

    def _sliding_bcs(self, fixed):
        d = derived_geometry(self.cfg)
        bcs = []
        dz = fem.locate_dofs_topological(self.V.sub(2), 2, fixed)
        bcs.append(fem.dirichletbc(PETSc.ScalarType(0), dz, self.V.sub(2)))
        X = self.V.tabulate_dof_coordinates()
        rf, z0 = d["r_f"], d["z_f0"]
        # ux = 0 at (0, +-rf), uy = 0 at (+-rf, 0): no net translation or rotation, radial growth free
        for comp, pts in [(0, [(0, rf, z0), (0, -rf, z0)]), (1, [(rf, 0, z0), (-rf, 0, z0)])]:
            nodes = [int(np.argmin(np.linalg.norm(X - np.array(p), axis=1))) for p in pts]
            dofs = np.array([3 * n + comp for n in nodes], dtype=np.int32)
            bcs.append(fem.dirichletbc(PETSc.ScalarType(0), dofs, self.V.sub(comp)))
        return bcs

    def solve(self, body_force=(0, 0, 0), traction_top=(0, 0, 0), dT=None):
        msh, v = self.msh, self.v
        b = fem.Constant(msh, np.array(body_force, dtype=PETSc.ScalarType))
        t = fem.Constant(msh, np.array(traction_top, dtype=PETSc.ScalarType))
        L = self.rho * ufl.dot(b, v) * self.dx + ufl.dot(t, v) * self.ds(TAGS["top_end"])
        if dT is not None:
            # thermal pre-stress term: C : (alpha dT I)
            lmbda = self.E * self.nu / ((1 + self.nu) * (1 - 2 * self.nu))
            mu = self.E / (2 * (1 + self.nu))
            L += (3 * lmbda + 2 * mu) * self.alpha * dT * ufl.div(v) * self.dx
        Lf = fem.form(L)
        rhs = assemble_vector(Lf)
        apply_lifting(rhs, [self.a], bcs=[self.bcs])
        rhs.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
        set_bc(rhs, self.bcs)
        uh = fem.Function(self.V)
        self.ksp.solve(rhs, uh.x.petsc_vec)
        uh.x.scatter_forward()
        return uh

    def von_mises(self, u, dT=None):
        s = sigma(u, self.E, self.nu, dT=dT, alpha=self.alpha if dT is not None else None)
        return project_dg1(von_mises(s), self.msh)

    def summarise(self, u, dT=None):
        cfg = self.cfg
        vm = self.von_mises(u, dT)
        rings = singular_rings(cfg)
        excl = cfg["postprocessing"]["corner_exclusion"]
        umax, uloc = nodal_max_displacement(u)
        out = {
            "max_displacement_m": umax,
            "max_displacement_location": uloc,
            "top_deflection_x_m": surface_average(u[0], self.ds(TAGS["top_end"]), self.msh),
            "top_displacement_z_m": surface_average(u[2], self.ds(TAGS["top_end"]), self.msh),
            "von_mises": {},
        }
        for body in BODIES:
            s_all, loc_all = body_max(vm, self.ct, body)
            s_ex, loc_ex = body_max(vm, self.ct, body, rings, excl)
            out["von_mises"][body] = {
                "max_Pa": s_all, "max_location": loc_all,
                "max_away_from_corners_Pa": s_ex, "away_location": loc_ex,
            }
        return out, vm


def save_fields(model, name, fields_u, fields_vm):
    """Save displacement (P2 grid) and von Mises (DG1 grid) for plotting with PyVista."""
    RESULTS_DIR.mkdir(exist_ok=True)
    grid = pv.UnstructuredGrid(*plot.vtk_mesh(model.V))
    for key, u in fields_u.items():
        grid.point_data[key] = u.x.array.reshape(-1, 3)
    grid.save(str(RESULTS_DIR / f"{name}_u.vtu"))
    if fields_vm:
        W = next(iter(fields_vm.values())).function_space
        g2 = pv.UnstructuredGrid(*plot.vtk_mesh(W))
        for key, f in fields_vm.items():
            g2.point_data[key] = f.x.array
        g2.cell_data["body"] = _cell_body(model)
        g2.save(str(RESULTS_DIR / f"{name}_vm.vtu"))


def _cell_body(model):
    n = model.msh.topology.index_map(3).size_local
    body = np.zeros(n, dtype=np.int32)
    for name in BODIES:
        body[model.ct.find(TAGS[name])] = TAGS[name]
    return body


def run(mesh_data, cfg, model=None, save=True):
    model = model or ElasticModel(mesh_data, cfg)
    d = derived_geometry(cfg)
    ld = cfg["loads"]
    A_end = np.pi * (d["r_co"] ** 2 - d["r_ci"] ** 2)
    g = ld["g"]
    cases = {
        "S1_cantilever": dict(traction_top=(ld["cantilever_force"] / A_end, 0, 0)),
        "S2_self_weight": dict(body_force=(0, 0, -g)),
        "S3_seismic_static_equivalent": dict(body_force=(ld["seismic_horizontal_g"] * g, 0, -g)),
    }
    results = {"n_dofs": model.ndofs}
    us, vms = {}, {}
    for name, kw in cases.items():
        u = model.solve(**kw)
        res, vm = model.summarise(u)
        results[name] = res
        us[name], vms[name] = u, vm
        print(f"{name}: top deflection x = {res['top_deflection_x_m'] * 1e3:.4f} mm, "
              f"max |u| = {res['max_displacement_m'] * 1e3:.4f} mm, "
              + ", ".join(f"{b} vm max {v['max_Pa'] / 1e6:.2f} MPa (away {v['max_away_from_corners_Pa'] / 1e6:.2f})"
                          for b, v in res["von_mises"].items()))
    # S1 deflection along the conductor mid-wall line in the neutral plane (x = 0), from the fixed plane up
    from common import eval_points
    r_mid = 0.5 * (d["r_ci"] + d["r_co"])
    z = np.linspace(d["z_f0"], d["z_c1"], cfg["postprocessing"]["axial_line_points"])
    pts = np.column_stack([np.zeros_like(z), np.full_like(z, r_mid), z])
    ux = eval_points(us["S1_cantilever"], pts)[:, 0]
    results["S1_cantilever"]["profile"] = {"z_m": z.tolist(), "ux_m": ux.tolist(), "line": "x = 0, y = mid-wall radius of conductor"}
    k_top = int(np.argmin(abs(z - d["L_core"])))
    results["S1_cantilever"]["ux_at_core_top_m"] = float(ux[k_top])
    # Total mass (check of self-weight)
    mass = fem.assemble_scalar(fem.form(model.rho * model.dx))
    results["total_mass_kg"] = float(mass)
    if save:
        save_fields(model, "static", us, vms)
        update_results("static", results)
    return model, us, results


if __name__ == "__main__":
    from geometry import create_mesh
    cfg = load_config()
    data, stats = create_mesh(cfg)
    run(data, cfg)
