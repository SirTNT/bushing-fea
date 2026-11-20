"""Stage 4: modal analysis. First natural frequencies and mode shapes with SLEPc (shift-invert about 0).

Same fixed support as the static cases (flange bottom face). Effective modal mass in x, y, z is
reported to identify each mode type.
"""
import numpy as np
import pyvista as pv
import ufl
from dolfinx import fem, plot
from dolfinx.fem.petsc import assemble_matrix
from petsc4py import PETSc
from slepc4py import SLEPc

from common import RESULTS_DIR, location, load_config, update_results


def solve_modes(model, n_modes):
    V, msh = model.V, model.msh
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    m_form = fem.form(model.rho * ufl.inner(u, v) * model.dx)
    M = assemble_matrix(m_form, bcs=model.bcs, diag=0.0)   # constrained dofs: K diag 1, M diag 0
    M.assemble()

    eps = SLEPc.EPS().create(msh.comm)
    eps.setOperators(model.A, M)
    eps.setProblemType(SLEPc.EPS.ProblemType.GHEP)
    eps.setType(SLEPc.EPS.Type.KRYLOVSCHUR)
    eps.setDimensions(nev=n_modes)
    eps.setTolerances(tol=1e-10, max_it=200)
    st = eps.getST()
    st.setType(SLEPc.ST.Type.SINVERT)
    st.setShift(0.0)
    eps.setTarget(0.0)
    eps.setWhichEigenpairs(SLEPc.EPS.Which.TARGET_MAGNITUDE)
    ksp = st.getKSP()
    ksp.setType("preonly")
    ksp.getPC().setType("cholesky")
    ksp.getPC().setFactorSolverType("mumps")
    eps.setFromOptions()
    eps.solve()
    nconv = eps.getConverged()
    assert nconv >= n_modes, f"only {nconv} modes converged"

    modes = []
    for i in range(nconv):
        lam = eps.getEigenvalue(i).real
        phi = fem.Function(V)
        eps.getEigenvector(i, phi.x.petsc_vec)
        phi.x.scatter_forward()
        err = eps.computeError(i, SLEPc.EPS.ErrorType.RELATIVE)
        modes.append((lam, phi, err))
    modes.sort(key=lambda m: m[0])
    return modes[:n_modes]


def align_degenerate_pairs(model, modes, rel_tol=1e-3):
    """Bending modes of an axisymmetric body come in pairs with (almost) equal frequency; the solver
    returns an arbitrary orthogonal combination. Rotate each pair so that the first mode moves in x
    and the second in y. Eigenvectors are mass-normalised, so the rotation keeps them orthonormal."""
    def part(phi, k):
        return fem.assemble_scalar(fem.form(model.rho * phi[k] * model.dx))
    i = 0
    while i < len(modes) - 1:
        (l1, p1, e1), (l2, p2, e2) = modes[i], modes[i + 1]
        if abs(l2 - l1) / l1 < 2 * rel_tol:
            a, b = part(p1, 0), part(p2, 0)
            n = np.hypot(a, b)
            if n > 0:
                c, s = a / n, b / n
                q1, q2 = fem.Function(model.V), fem.Function(model.V)
                q1.x.array[:] = c * p1.x.array + s * p2.x.array
                q2.x.array[:] = -s * p1.x.array + c * p2.x.array
                modes[i], modes[i + 1] = (l1, q1, e1), (l2, q2, e2)
            i += 2
        else:
            i += 1
    return modes


def mode_info(model, lam, phi, total_mass):
    msh = model.msh
    f = float(np.sqrt(lam) / (2 * np.pi))
    mod_mass = fem.assemble_scalar(fem.form(model.rho * ufl.inner(phi, phi) * model.dx))
    eff = {}
    for k, name in enumerate("xyz"):
        L = fem.assemble_scalar(fem.form(model.rho * phi[k] * model.dx))
        eff[name] = float(L ** 2 / mod_mass / total_mass)
    a = phi.x.array.reshape(-1, 3)
    X = model.V.tabulate_dof_coordinates()[: a.shape[0]]
    imax = int(np.argmax(np.linalg.norm(a, axis=1)))
    # crude classification from effective mass
    if max(eff.values()) < 0.01:
        kind = "torsion or local mode (low effective mass)"
    else:
        kdir = max(eff, key=eff.get)
        kind = f"bending in {kdir}-z plane" if kdir in "xy" else "axial"
    # mode below or above the flange: compare max amplitude location
    part = "upper part (above flange)" if X[imax, 2] > model_flange_z(model) else "lower part (below flange)"
    return {"frequency_Hz": f, "eff_mass_fraction": eff, "type": kind,
            "max_amplitude_location": location(X[imax]), "dominant_region": part,
            "modal_mass_kg": float(mod_mass)}


def model_flange_z(model):
    from common import derived_geometry
    return derived_geometry(model.cfg)["z_f0"]


def run(mesh_data, cfg, model=None, save=True):
    from static import ElasticModel
    model = model or ElasticModel(mesh_data, cfg)
    n = cfg["modal"]["n_modes"]
    modes = align_degenerate_pairs(model, solve_modes(model, n))
    total_mass = float(fem.assemble_scalar(fem.form(model.rho * model.dx)))
    out = {"n_dofs": model.ndofs, "total_mass_kg": total_mass, "modes": []}
    for i, (lam, phi, err) in enumerate(modes):
        info = mode_info(model, lam, phi, total_mass)
        info["relative_residual"] = float(err)
        out["modes"].append(info)
        print(f"mode {i + 1}: {info['frequency_Hz']:.2f} Hz  {info['type']:45s} {info['dominant_region']}  "
              f"eff mass x/y/z = {info['eff_mass_fraction']['x']:.3f}/{info['eff_mass_fraction']['y']:.3f}/"
              f"{info['eff_mass_fraction']['z']:.3f}")
    f1, f2 = out["modes"][0]["frequency_Hz"], out["modes"][1]["frequency_Hz"]
    out["f1_f2_relative_difference_pct"] = 100 * abs(f2 - f1) / f1
    if save:
        grid = pv.UnstructuredGrid(*plot.vtk_mesh(model.V))
        for i, (_, phi, _) in enumerate(modes):
            a = phi.x.array.reshape(-1, 3)
            grid.point_data[f"mode{i + 1}"] = a / np.abs(a).max()
        grid.save(str(RESULTS_DIR / "modal_u.vtu"))
        update_results("modal", out)
    return out


if __name__ == "__main__":
    from geometry import create_mesh
    cfg = load_config()
    data, _ = create_mesh(cfg)
    run(data, cfg)
