"""Shared helpers: configuration, tags, material fields, stress, point evaluation, results file."""
import json
from pathlib import Path

import numpy as np
import ufl
import yaml
from dolfinx import fem, geometry
from mpi4py import MPI

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = ROOT / "figures"
RESULTS_FILE = RESULTS_DIR / "results.json"

# Physical group tags (shared by geometry.py and the solvers)
TAGS = {
    "conductor": 1, "core": 2, "flange": 3,                    # volumes
    "fixed": 11,               # flange bottom face (fixed support; also convective in the thermal case)
    "top_end": 12,             # top conductor end face (terminal load)
    "convection": 13,          # outer surfaces of core and flange, except the flange bottom face
    "bore": 14,                # conductor bore (insulated)
    "conductor_exposed": 15,   # conductor outer surface outside the core and the bottom end face
    "if_conductor_core": 21,   # conductor to core interface
    "if_core_flange": 22,      # core to flange interface
}
BODIES = {"conductor": "copper", "core": "epoxy", "flange": "aluminium"}


def load_config():
    with open(ROOT / "config.yaml") as f:
        cfg = yaml.safe_load(f)
    for mat, props in cfg["materials"].items():
        for k, v in props.items():
            if not isinstance(v, (int, float)):
                raise TypeError(f"materials.{mat}.{k} is not a number: {v!r} (write exponents as e+9)")
    return cfg


def derived_geometry(cfg):
    """Radii and axial positions derived from config.yaml (z = 0 at the core bottom)."""
    g = cfg["geometry"]
    d = {
        "r_ci": g["conductor"]["inner_diameter"] / 2,
        "r_co": g["conductor"]["outer_diameter"] / 2,
        "r_k": g["core"]["outer_diameter"] / 2,
        "r_f": g["flange"]["outer_diameter"] / 2,
        "L_core": g["core"]["length"],
        "L_cond": g["conductor"]["length"],
        "z_f0": g["flange"]["z_bottom"],
        "t_f": g["flange"]["thickness"],
    }
    d["z_c0"] = d["L_core"] / 2 - d["L_cond"] / 2       # conductor bottom
    d["z_c1"] = d["L_core"] / 2 + d["L_cond"] / 2       # conductor top
    d["z_f1"] = d["z_f0"] + d["t_f"]
    return d


def singular_rings(cfg):
    """(r, z) of ring edges where the elastic stress field is singular (re-entrant or
    bi-material corners and the edges of the fixed face)."""
    d = derived_geometry(cfg)
    return [
        (d["r_k"], d["z_f0"]),   # flange bottom inner edge / fixed face inner edge
        (d["r_k"], d["z_f1"]),   # flange top inner edge (re-entrant)
        (d["r_f"], d["z_f0"]),   # fixed face outer edge
        (d["r_co"], d["L_core"]),  # conductor exits the core top
        (d["r_co"], 0.0),        # conductor exits the core bottom
    ]


def material_field(msh, cell_tags, cfg, prop):
    """Piecewise-constant (DG0) field of a material property, assigned by body."""
    Q = fem.functionspace(msh, ("DG", 0))
    f = fem.Function(Q)
    for body, mat in BODIES.items():
        cells = cell_tags.find(TAGS[body])
        f.x.array[cells] = cfg["materials"][mat][prop]
    f.x.scatter_forward()
    return f


def lame(E, nu):
    return E * nu / ((1 + nu) * (1 - 2 * nu)), E / (2 * (1 + nu))


def eps(u):
    return ufl.sym(ufl.grad(u))


def sigma(u, E, nu, dT=None, alpha=None):
    lmbda, mu = lame(E, nu)
    e = eps(u)
    if dT is not None:
        e = e - alpha * dT * ufl.Identity(3)
    return 2 * mu * e + lmbda * ufl.tr(e) * ufl.Identity(3)


def von_mises(s):
    dev = s - ufl.tr(s) / 3 * ufl.Identity(3)
    return ufl.sqrt(1.5 * ufl.inner(dev, dev))


def project_dg1(expr, msh):
    """Interpolate a scalar UFL expression into discontinuous P1 (per-body values are kept)."""
    W = fem.functionspace(msh, ("DG", 1))
    f = fem.Function(W)
    f.interpolate(fem.Expression(expr, W.element.interpolation_points))
    return f


def body_max(f, cell_tags, body, exclude_rings=None, exclusion=0.0):
    """Maximum of a DG field in one body and its location; optionally skip points close
    to singular ring edges."""
    W = f.function_space
    cells = cell_tags.find(TAGS[body])
    dofs = np.unique(W.dofmap.list[cells].flatten())
    X = W.tabulate_dof_coordinates()[dofs]
    vals = f.x.array[dofs]
    if exclude_rings:
        r = np.hypot(X[:, 0], X[:, 1])
        keep = np.ones(len(dofs), dtype=bool)
        for (r0, z0) in exclude_rings:
            keep &= np.hypot(r - r0, X[:, 2] - z0) > exclusion
        X, vals = X[keep], vals[keep]
    i = int(np.argmax(vals))
    return float(vals[i]), location(X[i])


def location(x):
    x = np.asarray(x, dtype=float)
    return {"x": round(float(x[0]), 4), "y": round(float(x[1]), 4), "z": round(float(x[2]), 4),
            "r": round(float(np.hypot(x[0], x[1])), 4)}


def eval_points(u, points):
    """Evaluate a Function at points (N, 3); NaN where a point is outside the mesh."""
    msh = u.function_space.mesh
    points = np.asarray(points, dtype=float)
    tree = geometry.bb_tree(msh, msh.topology.dim)
    cand = geometry.compute_collisions_points(tree, points)
    coll = geometry.compute_colliding_cells(msh, cand, points)
    out = np.full((len(points), int(np.prod(u.ufl_shape)) or 1), np.nan)
    idx, cells = [], []
    for i in range(len(points)):
        c = coll.links(i)
        if len(c) > 0:
            idx.append(i)
            cells.append(c[0])
    if idx:
        out[idx] = u.eval(points[idx], np.array(cells, dtype=np.int32))
    return out


def surface_average(expr, ds_tag, msh):
    num = fem.assemble_scalar(fem.form(expr * ds_tag))
    den = fem.assemble_scalar(fem.form(fem.Constant(msh, 1.0) * ds_tag))
    return float(msh.comm.allreduce(num, op=MPI.SUM) / msh.comm.allreduce(den, op=MPI.SUM))


def nodal_max_displacement(u):
    """Largest displacement magnitude at a node of the displacement field and its location."""
    V = u.function_space
    X = V.tabulate_dof_coordinates()
    a = u.x.array.reshape(-1, 3)
    n = a.shape[0]
    mag = np.linalg.norm(a, axis=1)
    i = int(np.argmax(mag))
    return float(mag[i]), location(X[:n][i])


def update_results(section, data):
    """Merge one section into results/results.json."""
    RESULTS_DIR.mkdir(exist_ok=True)
    res = json.loads(RESULTS_FILE.read_text()) if RESULTS_FILE.exists() else {}
    res[section] = data
    RESULTS_FILE.write_text(json.dumps(res, indent=2))


def read_results():
    return json.loads(RESULTS_FILE.read_text())


def pct_change(new, ref):
    return 100.0 * (new - ref) / ref
