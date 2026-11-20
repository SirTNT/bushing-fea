"""Stage 2: simplified dry-type condenser bushing in Gmsh (OpenCASCADE kernel).

Three coaxial bodies along z: copper conductor tube, epoxy core, aluminium flange.
The bodies are fragmented so the mesh is conforming at the bonded interfaces.
"""
import sys

import gmsh
import numpy as np
from mpi4py import MPI

from common import RESULTS_DIR, TAGS, derived_geometry, load_config, update_results

TOL = 1e-6


def _classify_surface(tag, d):
    """Return the physical group name of an exterior or interface surface from its bounding box."""
    x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(2, tag)
    r = max(abs(x0), abs(x1), abs(y0), abs(y1))
    flat = abs(z1 - z0) < TOL

    def near(a, b):
        return abs(a - b) < 1e-4

    if flat:
        z = z0
        if near(z, d["z_f0"]) and near(r, d["r_f"]):
            return "fixed"
        if near(z, d["z_c1"]):
            return "top_end"
        if near(z, d["z_c0"]):
            return "conductor_exposed"
        return "convection"          # core bottom/top annulus, flange top annulus
    if near(r, d["r_ci"]):
        return "bore"
    if near(r, d["r_co"]):
        if z0 > -TOL and z1 < d["L_core"] + TOL:
            return "if_conductor_core"
        return "conductor_exposed"
    if near(r, d["r_k"]):
        if z0 > d["z_f0"] - TOL and z1 < d["z_f1"] + TOL:
            return "if_core_flange"
        return "convection"
    if near(r, d["r_f"]):
        return "convection"
    raise RuntimeError(f"Unclassified surface {tag}: bbox {x0, y0, z0, x1, y1, z1}")


def build_model(cfg, h=None):
    """Create geometry, physical groups and a quadratic tetrahedral mesh in the active gmsh model.

    Returns a dictionary with mesh statistics.
    """
    d = derived_geometry(cfg)
    mc = cfg["mesh"]
    h = mc["h"] if h is None else h
    occ = gmsh.model.occ

    # Conductor tube
    c_out = occ.addCylinder(0, 0, d["z_c0"], 0, 0, d["L_cond"], d["r_co"])
    c_in = occ.addCylinder(0, 0, d["z_c0"], 0, 0, d["L_cond"], d["r_ci"])
    cond, _ = occ.cut([(3, c_out)], [(3, c_in)])
    # Epoxy core: solid around the conductor (annulus from conductor outer radius)
    k_out = occ.addCylinder(0, 0, 0, 0, 0, d["L_core"], d["r_k"])
    k_in = occ.addCylinder(0, 0, 0, 0, 0, d["L_core"], d["r_co"])
    core, _ = occ.cut([(3, k_out)], [(3, k_in)])
    # Flange ring
    f_out = occ.addCylinder(0, 0, d["z_f0"], 0, 0, d["t_f"], d["r_f"])
    f_in = occ.addCylinder(0, 0, d["z_f0"], 0, 0, d["t_f"], d["r_k"])
    flange, _ = occ.cut([(3, f_out)], [(3, f_in)])
    occ.fragment(cond + core + flange, [])
    occ.synchronize()

    # Identify volumes by their centre of mass radius
    vols = {}
    for _, v in gmsh.model.getEntities(3):
        x, y, z = gmsh.model.occ.getCenterOfMass(3, v)
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(3, v)
        rmax = max(abs(x0), abs(x1))
        name = ("conductor" if abs(rmax - d["r_co"]) < 1e-4 else
                "core" if abs(rmax - d["r_k"]) < 1e-4 else "flange")
        vols.setdefault(name, []).append(v)
    assert sorted(vols) == ["conductor", "core", "flange"], vols
    for name, vs in vols.items():
        gmsh.model.addPhysicalGroup(3, vs, TAGS[name], name)

    surfs = {}
    for _, s in gmsh.model.getEntities(2):
        surfs.setdefault(_classify_surface(s, d), []).append(s)
    for name, ss in surfs.items():
        gmsh.model.addPhysicalGroup(2, ss, TAGS[name], name)

    # Mesh size: global h, refined at the conductor surface and in the flange zone
    cond_surfs = surfs["if_conductor_core"] + [s for s in surfs["conductor_exposed"]
                                               if abs(gmsh.model.getBoundingBox(2, s)[2] -
                                                      gmsh.model.getBoundingBox(2, s)[5]) > TOL]
    fd = gmsh.model.mesh.field
    dist = fd.add("Distance")
    fd.setNumbers(dist, "SurfacesList", cond_surfs)
    fd.setNumber(dist, "Sampling", 200)
    thr = fd.add("Threshold")
    fd.setNumber(thr, "InField", dist)
    fd.setNumber(thr, "SizeMin", mc["conductor_size_ratio"] * h)
    fd.setNumber(thr, "SizeMax", h)
    fd.setNumber(thr, "DistMin", 0.0)
    fd.setNumber(thr, "DistMax", mc["refine_distance_ratio"] * h)
    box = fd.add("Box")
    hf = mc["flange_size_ratio"] * h
    fd.setNumber(box, "VIn", hf)
    fd.setNumber(box, "VOut", h)
    for key, val in [("XMin", -d["r_f"]), ("XMax", d["r_f"]), ("YMin", -d["r_f"]), ("YMax", d["r_f"]),
                     ("ZMin", d["z_f0"] - 2 * hf), ("ZMax", d["z_f1"] + 2 * hf)]:
        fd.setNumber(box, key, val)
    fd.setNumber(box, "Thickness", 2 * h)
    fmin = fd.add("Min")
    fd.setNumbers(fmin, "FieldsList", [thr, box])
    fd.setAsBackgroundMesh(fmin)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)  # HXT
    gmsh.option.setNumber("Mesh.Optimize", 1)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(mc["element_order"])
    gmsh.model.mesh.optimize("HighOrder")  # untangle curved elements if needed

    # Statistics
    node_tags, _, _ = gmsh.model.mesh.getNodes()
    etypes, etags, _ = gmsh.model.mesh.getElements(3)
    tets = np.concatenate(etags)
    q = gmsh.model.mesh.getElementQualities(tets, "minSICN")
    gamma = gmsh.model.mesh.getElementQualities(tets, "gamma")
    stats = {
        "h_global_m": h,
        "h_conductor_m": mc["conductor_size_ratio"] * h,
        "h_flange_m": hf,
        "element_type": gmsh.model.mesh.getElementProperties(int(etypes[0]))[0],
        "n_nodes": int(len(node_tags)),
        "n_elements": int(len(tets)),
        "quality_minSICN_min": float(q.min()),
        "quality_minSICN_mean": float(q.mean()),
        "quality_gamma_min": float(gamma.min()),
        "quality_gamma_mean": float(gamma.mean()),
        "share_minSICN_below_0p2": float(np.mean(q < 0.2)),
    }
    return stats


def create_mesh(cfg, h=None, comm=MPI.COMM_WORLD):
    """Build the model and return (dolfinx MeshData, stats)."""
    from dolfinx.io import gmsh as dgmsh

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("bushing")
    stats = build_model(cfg, h)
    data = dgmsh.model_to_mesh(gmsh.model, comm, 0, gdim=3)
    gmsh.finalize()
    return data, stats


if __name__ == "__main__":
    cfg = load_config()
    gmsh.initialize(sys.argv)
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("bushing")
    stats = build_model(cfg)
    RESULTS_DIR.mkdir(exist_ok=True)
    gmsh.write(str(RESULTS_DIR / "bushing.msh"))
    gmsh.finalize()
    d = derived_geometry(cfg)
    update_results("mesh", stats)
    update_results("geometry", {k: float(v) for k, v in d.items()})
    for k, v in stats.items():
        print(f"{k:28s} {v}")
