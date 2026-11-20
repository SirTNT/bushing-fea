"""Run the whole study: mesh, static, modal, thermal, thermal-structural, hand calculations,
verification, mesh convergence, figures, report tables and documents.

Usage:  python src/run_all.py                 (everything)
        python src/run_all.py --skip-convergence   (keeps the convergence results already stored)
"""
import argparse
import gc
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import dolfinx  # noqa: E402

from common import RESULTS_DIR, RESULTS_FILE, load_config, update_results  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-convergence", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    cfg = load_config()
    RESULTS_DIR.mkdir(exist_ok=True)
    old = json.loads(RESULTS_FILE.read_text()) if RESULTS_FILE.exists() else {}
    RESULTS_FILE.write_text("{}")
    if args.skip_convergence and "convergence" in old:
        update_results("convergence", old["convergence"])

    import gmsh
    import numpy as np
    from petsc4py import PETSc
    update_results("software", {"dolfinx": dolfinx.__version__, "gmsh": gmsh.__version__,
                                "petsc": ".".join(map(str, PETSc.Sys.getVersion())), "numpy": np.__version__})
    update_results("inputs", cfg)

    import geometry, hand_calcs, modal, static, thermal, thermo_mech
    from common import derived_geometry

    print("== Stage 2: geometry and mesh", flush=True)
    data, stats = geometry.create_mesh(cfg)
    update_results("mesh", stats)
    update_results("geometry", {k: float(v) for k, v in derived_geometry(cfg).items()})
    print(stats, flush=True)

    print("== Hand calculations", flush=True)
    hand_calcs.run(cfg)

    print("== Stage 3: static", flush=True)
    model, us, _ = static.run(data, cfg)

    print("== Stage 4: modal", flush=True)
    modal.run(data, cfg, model=model)

    print("== Stage 5: thermal", flush=True)
    T, _ = thermal.run(data, cfg)

    print("== Stage 6: thermal-structural", flush=True)
    thermo_mech.run(data, cfg, T, model=model, u_s1=us["S1_cantilever"])

    print("== Verification", flush=True)
    for row in hand_calcs.verify():
        print(f"  {row['check']:85s} analytical {row['analytical']:9.4f}  FEA {row['fea']:9.4f} {row['unit']:3s} "
              f"diff {row['difference_pct']:+7.2f} %")

    model.ksp.destroy()
    del model, us, T, data
    gc.collect()

    if not args.skip_convergence:
        print("== Stage 7: mesh convergence", flush=True)
        import convergence
        convergence.run(cfg)

    print("== Stage 8: figures", flush=True)
    import make_figures
    make_figures.main()

    print("== Stage 9/10: report tables and documents", flush=True)
    update_results("run_time_s", time.time() - t0)
    import make_tables
    make_tables.main()
    print(f"done in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
