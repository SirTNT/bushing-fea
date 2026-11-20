"""Stage 7: mesh convergence study for S1 (top deflection, maximum stress) and f1, on three meshes."""
import gc

from common import load_config, update_results
from geometry import create_mesh
from hand_calcs import sections
from common import pct_change


def run(cfg):
    from modal import solve_modes
    from static import ElasticModel

    A_end = sections(cfg)["A_c"]
    rows = []
    for h in cfg["mesh"]["convergence_h"]:
        data, stats = create_mesh(cfg, h=h)
        model = ElasticModel(data, cfg)
        u = model.solve(traction_top=(cfg["loads"]["cantilever_force"] / A_end, 0, 0))
        s1, _ = model.summarise(u)
        model.ksp.destroy()          # free the static factorisation before the modal solve (memory)
        gc.collect()
        modes = solve_modes(model, 2)
        f1 = float(modes[0][0] ** 0.5 / (2 * 3.141592653589793))
        vm = s1["von_mises"]
        row = {"h_m": h, "n_elements": stats["n_elements"], "n_nodes": stats["n_nodes"], "n_dofs": model.ndofs,
               "quality_minSICN_min": stats["quality_minSICN_min"],
               "top_deflection_m": s1["top_deflection_x_m"],
               "max_von_mises_Pa": max(v["max_Pa"] for v in vm.values()),
               "max_von_mises_body": max(vm, key=lambda b: vm[b]["max_Pa"]),
               "max_von_mises_away_from_corners_Pa": max(v["max_away_from_corners_Pa"] for v in vm.values()),
               "von_mises_per_body": {b: {"max_Pa": v["max_Pa"], "away_Pa": v["max_away_from_corners_Pa"]}
                                      for b, v in vm.items()},
               "f1_Hz": f1}
        rows.append(row)
        print(f"h = {h * 1e3:.1f} mm: {row['n_elements']} elements, {row['n_dofs']} dofs, "
              f"top deflection {row['top_deflection_m'] * 1e3:.4f} mm, max vm {row['max_von_mises_Pa'] / 1e6:.2f} MPa "
              f"(away {row['max_von_mises_away_from_corners_Pa'] / 1e6:.2f}), f1 {f1:.3f} Hz", flush=True)
        del model, modes, u, data
        gc.collect()
    m, f = rows[-2], rows[-1]
    change = {
        "top_deflection_pct": pct_change(f["top_deflection_m"], m["top_deflection_m"]),
        "max_von_mises_pct": pct_change(f["max_von_mises_Pa"], m["max_von_mises_Pa"]),
        "max_von_mises_away_from_corners_pct": pct_change(f["max_von_mises_away_from_corners_Pa"],
                                                           m["max_von_mises_away_from_corners_Pa"]),
        "f1_pct": pct_change(f["f1_Hz"], m["f1_Hz"]),
    }
    per_body = {b: {"max_pct": pct_change(f["von_mises_per_body"][b]["max_Pa"], m["von_mises_per_body"][b]["max_Pa"]),
                    "away_pct": pct_change(f["von_mises_per_body"][b]["away_Pa"], m["von_mises_per_body"][b]["away_Pa"])}
                for b in f["von_mises_per_body"]}
    out = {"rows": rows, "change_medium_to_fine": change, "change_medium_to_fine_per_body": per_body}
    update_results("convergence", out)
    print("medium to fine change:", {k: round(v, 3) for k, v in change.items()})
    return out


if __name__ == "__main__":
    run(load_config())
