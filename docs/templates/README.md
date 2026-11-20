# Transformer Bushing FEA: Static, Modal and Thermal-Structural

A self-directed, simplified finite element study of a dry-type (resin-impregnated) transformer bushing, built only with open-source Python tools (Gmsh, FEniCSx/dolfinx, SLEPc, PyVista, Matplotlib). The model has three bonded bodies (copper conductor tube, epoxy core, aluminium flange) and covers a static terminal (cantilever) load, self-weight, a static-equivalent horizontal acceleration, modal analysis, steady-state Joule heating and the resulting thermal stress. Every FEA result is checked against a hand calculation. The load cases are *based on* the type of checks in IEC 60137 and IEEE 693, but the study does **not** claim compliance with these standards, qualification or a rated design. All material values are typical literature values, assumed.

![Geometry and mesh](figures/fig1_geometry_mesh.png)

## Key results (mesh h = {{mesh-h-mm}} mm, {{mesh-nel}} quadratic tetrahedra, {{ndofs}} displacement DOFs)

| Quantity | Value |
|---|---:|
| S1 top deflection under F = {{F-N}} N | {{S1-top-mm}} mm |
| First natural frequency (bending, two near-equal directions) | {{f1}} Hz and {{f2}} Hz |
| Hot spot at I = {{I-A}} A, ambient {{T-amb}} C | {{Tmax}} C |
| Thermal von Mises stress, max in conductor / core / flange (fixed flange face) | {{TH-conductor-max}} / {{TH-core-max}} / {{TH-flange-max}} MPa |

## Verification against hand calculations

{{verification-table-md}}

The two-segment beam and 1D conduction checks are the meaningful ones; the uniform-beam and free-copper rows are kept to show why the simpler formulas do not fit this geometry. See `report/report.pdf` for the explanations.

## How to run

```bash
conda env create -f environment.yml      # or: micromamba create -f environment.yml
conda activate bushing-fea
python src/run_all.py                    # all stages ({{run-time-min}} minutes on the test machine)
python src/run_all.py --skip-convergence # keeps the stored convergence results, much faster
cd report && pdflatex report.tex && pdflatex report.tex
```

`run_all.py` regenerates `results/results.json`, every figure in `figures/`, the tables in `report/generated/` and the numbers in this README, `docs/EXPLAINER.md` and `docs/CV_ENTRY.md` (from `docs/templates/`). All parameters are in `config.yaml`. Tested with dolfinx {{ver-dolfinx}}, Gmsh {{ver-gmsh}}, PETSc {{ver-petsc}}. The fine mesh of the convergence study needs a lot of memory; reduce `convergence_h` in `config.yaml` on a smaller machine.

## Repository layout

```
config.yaml      all geometry, material, load and mesh parameters (SI units)
src/             geometry.py, static.py, modal.py, thermal.py, thermo_mech.py, hand_calcs.py,
                 convergence.py, make_figures.py, make_tables.py, run_all.py
results/         results.json and field files (.vtu, open in ParaView)
figures/         PNG figures, 300 dpi
report/          report.tex (simulation report) and generated tables
docs/            EXPLAINER.md (plain-language walkthrough), CV_ENTRY.md
```

## Limitations

- Simplified geometry: no silicone sheds or housing, no condenser foils, no oil-side details, no draw lead in the bore.
- Perfect bonding between bodies, linear elastic isotropic materials, properties at a single temperature.
- The flange bottom face is fully fixed. This idealised mount over-predicts thermal stress near the flange; a sliding-support sensitivity case is included to show how much.
- Thermal stress assumes a stress-free state at {{T-ref}} C. A real resin-impregnated core is cured at elevated temperature and has residual stress, so the thermal stress values are indicative only.
- Static-equivalent acceleration only; no response spectrum or time-history seismic analysis to IEEE 693.
- Resistivity at 20 C (no temperature dependence), no short-circuit electromagnetic forces, no dielectric losses.
