"""Build a single Google Colab notebook (bushing_fea_colab.ipynb) from the repository files.

Each source file is written by a %%writefile cell, so the notebook runs exactly the code in this
repository. Run:  python tools/make_notebook.py
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
cells = []


def md(text):
    cells.append({"cell_type": "markdown", "metadata": {}, "source": text.strip("\n").splitlines(keepends=True)})


def code(text):
    cells.append({"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [],
                  "source": text.strip("\n").splitlines(keepends=True)})


def writefile(rel):
    content = (ROOT / rel).read_text()
    code(f"%%writefile {rel}\n{content}")


md("""
# Transformer Bushing FEA: Static, Modal and Thermal-Structural

Simplified dry-type bushing (copper conductor, epoxy core, aluminium flange) analysed with open-source tools:
Gmsh for geometry and mesh, FEniCSx (dolfinx 0.11) for the finite element analysis, SLEPc for the modal analysis,
PyVista and Matplotlib for figures. Every number in the generated report and documents comes from
`results/results.json`.

**How to use in Google Colab**
1. Run the first code cell. It installs a conda base (condacolab) and **restarts the runtime once**. This is expected.
2. After the restart, run all cells from the top again (*Runtime > Run all*). The install cell then installs
   dolfinx and the other packages from conda-forge (a few minutes) and is skipped on later runs.
3. The full study takes roughly 15 to 30 minutes. The fine mesh of the convergence study is chosen from the
   available RAM (see the run cell).

The load cases are based on the type of checks in IEC 60137 and IEEE 693; the study does not claim compliance,
qualification or a rated design. All material values are typical literature values, assumed.
""")

md("## 1. Install (condacolab, then conda-forge packages)")
code("""
import importlib.util, sys
if importlib.util.find_spec("dolfinx") is None:
    if "google.colab" in sys.modules:
        if importlib.util.find_spec("condacolab") is None:
            !pip install -q condacolab
        import condacolab
        if not condacolab.check():
            condacolab.install()          # restarts the runtime; run all cells again afterwards
    else:
        raise RuntimeError("dolfinx not found: create the conda environment from environment.yml first")
else:
    import dolfinx
    print("dolfinx", dolfinx.__version__, "already installed")
""")
code("""
import importlib.util
if importlib.util.find_spec("dolfinx") is None:
    !mamba install -y -q -c conda-forge "fenics-dolfinx=0.11" mpich petsc4py slepc4py gmsh python-gmsh pyvista matplotlib numpy scipy pyyaml pandas
import dolfinx, gmsh
from petsc4py import PETSc
print("dolfinx", dolfinx.__version__, "| gmsh", gmsh.__version__, "| PETSc", PETSc.Sys.getVersion())
""")

md("## 2. Project folder")
code("""
import os
os.makedirs("bushing-fea/src", exist_ok=True)
os.makedirs("bushing-fea/docs/templates", exist_ok=True)
os.makedirs("bushing-fea/report", exist_ok=True)
os.makedirs("bushing-fea/tools", exist_ok=True)
%cd bushing-fea
""")

md("## 3. Parameters (`config.yaml`) and environment file")
writefile("config.yaml")
writefile("environment.yml")

md("""
## 4. Source code
`common.py` holds shared helpers (tags, material fields, stress, point evaluation, results file).
""")
for name, title in [
    ("common.py", None),
    ("smoke_test.py", "Stage 0 smoke test: Poisson problem on the unit cube"),
    ("geometry.py", "Stage 2: geometry, physical groups and quadratic tetrahedral mesh (Gmsh)"),
    ("hand_calcs.py", "Analytical checks: beam theory, 1D radial conduction, verification table"),
    ("static.py", "Stage 3: static load cases S1, S2, S3"),
    ("modal.py", "Stage 4: modal analysis with SLEPc (shift-invert)"),
    ("thermal.py", "Stage 5: steady-state Joule heating"),
    ("thermo_mech.py", "Stage 6: thermal stress, combination with S1, sliding-support sensitivity"),
    ("convergence.py", "Stage 7: mesh convergence"),
    ("make_figures.py", "Stage 8: figures"),
    ("make_tables.py", "Stages 9 and 10: report tables and generated documents"),
    ("run_all.py", "Runs everything"),
]:
    if title:
        md(f"### `{name}`: {title}")
    writefile(f"src/{name}")

md("""
## 5. Templates for the generated documents and the LaTeX report
These contain `{{placeholders}}` (Markdown) and `\\V{key}` macros (LaTeX) that are filled from `results.json`.
""")
for rel in ["docs/templates/README.md", "docs/templates/EXPLAINER.md", "docs/templates/CV_ENTRY.md",
            "report/report.tex"]:
    writefile(rel)

md("## 6. Smoke test")
code("import sys\n!{sys.executable} src/smoke_test.py")

md("""
## 7. Run the full study
The fine convergence mesh (h = 13.5 mm, about 1 million degrees of freedom) needs about 13 GB of RAM with the
direct solver. On a machine with less memory the fine size is set to 16 mm; this change is written into
`config.yaml` and therefore recorded in `results.json`.
""")
code("""
import os, re
ram_gb = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / 1e9
print(f"RAM: {ram_gb:.1f} GB")
if ram_gb < 20:
    txt = open("config.yaml").read()
    txt = re.sub(r"convergence_h: \\[.*?\\]", "convergence_h: [0.030, 0.020, 0.016]", txt)
    open("config.yaml", "w").write(txt)
    print("fine convergence mesh set to h = 16 mm for this machine")
""")
code("""
import subprocess, sys
ret = subprocess.run([sys.executable, "src/run_all.py"]).returncode
assert ret == 0, f"run_all.py failed with exit code {ret}"
""")

md("## 8. Results")
code("""
import json
r = json.load(open("results/results.json"))
print(f"S1 top deflection: {r['static']['S1_cantilever']['top_deflection_x_m'] * 1e3:.3f} mm")
print("Natural frequencies (Hz):", [round(m['frequency_Hz'], 2) for m in r['modal']['modes']])
print(f"Hot spot: {r['thermal']['T_max_C']:.1f} C")
print()
print(f"{'Check':95s} {'Analytical':>11s} {'FEA':>10s} {'Diff %':>8s}")
for row in r["verification"]:
    print(f"{row['check'] + ' (' + row['unit'] + ')':95s} {row['analytical']:11.4f} {row['fea']:10.4f} {row['difference_pct']:+8.2f}")
""")
code("""
from IPython.display import Image, display
import glob
for f in sorted(glob.glob("figures/*.png")):
    print(f)
    display(Image(filename=f, width=750))
""")
code("""
from IPython.display import Markdown
display(Markdown(open("docs/CV_ENTRY.md").read()))
""")

md("""
## 9. Optional: compile the report to PDF
Installs a LaTeX distribution with apt (about 5 minutes on Colab).
""")
code("""
!apt-get -qq update && apt-get -qq install -y --no-install-recommends texlive-latex-recommended texlive-latex-extra texlive-science texlive-fonts-recommended > /dev/null
%cd report
!pdflatex -interaction=nonstopmode report.tex > /dev/null; pdflatex -interaction=nonstopmode report.tex | tail -2
%cd ..
""")

md("## 10. Download everything as a zip")
code("""
!zip -qr bushing-fea-results.zip results/results.json figures report docs README.md config.yaml src
try:
    from google.colab import files
    files.download("bushing-fea-results.zip")
except ImportError:
    print("written bushing-fea-results.zip")
""")

nb = {"cells": cells, "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                                   "language_info": {"name": "python"}, "colab": {"provenance": []}},
      "nbformat": 4, "nbformat_minor": 5}
for i, c in enumerate(nb["cells"]):
    c["id"] = f"cell-{i:02d}"
out = ROOT / "bushing_fea_colab.ipynb"
out.write_text(json.dumps(nb, indent=1))
print("written", out, "with", len(cells), "cells")
