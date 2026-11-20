# CV entry

Every number below is filled from `results/results.json` by `src/make_tables.py`.

```latex
\projline{Transformer Bushing FEA: Static, Modal and Thermal-Structural}{FEniCSx, Gmsh, Python, self-directed, {{date}}}
\begin{itemize}
  \item Built an open-source 3D finite element model of a simplified transformer bushing in Python: static terminal load, modal analysis ($f_1$ = {{f1}} Hz) and Joule-heating thermal-structural coupling.
  \item Verified against beam theory (deflection within {{ver2-absdiff}}\,\%, $f_1$ within {{ver4-absdiff}}\,\%) and analytical radial heat conduction (within {{ver8-absdiff}}\,\%); mesh convergence study; github.com/SirTNT/bushing-fea.
\end{itemize}
```

## Skills section, one line

```latex
Open-source FEA in Python (FEniCSx/dolfinx, Gmsh, SLEPc): static, modal and coupled thermal-structural analysis with analytical verification.
```

## Notes for using these lines

- "Analytical radial heat conduction (within {{ver8-absdiff}} %)" refers to the temperature drop across the epoxy at core mid-height, with the 1D formula fed the local radial heat flow taken from the FEA. With the locally generated heat instead, the difference is {{ver7-absdiff}} %, because heat from the insulated conductor ends flows along the copper into the core region. Be ready to explain both numbers.
- "Deflection within {{ver2-absdiff}} %" is against the two-segment beam (core plus conductor, then bare conductor). The uniform composite beam differs by {{ver1-absdiff}} % and is not a fair comparison for this geometry.
- "$f_1$ within {{ver4-absdiff}} %" is against a 1D two-segment beam model. The uniform cantilever formula with composite properties differs by {{ver3-absdiff}} %, because it places the full composite mass on the bare conductor at the tip.
