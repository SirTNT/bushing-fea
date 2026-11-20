# Explainer: the bushing FEA project in plain language

This note is for interview preparation. It explains what each part of the study does, why a bushing designer would care, how the inputs were chosen, what the checks show and where the model is weak. All numbers are filled in from `results/results.json`.

## 1. What is modelled

A transformer bushing carries the high-voltage conductor through the grounded tank cover. In a resin-impregnated (dry) bushing, the conductor is surrounded by an epoxy core, and a metal flange bolts the bushing to the tank. The model keeps only these three parts:

- copper conductor tube, {{Dco-mm}} mm outer and {{Dci-mm}} mm inner diameter, {{Lcond-m}} m long, with an empty bore,
- epoxy core, {{Dk-mm}} mm diameter, {{Lcore-m}} m long, bonded to the conductor,
- aluminium flange ring, {{Df-mm}} mm outer diameter, {{tf-mm}} mm thick, bonded to the core with its bottom face {{zf0-m}} m above the core bottom.

Above the flange the bushing is a cantilever {{L-beam-m}} m long: {{a-beam-m}} m of core with conductor, then {{b-beam-m}} m of bare conductor (the terminal stud). Below the flange hangs a shorter part that would sit in the transformer.

Left out on purpose: sheds and housing, condenser foils, oil-side details, the draw lead in the bore, bolts and the tank. All contacts are perfectly bonded, and all materials are linear elastic and isotropic with typical literature values (assumed, not from a datasheet).

The mesh has {{mesh-nel}} quadratic (10-node) tetrahedra with a global size of {{mesh-h-mm}} mm, refined to {{mesh-hc-mm}} mm at the conductor and {{mesh-hf-mm}} mm at the flange. The worst element has a scaled Jacobian of {{mesh-qmin}} (mean {{mesh-qmean}}), which is acceptable.

## 2. The analyses and why a designer cares

**Static, S1 cantilever.** A line connected to the terminal pulls sideways on the bushing (conductor tension, wind, ice). IEC 60137 has a cantilever withstand test for this. I applied an illustrative {{F-N}} N horizontally on the top end face; the real test load depends on the bushing rating and must be taken from the standard. Result: top deflection {{S1-top-mm}} mm. The highest stress is in the copper where it leaves the core top ({{S1-conductor-max}} MPa at the corner, {{S1-conductor-away}} MPa at {{excl-mm}} mm from it), because the bare stud carries the full bending moment on a small section. The flange sees only {{S1-flange-away}} MPa away from corners.

**Static, S2 self-weight.** A sanity case. The bushing model weighs {{mass-kg}} kg, and the largest displacement is only {{S2-umax-um}} micrometres, so self-weight matters little for this short, stiff model.

**Static, S3 static-equivalent seismic.** A horizontal acceleration of {{seis-g}} g on all mass, plus gravity. The top moves {{S3-top-um}} micrometres. This is a quick screening check, not a seismic qualification. IEEE 693 uses response spectra or shake-table tests, and the dynamic amplification depends on the natural frequency (next point).

**Modal.** Earthquake ground motion has most of its energy at low frequencies, so a bushing whose first natural frequency is well above that responds almost as a rigid body, while a low frequency means the ground motion is amplified. The first two modes are bending in two perpendicular directions at {{f1}} Hz and {{f2}} Hz. For an axisymmetric body they should be equal; the {{f12-diff-pct}} % difference comes only from the unstructured mesh, so this is a useful symmetry check. Modes 3 and 4 ({{f3}} Hz) are bending of the part below the flange, and modes 5 and 6 ({{f5}} Hz) are second bending of the upper part. In a real bushing the long porcelain or composite housing, the oil and a flexible tank cover lower these frequencies a lot.

**Thermal.** Load current heats the conductor (I squared R). The heat must cross the epoxy, which conducts poorly, and leave by natural convection. The insulation temperature limits the bushing rating (IEC 60137 sets temperature-rise limits for the hot spot). At {{I-A}} A the conductor generates {{Qp-Wm}} W per metre, {{Qtot-W}} W in total. With convection h = {{h-conv}} W/(m2 K) to {{T-amb}} C air, the hot spot is {{Tmax}} C, at the top end of the conductor. This is because I treated the exposed conductor ends as insulated, so their heat has to flow along the copper into the core. In a real bushing the terminal connection would carry heat away (or bring heat in from a hot line), so the hot spot location depends strongly on that boundary condition.

**Thermal-structural.** Copper expands about half as much as epoxy per kelvin ({{alpha-cu}} versus {{alpha-ep}} micrometres per metre per kelvin). When the bonded bushing warms up, the two materials fight each other and stress builds up at the interfaces. With a stress-free state at {{T-ref}} C and the flange face fully fixed, the maximum von Mises stress is {{TH-conductor-max}} MPa in the conductor, {{TH-core-max}} MPa in the core and {{TH-flange-max}} MPa in the flange, all at the flange level. These numbers are much larger than the {{F-N}} N cantilever stresses, and adding S1 changes them very little.

I checked how much of this comes from the idealised fixed support by repeating the case with a sliding support (axial restraint only, free radial growth). The maxima drop to {{THSL-conductor-max}}, {{THSL-core-max}} and {{THSL-flange-max}} MPa. So a large part of the fixed-support result is caused by the rigid mount, and a real assessment needs the tank cover and bolt stiffness. A second strong effect is the stress-free temperature: a real resin core is cured hot on the conductor, so its residual stress state is not zero at {{T-ref}} C.

## 3. How the boundary conditions and loads were chosen

| Item | Choice in the model | What a real design would do |
|---|---|---|
| Support | Flange bottom face fully fixed | Bolted joint on a flexible tank cover, possibly with contact and bolt preload |
| Terminal load | {{F-N}} N uniform traction on the top end face, illustrative | Rated cantilever test load from IEC 60137, applied at the terminal point |
| Seismic | {{seis-g}} g static horizontal plus gravity | Response spectrum analysis to IEEE 693, with the tank flexibility (amplification factors) |
| Heat source | I squared R with resistivity at 20 C | Temperature-dependent resistivity, skin effect, dielectric losses, oil and tank temperatures |
| Cooling | h = {{h-conv}} W/(m2 K) to {{T-amb}} C on core and flange | Correlations for natural convection and radiation on the real housing, oil-side heat transfer |
| Thermal stress | Stress-free at {{T-ref}} C, bonded | Cure temperature and residual stress of the resin, temperature-dependent epoxy modulus |

## 4. What the verification shows

{{verification-table-md}}

**Static deflection.** The uniform composite beam ({{ver1-an}} mm) is not a fair check, because it gives the bare {{b-beam-m}} m stud the bending stiffness of the whole core. The two-segment beam ({{ver2-an}} mm) is the right check, and the FEA is {{ver2-diff}} % more flexible. About {{share-overhang-pct}} % of the beam-theory deflection comes from the stud (its own bending plus the rotation of the core top). The beam model clamps the core rigidly at the flange and assumes that the conductor is rigidly built in at the core top; in 3D the epoxy around the conductor and at the flange deforms locally, which adds flexibility. Shear deformation of the core, which beam theory ignores, adds a little more. Figure 8 compares the deflection lines.

**First natural frequency.** The uniform formula gives {{f1-uniform}} Hz against {{f1}} Hz from FEA. The error is large because the formula assumes the composite mass ({{meff}} kg/m) all the way to the tip, while the stud only has {{mc}} kg/m, and tip mass controls the first bending frequency. A 1D beam model with the two real segments gives {{f1-twoseg}} Hz, {{ver4-diff}} % from FEA. The FEA is slightly higher, probably because the flange stiffens the core over its {{tf-mm}} mm height, so the effective clamp is a little higher than the flange bottom face.

**Radial temperature.** At core mid-height the 1D long-cylinder solution with the local generation ({{Qp-Wm}} W/m) predicts {{mid-Tco-an}} C at the conductor surface; the FEA gives {{mid-Tco-fea}} C. The FEA is higher because the heat generated in the insulated conductor ends travels along the copper into the core section, so the radial heat flow at mid-height is {{mid-Qlocal}} W/m ({{mid-Qlocal-excess-pct}} % more than generated locally). If the 1D formula uses that local heat flow, the temperature drop across the epoxy matches the FEA within {{ver8-absdiff}} %. So the conduction physics and the mesh are right; the difference is a real 3D effect of a short bushing with insulated ends. This is the same kind of check as the analytical comparison in my thesis on multilayer cable insulation. The heat balance also closes: {{Q-gen}} W generated and {{Q-out}} W convected.

**Thermal elongation.** Free copper expansion over the {{L-beam-m}} m above the flange gives {{el-free-mm}} mm, but the FEA top displacement is {{el-fea-mm}} mm. The epoxy is {{EA-ratio}} times stiffer than the copper in the axial direction (E times A) and expands twice as much, so the bonded core pulls the conductor with it. An EA-weighted composite estimate gives {{el-comp-mm}} mm, within {{ver10-absdiff}} % of the FEA.

## 5. Mesh convergence

| Mesh | h (mm) | Elements | Top deflection (mm) | Max stress (MPa) | Max stress away from corners (MPa) | f1 (Hz) |
|---|---:|---:|---:|---:|---:|---:|
| Coarse | {{cv-coarse-h}} | {{cv-coarse-nel}} | {{cv-coarse-top}} | {{cv-coarse-vm}} | {{cv-coarse-vmaway}} | {{cv-coarse-f1}} |
| Medium | {{cv-medium-h}} | {{cv-medium-nel}} | {{cv-medium-top}} | {{cv-medium-vm}} | {{cv-medium-vmaway}} | {{cv-medium-f1}} |
| Fine | {{cv-fine-h}} | {{cv-fine-nel}} | {{cv-fine-top}} | {{cv-fine-vm}} | {{cv-fine-vmaway}} | {{cv-fine-f1}} |

From medium to fine, the top deflection changes by {{cv-chg-topdeflectionpct}} % and f1 by {{cv-chg-f1pct}} %, so these are converged.

The stresses are a different story, and it is worth explaining honestly:

- The **peak** stresses rise with every refinement (epoxy {{cv-coarse-core-max}}, {{cv-medium-core-max}}, {{cv-fine-core-max}} MPa; flange {{cv-coarse-flange-max}}, {{cv-medium-flange-max}}, {{cv-fine-flange-max}} MPa). They sit at sharp re-entrant or bi-material corners, where linear elasticity gives an infinite stress, so a finer mesh only shows a higher number.
- The maxima at **{{excl-mm}} mm from those corners** do not keep rising. In the copper they are converged ({{cv-chg-conductor-away}} %). In the epoxy and flange they still move by {{cv-chg-core-away}} % and {{cv-chg-flange-away}} % (epoxy {{cv-coarse-core-away}}, {{cv-medium-core-away}}, {{cv-fine-core-away}} MPa), because {{excl-mm}} mm is only one or two elements near the flange and the point sits in a steep gradient.

So the global quantities are converged, the copper stress is converged, and the local stresses at the flange are not. To get a usable local stress I would model the real fillet radius with a local mesh refinement, or use a linearised or hot-spot stress method.

## 6. Ten likely interview questions

**1. Why FEniCSx and not Ansys?**
I work on a MacBook, where Ansys does not run, and I wanted a project anyone can rerun from GitHub. FEniCSx is a research-grade open-source FE library; writing the weak forms myself made me check every boundary condition and material assignment. The physics and verification steps are the same as I used with Ansys in my thesis, and I would expect to switch back to Ansys Mechanical at work.

**2. Why does the stress not converge at the corner?**
At a sharp re-entrant corner or where two materials meet at an edge, the elastic solution has a singularity: stress grows without bound as you approach the point. Each mesh refinement samples closer to it, so the peak keeps rising. I report the peak and the maximum at {{excl-mm}} mm from the singular edges. The second value does not keep growing, but in the epoxy it still moves by {{cv-chg-core-away}} % between my medium and fine meshes, so I would not call it converged either. In a real design I would model the actual fillet radius or use a hot-spot or linearised stress method.

**3. Your deflection is {{ver2-absdiff}} % off the beam formula. Is the model wrong?**
No. The beam model assumes a perfect clamp at the flange and a rigid joint where the stud leaves the core. Figure 8 shows the deflection line, so I can see where the extra flexibility enters. A few percent extra flexibility from local deformation is expected for a composite with a soft core. The mesh study shows the FEA deflection is converged to {{cv-chg-topdeflectionpct}} %.

**4. How would you do a proper seismic analysis?**
Follow IEEE 693: get the bushing natural frequencies with the real housing and the flexible tank cover (often with an amplification factor for the tank), then run a response spectrum analysis with the required spectrum and damping, combine modes (SRSS or CQC) and directions, and compare the stresses at the flange, the porcelain or composite housing and the gaskets with allowables. For high qualification levels it ends in a shake-table test. The effective modal masses I computed ({{f1-effx}} % in mode 1) are the input to such a modal combination.

**5. What would short-circuit forces add?**
During a short circuit, current many times the rated value flows for a short time. The magnetic field from neighbouring phases gives dynamic electromagnetic forces on the conductor and terminal, and the extra I squared R heats the conductor quickly (adiabatically). Both would be new load cases: a transient structural case with the force time history (or an equivalent static load from IEC 60865-type calculations) and a short transient thermal case.

**6. Why is the hot spot at the top end and not in the middle?**
Because I made the exposed conductor ends insulated, so heat generated there can only escape by flowing along the copper into the core. In a real installation the terminal connection is a heat path, either out or in, so the hot spot location is a modelling choice I would discuss with the thermal team.

**7. The thermal stresses are high. Is the bushing failing?**
I would not conclude that. The values depend strongly on two assumptions: the fully fixed flange face (the sliding support roughly halves the peaks) and a stress-free state at {{T-ref}} C. The real cure temperature, residual stress, the temperature dependence of epoxy properties and the real mounting stiffness all need to be included before comparing with allowables. What the study does show is where the stress concentrates (flange level and conductor exit) and that thermal mismatch dominates over the cantilever load in this model.

**8. Why can you add the thermal and S1 stresses?**
The model is linear (small strains, linear elastic, bonded contacts, no friction). For a linear model the displacement fields superpose. I added the displacement fields and recalculated the von Mises stress, because von Mises itself is not additive.

**9. How did you check the modal results?**
Three ways: the two bending modes have almost equal frequencies ({{f12-diff-pct}} % apart) as symmetry requires; the effective masses show each mode's direction; and a 1D two-segment beam model gives {{f1-twoseg}} Hz against {{f1}} Hz. f1 also changes by only {{cv-chg-f1pct}} % between the medium and fine meshes.

**10. What would you add next?**
The real shed and housing geometry, the oil side and tank cover stiffness, contact instead of bonding at the flange, temperature-dependent resistivity and material properties, a response spectrum seismic case and a short-circuit force case. I would also compare against test data if available, because verification against hand calculations only shows that the model solves the equations correctly, not that the equations represent the real bushing.
