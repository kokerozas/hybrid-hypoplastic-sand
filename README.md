# hybrid-hypoplastic-sand

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23439637.svg)](https://doi.org/10.5281/zenodo.23439637)

MATLAB implementation of a hybrid hypoplastic constitutive model for sands under monotonic and cyclic
loading, up to and beyond initial liquefaction. It combines the hypoplasticity of von Wolffersdorff
(1996), the intergranular strain of Niemunis & Herle (1997), the fabric tensor of Liao et al. (2022),
the semifluidized state of Barrero et al. (2020) in the form of Liao et al. (2022), and two scalar
strain memories (cyclic and monotonic) with a ramp of the mobilized friction angle. It is calibrated
on Karlsruhe fine sand (KFS) and Bío-Bío sand (BBS).

Supporting code of the MSc thesis of Jorge Rozas P. (Universidad Católica de la Santísima Concepción)
and of a manuscript in preparation for *Computers and Geotechnics*.

Requires MATLAB R2020a or later, no toolboxes.

## Run an element test

The workflow is the one of the original constitutive driver of C. Tamagnini (2012), with the
modifications for hypoplasticity of D. Mašín (ALERT Olek Zienkiewicz School 2012, Dresden):

```
driver_HYP_istr.m      the script you run
  input_data.m         what to run, material constants (parms) and loading path (nspb, path_info)
  init_state.m         initial state vector (y0)
  updateModel          integration over the loading path (model/)
  xyplot_vmh.m         plots: model in red, measured curves in blue
```

Open `input_data.m`, edit its first block and run `driver_HYP_istr`:

```matlab
SAND = 'KFS';   % 'KFS' | 'BBS'
TEST = 'TCU';   % 'TCU' undrained cyclic | 'TMU' / 'TMU2' undrained monotonic | 'TMD' drained monotonic | 'OED' oedometer
ID   = 10;      % a published test, or [] to use your own initial state and loading
```

With these defaults (KFS, TCU, test 10) the driver integrates the undrained cyclic test TCUI10, prints
`N_liq = 14.008 cycles (first p'/p'0 <= 0.05; counting: halfwave)   test: 15` and plots the model next
to the measured curves.

With `ID = []` the test starts from the values below in the same block: `p0` and `e0`, the deviator
amplitude `q_amp` and the number of cycles `NCYC` (TCU), the final axial strain `EPS_A` (TMU, TMD), or
the initial stress and strain increments of the oedometer. The loading path follows the convention of
the original driver: `path_info = [icode; nstep; DX]`, one column per branch, with the control codes of
`model/constraints.m` (1 undrained TX, 2 oedometer, 4 drained TX with constant radial stress, 7 undrained
simple shear, ...). In TCU every branch ends when the deviator stress reaches `+q_amp` or `-q_amp`.

Published tests:

| sand | test | ID | initial state from |
|---|---|---|---|
| KFS | TCU | 8, 9, 10, 18, 19, 20 | `params/params_KFS.m` |
| KFS | TMD | 1–25 | `data/KFS/TMD/TMD_index.txt` |
| KFS | TMU, TMU2 | 1–6, 1–8 | `data/KFS/TMU/TMU_index.txt`, `data/KFS/TMU2/TMU2_index.txt` |
| KFS | OED | 1–12 | the measured curve at the initial stress; loaded up to the maximum measured stress |
| BBS | TCU | 6, 7, 8, 9 (Tx06–Tx09) | `params/params_BBS.m` |
| BBS | TMD | 1–4 (thesis Fig. 11) | `params/params_BBS.m` |

The thesis calibrated the oedometric indices of KFS with a per-specimen loading recipe that the driver
does not repeat, so the OED curves are an illustration, not a reproduction.

For TCU the driver prints the number of cycles to initial liquefaction, the first state with
p'/p'₀ ≤ 0.05.

## Reproduce the main numbers of the thesis

```matlab
reproduce_thesis
```

reruns 49 element tests and compares 113 quantities with `reference/reference_values.csv`: cycles to
liquefaction of the six KFS and four BBS cyclic tests (thesis Table 4), peak and final deviator stress
and minimum volumetric strain of the 25 KFS drained tests, and minimum p' and final q of the 14 KFS
undrained runs (12 distinct tests: TMU2 2 repeats TMU 1 and TMU2 6 repeats TMU 4). It does not cover the
oedometer calibration of KFS, the drained tests of BBS (Fig. 11) or the tests with constant strain
amplitude. It prints its running time at the end.

The reference values come from the thesis run of 2026-09-23 (MATLAB R2020a, this same constitutive
core). Each check is reported as EXACT (relative difference ≤ 1e-12), MATCH (cycles to liquefaction
within 0.01 cycles, monotonic quantities within a relative 1e-3) or FAIL; both differences are absolute
when the reference value is 0. On the platform of the reference run all
checks are expected to be EXACT; on another release or processor a difference of one integration step
is possible and shows up as MATCH. Verified on 2026-10-16 with MATLAB R2020a: 113 checks, 113 EXACT,
largest difference 0.0e+00, in 10.1 min.

Model against test, cycles to initial liquefaction:

| test | N model | N test |
|---|---|---|
| KFS TCUI8 / 9 / 10 | 299.0 / 76.0 / 14.0 | 249 / 100 / 15 |
| KFS TCUI18 / 19 / 20 | 51.0 / 10.0 / 4.5 | 54 / 15 / 6 |
| BBS Tx06 / 07 / 08 / 09 | 9.0 / 11.0 / 2.0 / 19.0 | 11 / 9 / 2 / 20 |

KFS counts the cycles with the half-wave rule (N advances 0.5 at every reversal of q) and BBS with the
position inside the cycle, as in the thesis runs; both give the cycle of the same state
(`tools/cyclic_test.m`).

## Layout

```
driver_HYP_istr.m   input_data.m   init_state.m   xyplot_vmh.m   reproduce_thesis.m
model/        the constitutive core: 12 files, byte-identical to the calibrated core (docs/MD5.csv;
              updateModel.m is the copy with the progress messages commented out)
params/       params_KFS.m, params_BBS.m: material constants, initial internal state, published tests
tools/        cyclic_test.m (TCU cycle by cycle), cycles_halfwave.m, read_index.m
data/KFS/     test data of Wichtmann & Triantafyllidis (2016) used in the thesis (TCU, TMD, TMU, TMU2, OED)
data/BBS/     README only: the Bío-Bío data are not redistributed
reference/    reference values for reproduce_thesis.m
docs/         parameter and state slots (NOMENCLATURE.md), sign convention, MD5 of the core
CHANGELOG.md  release history
```

The two sets of each sand are in `params/`: KFS `'cyclic'` (TCU) and `'monotonic'` (TMD, TMU, OED),
which differ only in the error tolerance and in four slots of the semifluidized state; BBS `'cyclic'`
(11 m layer, Saldaña Sotelo 2023) and `'monotonic'` (the drained tests of Cuitiño et al. 2015 and
Brown et al. 2024 in Fig. 11). `parms(34)` is **η**, the rate of the cyclic memory, and `parms(35)` is **κ**, the rate of the
monotonic memory; earlier versions named them the other way round, and the state slots keep the
historical names `Strain_total` (y(28) = ε_cyc) and `Strain_total2` (y(30) = ε_mon)
(`docs/NOMENCLATURE.md`). The comments inside `model/` are in Spanish and refer to the development
record of the thesis; they do not affect the results. Some of them quote superseded values or
conclusions (for example `model/M_istr.m` L95–103 and L126–132): the values that run are those in
`params/`, and where a comment disagrees with the manuscript, the manuscript governs.

Sign convention: the constitutive equations work with compression negative; the state `y0`, the
loading path, the outputs and the data use compression positive. `model/updateModel*.m` cross the
boundary with `model/espejo_estado.m`.

## Cyclic memory and the detector y_δ

The cyclic memory ε_cyc ∈ [0, 1] (state `y(28)`) extends Poblete et al. (2016) and Fuentes et al.
(2020): it grows under repeated cycles and raises the exponent χ of the intergranular strain. It
evolves with the dimensionless detector y_δ, which is large while the strain rate keeps loading along
the intergranular strain δ and zero right after a reversal:

```
d(eps_cyc)/dt = eta * (1 - y_delta - eps_cyc) * ||d(eps)/dt||
y_delta       = rho^chi * < delta_hat : d(eps)/dt > / ||d(eps)/dt||,   0 <= y_delta <= 1
chi           = chi_0 + eps_cyc * (chi_max - chi_0)
```

with ρ = ‖δ‖/R ≤ 1, δ̂ = δ/‖δ‖ and ⟨·⟩ the Macaulay bracket. The detector is `yh` in
`model/M_istr.m` (L86–92) and is selected by `parms(36) = 1`, which both parameter files enforce; the
map from each equation to its line of code is in `docs/NOMENCLATURE.md`.

The form of the detector y_δ was searched for and tested with the assistance of Claude (Anthropic);
the author reviewed it and takes responsibility for it.

## Data

- **Karlsruhe fine sand**: the stress-controlled cyclic, drained, undrained and oedometric tests used in
  the thesis, from Part I of the database published by Wichtmann & Triantafyllidis (2016); they remain
  the work of their authors. The tests with strain cycles (Part II) are not included.
  See `data/KFS/README.md`.
- **Bío-Bío sand**: the measured curves are **not** redistributed. The published tests run from the
  initial states in `params/params_BBS.m`; to plot a measured curve, place the file as described in
  `data/BBS/README.md`.

## References

- Bardet, J.P., Choucair, W. (1991). A linearized integration technique for incremental constitutive
  equations. *Int. J. Numer. Anal. Methods Geomech.* 15(1), 1–19.
- Barrero, A.R., Taiebat, M., Dafalias, Y.F. (2020). Modeling cyclic shearing of sands in the
  semifluidized state. *Int. J. Numer. Anal. Methods Geomech.* 44(3), 371–388. doi:10.1002/nag.3007
- Brown, F., Villalobos, F.A., Prada-Sarmiento, L.F. (2024). Uplift-hypoplastic macroelement model for
  shallow foundations in sand under cyclic combined loading. *Soil Dyn. Earthq. Eng.* 185, 108848.
- Cuitiño, M., Poblete, M., Wichtmann, T., Triantafyllidis, T. (2015). Prediction of permanent
  settlements of foundations on Bío Bío sand subjected to high-cyclic loading. In *From Fundamentals
  to Applications in Geotechnics*, 1285–1292.
- Fuentes, W., Wichtmann, T., Gil, M., Lascarro, C. (2020). ISA-Hypoplasticity accounting for cyclic
  mobility effects for liquefaction analysis. *Acta Geotech.* 15, 1513–1531.
  doi:10.1007/s11440-019-00846-2
- Liao, D., Yang, Z., Wang, S., Wu, W. (2022). Hypoplastic model with fabric change effect and
  semifluidized state for post-liquefaction cyclic behavior of sand. *Int. J. Numer. Anal. Methods
  Geomech.* 46(17), 3154–3177. doi:10.1002/nag.3444
- Niemunis, A., Herle, I. (1997). Hypoplastic model for cohesionless soils with elastic strain range.
  *Mech. Cohes.-Frict. Mater.* 2(4), 279–299. doi:10.1002/(SICI)1099-1484(199710)2:4<279::AID-CFM29>3.0.CO;2-8
- Poblete, M., Fuentes, W., Triantafyllidis, T. (2016). On the simulation of multidimensional cyclic
  loading with intergranular strain. *Acta Geotech.* 11(6), 1263–1285. doi:10.1007/s11440-016-0492-2
- Saldaña Sotelo, H.A. (2023). *Modelamiento numérico de asentamientos inducidos por licuación en la
  zona subductiva chilena*. MSc thesis, Universidad de Concepción.
- Saldaña-Sotelo, H., Montalva, G., Escribano, D., Núñez-Jara, S., Tiznado, J.C. (2024).
  Two-dimensional nonlinear dynamic analysis of a liquefaction case history considering spatial
  variability and long-duration megathrust earthquakes. *Soil Dyn. Earthq. Eng.* 179, 108511.
  doi:10.1016/j.soildyn.2024.108511
- von Wolffersdorff, P.-A. (1996). A hypoplastic relation for granular materials with a predefined
  limit state surface. *Mech. Cohes.-Frict. Mater.* 1(3), 251–271.
  doi:10.1002/(SICI)1099-1484(199607)1:3<251::AID-CFM13>3.0.CO;2-3
- Wichtmann, T., Triantafyllidis, T. (2016). An experimental database for the development, calibration
  and verification of constitutive models for sand with focus to cyclic loading. Part I: tests with
  monotonic loading and stress cycles. *Acta Geotech.* 11, 739–761. doi:10.1007/s11440-015-0402-z.
  Part II: tests with strain cycles and combined loading. *Acta Geotech.* 11, 763–774.
  doi:10.1007/s11440-015-0412-x

## License and credits

The MIT License (`LICENSE`) covers the work of J. Rozas P.: the constitutive equations of the hybrid
model and the files written for this repository. The integration driver and the structure of the workflow derive from the
constitutive driver of C. Tamagnini (2012) with the modifications for hypoplasticity of D. Mašín,
freely available from https://soilmodels.com without an explicit license; their original portions are
not covered by the MIT License and remain the work of their authors, as do the KFS data. `NOTICE.md` gives the status of every file.
Please cite this software as in `CITATION.cff`; release v2026.10 is archived at Zenodo,
[doi:10.5281/zenodo.23439637](https://doi.org/10.5281/zenodo.23439637).

Questions and problems: GitHub Issues of this repository.
