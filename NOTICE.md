# Notices

The MIT License (`LICENSE`) covers the work of J. Rozas P. The files derived from the constitutive
driver of C. Tamagnini and D. Mašín contain, besides that work, original portions of the driver, which
keep the terms of its distribution. The test data in `data/KFS` remain the work of their authors.

## Constitutive driver of C. Tamagnini and D. Mašín

The workflow (`driver_HYP_istr.m`, `input_data.m`, `init_state.m`, `xyplot_vmh.m`) and the following
files of `model/` derive from the constitutive driver for hypoplasticity with intergranular strain
written by C. Tamagnini (2012), with modifications for hypoplasticity by D. Mašín, distributed as
course material of the ALERT Olek Zienkiewicz School 2012 "Constitutive Modeling of Soils" (Dresden).
The driver integrates the equations with an explicit adaptive Runge–Kutta scheme and the
mixed-control algorithm of Bardet & Choucair (1991).

| file | relation to the original driver |
|---|---|
| `constraints.m`, `stress_mult.m`, `stress_trace.m` | unchanged (line endings aside) |
| `updateModel.m`, `hypo_update.m` | adapted: sign mapping, failure code, error norm |
| `updateModel_stress.m` | adapted from `updateModel.m` for branches that end at a stress amplitude |
| `f_hyp.m`, `hyp_stif.m`, `M_istr.m`, `inv_s.m`, `inv_e.m` | started from the original files; the constitutive equations of the hybrid model are by J. Rozas P. |

The original driver is distributed freely by its authors at https://soilmodels.com, without an
explicit license. In the files above, the original portions (C. Tamagnini and D. Mašín, 2012) are not
covered by the MIT License and remain the work of their authors; their original header is kept in
`model/updateModel.m` and `model/updateModel_stress.m` ("written by C. Tamagnini - jan. 2012") and
quoted in `driver_HYP_istr.m`. The modifications and the constitutive equations of the hybrid model
(copyright 2026 J. Rozas P.) are under the MIT License. `espejo_estado.m`, `params/`, `tools/`,
`reproduce_thesis.m`, `reference/`, the documentation (`*.md`, `docs/`) and the metadata
(`CITATION.cff`, `.gitignore`, `.gitattributes`) are entirely by J. Rozas P. and are under the MIT
License.

The form of the detector y_δ of the cyclic memory (`model/M_istr.m` L86–92) was searched for and
tested with the assistance of Claude (Anthropic); the author reviewed it and takes responsibility
for it.

Bardet, J.P., Choucair, W. (1991). A linearized integration technique for incremental constitutive
equations. *International Journal for Numerical and Analytical Methods in Geomechanics* 15(1), 1–19.

## Experimental data

`data/KFS/` contains test data on Karlsruhe fine sand from the database published by
T. Wichtmann and T. Triantafyllidis. All of them (oedometric, monotonic drained and undrained
triaxial, and undrained triaxial tests with stress cycles) belong to Part I:

- Wichtmann, T., Triantafyllidis, T. (2016). An experimental database for the development,
  calibration and verification of constitutive models for sand with focus to cyclic loading.
  Part I: tests with monotonic loading and stress cycles. *Acta Geotechnica* 11(4), 739–761.

The tests with strain cycles of Part II (*Acta Geotechnica* 11(4), 763–774), also used in the thesis,
are not included. The data remain the work of their authors; cite the paper above when you use them. The index files
`data/KFS/*/*_index.txt` were compiled by J. Rozas P. from those data.

The Bío-Bío sand data are not included: Saldaña Sotelo, H.A. (2023), MSc thesis, Universidad de
Concepción; Saldaña et al. (2024), *Soil Dynamics and Earthquake Engineering* 179, 108511;
Cuitiño et al. (2015), in *From Fundamentals to Applications in Geotechnics*, 1285–1292;
Brown et al. (2024), *Soil Dynamics and Earthquake Engineering* 185, 108848.
