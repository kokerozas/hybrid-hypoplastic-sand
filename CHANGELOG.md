# Changelog

## v2026.10

First public release. Supporting code of the MSc thesis of J. Rozas P. (Universidad Católica de la
Santísima Concepción) and of a manuscript in preparation for *Computers and Geotechnics*.

- Constitutive core (`model/`), byte-identical to the calibrated core of the thesis (`docs/MD5.csv`):
  equations in the compression-negative convention, laboratory boundary in `model/espejo_estado.m`.
- Element-test workflow of the original constitutive driver of C. Tamagnini and D. Mašín:
  `driver_HYP_istr.m`, `input_data.m`, `init_state.m`, `xyplot_vmh.m`.
- Parameter sets of Karlsruhe fine sand and Bío-Bío sand, cyclic and monotonic (`params/`).
- Undrained cyclic tests cycle by cycle with the two cycle counts of the thesis (`tools/cyclic_test.m`).
- `reproduce_thesis.m`: 49 element tests and 113 reference values of the thesis run of 2026-09-23;
  verified with MATLAB R2020a on 2026-10-16 (113 EXACT, largest difference 0.0e+00, 10.1 min).
- Test data of Karlsruhe fine sand (Part I of Wichtmann & Triantafyllidis 2016); the Bío-Bío data are
  not redistributed.
