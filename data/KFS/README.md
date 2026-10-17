# Karlsruhe fine sand

Test data from Part I of the database of Wichtmann & Triantafyllidis (2016), *Acta Geotechnica* 11(4),
739–761 (see `NOTICE.md`). Only the tests used in the thesis are included, one
folder per test type, each with an index file written by the author. The TCU index lists the 22
cyclic tests of the database; only the six used are included. The tests with strain cycles used in
the thesis are not included.

| folder | index (one row per test) | files | columns used |
|---|---|---|---|
| `TCU/` | `TCU_index.txt`: e0; p0; q_amp; CSR (not read: the tests are tabulated in `params/params_KFS.m`) | `TCUI8, 9, 10, 18, 19, 20.dat` | 7 p' [kPa], 8 q [kPa], 9 ε₁ [%] |
| `TMD/` | `TMD_index.txt`: ε₁ at the end; p' at the first reading (p0); e0 | `TMD1.dat` … `TMD25.dat` | 1 ε₁ [%], 2 ε_v [%], 6 q, 7 p' |
| `TMU/` | `TMU_index.txt`: – ; p0; e0 | `TMU1.dat` … `TMU6.dat` | 1 ε₁ [%], 7 p', 8 q |
| `TMU2/` | `TMU2_index.txt`: – ; p0; e0 | `TMU1.dat` … `TMU8.dat` | 1 ε₁ [%], 7 p', 8 q |
| `OED/` | `OED_index.txt`: test; I_D0; e at the first reading (not read) | `OE1.dat` … `OE12.dat` | 1 σ₁ [kPa], 3 e |

Strains are in % in the files; the model works with absolute strains. The first column of the TMU and
TMU2 indices has no meaning (it repeats the TMD index) and is not read. TMU2 2 is the same test as
TMU 1, and TMU2 6 the same as TMU 4. `TMD10.dat` has a header
without a line of units; `importdata` reads it like the others.
