# Bío-Bío sand: data not included

The measured curves of Bío-Bío sand are not redistributed. The model does not need them: the published
tests run from the initial states tabulated in `params/params_BBS.m`.

| test | source | initial state in `params_BBS.m` |
|---|---|---|
| undrained cyclic Tx06–Tx09 (ID 6–9) | Saldaña Sotelo (2023), MSc thesis, Universidad de Concepción, Annex 4.1, Table A.4.1; Saldaña et al. (2024), *Soil Dynamics and Earthquake Engineering* 179:108511 | `TESTS`: e0 and D_r of the laboratory index, CSR = q_amp/p0 (twice the printed ratio), p0 (first measured p'), N from Table A.4.1 |
| drained monotonic DR1–DR4 (ID 1–4, thesis Fig. 11) | Cuitiño et al. (2015), *From Fundamentals to Applications in Geotechnics*, 1285–1292; Brown et al. (2024), *Soil Dynamics and Earthquake Engineering* 185:108848 | `TESTS.TMD`: p0 = 100 kPa, e0; set `'monotonic'` |

To plot a measured curve next to the model, place your copy of the files here:

| folder | file | columns (no header) |
|---|---|---|
| `TCU - Hector Saldana/` | `CSR1.dat` … `CSR4.dat` (= Tx06 … Tx09) | 1 ε₁ [%], 2 p' [kPa], 3 q [kPa] |
| `TMD_biobio_cuitino/` | `DR1.dat` … `DR4.dat` | 1 ε₁ [%], 2 q [kPa], 4 ε_v [%] |

`input_data.m` reads them only when they exist.
