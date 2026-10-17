# Parameter vector `parms` and state vector `y`

Read this before changing the code: two slots had inverted names in earlier versions, and two state
variables keep historical names.

## Names to watch

| slot | name in model/M_istr.m | earlier name | **what it is** | KFS | BBS |
|---|---|---|---|---|---|
| `parms(34)` | `eta_rate_cyc` | `kappa_rate` | **η**, rate of the cyclic memory `y(28)` | 20 | 0.35 |
| `parms(35)` | `kappa_rate_mon` | `eta_rate` | **κ**, rate of the monotonic memory `y(30)` | 9.9575 | 15 |
| `y(28)` | `Strain_total` | | **ε_cyc**, cyclic strain memory | | |
| `y(30)` | `Strain_total2` | | **ε_mon**, monotonic strain memory | | |

## `parms`, slot by slot

Values of the sets in `params/params_KFS.m` (`'cyclic'`) and `params/params_BBS.m` (`'cyclic'`).
"Read by" names the file of `model/` that uses the slot.

| # | name | meaning | read by | KFS | BBS |
|---|---|---|---|---|---|
| 1 | `max_ksub` | maximum substeps per step | hypo_update | 1000 | 1000 |
| 2 | `err_tol` | error tolerance of the adaptive embedded Runge–Kutta 3(2) | hypo_update | 1e-3 (cyclic) / 1e-4 (monotonic) | 1e-4 |
| 3 | `tol_NR` | tolerance of the Newton iteration (mixed or stress control) | f_hyp | 1e-10 | 1e-10 |
| 4 | `phi_c` [°] | critical friction angle | hyp_stif | 33.11 | 36.0 |
| 5 | `h_s` [kPa] | granular hardness (Bauer) | hyp_stif | 4e6 | 6.08e6 |
| 6 | `n` | exponent of Bauer's compression law | hyp_stif | 0.27 | 0.3020 |
| 7 | `e_d0` | minimum void ratio at p = 0 | hyp_stif | 0.677 | 0.596 |
| 8 | `e_c0` | critical void ratio at p = 0 | hyp_stif | 1.054 | 0.958 |
| 9 | `e_i0` | maximum void ratio at p = 0 | hyp_stif | 1.212 | 1.1017 |
| 10 | `alpha` | pycnotropy exponent (factor f_d) | hyp_stif | 0.8311 | 0.48 |
| 11 | `b_phi` | friction ramp: φ_c ← φ_c[(1−b_φ) + b_φ·ε_mon] (not the β of von Wolffersdorff) | hyp_stif | 0.1114 | 0.1 |
| 12 | `m_R` | stiffness multiplier after a 180° reversal | M_istr, f_hyp | 2.2 | 1.62 |
| 13 | `m_T` | stiffness multiplier after a 90° change | M_istr | 1.1 | 1.10 |
| 14 | `R` | saturation norm of the intergranular strain (ρ = ‖δ‖/R) | M_istr | 1e-4 | 1e-4 |
| 15 | `beta_r` | evolution of the intergranular strain | M_istr | 0.15 | 0.15 |
| 16 | `chi_0` | χ₀, exponent χ at zero cyclic memory | M_istr | 5.5 | 5.5 |
| 17 | `c` | cohesion (off) | hyp_stif | 0 | 0 |
| 18 | `r` | v_i = K_i/G_i | hyp_stif | 1.0902 | 1.13 |
| 19 | `q_br` | not read | — | 9.9575 | 2.0 |
| 20 | `c_z` | fabric: rate | hyp_stif | 500 | 300 |
| 21 | `z_max` | fabric: saturation | hyp_stif | 9.1084 | 12.0 |
| 22 | `c_l` | semifluidized state: rate | hyp_stif | 24 | 40 |
| 23 | `n_l` | semifluidized state: exponent | hyp_stif | 16 | 8 |
| 24 | `p_th` [kPa] | semifluidized state: threshold mean stress | hyp_stif | 10 | 10 |
| 25 | `c_r` | semifluidized state: volumetric term (inactive in undrained tests) | hyp_stif | 1000 | 1000 |
| 26 | `lambda_1` | semifluidized state: magnitude | hyp_stif | 0.627986 | 0.07 |
| 27 | `lambda_2` | semifluidized state: density sensitivity | hyp_stif | 1.5 | 3.0 |
| 28 | `e0` | documentary, not read (the void ratio enters through `y(13)`) | — | 0.798 | 0.7427 |
| 29 | `chi_max` | χ_max, exponent χ at saturated cyclic memory | M_istr | 10 | 10 |
| 30 | `Omega` | shape of the failure surface in the π-plane: 0 Drucker–Prager, 1 Matsuoka–Nakai | hyp_stif | 0 | 0 |
| 31 | `a` | pressure factor (p_in/p_inr)^a of Barrero et al. (2020); **must stay 0** in this release (Liao et al. 2022 Eq. 10) | hyp_stif | 0 | 0 |
| 32 | `p_inr` | only with a ≠ 0 | hyp_stif | 0 | 0 |
| 33 | `p_in` | only with a ≠ 0 | hyp_stif | 0 | 0 |
| 34 | **η** (`eta_rate_cyc`) | rate of ε_cyc | M_istr | 20 | 0.35 |
| 35 | **κ** (`kappa_rate_mon`) | rate of ε_mon | M_istr | 9.9575 | 15 |
| 36 | `y_h` | form of the cyclic detector: **must be 1**, the dimensionless y_δ of the thesis and the manuscript (code variable `yh`); any other value, or a vector shorter than 36 slots, runs an earlier projection −ρ^χ N:Δε that is not the published model | M_istr | 1 | 1 |
| 37 | `c_h` | retired (h_pt = 1) | hyp_stif | 0 | 0 |
| 38 | `h0` | inactive with c_h = 0 | hyp_stif | 0.25 | 0.25 |
| 39 | `g` | experimental floor on ρ^χ in the N term, not part of the published model: **must be 0** | M_istr | 0 | 0 |
| 40 | `k_chi` | experimental state dependence of χ_max, not part of the published model: **must be 0** | f_hyp | 0 | 0 |
| 41 | `psi_ref` | inactive with k_chi = 0 | f_hyp | −0.218 | −0.218 |
| 42 | `tent` | inactive with c_h = 0 | hyp_stif | 2 | 2 |
| 43 | `f_l` | optional; when absent or 0 the core uses f_l = 0.01 (Liao et al. 2022) | hyp_stif | — | — |

The monotonic set of KFS (`params_KFS('monotonic')`) differs in slots 2 (1e-4), 22 (7.4), 23 (8),
26 (0.06369) and 27 (4.5). The `'monotonic'` set of BBS differs in slots 7–9 (0.610, 0.990, 1.1385),
12 (2.2), 21 (8.0), 26 (0.15), 34 (20), the documentary slot 28 and z(0) = 0.5 (−1,−1,2,0,0,0)/√6.

## The state vector `y` (31)

| index | content |
|---|---|
| 1:6 | strain ε (Voigt order 11 22 33 12 23 13; shear strains engineering) |
| 7:12 | effective stress σ [kPa] |
| 13 | void ratio e |
| 14:19 | intergranular strain δ |
| 20:25 | fabric tensor z |
| 26 | variable l of the semifluidized state |
| 27 | not used |
| 28 | **ε_cyc** (`Strain_total`) |
| 29 | not used |
| 30 | **ε_mon** (`Strain_total2`) |
| 31 | not used |

Outputs of `updateModel`, `updateModel_stress` and `tools/cyclic_test.m`: `HARD(:,k) = y(12+k)` in the
laboratory convention, so e = `HARD(:,1)`, ε_cyc = `HARD(:,16)` and ε_mon = `HARD(:,18)`. y_δ is not an
output.

## Cyclic memory: equation to code

The equations of the cyclic memory (README, section *Cyclic memory and the detector y_δ*) and the lines
of `model/` that implement them. `model/` is fixed by `docs/MD5.csv`, so the line numbers are stable.

| quantity | code | where |
|---|---|---|
| metric for Voigt vectors with engineering shear, M2 = diag(1,1,1,½,½,½) | `M2` | `M_istr.m` L7–12 |
| χ = χ₀ + ε_cyc (χ_max − χ₀) | `chi`; χ₀ = `parms(16)`, χ_max = `parms(29)`, ε_cyc = `Strain_total` = y(28) | `M_istr.m` L25, L28 |
| ρ = ‖δ‖/R, at most 1 | `istrain_norm`, `rho` | `M_istr.m` L32, L41–45 |
| δ̂ = δ/‖δ‖ | `istrain_dir` | `M_istr.m` L36–38 |
| δ̂ : ε̇ | `istraindir_D` | `M_istr.m` L40 |
| ‖ε̇‖ | `dev_strain_rate_norm` (the full norm, despite the name), `norm_deps_M2` | `M_istr.m` L73, L87 |
| y_δ = ρ^χ ⟨δ̂ : ε̇⟩ / ‖ε̇‖ (parms(36) = 1) | `yh` | `M_istr.m` L86–92 |
| η | `eta_rate_cyc` = `parms(34)` | `M_istr.m` L104–108 |
| dε_cyc/dt = η (1 − y_δ − ε_cyc) ‖ε̇‖ | `d_strain` | `M_istr.m` L117 |
| ρ^χ in the stiffness (intergranular strain) | `Mep` | `M_istr.m` L51, L61–63, L67 |
| integration of ε_cyc | `F(28,1) = dStrain` | `f_hyp.m` L74, L83, L119, L138 |
| ε_cyc in the error norm of the integrator | `active_idx` | `hypo_update.m` L84 |

ε_cyc is not mirrored at the boundary (`espejo_estado.m` L15 changes slots 1:12 and 14:25 only), and
δ̂ : ε̇ and ‖ε̇‖ are the same in both sign conventions, so y_δ is too.

## Sign convention

The constitutive equations in `model/` (`hyp_stif.m`, `M_istr.m`, `f_hyp.m`, `hypo_update.m`,
`inv_s.m`, `inv_e.m`) work with compression **negative**, as von Wolffersdorff (1996), Niemunis &
Herle (1997) and Liao et al. (2022). Everything outside (the state `y0`, the loading path, the output
arrays and the data) uses compression **positive**. The drivers `updateModel.m` and
`updateModel_stress.m` cross the boundary with `espejo_estado.m`, which changes the sign of the
strain, stress, intergranular strain and fabric slots (an exact involution).
Initial states in the laboratory convention: KFS δ(0) = −0.82 R/√3 (1,1,1,0,0,0) and
z(0) = 0.5 (−1,−1,2,0,0,0)/√6; BBS δ(0) = −0.8 R/√3 (1,1,1,0,0,0) and z(0) = 0.

`inv_s.m` returns p = −tr σ/3 > 0, q and sin 3θ of the model stress; `inv_e.m` returns the
volumetric strain (positive in contraction) and the shear strain.
