# SHOC "in-and-out" standalone driver and ERF integration notes

`shoc_in_and_out.cpp` runs a single-column, SHOC-only DYCOMS RF01 simulation
from the AMR team's ASCII initial-condition files, mirroring EAM's
`turb_standalone` ("in-and-out") mode. It is the reference for making a host
model that embeds the EAMxx C++ SHOC (e.g. ERF) reproduce the E3SM Fortran
in-and-out runs.

## Usage

```
shoc_in_and_out -d <dir with ShocInOut_IC_*.txt> [-f] [-dt 60] [-s 360]
                [-x exp1|exp2] [-dx 100000] [-o prefix]
```

- default engine: EAMxx C++ `shoc_main`; `-f`: the reference Fortran
  `shoc.F90` compiled into EAMxx (with the thv fix).
- `-x exp2` (default): `-s` calls to `shoc_main` with `nadv=1`, per-substep
  text output after each call (matches EAM `l_shoc_outer_loop=.true.` and the
  natural ERF call pattern of one call per host timestep).
- `-x exp1`: one call with `nadv=<s>`; only the final state is written.

Output format is byte-compatible with the Fortran writer
(`(I5,1X,I4,5(1X,F17.14),1X,F16.12,1X,F17.10,3(1X,F12.4))`), columns:
time, k (0 = top), u, v, tke, qv, qc, T, P, lscale×3 (zeros here — the
Larson length scale is an EAM-side diagnostic).

## What a host (ERF) must do for an apples-to-apples run

1. **Call `shoc_main` directly.** Do not go through the EAMxx atmosphere
   process (`eamxx_shoc_process_interface`): its `SHOCPreprocess` recomputes
   thv/dz/grids/fluxes from host fields, applies TMS
   (`apply_tms`), flux-consistency checks and post-process clamping — none of
   which exist in the EAM Fortran in-and-out path.
2. **Feed the ASCII IC values verbatim** (the same
   `ShocInOut_IC_{surface_vars,zi_grid,zt_grid}.txt` files, 17-significant-
   digit v4+ format). Files are surface-first; `shoc_main` in both languages
   is top-first (k=0 = model top) — flip rows.
3. **Fixed across all substeps/calls:** zt/zi grids, `pres`/`presi`/`pdel`,
   `inv_exner`, `w_field` (wm_zt), the four surface fluxes, `dx`/`dy`
   (1e5 m for these runs), `phis = 0`, `wtracer_sfc = 0`.
4. **Prognostic (carry between calls, never reset):** `tke` (clip at
   `mintke = 4e-4` before the first call only), `thetal`, `qw`, `u_wind`,
   `v_wind`, `shoc_ql`, `shoc_cldfrac`, `tk`, `tkh`, `wthv_sec`
   (initialize to 0), `host_dse` (initialize as `cp*T + g*z` with
   `T = thetal/inv_exner + (Lv/cp)*shoc_ql`; it does not feed back into the
   prognostics — only the energy fixer adjusts it).
5. **Runtime options:** EAM defaults, identical to what
   `shoc_functions_f90.cpp::shoc_main_f` hardcodes:
   `lambda_low=0.001, lambda_high=0.04, lambda_slope=2.65,
   lambda_thresh=0.02, thl2tune=qw2tune=qwthl2tune=w2tune=1.0,
   length_fac=0.5, c_diag_3rd_mom=7.0, Ckh=Ckm=0.1`. Physics constants in
   `physics_constants.hpp` already match EAM's `physconst` values.
   `npbl = nlev` (all levels of these 3-km-domain grids are below the
   400 hPa PBL pressure cap).
6. **Carry the thv-staleness fix.** `shoc_main_impl.hpp` on this branch
   recomputes `shoc_thv = shoc_tabs*inv_exner*(1 + zvir*shoc_qv - shoc_ql)`
   each substep and uses it in `shoc_length` instead of the (stale) `thv`
   input. An ERF copy of the EAMxx SHOC taken from upstream lacks this —
   apply the same change (commit 63e246b4ee).
7. **Precision:** build with `SCREAM_DOUBLE_PRECISION=ON` (default).

## Validated agreement (DYCOMS RF01, dz = 10/20/50/100 m, dt = 60 s, 6 h)

- C++ vs Fortran engine, same driver, same ICs: max relative differences
  1e-14 .. 4e-11 over the full run, all grids, exp1 and exp2; pressure
  bit-identical. exp1 and exp2 produce identical prognostics (the energy
  fixer only touches host_dse).
- Full E3SM SCM (Intel, `l_turb_standalone=.true.`, this branch) vs this
  driver (GNU): bounded ~1e-7 .. 2e-5 relative over 6 h — pure
  compiler-transcendental noise (at step 1, u/v/tke are bit-identical; only
  the PDF-diagnosed qc and its derivatives differ).

## E3SM side reminders

- The EAM run scripts must set `l_turb_standalone = .true.` in `user_nl_eam`
  on branches where it is a namelist flag (older branches had the mode
  hardwired; without the flag the case silently runs as a normal SCM step).
- exp1 vs exp2 in EAM: `l_shoc_outer_loop = .false./.true.` with
  `shoc_timestep = 60` and one 21600 s host step.
