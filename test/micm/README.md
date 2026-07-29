# MICM shim known-answer harness

Compiles the real `src/chemistry/mozart/mo_micm.F90` against minimal stubs
of the CAM modules it uses (configured for the `trop_mam4` mechanism) and
runs it on the generated MICM configuration
(`src/chemistry/pp_trop_mam4/micm/`, produced by
`tools/micm/gen_micm_config.py`). With rates held fixed over the step —
exactly what the strategy-A rate injection does — trop_mam4 is linear in
the solution species, so the result of one solve is checked against exact
closed-form solutions per species (plus a total-sulfur invariant), for
every cell of a multi-cell chunk with per-cell-distinct rates, on both the
full-block and padded-lanes code paths.

## Build

Requires a MUSICA install built with the same options as
`cime_config/buildlib` `_build_musica` (MICM only, Fortran interface,
vendored deps from `libraries/`):

```
cmake <CAM>/libraries/musica -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=$PWD/install -DCMAKE_INSTALL_LIBDIR=lib \
  -DMUSICA_BUILD_FORTRAN_INTERFACE=ON -DMUSICA_ENABLE_TUVX=OFF \
  -DMUSICA_ENABLE_CARMA=OFF -DMUSICA_ENABLE_MIAM=OFF -DMUSICA_ENABLE_MIEM=OFF \
  -DMUSICA_ENABLE_TESTS=OFF -DMUSICA_SET_MICM_DEFAULT_VECTOR_SIZE=128 \
  -DFETCHCONTENT_TRY_FIND_PACKAGE_MODE=NEVER \
  -DFETCHCONTENT_SOURCE_DIR_MICM=<CAM>/libraries/micm \
  -DFETCHCONTENT_SOURCE_DIR_MECHANISM_CONFIGURATION=<CAM>/libraries/mechanism-configuration \
  -DFETCHCONTENT_SOURCE_DIR_FMT=<CAM>/libraries/fmt \
  -D"FETCHCONTENT_SOURCE_DIR_YAML-CPP=<CAM>/libraries/yaml-cpp"
cmake --build . --target install -j 8

gfortran -DMICM -O1 -fallow-argument-mismatch -o harness_micm \
  cam_stubs.F90 mpi_stub.F90 ../../src/chemistry/mozart/mo_micm.F90 harness_micm.F90 \
  -Iinstall/include/musica/fortran -Linstall/lib \
  -lmusica-fortran -lmusica -lmechanism_configuration -lyaml-cpp -lfmt -lstdc++
```

## Run

Write a `harness.nml`:

```
&micm_opts
 micm_active = .true.
 micm_config_path = '<CAM>/src/chemistry/pp_trop_mam4/micm'
 micm_rxt_map_path = '<CAM>/src/chemistry/pp_trop_mam4/micm/rxt_map.txt'
 micm_solver_type = 'rosenbrock'
/
```

then `./harness_micm harness.nml` — prints `PASS` or per-cell mismatches.
`rosenbrock` and `rosenbrock_standard` pass at the 2e-3 relative
tolerance; `backward_euler` converges but carries the expected first-order
accuracy error (~0.3–0.8% at an 1800 s step), so it fails the Rosenbrock-
calibrated tolerance — informational, not a defect.

## Negative controls

The test fails (verified) when the reaction map is corrupted, e.g.:

- reactant-count exponent wrong (`8 1 USER.SO2_OH_M` -> `8 2 ...`):
  SO2/S-total mismatches
- EMISSION yield wrong (`5 0 EMIS.usr_HO2_HO2 1.0` -> `... 2.0`):
  H2O2 mismatches
- unknown rate parameter name: init aborts
- header count mismatch: init aborts

# Strategy A vs B equivalence harness (harness_ab_t1s1.F90)

Validates the `--native` generator mode for `trop_strat_mam5_t1s1`:
compiles the real generated t1s1 rate code (chem_mods, mo_sim_dat,
mo_setrxt, mo_adjrxt, mo_phtadj, mo_jpl, mo_tracname) to produce CAM-truth
post-adjrxt rate constants, then solves identical initial states through
mo_micm with the all-injected config (`micm/`) and the native-rate-law
config (`micm_native/`, 368 Arrhenius/Troe/photolysis primitives). The
pass criterion is per-cell tendency agreement at a small step (canonical
`delt = 1e-4 s`), which isolates the rate representation from
stiff-trajectory sensitivity.

Build: as above, but with `cam_stubs_ab.F90` (no chem_mods/mo_tracname
stubs) plus the real t1s1 files listed above, then

```
./harness_ab_t1s1 A.nml B.nml <pp>/micm/rxt_map.txt 1.e-4
```

where A.nml/B.nml are micm_opts namelists pointing at `micm/` and
`micm_native/` respectively (`micm_abort_on_nonconvergence = .false.`).
Verified results: clean configs agree to 2.75e-7 in concentrations with
all tendencies within 1e-3; a +10% perturbation of a single native
Arrhenius coefficient produces a 1.1e-1 tendency divergence (FAIL);
365/368 native coefficients independently match MUSICA's own v0 TS1
configuration (the 3 differences are newer C3H7O2 JPL constants in CAM's
current mechanism).
