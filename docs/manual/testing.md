# Testing

The repository includes regression tests under `tests/`.

## Available Test Suites

- `pert_test`
- `rad_test`
- `volc_test`

Run all test suites:

```bash
make -C src check
```

Run a single suite:

```bash
make -C src pert_test
make -C src rad_test
make -C src volc_test
```

## What The Tests Do

The test scripts:

- set `LD_LIBRARY_PATH` to the locally built libraries in `libs/build/lib`
- set `OMP_NUM_THREADS=4`
- reconstruct the bundled split Level-1B test archive in `tests/data/`
- create fresh output directories inside each suite
- compare generated output against `data.ref/`

Because the tests regenerate files locally, running them changes contents under
the corresponding `tests/*/data/` directories.

## Test Data

The test input archive is stored in:

- `tests/data/cris_l1b.zip`
- `tests/data/cris_l1b.z01`

The suite scripts combine and unpack the archive before executing the binaries.

## When To Run Which Test

- use `rad_test` after changes related to `spec2tab` or `map_rad`
- use `pert_test` after changes related to perturbation, noise, mapping, or
  variance logic
- use `volc_test` after changes related to volcanic index calculations
- run `make -C src check` after shared-library changes such as edits in
  `libcris.c`
