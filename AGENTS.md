# AGENTS.md

## Purpose

This repository contains the CrIS Code Collection, a C-based toolkit for processing and analyzing Cross-track Infrared Sounder satellite data.

Agents working in this repository should prefer small, verifiable changes that preserve existing command-line behavior and numerical outputs.

## Repository Layout

- `src/`: C sources, headers, and the main `Makefile`. Binaries are built here.
- `tests/`: regression-style test suites and input/reference data.
- `libs/`: bundled source archives plus `build.sh` for local dependency builds.
- `docs/`: MkDocs and Doxygen configuration plus manual content.
- `README.md`: installation and build overview.

## Build Workflow

The main build happens in `src/Makefile`.

- Build executables: `make -C src`
- Run all tests: `make -C src check`
- Run one test suite: `make -C src pert_test`, `make -C src rad_test`, or `make -C src volc_test`
- Clean build artifacts: `make -C src clean`
- Build docs: `make -C src doxygen` or `make -C src mkdocs`

Run `make -C src clean` before rebuilding with changed compiler options; Make does not track command-line flag changes. Run `make -C src check` without `-j`.

## Dependencies

The code expects GSL and netCDF headers and libraries under `../libs/build` by default.

- Include path is configured with `INCDIR += -I ../libs/build/include`
- Library path is configured with `LIBDIR += -L ../libs/build/lib`
- If dependencies are missing, check whether they already exist on the system before changing the `Makefile`
- Bundled library sources live in `libs/`, and `libs/build.sh` is the project-provided way to build local copies

Do not casually change compiler flags, include paths, or linker settings in `src/Makefile` unless the task is specifically about build configuration.

## Test Behavior

The regression tests in `tests/*/run.sh` do the following:

- set `LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH`
- set `OMP_NUM_THREADS=4`
- reconstruct `tests/data/cris_l1b.zip` into a full test input file
- regenerate outputs into a local `data/` directory
- compare outputs against `data.ref/`

Because these tests rewrite their local `data/` directories and unpack inputs, expect filesystem changes under `tests/` when running them.

## Coding Conventions

Match the existing style in `src/`:

- Keep code in C, not C++
- Preserve the current brace and indentation style already present in each file
- Prefer existing naming patterns such as `snake_case` for functions and short, domain-specific variable names
- Keep comments brief and functional; many files use block separators and short inline comments ending with `...`
- Avoid introducing new dependencies or large abstractions in this codebase
- Before changing behavior, check existing callers, tests, and documented expectations; question assumptions and explain non-obvious choices
- Reuse existing data structures where possible; check memory and affected data paths before adding a `struct`
- Avoid new `static` helper functions; use a macro for a small operation when its arguments can be evaluated safely, or a regular function declared and documented in `src/libcris.h`
- Keep public declarations synchronized between `.h` and `.c` files
- Format C sources and headers from `src/` with `make indent` after a successful build; review the diff because the target rewrites all `*.c` and `*.h` files, including copied JURASSIC sources

This project builds with strict warnings and `-Werror`, so any change should compile cleanly without new warnings.

## Change Guidelines

- Prefer minimal diffs in numerics-heavy code
- Investigate test failures before changing reference files; update only affected baselines after verifying and explaining the new results
- Do not remove bundled libraries or large test assets
- Be careful with memory usage: several structures allocate large static arrays
- Preserve CLI compatibility for existing executables in `src/`

## Documentation

Keep the MkDocs manual in `docs/manual/` and Doxygen API documentation current with relevant code, build, and workflow changes. Build affected docs with `make -C src mkdocs` and/or `make -C src doxygen` when the tools are available. For new or revised API Doxygen comments, use a substantive description beyond `@brief`, document parameters and behavior where relevant, and include `@author` with the responsible author's name. Follow the detailed `/** ... */` examples in `src/jurassic.h` without editing that copied dependency for style alone; do not invent attribution.

## Verification Expectations

When making code changes, verify as close to the modified area as possible:

- build with `make -C src`
- run the most relevant regression suite
- if a change affects shared code such as `libcris.c`, prefer `make -C src check`

If verification cannot be run because dependencies are unavailable or the environment blocks it, state that explicitly.
