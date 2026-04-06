# Installation

## Requirements

The project is developed for Linux and expects a standard C toolchain, including:

- `gcc`
- `make`

The code also depends on:

- GNU Scientific Library (GSL)
- netCDF-C

Bundled source archives for the required libraries are stored in `libs/`.

## Clone The Repository

```bash
git clone https://github.com/slcs-jsc/cris.git
cd cris
```

## Build Local Dependencies

If compatible system libraries are not already installed, build the bundled
copies:

```bash
cd libs
./build.sh
```

By default the source `Makefile` expects headers and libraries under
`../libs/build`.

## Compile The Code

Build the executables from the `src/` directory:

```bash
make -C src
```

The resulting binaries remain in `src/`.

## Build Notes

- The default build is strict and uses `-Werror`, so warnings stop the build.
- Static linking can be enabled with `STATIC=1`.
- If static linking causes problems on a target system, rebuild without it.
- Include and library paths are controlled in `src/Makefile`.

## Build Documentation

Build the MkDocs manual:

```bash
make -C src mkdocs
```

Build the Doxygen documentation:

```bash
make -C src doxygen
```
