# CrIS Code Collection

The CrIS Code Collection provides command-line tools for processing and analyzing
remote sensing observations from the Cross-track Infrared Sounder (CrIS)
instruments.

The repository is centered around a set of C executables in `src/` that read
CrIS Level-1B granules, derive tabular products, and generate perturbation and
variance fields for downstream analysis.

## What This Project Contains

The current codebase provides tools for:

- converting between calendar and mission time representations
- extracting spectra or geolocated radiance maps from CrIS Level-1B files
- generating perturbation products in netCDF format
- estimating perturbation noise and variance products
- deriving volcanic ash and SO2 indicator tables

## Repository Structure

- `src/`: source files, headers, and the main build targets
- `tests/`: regression tests and reference outputs
- `libs/`: bundled third-party library sources and local build helper
- `docs/`: MkDocs manual and Doxygen configuration

## Typical Workflow

1. Build the bundled libraries if system libraries are not already available.
2. Compile the executables in `src/`.
3. Run the regression tests in `tests/`.
4. Use the command-line tools on one or more CrIS Level-1B input granules.

## Key Tools

The most commonly used programs are:

- `spec2tab`: extract a single spectrum for one footprint
- `map_rad`: extract a geolocated map for one wavenumber
- `perturbation`: generate perturbation data products in netCDF format
- `map_pert`: extract perturbation map data from a perturbation file
- `noise_pert`: estimate perturbation noise
- `variance`: compute gridded perturbation variance products
- `volcano`: derive volcanic ash and SO2 indicator tables

See the [Tools](tools.md) page for details and examples.

## Documentation

- Use this MkDocs manual for installation, workflow, and tool overviews.
- Use the Doxygen documentation for API-level details of the C code and data
  structures.

## Contact

Research users who need support or want to discuss the software can contact:

Dr. Lars Hoffmann, <l.hoffmann@fz-juelich.de>  
Jülich Supercomputing Centre, Forschungszentrum Jülich
