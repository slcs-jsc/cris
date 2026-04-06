# Tools

This page summarizes the main command-line programs built from `src/`.

Each executable expects a control-file style first argument named `<ctl>`.
In practice, the programs parse named parameters directly from the command line
via tokens such as `APO 1` or `NU 962.5`.

## Time Conversion Utilities

These small tools convert between date and mission-time representations:

- `day2doy`: convert `year month day` to day-of-year
- `doy2day`: convert `year doy` to calendar date
- `time2jsec`: convert date/time fields to mission seconds
- `jsec2time`: convert mission seconds back to calendar time

## Radiance Extraction

### `spec2tab`

Extracts a single CrIS spectrum for one footprint from a Level-1B granule.

Required arguments:

```text
spec2tab <ctl> <l1b_file> <spec.tab>
```

Supported selection controls visible in the source include:

- `TRACK`
- `XTRACK`
- `IFOV`
- `LON`
- `LAT`
- `APO`

If `LON` and `LAT` are provided, the program searches for the nearest covered
footprint in the granule.

Example:

```bash
src/spec2tab - tests/data/input.nc out/spec.tab TRACK 0 XTRACK 0 IFOV 0
```

### `map_rad`

Extracts all footprints for one selected wavenumber from one or more Level-1B
granules.

Required arguments:

```text
map_rad <ctl> <map.tab> <l1b_file1> [<l1b_file2> ...]
```

Supported controls visible in the source include:

- `NU`
- `APO`

Example:

```bash
src/map_rad - out/map.tab tests/data/input.nc NU 962.5
```

## Perturbation Products

### `perturbation`

Builds perturbation fields from one or more Level-1B files and writes a netCDF
output product.

Required arguments:

```text
perturbation <ctl> <pert.nc> <l1b_file1> [<l1b_file2> ...]
```

Supported controls visible in the source include:

- `APO`
- `BIAS`

This program computes multiple derived brightness-temperature products,
including 4 micron and 15 micron perturbation fields.

### `map_pert`

Extracts tabular map output from a perturbation netCDF file.

The regression tests call it in this form:

```bash
src/map_pert - out/map_4mu.tab out/pert.nc PERTNAME 4mu
```

### `noise_pert`

Estimates perturbation noise from a perturbation netCDF file.

The regression tests call it in this form:

```bash
src/noise_pert - out/pert.nc out/noise_4mu.tab PERTNAME 4mu
```

### `variance`

Computes gridded perturbation variance products from one or more perturbation
files.

Required arguments:

```text
variance <ctl> <var.tab> <pert1.nc> [<pert2.nc> ...]
```

Supported controls visible in the source include:

- `PERTNAME`
- `NX`, `NY`
- `LON0`, `LON1`
- `LAT0`, `LAT1`
- `THRESH_GW`
- `THRESH_DC`
- `DT_TROP`
- `DT230`
- `NU`
- `DC`
- `OUTPUT`

Example used by the regression tests:

```bash
src/variance - out/var_4mu.tab out/pert.nc PERTNAME 4mu NX 60 NY 30
```

## Volcanic Emission Diagnostics

### `volcano`

Writes tabular diagnostics for cloud, ash, and SO2 indicator channels from one
or more Level-1B granules.

Required arguments:

```text
volcano <ctl> <out.tab> <l1b_file1> [<l1b_file2> ...]
```

Supported controls visible in the source include:

- `APO`

Example:

```bash
src/volcano - out/volcano.tab tests/data/input.nc APO 1
```

## Output Formats

Most tools write either:

- plain text tables with self-describing column headers
- netCDF files for structured perturbation products

The text-table tools write column descriptions at the top of each output file,
which is the quickest way to inspect exact field definitions.
