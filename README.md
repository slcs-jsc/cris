# CrIS Code Collection

The CrIS Code Collection enables data processing and analysis for remote sensing observations captured by the Cross-track Infrared Sounder satellite instruments.

[![release (latest by date)](https://img.shields.io/github/v/release/slcs-jsc/cris)](https://github.com/slcs-jsc/cris/releases)
[![commits since latest release (by SemVer)](https://img.shields.io/github/commits-since/slcs-jsc/cris/latest)](https://github.com/slcs-jsc/cris/commits/master)
[![last commit](https://img.shields.io/github/last-commit/slcs-jsc/cris.svg)](https://github.com/slcs-jsc/cris/commits/master)
[![top language](https://img.shields.io/github/languages/top/slcs-jsc/cris.svg)](https://github.com/slcs-jsc/cris/tree/master/src)
[![code size](https://img.shields.io/github/languages/code-size/slcs-jsc/cris.svg)](https://github.com/slcs-jsc/cris/tree/master/src)
[![repo size](https://img.shields.io/github/repo-size/slcs-jsc/cris.svg)](https://github.com/slcs-jsc/cris/tree/master/src)
[![codacy](https://api.codacy.com/project/badge/Grade/a9de7b2239f843b884d2a4eb583726c9)](https://app.codacy.com/gh/slcs-jsc/cris?utm_source=github.com&utm_medium=referral&utm_content=slcs-jsc/cris&utm_campaign=Badge_Grade_Settings)
[![codecov](https://codecov.io/gh/slcs-jsc/cris/branch/master/graph/badge.svg?token=4X6IEHWUBJ)](https://codecov.io/gh/slcs-jsc/cris)
[![tests](https://img.shields.io/github/actions/workflow/status/slcs-jsc/cris/tests.yml?branch=master&label=tests)](https://github.com/slcs-jsc/cris/actions)
[![docs](https://img.shields.io/github/actions/workflow/status/slcs-jsc/cris/docs.yml?branch=master&label=docs)](https://slcs-jsc.github.io/cris)
[![license](https://img.shields.io/github/license/slcs-jsc/cris.svg)](https://github.com/slcs-jsc/cris/blob/master/COPYING)
[![doi](https://zenodo.org/badge/DOI/10.5281/zenodo.14790896.svg)](https://doi.org/10.5281/zenodo.14790896)

## Installation

This documentation describes the installation on a Linux system.
A number of standard tools such as the GNU Compiler Collection (gcc)
and 'make' are required for installation.

Start by downloading the source code from the git repository:

    git clone https://github.com/slcs-jsc/cris.git

Change to the directory cris/ which holds source codes,
libraries, documentation, etc:

    cd cris

The GNU Scientific Library (https://www.gnu.org/software/gsl) is
required for numerical calculations and the Unidata netCDF library
(http://www.unidata.ucar.edu/software/netcdf) is needed for file-I/O.
Copies of these libraries can be found in the repository, if they are
not available on your system. A script is provided to build the
libraries:

    cd lib
    ./build.sh

Next, change to the source directory and edit the Makefile according
to your needs. In particular, check the paths to the libraries (INCDIR
and LIBDIR). Then try to compile the code:

    cd ../src
    emacs Makefile
    make

The binaries will be linked statically, i.e., they can be copied to
other machines. Sometimes static compilations causes problems, in
particular in combination with MPI. In this case remove the '-static'
flag from the CFLAGS in the Makefile and compile again.

By default we use rather strict compiler warnings.  All warning
messages will be turned into errors and no binaries will be
produced. This behavior is enforced by the flag '-Werror'.

The binaries will remain in the src/ directory.

## License

The CrIS Code Collection is distributed under the GNU GPL v3.
Software libraries distributed along with this software package may
have their own licenses and copyrights, please see corresponding
documentation.

## Contact

We are interested in sharing the CrIS Code Collection for research applications.

Please do not hesitate to contact us if you have any further questions:

Dr. Lars Hoffmann

Jülich Supercomputing Centre, Forschungszentrum Jülich

e-mail: <l.hoffmann@fz-juelich.de>
