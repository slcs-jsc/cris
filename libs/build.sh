#! /bin/bash

# Setup...
#rm -rf build
target=$(mkdir -p build && cd build && pwd)

# Prepare directories...
mkdir -p $target/src $target/bin $target/lib $target/man/man1 \
    && cp *tar.gz $target/src \
    && cd $target/src \
    && for f in $(ls *tar.gz) ; do tar xvzf $f ; done \
    || exit

# GSL...
dir=gsl-2.7
cd $target/src/$dir \
    && ./configure --prefix=$target \
    && make && make check && make install && make clean \
	|| exit

# zlib...
dir=zlib-1.2.12
cd $target/src/$dir \
    && ./configure --prefix=$target \
    && make && make check && make install && make clean \
	|| exit

# HDF5...
dir=hdf5-1.12.1
cd $target/src/$dir \
    && ./configure --prefix=$target --with-zlib=$target --enable-hl \
    && make && make check && make install && make clean \
	|| exit

# netCDF...
dir=netcdf-c-4.8.1
cd $target/src/$dir \
    && CPPFLAGS=-I$target/include LDFLAGS=-L$target/lib ./configure --prefix=$target --disable-dap --disable-nczarr \
    && make && make install && make clean \
	|| exit
