#! /bin/bash

# Set environment...
export LD_LIBRARY_PATH=../../libs/build/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=4

# Setup...
cris=../../src

# Uncompress test data...
cd ../data \
    && zip -s 0 cris_l1b.zip --out cris_l1b_full.zip \
    && unzip -o cris_l1b_full.zip \
    && cd - || exit

# Create directory...
rm -rf data && mkdir -p data

# Extract spectrum...
$cris/spec2tab - ../data/SNDR.SNPP.CRIS.20220115T0000.m06.g001.L1B.std.v03_08.G.220115064421.nc data/spec_noapo.tab
$cris/spec2tab - ../data/SNDR.SNPP.CRIS.20220115T0000.m06.g001.L1B.std.v03_08.G.220115064421.nc data/spec_apo.tab APO 1

# Extract map...
$cris/map_rad - data/map_noapo.tab ../data/SNDR.SNPP.CRIS.20220115T0000.m06.g001.L1B.std.v03_08.G.220115064421.nc NU 962.5
$cris/map_rad - data/map_apo.tab ../data/SNDR.SNPP.CRIS.20220115T0000.m06.g001.L1B.std.v03_08.G.220115064421.nc NU 962.5 APO 1

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.tab) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
