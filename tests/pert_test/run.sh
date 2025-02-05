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
rm -rf data && mkdir -p data || exit

# Create perturbation file...
$cris/perturbation - data/pert.nc \
		   ../data/SNDR.SNPP.CRIS.20220115T0000.m06.g001.L1B.std.v03_08.G.220115064421.nc

# Loop over channel sets...
for pert in 4mu 15mu_low 15mu_high ; do

    # Get map data...
    $cris/map_pert - data/pert.nc data/map_$pert.tab PERTNAME $pert
    
    # Estimate noise...
    $cris/noise_pert - data/pert.nc data/noise_$pert.tab PERTNAME $pert
    
    # Get variance...
    $cris/variance - data/var_$pert.tab data/pert.nc PERTNAME $pert NX 60 NY 30
    
done

# Compare files...
echo -e "\nCompare results..."
error=0
for f in $(ls data.ref/*.nc data.ref/*.tab) ; do
    diff -q -s data/"$(basename "$f")" "$f" || error=1
done
exit $error
