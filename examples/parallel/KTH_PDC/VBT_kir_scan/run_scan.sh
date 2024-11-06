#!/bin/bash

python setup_kir_scan.py networks/wt_kir_scan_1

export sim_time=0.2

for kir_factor in $(seq 1.0 -0.2 0.0); do
    echo "Simulating kir_factor $kir_factor"
    python simulate_kir_scan.py networks/wt_kir_scan_1 --kir_factor $kir_factor --time $sim_time --output networks/wt_kir_scan_1/simulation/output_kir_$kir_factor.hdf5
done



