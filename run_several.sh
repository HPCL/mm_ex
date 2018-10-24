#!/bin/bash

module load intel/19 caliper papi

export OMP_NUM_THREADS=8
export CALI_CONFIG_FILE=/home/users/gravelle/soft/src/Caliper/examples/configs/papi_cycles.conf

export CALI_REPORT_FILENAME=original_cyc.json
export CALI_PAPI_COUNTERS=PAPI_TOT_CYC
./original



