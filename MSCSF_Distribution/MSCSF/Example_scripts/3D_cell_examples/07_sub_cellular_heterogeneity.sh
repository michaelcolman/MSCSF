#!/bin/sh

# example using one map file and setting SERCA and LTCC to this map, while setting RyR to random

./model_single_3D SERCA_het On SERCA_map_file GRF_example_map.dat RyR_het random LTCC_het map LTCC_het_map_file GRF_example_map.dat Reference Sub_cellular_het_test
