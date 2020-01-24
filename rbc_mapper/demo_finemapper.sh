#!/bin/sh

#1. run without optimization
python rbc_finemapper.py -f ./../../data/rbci-files/pools_only_on_9311_6MM_PCRpositions_filter.rbci -o ./out/pools_only_on_9311_6MM_PCRpositions_filter.rbci --n_res 178 --n_sus 280 --dp-max-ratio 1.5 --filter-flags PASS --min-qual 200 --chromosome 9311chr07 --resolution 1000 --window_size 1E8 --phenoNoise 1E-10


#1. run with recombination optimization
python rbc_finemapper.py -f ./../../data/rbci-files/pools_only_on_9311_6MM_PCRpositions_filter.rbci -o ./out/pools_only_on_9311_6MM_PCRpositions_filter.rbci.OPTRECOMB --n_res 178 --n_sus 280 --dp-max-ratio 1.5  --filter-flags PASS --min-qual 200 --chromosome 9311chr07 --resolution 1000 --window_size 1E8 --phenoNoise 1E-10 --opt_recombination

#1. run with recombination optimization
python rbc_finemapper.py -f ./../../data/rbci-files/pools_only_on_9311_6MM_PCRpositions_filter.rbci -o ./out/pools_only_on_9311_6MM_PCRpositions_filter.rbci.OPTRECMOB.OPTNOISE --n_res 178 --n_sus 280 --dp-max-ratio 1.5 --filter-flags PASS --min-qual 200 --chromosome 9311chr07 --resolution 1000 --window_size 1E8 --phenoNoise 1E-10 --opt_recombination --opt_eps
