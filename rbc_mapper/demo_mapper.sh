#!/bin/sh

#testcode used January 2020
python rbc_mapper_vienna.py -f ../mapper_test.rbci -o ../mapper_test --n_res 178 --n_sus 280 --chromosome 9311chr07 --max-hs 0.1 --dp-max-ratio 1.5 --trust-gaTK --filter-flags ExactAN,NoSBtag --min-qual 200 --resolution 100000 --phenoNoise 1E-10

#old
#1. run without optimization
python rbc_mapper.py -f ./../../data/rbci-files/pools_only_on_9311_6MM.SW.dm.BOTH.phased.NewFilter.rbci -o ./out/pools_only_on_9311_6MM.SW.dm.BOTH.phased.NewFilter.rbc --n_res 178 --n_sus 280 --max-hs 0.1 --dp-max-ratio 1.5 --trust-gaTK --filter-flags PASS --min-qual 200 --chromosome 9311chr07 --resolution 100000 --phenoNoise 1E-10
#2. run with optimization of recombination
#python rbc_mapper.py -f ./../../data/rbci-files/pools_only_on_9311_6MM.SW.dm.BOTH.phased.NewFilter.rbci -o ./out/pools_only_on_9311_6MM.SW.dm.BOTH.phased.NewFilter.rbc.OPTRECOMB --n_res 178 --n_sus 280 --max-hs 0.1 --dp-max-ratio 1.5 --trust-gaTK --filter-flags PASS --min-qual 200 --chromosome 9311chr07 --resolution 100000 --phenoNoise 1E-10 --opt_recombination
#3. run with optimization of recombination and error varaince
#python rbc_mapper.py -f ./../../data/rbci-files/pools_only_on_9311_6MM.SW.dm.BOTH.phased.NewFilter.rbci -o ./out/pools_only_on_9311_6MM.SW.dm.BOTH.phased.NewFilter.rbc.OPTRECOMB.OPTNOISE --n_res 178 --n_sus 280 --max-hs 0.1 --dp-max-ratio 1.5 --trust-gaTK --filter-flags PASS --min-qual 200 --chromosome 9311chr07 --resolution 100000 --phenoNoise 1E-10 --opt_recombination --opt_eps
