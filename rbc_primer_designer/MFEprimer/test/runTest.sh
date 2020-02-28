#!/bin/sh
#../MFEprimer.py -i primers.fa -d ./database.seq -T F -e 10000 -W 11 -a 2 -s 0.2
#../MFEprimer.py -i primers.fa -d ./database.seq -T F -e 10000 -W 11 -a 2 --ppc_cutoff=0.2 --seq=primer.mfe.fa -o primer.mfe
../MFEprimer.py -i primers.fa -d ./database.seq -T F -e 10000 -W 4 -a 2 --ppc_cutoff=0.2 --seq=primer.mfe.fa -o primer.mfe
