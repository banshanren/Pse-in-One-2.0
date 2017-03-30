#!/bin/bash

python pse.py temp/data_pos.fasta temp/pse_pos_svm.txt DNA PseDNC -lamada 3 -w 0.05 -f svm -l +1
python pse.py temp/data_neg.fasta temp/pse_neg_svm.txt DNA PseDNC -lamada 3 -w 0.05 -f svm -l -1
python feature_combine.py temp/pse_pos_svm.txt temp/
# python train.py temp/pos_svm.txt temp/neg_svm.txt -m dna.model -c 1 9 2 -g -3 1 2 -v 10
