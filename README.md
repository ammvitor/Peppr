# Peppr
Peppr - Automated pipeline for peptide generation, docking and scoring

For test
 cd ./test
 rm -rf *
 cp ../Peppr_Stage2/test/test.trg .
 python ../Peppr_Stage1/Peppr_main.py -c 2 -i TLXDYL -t test.trg -p 50 -N 2 -n 800
 python ../Peppr_Stage2/__main__sulphotyrosine_2.py
