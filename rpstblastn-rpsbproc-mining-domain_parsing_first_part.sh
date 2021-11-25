#!/bin/bash
#export CDDATA="/way/to/cddblast/ncbi-blast-N.N.N+-src/c++/ReleaseMT/bin/data"
source ~/.profile
source ~/.bashrc

#STEP 1 - finding aRH and its neighborhood
rpstblastn -query $1 -db $DARTS/customCDD/Step1 -outfmt 11 -out $2.step1.rpstblastn.out -num_threads 4

rpsbproc -i $2.step1.rpstblastn.out -o $2.step1.rpsbproc.result -e 0.0001 -d $CDDATA -m rep -t doms

python $DARTS/new_1step_rpsbalstn-rpsbproc_parser.py $1 $2.step1.rpsbproc.result $2 

rm $2.step1.rpstblastn.out

#STEP 2 - mining all domains in one file, elements in another 4 files with its Scores
rpstblastn -query aRNH_and_approximates_$2_genome.fa -db $DARTS/customCDD/Step3 -outfmt 11 -out $2.step3.rpstblastn.out -num_threads 4

rpsbproc -i $2.step3.rpstblastn.out -o $2.step3.rpsbproc.result -e 0.0001 -d $CDDATA -m rep -t doms

python $DARTS/new_fix_mining.py aRNH_and_approximates_$2_genome.fa $2.step3.rpsbproc.result $2 

#STEP - clustering and choosing delegates
python $DARTS/parse_prot_seq_domains.py prot_domains.fa $2_elements 

python $DARTS/replacing_star_to_X.py $2 'gRT' 
python $DARTS/replacing_star_to_X.py $2 'gRH' 
python $DARTS/replacing_star_to_X.py $2 'aRH' 
python $DARTS/replacing_star_to_X.py $2 'GAG' 
python $DARTS/replacing_star_to_X.py $2 'INT' 
python $DARTS/replacing_star_to_X.py $2 'PRo' 

