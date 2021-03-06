#!/bin/bash
#export CDDATA="/way/to/cddblast/ncbi-blast-N.N.N+-src/c++/ReleaseMT/bin/data"
# variables: $1 - genome assembly file or its split, $2 - project name,
source ~/.profile
source ~/.bashrc

#STEP 1 - finding aRH and its neighborhood
rpstblastn -query $1 -db $DARTS/customCDD/Step1 -outfmt 11 -out $2.step1.rpstblastn.out -num_threads 15

rpsbproc -i $2.step1.rpstblastn.out -o $2.step1.rpsbproc.result -e 0.0001 -d $CDDATA -m rep -t doms

python $DARTS/intermediate_scripts/1step_parser.py $1 $2.step1.rpsbproc.result $2 

rm $2.step1.rpstblastn.out

#STEP 2 - mining all domains in one file, elements in another 4 files with its Scores
rpstblastn -query aRNH_and_approximates_$2_genome.fa -db $DARTS/customCDD/Step3 -outfmt 11 -out $2.step3.rpstblastn.out -num_threads 15

rpsbproc -i $2.step3.rpstblastn.out -o $2.step3.rpsbproc.result -e 0.0001 -d $CDDATA -m rep -t doms

python $DARTS/intermediate_scripts/mining.py aRNH_and_approximates_$2_genome.fa $2.step3.rpsbproc.result $2 

#STEP - clustering and choosing delegates
python $DARTS/intermediate_scripts/parse_prot_seq_domains.py prot_domains.fa $2_elements 

python $DARTS/intermediate_scripts/replacing_star_to_X.py $2 'gRT' 
python $DARTS/intermediate_scripts/replacing_star_to_X.py $2 'gRH' 
python $DARTS/intermediate_scripts/replacing_star_to_X.py $2 'aRH' 
python $DARTS/intermediate_scripts/replacing_star_to_X.py $2 'GAG' 
python $DARTS/intermediate_scripts/replacing_star_to_X.py $2 'INT' 
python $DARTS/intermediate_scripts/replacing_star_to_X.py $2 'PRo' 


mmseqs easy-cluster gRT_$2_elements.faa gRT_$2_elements.mmseqs80-80 tmp --min-seq-id 0.8 -c 0.8 #--cov-mode 1
python $DARTS/intermediate_scripts/mod_mmseq_outfile_parser.py gRT_$2_elements.mmseqs80-80_cluster.tsv gRT_$2_elements.faa
python $DARTS/intermediate_scripts/del_cluster_amounts_from_name.py gRT_$2_elements.faa mmseq_delegates_gRT_$2_elements.faa
python $DARTS/intermediate_scripts/find_another_elements_delegate.py $2 rewrited_mmseq_delegates_gRT_$2_elements.faa ### OR gRT_$2_elements_after_cdhit80_new_delegates.faa
for i in {'gRH','aRH','GAG','INT','PRo'}; do python $DARTS/intermediate_scripts/turn_back_cluster_amounts_to_name.py mmseq_delegates_gRT_$2_elements.faa $i-domain_$2_elements_delegates.faa; done
mkdir mmseq2-80
cp mmseq_delegates_gRT_$2_elements.faa ./mmseq2-80/gRT-$2_delegates.faa
for i in {'gRH','aRH','GAG','INT','PRo'}; do cp rewrited_$i-domain_$2_elements_delegates.faa ./mmseq2-80/$i-$2_delegates.faa; done
python $DARTS/intermediate_scripts/adapt_coordinates.py $2_elements Coords_1step qwerty_table


