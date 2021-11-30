#!/bin/bash
#export CDDATA="/way/to/cddblast/ncbi-blast-N.N.N+-src/c++/ReleaseMT/bin/data"
# variables: $1 - genome assembly file or its split, $2 - project name, $3 - "GYPSY" or all RT search
source ~/.profile
source ~/.bashrc

#STEP 1 - finding RT and its neighborhood
if [ "$3" = "GYPSY" ]; then
	rpstblastn -query $1 -db $DARTS/customCDD/Step2GYPSY -outfmt 11 -out $2.step2g.rpstblastn.out -num_threads 4
	rpsbproc -i $2.step2g.rpstblastn.out -o $2.step2g.rpsbproc.result -e 0.0001 -d $CDDATA -m rep -t doms
	python3 $DARTS/intermediate_scripts/1step_parser_RT.py $1 $2.step2g.rpsbproc.result $2
else
	rpstblastn -query $1 -db $DARTS/customCDD/Step2 -outfmt 11 -out $2.step2.rpstblastn.out -num_threads 4
	rpsbproc -i $2.step2.rpstblastn.out -o $2.step2.rpsbproc.result -e 0.0001 -d $CDDATA -m rep -t doms
	python3 $DARTS/intermediate_scripts/1step_parser_RT.py $1 $2.step2.rpsbproc.result $2
fi

rm $2.step2g.rpstblastn.out
rm $2.step2.rpstblastn.out

#STEP 2 - mining all domains in one file, elements in another 4 files with its Scores
if [ "$3" = "GYPSY" ]; then
	rpstblastn -query RT_and_approximates_$2_genome.fa -db $DARTS/customCDD/Step3 -outfmt 11 -out $2.step3.rpstblastn.out -num_threads 4
	rpsbproc -i $2.step3.rpstblastn.out -o $2.step3.rpsbproc.result -e 0.0001 -d $CDDATA -m rep -t doms
	python3 $DARTS/intermediate_scripts/mining_RT.py RT_and_approximates_$2_genome.fa $2.step3.rpsbproc.result $2 $3
else
	rpstblastn -query RT_and_approximates_$2_genome.fa -db $DARTS/customCDD/StepX4 -outfmt 11 -out $2.stepX4.rpstblastn.out -num_threads 4
	rpsbproc -i $2.stepX4.rpstblastn.out -o $2.stepX4.rpsbproc.result -e 0.0001 -d $CDDATA -m rep -t doms
	python3 $DARTS/intermediate_scripts/mining_RT.py RT_and_approximates_$2_genome.fa $2.stepX4.rpsbproc.result $2
fi

rm $2.step3.rpstblastn.out
rm $2.stepX4.rpstblastn.out

#STEP - clustering and choosing delegates
python3 $DARTS/intermediate_scripts/parse_prot_seq_domains_RT.py prot_domains.fa $2_elements 

python3 $DARTS/intermediate_scripts/replacing_star_to_X.py $2 'gRT' 
python3 $DARTS/intermediate_scripts/replacing_star_to_X.py $2 'gRH' 
python3 $DARTS/intermediate_scripts/replacing_star_to_X.py $2 'aRH' 
python3 $DARTS/intermediate_scripts/replacing_star_to_X.py $2 'GAG' 
python3 $DARTS/intermediate_scripts/replacing_star_to_X.py $2 'INT' 
python3 $DARTS/intermediate_scripts/replacing_star_to_X.py $2 'PRo' 


mmseqs easy-cluster gRT_$2_elements.faa gRT_$2_elements.mmseqs80-80 tmp --min-seq-id 0.8 -c 0.8 #--cov-mode 1
python3 $DARTS/intermediate_scripts/mod_mmseq_outfile_parser.py gRT_$2_elements.mmseqs80-80_cluster.tsv gRT_$2_elements.faa
python3 $DARTS/intermediate_scripts/del_cluster_amounts_from_name.py gRT_$2_elements.faa mmseq_delegates_gRT_$2_elements.faa
python3 $DARTS/intermediate_scripts/find_another_elements_delegate.py $2 rewrited_mmseq_delegates_gRT_$2_elements.faa ### OR gRT_$2_elements_after_cdhit80_new_delegates.faa
for i in {'gRH','INT','GAG','PRo'}; do python $DARTS/intermediate_scripts/turn_back_cluster_amounts_to_name.py mmseq_delegates_gRT_$2_elements.faa $i-domain_$2_elements_delegates.faa; done
mkdir mmseq2-80
cp mmseq_delegates_gRT_$2_elements.faa ./mmseq2-80/gRT-$2_delegates.faa
for i in {'gRH','INT','GAG','PRo'}; do cp rewrited_$i-domain_$2_elements_delegates.faa ./mmseq2-80/$i-$2_delegates.faa; done
python $DARTS/intermediate_scripts/adapt_coordinates.py $2_elements Coords_1step qwerty_table

