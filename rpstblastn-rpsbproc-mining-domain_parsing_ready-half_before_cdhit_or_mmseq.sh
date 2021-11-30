#!/bin/bash
#export CDDATA="/way/to/cddblast/ncbi-blast-N.N.N+-src/c++/ReleaseMT/bin/data"
source ~/.profile
source ~/.bashrc

mmseqs easy-cluster gRT_$2_elements.faa gRT_$2_elements.mmseqs80-80 tmp --min-seq-id 0.8 -c 0.8 #--cov-mode 1
python $DARTS/intermediate_scripts/mod_mmseq_outfile_farser.py gRT_$2_elements.mmseqs80-80_cluster.tsv gRT_$2_elements.faa
python $DARTS/intermediate_scripts/del_cluster_amounts_from_name.py gRT_$2_elements.faa mmseq_delegates_gRT_$2_elements.faa
python $DARTS/intermediate_scripts/find_another_elements_delegate.py $2 rewrited_mmseq_delegates_gRT_$2_elements.faa ### OR gRT_$2_elements_after_cdhit80_new_delegates.faa
for i in {'gRH','aRH','GAG','INT','PRo'}; do python $DARTS/intermediate_scripts/turn_back_cluster_amounts_to_name.py mmseq_delegates_gRT_$2_elements.faa $i-domain_$2_elements_delegates.faa; done
mkdir mmseq2-80
cp mmseq_delegates_gRT_$2_elements.faa ./mmseq2-80/gRT-$2_delegates.faa
for i in {'gRH','aRH','GAG','INT','PRo'}; do cp rewrited_$i-domain_$2_elements_delegates.faa ./mmseq2-80/$i-$2_delegates.faa; done


