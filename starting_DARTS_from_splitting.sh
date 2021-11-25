#!/bin/bash

# $1 - start domain; $2 - genome; $3 - project name; $4 GYPSY ?

Test=$(python $DARTS/is_splitting_needed.py $2)
if [ "$Test" = "1" ]; then
	python $DARTS/genome_splitter_choose_file.py $2 $3
	gzip $2
	for j in $3*.fa
	do
		Test2=$(python $DARTS/is_splitting_needed.py $j)
		if [ "$Test2" = "2" ]; then
			python $DARTS/alternative_sequence_splitter_choose_file.py $j $j
			rm $j
		fi
	done
elif [ "$Test" = "2" ]
then
	python $DARTS/alternative_sequence_splitter_choose_file.py $2 $3
	gzip $2
else
	mv $2 $3.fa
fi

for i in $3*.fa
do
	mkdir $i-folder
	cd $i-folder
	ln -s ../$i
	if [ "$1" = "aRH" ]; then
		sh $DARTS/rpstblastn-rpsbproc-mining-domain_parsing_change_cdhit_or_mmseq_mode.sh $i $i lol $4 lol #last two have no meaning
	else
		sh $DARTS/rpstblastn-rpsbproc-mining-domain_parsing_no_aRH_change_cdhit_or_mmseq_mode.sh $i $i lol $4 lol
	fi
	cd ..
done


if [ "$Test" != "3" ]
then
	mkdir $2_all
	for i in $3*-folder
	do
		cd $i-folder/
		cp $i_elements ../$2_all
		for domains in {'gRH','aRH','GAG','INT','PRo','gRT'}; do cp $domains_$2_elements.faa ../$2_all; done
		cd ..
	done
	cd $2_all
	for domains in {'gRH','aRH','GAG','INT','PRo','gRT'}; do cat $domains_* $domains_$2_all_elements.faa; done
	cat *_elements $2_all_elements
	rm *.fa_*
	sh $DARTS/rpstblastn-rpsbproc-mining-domain_parsing_ready-half_before_cdhit_or_mmseq-50.sh lol $2_all lol $4 lol
	cd ..

fi


