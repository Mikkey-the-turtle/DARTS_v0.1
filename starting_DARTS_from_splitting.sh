#!/bin/bash

# $1 - start domain; $2 - genome; $3 - project name; $4 GYPSY ?

Test=$(python $DARTS/intermediate_scripts/is_splitting_needed.py $2)
if [ "$Test" = "1" ]; then
	python $DARTS/intermediate_scripts/genome_splitter.py $2 $3
	gzip $2
	for j in $3*.fa
	do
		Test2=$(python $DARTS/intermediate_scripts/is_splitting_needed.py $j)
		if [ "$Test2" = "2" ]; then
			python $DARTS/intermediate_scripts/sequence_splitter.py $j $j
			rm $j
		fi
	done
elif [ "$Test" = "2" ]
then
	python $DARTS/intermediate_scripts/sequence_splitter.py $2 $3
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
		sh $DARTS/core_DARTS_aRH_search.sh $i $i $4 
	else
		sh $DARTS/core_DARTS_RT_search.sh $i $i $4 lol #last has no meaning
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
	sh $DARTS/core_DARTS_clustering50.sh lol $2_all $4 lol #last has no meaning
	cd ..

fi


