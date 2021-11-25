#!/bin/bash
if [ "$1" = "" ]; then
	echo "\$1 - basic domain alignment, \$2 - current domain representative\`s file from studied domain, \$3 -mproject name, \$4 iqtree - start iqtree after fasttree if necessary"
	echo "(Needs Mafft, FastTree & IQTree tools)"
else
	mafft $1 > $3-step1.fa
	mafft --add $2 --reorder $3-step1.fa > $3-step2.fa
	fasttree -lg $3-step2.fa  > $3-step3.fasttree.treefile
fi

Tree="iqtree"
if [ "$4" = "$Tree" ]; then
	iqtree -s $3-step2.fa -nt 4 -bb 1000 -alrt 1000 -m TEST
fi

