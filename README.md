# DARTS v0.1      README file, November 2021
================================================================================

CONTENTS: 

1. SUMMARY
2. DEPENDENCIES AND REQUIREMENTS
3. INSTALLATION OF THE PIPELINE PACKAGE
5. RUN THE PROGRAMS
6. SAMPLE OUTPUT

--------------------------------------------------------------------------------

SUMMARY

DARTS, Domain-Assosiated RetroTransposon Search, is a pipeline developed for 
searching, mining and annotation of LTR retrotransposons.
DARTS pipeline is based om standalone version of Reverse Position-Specific BLAST 
(RPS-BLAST), also known as CD-Search (Conserved Domain Search).
The pipeline starts with two rounds of RPS-BLAST domain search. First round is 
about finding the central domain (here, it may be aRNH or RT of either Ty3/Gypsy
LTR retrotransposons or including Ty1/Copia, Bep/Pao, Retroviruses). Then, all
hits are extracted with their flanking regions (about +-7500 bp) for the second
RPS-BLAST run, which is searching for 6 domains: reverse transcriptase, protease, 
gag polyprotein, integrase, ribonuclease H and additional/achaeal ribonuclease H,
annotated by DARTS as gRT, PRo, GAG, INT, gRH and aRH subsequently.
Its outputs help to annotate genome assemblies for repetative elements, produce 
clusters with representative elements for phylogenetic analysis and generate
amoni acid sequences of each studied domain for alignments, and full-length 
nucleotide element sequences.
For each element the pipeline uses its own score of wholeness. Choice of representative
element of a cluster is based on that score.

--------------------------------------------------------------------------------

DEPENDENCIES AND REQUIREMENTS

Before starting work with DARTS, a few dependencies must be installed.
First, it is necessary to upload and configure rpstblastn and rpsbproc commands 
from blast+ package. The way how to do that is described in the following manual:
https://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/README

Almost all DARTS scripts are based on Python programming language (developing 
version was python3.6). The Biopython package is necessary too. Information about
installation of Python and Biopython package is shown at 
https://www.python.org/downloads/ and https://biopython.org/wiki/Download subsequently.

Next, for clusterisation step, the MMseqs2 tool is needed. The way to installation 
described in the following manual:
https://mmseqs.com/latest/userguide.pdf

--------------------------------------------------------------------------------

INSTALLATION OF THE PIPELINE PACKAGE

To use the DARTS pipeline, user need to provoke the following tips:

1. Installation of dependencies (see DEPENDENCIES AND REQUIREMENTS).
2. Download DARTS files in a separate folder.
3. Unzip customCDD.zip in the same folder. This container carries local databases which DARTS is needed.
4. Write in /.bashrc or /.profile ways to both cddblast data folder and to DARTS scripts two following strings:

export DARTS="/way/to/DARTS/scripts" 
export CDDATA="/way/to/cddblast/ncbi-blast-N.N.N+-src/c++/ReleaseMT/bin/data"

The cddblast data folder is normally presented by the following files :
    bitscore_specific.txt  cddannot_generic.dat  cdtrack.txt
    cddannot.dat           cddid.tbl             family_superfamily_links

Now DARTS is ready to use.

--------------------------------------------------------------------------------

RUN THE PROGRAMS

After adding variable DARTS(export DARTS="/way/to/DARTS/scripts") showing the path
to scripts folder, DARTS can be run in the command line by:
sh DARTS/starting_DARTS_from_splitting.sh
The same situation will be for all scripts in the same folder - they can be started 
using "sh/python DARTS/%script_name.sh/py".

Here are the main DARTS scripts to deal with and how to run them.

1. starting_DARTS_from_splitting.sh

The main DARTS script produces all the pipelinedescribed in the issue. It provokes 
the test of genome assembly splitting being necessary,splits, if needed, then starts
core part of the pipeline, and, if the splitting took place, reunite output of 
each batch to produce folder of all resuts and new tip of clustering.
The script invocation is:

sh DARTS/starting_DARTS_from_splitting.sh $1 $2 $3 $4

First variable the current script needs is a domain of the first step of RPS-BLAST.
It may be either 'aRH' or 'RT'. If there will be something different, script will 
interpret it as 'RT'.

Second variable is a genome assembly file. 
! Warning! All script are sensitive with the ways. Try to run them  in the current folder, 
where genome assembly is put. Or make a "ln -s way/to/genome_assemby.fa" in the directory
you want DARTS to be run.

Third variable is a %project_name. Most of output files will be called partly with that name.
Fourth varaible needs only if user want to search only LTR retrotransposons of Ty3/Gypsy.

THIS PART WILL BE UPDATED.


--------------------------------------------------------------------------------

SAMPLE OUTPUT

Here is a list of file types user may found in DARTS outputs:
1. 
2. "(%genome_name).fa_elements" (example - Athal.fa_elements) - FASTA file of nucleotide sequences of found elements.
3. "(%domain_type)_(%genome_name).fa_elements.faa" (GAG_Athal.fa_elements.faa) - FASTA file of amino acid sequences of %domain_type found in all elements from previous file (1.). Altogether, there are 6 domains - GAG (gag polyprotein), PRo (Protease), RT (reverse transcriptase, annotated as gRT), RNH (ribonuclease H, annotated as gRH), INT (integrase) and aRNH (additional/archaeal RNH, annotated as aRH). 
4. folder "mmseq2-80" with the same files as it was in 2., but only from representative elements from each cluster after clusterisation.
5. "Natural_table" - processed mmseq2 output table of clusters (all chimeric elements excluded, better representative elements chosen).
6. "coordinates_table_of_(%genome_name).fa_elements" - BAS-like file, where for each element found the following information is presented: 1) contig from assembly, 2-3) start and end coordinates for 4) LTR1/LTR2/internal (domain-containing part), 5)element annotation and 6) stand (+/-).


