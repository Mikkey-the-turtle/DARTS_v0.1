DARTS v0.1 README file, November 2021
================================================================================

CONTENTS:
SUMMARY
DEPENDENCIES AND REQUIREMENTS
INSTALLATION OF THE PIPELINE PACKAGES
RUNNING THE SCRIPTS
SAMPLE OUTPUT

--------------------------------------------------------------------------------

SUMMARY

DARTS, Domain-Associated RetroTransposon Search, is a pipeline developed for searching, mining and annotation of LTR retrotransposons. DARTS pipeline is based on a standalone version of Reverse Position-Specific BLAST (RPS-BLAST), implemented online as CD-Search (Conserved Domain Search, https://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi). 

The pipeline starts with two rounds of RPS-BLAST domain search. First round is to find the central domain (here, it may be aRNH or RT of either Ty3/Gypsy LTR retrotransposons or including RT of Ty1/Copia, Bel/Pao, and Retroviruses). Then, all hits are extracted with their flanking regions (about +/-7500 bp) for the second RPS-BLAST run, which searches for 6 additional domains: reverse transcriptase, protease, gag polyprotein, integrase, ribonuclease H, and additional/archaeal ribonuclease H, annotated by DARTS as gRT, PRo, GAG, INT, gRH and aRH, respectively. 

DARTS is aimed to help researchers in annotation of the genome assemblies for repetitive elements, to produce clusters with representative elements for phylogenetic analysis, generate amino acid sequences of each studied domain for alignments, and extract full-length nucleotide element sequences with annotated structure.

--------------------------------------------------------------------------------

DEPENDENCIES AND REQUIREMENTS

The pipeline was tested and worked on a 64-bit Linux machine with installed Ubuntu v18.0 and v20.0. The PC had 19Gb RAM, although DARTS can analyse small genomes on PCs with 4 Gb RAM. Processor - Intel(R) core(TM) i5-3470 CPU @ 3.20GHz. 

! Warning: the pipeline outputs take a lot of space on disk (compared to genome assembly sizes); for its intermediate results DARTS needs at least 2-3 times more free space than the unpacked genome assembly takes.

Before starting working with DARTS, a few dependencies must be installed. First, it is necessary to upload and configure the rpstblastn and rpsbproc tools from the blast+ package. To do this, following the manual: https://ftp.ncbi.nih.gov/pub/mmdb/cdd/rpsbproc/README

Almost all DARTS scripts are based on Python programming language (developing version was python 3.6). The Biopython library is necessary too. Information about installation of the Python and Biopython packages is shown at https://www.python.org/downloads/ and https://biopython.org/wiki/Download respectively.

Next, for clustering, the MMseqs2 tool is required. The installation is described in the following manual: https://mmseqs.com/latest/userguide.pdf

For the following analysis there is a small script (not a part of DARTS), that launches alignments by MAFFT and subsequent phylogenetic tree building by FastTree and IQTree tools. Information on installation of these tools is provided at the following links: https://mafft.cbrc.jp/alignment/software/, http://www.microbesonline.org/fasttree/, http://www.iqtree.org/doc/iqtree-doc.pdf

--------------------------------------------------------------------------------

INSTALLATION OF THE PIPELINE PACKAGES

To use the DARTS pipeline, a user needs to perform the following steps:
1. Install the dependencies (see DEPENDENCIES AND REQUIREMENTS).
2. Download DARTS files in a separate folder.
3. Unzip customCDD.zip in the same folder. This container carries local databases which are required by DARTS.
4. Write the paths to rpsblast, rpsbproc and DARTS scripts to the system PATH. For example put them to ~/.bashrc or ~/.profile:
   export DARTS="/way/to/DARTS/scripts" 
   export CDDATA="/way/to/cddblast/ncbi-blast-N.N.N+-src/c++/ReleaseMT/bin/data" (N.N.N - the version of the cddblast)
The cddblast data folder should have the following files: bitscore_specific.txt, cddannot_generic.dat, cdtrack.txt, cddannot.dat, cddid.tbl, family_superfamily_links
Now DARTS is ready to use.

--------------------------------------------------------------------------------

RUNNING THE SCRIPTS

After adding the DARTS scripts to the path, DARTS can be run in the command line by executin:
  sh DARTS/main_DARTS.sh
Other DARTS scripts can launched similarly: 
  sh DARTS/%script_name.sh
or, for untermediate scripts:
  python DARTS/intermediate_scripts/%script_name.py 

Here are the main DARTS scripts and their usage:
1. main_DARTS.sh
The main DARTS script launches a full analysis pipeline as described in (Biryukov and Ustyantsev 2021; submitted to the Genes MDPI journal).  The script initiates the test of genome assembly size, and splits it into several batches if necessary, then DARTS launches the core part of the pipeline, and, if the splitting took place, joins the output from each batch, and creates a directory with all the results. The script invocation is:
  sh DARTS/main_DARTS.sh $1 $2 $3 $4
Where, $1 - name of the protein domain for the first step of RPS-BLAST search. For now, it may be either 'aRH' or 'RT'. RT will be run by default.
$2 - a genome assembly file. ! Warning ! The script must be launched from the same folder where the genome file is. Alternatively, you can make a symbolic link to the file and put it in a different folder: "ln -s way/to/genome_assemby.fa".
$3 - a %project_name. Most of the output file names will contain the %project_name.
$4 - a choice between all LTR retrotransposons or only Ty3/Gypsy members.

2. core_DARTS_aRH_search.sh
The DARTS pipeline without splitting, searches the aRNH domain at the first RPS-BLAST step. To run type:
  sh DARTS/core_DARTS_aRH_search.sh $1 $2
The arguments are: $1 - genome assembly file or file from its split batch. $2 - %project_name.

3. core_DARTS_RT_search.sh
The DARTS pipeline without splitting, searches the RT domain at the first RPS-BLAST step. To run type:
  sh DARTS/core_DARTS_RT_search.sh $1 $2 $3
The arguments are: $1 - genome assembly file or file from its split batch. $2 - %project_name. $3 - "GYPSY" mode (only search for Ty3/Gypsy elements) or "all" mode (+Ty1/Copia, Bel/Pao, retroviruses). The "all" mode is a default option.

4. clustering50_DARTS.sh
If user has ready "(%project_name).fa_elements" and "(%domain_type)_(%project_name).fa_elements.faa" files, the script will run all the DARTS steps from clustering with MMseqs2 at 50% identity level. To run type:
  sh DARTS/clustering50_DARTS.sh $1
The only variable the script needs is the %project_name. ! This variable must be the same as in the parts of input file names.
This script is used for final clustering of split batches in the folder "(@project_name)_all".

5. clustering80_DARTS.sh
The script is the same as the previous, but with another identity level of clustering - 80%. To run type:
  sh DARTS/clustering80_DARTS.sh $1
The only variable the script needs is the %project_name. ! This variable must be the same as in the parts of input file names. This script will be launched automatically during execution of the "core_DARTS_RT_search.sh" and "core_DARTS_aRH_search.sh" scripts.

6. mafft_2_steps.sh
This is a secondary script for alignment and phylogenetic tree building. It uses MAFFT for alignment of amino acid sequences of protein domains and uses FastTree and IQTree for phylogeny. To run type:
  sh DARTS/mafft_2_steps.sh $1 $2 $3 $4
The arguments are: $1 - a basal (profile) alignment to which a user wants to add the results. $2 - a *.fasta file with sequences to add. $3 - %project_name. $4 - "iqtree" or nothing - to use or not use IQTree.
The script will produce: $3-step1.fa (aligned $1 file) and  $3-step2.fa (alignment of sequences from $1 and $2 (mafft uses '--add' function)). $3-step3.fasttree.treefile - treefile made by FastTree tool. $3-step2.fa.X - X IQTree output files, including $3-step2.fa.treefile

--------------------------------------------------------------------------------

SAMPLE OUTPUT

Here is a list of file types a user may found in the DARTS outputs:
1. "(%genome_name).fa.step?.rpsbproc.result" - result tables of rpsbproc steps. Output files of rpstblastn are removed after the rpsbproc tables are produced.
2. "(%genome_name).fa_elements" (for example - Athal.fa_elements) - FASTA file with the nucleotide sequences of all found elements.
3. "(%genome_name).fa_%quality_elements" - previous file splitted into score gradations.
4. "(%domain_type)_(%genome_name).fa_elements.faa" (for example GAG_Athal.fa_elements.faa) - FASTA file with amino acid sequences of %domain_type found in all elements from the    previous file (1.). Altogether, there will be 6 domains - GAG (gag polyprotein), PRo (Protease), RT (reverse transcriptase, annotated as gRT), RNH (ribonuclease H, annotated    as gRH), INT (integrase), and aRNH (additional/archaeal RNH, annotated as aRH).
4. folder "mmseq2-80" contains the same files as in (2.), but the sequences are from representative elements from each lineage after clustering.
5. "Natural_table" - processed mmseq2 output table of clusters (all chimeric elements are excluded, the best representative elements are chosen).
6. "coordinates_table_of_(%genome_name).fa_elements" - BED-like formated file, where for each element found the following information is presented: 1) contig from assembly, 2-3)    start and end coordinates for 4) LTR1/LTR2/internal (domain-containing part), 5) the element structure annotation and 6) strand (+/-).
7. "Coords_1step" - intermediate file required to produce the previous file.
8. "prot_domains.fa" - all amino acid sequences of the extracted domains. From this file other following files are produced "(%domain_type)_(%genome_name).fa_elements.faa"
9. "aRNH/RT_and_approximates_Athal.fa_genome.fa" - nucleotide sequences of the identified hypothetical elements (first domain being found and its flanking regions). The results    of the first step of the RPS-BLAST search.

For most applications a user may need the following files: 2. "(%genome_name).fa_elements" file - to work with full-length element nucleotide sequences. 4. "(%domain_type)(%genome_name).fa_elements.faa" files - to work with all the domains of a given type that were found. 5. Folder "mmseq2-80" - for alignment and phylogeny. 6. "coordinates_table_of(%genome_name).fa_elements" - for genome annotation.
