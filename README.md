# ig_motif

###This repository contains scripts for the computational analysis as described in the manuscript:

Bowen A, Wear MP, Cordero RJB, Oscarson S, Casadevall A. A monoclonal antibody to *Cryptococcus neoformans* glucuronoxylomannan manifests hydrolytic activity for both peptides and polysaccharides. *J Biol Chem*. 2016. In Press.

###Three scripts are included in this repository:

1. PDB analysis.R
This script can be used to download all PDB files from the Protein Data Bank or selected sets from lists in the PDBid_lists.xlsx file. This script also analyzes all PDB files for duplicates, and identifies any files that contain antibody V region sequences. This analysis requires the comparison of PDB file protein chains to sequences in the "IMGT Gene-DB AA V gaps" directory.

2. Cab_seq_tree1-0.R
This script is used to generate circular dendrogram plots of heavy or light chain sequences analyzed by scripts 1 and 3. Examples of the plots are shown in the "figs" directory and were used to create Figure 10 in the JBC manuscript. This script was also used to analyze the amino acid conservation of V region sequences and to generate Figure 9 in the JBC manuscript. In order to generate the sequence conservation figures, fasta files of the sequences are created in the output folder. These are then aligned using the online Clustal Omega tool (http://www.ebi.ac.uk/Tools/msa/clustalo/) with default parameters. Aligned sequences are saved as fasta files in the alignments directory. These aligned files are used by the script to generate the conservation figure.

3. PDB templates1-3.R
This script contains the main IgMotif algorithm where a three-residue seed motif can be specified and compared to possible motifs from other protein structures in the PDB. Script 1 must be run first to download and analyze the appropriate PDB files.

###Directories and sub-directories that need to be pointed to in the beginning of these scripts are listed below:

* IgMotif
 * IMGT Gene-DB AA V gaps
 * alignments
 * figs
 * output
 * mean temps
  * seeds
* R Working Directory
  * IgMotif_temp
    * PDB_DB
      * PDB_files

###Files needed to run these scripts are listed below and present in the repository:

* PDBid_lists.xlsx (in the parent directory)
* all ".fasta" files listed in the "IMGT Gene-DB AA V gaps" directory

###Contact Anthony Bowen with questions or difficulties running these scripts
* anthony.bowen at med.einstein.yu.edu
