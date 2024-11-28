This repository is made up of scripts written in bash, python and R.

#BASH SCRIPTS
The bash scripts involve the manipulation of sequence data for the purposes of bioinformatics analysis. This analysis ranges from quality control, assembly, taxonomy identification, genome annotation and phylogenetic analysis. 
The bash script has prerequisite tools that need to be installed in order to be of use. These	are;
1. Trimmomatic
2. BBtools
3. Qiime2
3. Dorado (base calling nanopore reads into fast)
4. fastqc (quality checks)
5. Parsnp
6. Roary
7. Mafft
8. IQtree
9. Kraken2

#Python Scripts
The Python scrips are designed for file manipulation including but not limited to sequence files. These operations include, splitting genomes from downloaded fasta files, combining kraken2 reports and generating readmapping visualizations. Python needs to be installed either through Anaconda or Miniconda for these to work.

#R Scripts
These scripts are designed for metagenomic diversity analysis through Alpha diversity (Shannon and Simpson) and Beta diversity (Bray-curtis dissimilarity into PERMANOVA). The scripts are designed to work with the OTU table output gotten from Qiime2 and Kraken2.
  
 
