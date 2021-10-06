# WasteWaterVariantCalling

## WasteWaterVariantCalling.wdl

Performs variant calling on SARS-CoV-2 in waster water samples and identifies mutations in the Spike gene associated with known VOCs and VUIs on the Google Cloud Terra Platform.
Takes as input an Array[File] of coordinate sorted bams from the following workflows: illumina-preprocessing-assembly (for paired-end illumina data), nanopore-preprocessing-assembly (for nanopore data), and COVIDSeqSEassembly (for single-end sequencing using the COVIDSeq protocol).

inputs needed: Array[File] of coordinate sorted bams to be analyzed, Array[File] of sample ids, fasta reference genome for SARS-CoV-2 saved as workspace data, a bed file containing 1-based reference positions for mutations of interest (Spike VOC and VUI mutations) saved as workpace data, a 3-column tsv with annotations for them muations of interest (AA_change, Nucl_change, Lineage), the google bucket path to where outputs should be transferred (written as a String in double quotes)

1. performs variant calling and variant filtration using Freebayes
2. uses bcftools view to pull out VOC/VUI-associated spike mutations
3. generates a summary tsv of the allele frequencies of spike mutations in the bam files for all samples in both wide and long formats, and a long format summary of the depth of coverage for spike mutations
4. transfers outputs to the user's chosen google bucket location

citation for bcftools:
Twelve years of SAMtools and BCFtools
Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li
GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008

citation for freebayes:
Garrison E, Marth G. Haplotype-based variant detection from short-read sequencing. arXiv preprint arXiv:1207.3907 [q-bio.GN] 2012

Docker images used:
bcftools: quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0
samtools: staphb/samtools:1.10
Theiagen general utilities: theiagen/utility:1.0
Datamash: rapatsky/debian
Freebayes: wgspipeline/freebayes:v0.0.1
