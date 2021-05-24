# WasteWaterVariantCalling

## WasteWaterVariantCalling.wdl

Performs variant calling on SARS-CoV-2 in waster water samples and identifies mutations in the Spike gene associated with known VOCs and VUIs.
Takes as input an Array[File] of coordinate sorted bams from the following workflowsP illumina-preprocessing-assembly (for paired-end illumina data), nanopore-preprocessing-assembly (for nanopore data), and COVIDSeqSEassembly (for single-end sequencing using the COVIDSeq protocol).
