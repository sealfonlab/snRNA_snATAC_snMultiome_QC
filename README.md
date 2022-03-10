# snRNA_snATAC_snMultiome_QC

These R scripts can be used for calculating quality control metrics from the assay specific CellRanger count outputs.

For snRNAseq, use snrna.R to calculate the fraction of cells with mitochondrial reads >= 20%, the fraction of cells with XIST expression and the sex determination based on XIST expression. 

For snATACseq, use snatac.R to determine the presence of the nucleosome free region and mononucleosome peak. 

For snMultiome (ATAC + Gene Expression), use snmultiome.R which combines both calculations from snrna.R and snatac.R.
