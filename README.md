# snRNA_snATAC_snMultiome_QC

These R scripts can be used for calculating quanlity control metrics from the assay specific CellRanger count outputs.

For snRNAseq, use snrna.R to calculate the fraction of cells with mitochondrial reads >= 20%, the fraction of cells with XIST expression and the gender determination based of XIST expression. 

For snATACseq, use snatac.R to determine the present of the nucleosome free region peak and mononucleosome peak. 

For snMultiome (ATAC + Gene Expression), use snmultiome.R which combined both calculations from snrna.R and snatac.R.
