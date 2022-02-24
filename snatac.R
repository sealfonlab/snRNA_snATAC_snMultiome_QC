#!/usr/local/R/4.0.3/bin/R

library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(ChIPseeker)
library(ChIPpeakAnno)
library(Seurat)
library(Signac)
library(GenomicRanges)
library(GenomeInfoDb)
library(ggplot2)
library(stringr)


path <- "/path/to/cellranger-count/output/folder/outs/"
sample.type <- "human" | "mouse"

#Read raw data
data <- Read10X_h5(filename = paste0(path,"/filtered_peak_bc_matrix.h5"))
metadata <- read.csv(file = paste0(path,"/singlecell.csv"), header = TRUE, row.names = 1)
fragment.path <- paste0(path,"/fragments.tsv.gz")

# extract gene annotations from EnsDb
if (sample.type == "human"){
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotations),sep=""), pattern="chrMT", replacement="chrM")
    seqlevels(annotations) <- ucsc.levels
    genome(annotations) <- "hg38"
} else if (sample.type == "mouse"){
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
    ucsc.levels <- str_replace(string=paste("chr",seqlevels(annotations),sep=""), pattern="chrMT", replacement="chrM")
    seqlevels(annotations) <- ucsc.levels
    genome(annotations) <- "mm10"
}

#Create Seurat Object
chrom_assay <- CreateChromatinAssay(counts = data, sep = c(":", "-"), fragments = fragment.path, annotation = annotations)
atac.data <- CreateSeuratObject(counts = chrom_assay, assay = "peaks", min.cells = 1, meta.data = metadata)

#Prepare fragment length density plot
fragments.use <- FragmentHistogram(object = atac.data, region = 'chr1-1-200000000')
fragments.use <- ggplot_build(fragments.use)
fragments.use <- fragments.use$plot$data

NFR.UPPER.LIMIT <- 150
MONO.NUC.LOWER.LIMIT <- 150
MONO.NUC.UPPER.LIMIT <- 300

#Narrowing the fragment length to focus on
df.fragments.use <- as.data.frame(x = fragments.use[which(fragments.use[,"length"] < 600), c('cell', 'length')])
density.fragments.use <- density(df.fragments.use[,"length"], adjust = 2)

#Find peaks
peaks.index <- which(diff(sign(diff(density.fragments.use$y))) == -2)
peaks.rank <- rank(-density.fragments.use$y[peaks.index+1])

#Indentify peaks
mono.nuc <- FALSE
nfr <- FALSE
for (i in 1:length(peaks.index)){
  density.index <- peaks.index[i]+1
  if (density.fragments.use$x[density.index] <= NFR.UPPER.LIMIT && peaks.rank[i] <= 2){
    nfr <- TRUE
  } else if (density.fragments.use$x[density.index] >= MONO.NUC.LOWER.LIMIT && density.fragments.use$x[density.index] <= MONO.NUC.UPPER.LIMIT && peaks.rank[i] <= 2){
    mono.nuc <- TRUE
  }
}


#Return all result
cat(paste0(c(mono.nuc, nfr), collapse = ","))
