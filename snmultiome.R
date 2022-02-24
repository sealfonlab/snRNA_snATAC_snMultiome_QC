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
library(biomaRt)

path <- "/path/to/cellranger-count/output/folder/outs/"
sample.type <- "human" | "mouse"

data <- Read10X_h5(paste0(path,"/filtered_feature_bc_matrix.h5"))

fragment.path <- paste0(path, "/atac_fragments.tsv.gz")

metadata <- read.csv(
  file = paste0(path, "/per_barcode_metrics.csv"),
  header = TRUE,
  row.names = 1
)

#Read GEX side of the data
multiome <- CreateSeuratObject(
  counts = data$`Gene Expression`,
  assay = "RNA"
)

multiome <- AddMetaData(
  object = multiome,
  metadata = metadata
)

#Load ensembl data
if (sample.type == "human"){
    mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
} else if (sample.type == "mouse"){
    mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
}


#Calculate percentage of mitochondrial genes
if (sample.type == "human"){
    multiome[["percent.mt"]] <- PercentageFeatureSet(multiome, pattern = "^MT-")
} else if (sample.type == "mouse"){
    multiome[["percent.mt"]] <- PercentageFeatureSet(multiome, pattern = "^mt-")
}

percent.mt.20 <- length(which(multiome$percent.mt >= 20))/length(multiome$percent.mt)*100

#Calculate percentage of cell with XIST expression
exp.matrix <-  as.matrix(GetAssayData(multiome, slot = "data"))

xist.index <- match("XIST", toupper(rownames(exp.matrix)))
if (!is.na(xist.index)){
    percent.xist <- length(which(exp.matrix [xist.index,] > 0))/ncol(exp.matrix )
}else{
    percent.xist <- 0
}

#Determine gender
if(percent.xist > 0.1){
    rna.gender <- "Female"
}else{
    rna.gender <- "Male"
}

# Prepared to read ATAC side of the data
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


# create ATAC assay and add it to the object
multiome[["ATAC"]] <- CreateChromatinAssay(
  counts = data$Peaks,
  sep = c(":", "-"),
  fragments = fragment.path,
  annotation = annotations
)

DefaultAssay(multiome) <- "ATAC"

#Prepare fragment length density plot
fragments.use <- FragmentHistogram(object = multiome, region = 'chr1-1-200000000')
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
cat(paste0(c(percent.mt.20, percent.xist, rna.gender, mono.nuc, nfr), collapse = ","))
