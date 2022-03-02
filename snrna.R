library(Seurat)
library(biomaRt)

path <- "/path/to/cellranger-count/output/folder/outs/"
sample.type <- "human" | "mouse"

#Load ensembl rna.data
if (sample.type == "human"){
    mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
} else if (sample.type == "mouse"){
    mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")
}


#Read raw rna.dataset
h5.data <- Read10X_h5(filename = paste0(path, "/filtered_feature_bc_matrix.h5"))

#Create the Seurat object with the raw rna.data
rna.data <- CreateSeuratObject(counts = h5.data, min.cells = 3, min.features = 200)
#Calculate percentage of mitochondrial genes
if (sample.type == "human"){
    rna.data[["percent.mt"]] <- PercentageFeatureSet(rna.data, pattern = "^MT-")
} else if (sample.type == "mouse"){
    rna.data[["percent.mt"]] <- PercentageFeatureSet(rna.data, pattern = "^mt-")
}

percent.mt.20 <- length(which(rna.data$percent.mt >= 20))/length(rna.data$percent.mt)*100

#Calculate percentage of cell with XIST expression
exp.matrix <-  as.matrix(GetAssayData(rna.data, slot = "data"))

xist.index <- match("XIST", toupper(rownames(exp.matrix)))
if (!is.na(xist.index)){
    percent.xist <- length(which(exp.matrix [xist.index,] > 0))/ncol(exp.matrix )
}else{
    percent.xist <- 0
}

#Determine gender
if(percent.xist > 0.1){
    gender <- "Female"
}else{
    gender <- "Male"
}


#Return all result
cat(paste0(c(percent.mt.20, percent.xist, gender), collapse = ","))
