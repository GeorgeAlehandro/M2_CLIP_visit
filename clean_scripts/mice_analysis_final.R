library(rVictis)
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(SeuratObject)
library(ggplot2)
library(gplots)
library(patchwork)
library(SingleR)
library(Rmagic)
library(SingleCellExperiment)
library(SummarizedExperiment)
# set directory to where the data lies
# setwd('/home/rstudio/Ahmad_workdir/')
# DATAPATH = 'cellranger/Mice/'

#' \code{preproc.cellranger}
#' @param folder character. Folder where the 10XGenomics data to read lies.
#' @details Reads the gene data and creates a seurat object, then applies some pre-processing on the cells, genes, and row/column names
#' @return Returns the pre-processed seurat object.
#' @export

preproc.cellranger <- function(folder){
  data<-Read10X(data.dir = paste(DATAPATH, folder, sep=''))
  seurat.obj<-CreateSeuratObject(counts=data,min.features = 200, min.cells = 3,project=folder)
  # Filters out cells that express human genes, only keeps cells that express murine genes
  rna <- seurat.obj@assays$RNA
  tmp.genes<-gsub("-.*","",rownames(rna))
  hum <- names(table(tmp.genes))[1]
  mur <- names(table(tmp.genes))[2]
  hcts <- colSums(rna[which(tmp.genes==hum),])
  mcts <- colSums(rna[which(tmp.genes==mur),])
  ratio <- mcts/hcts
  plot(sort(ratio),yaxp  = c(0, 2500, 50))
  grid()
  cutoff <- 50 # can change this cutoff based on what the above plot shows
  rna <- rna[,which(ratio>cutoff)]
  rna@Dimnames[[1]] <- gsub(".*-","",rna@Dimnames[[1]]) # remove gene prefix
  rna@Dimnames[[2]] <- sapply(rna@Dimnames[[2]], function(x){substr(x, 1, nchar(x)-2)}, USE.NAMES = F) # remove -1 at end of cell name

  seurat.obj <- CreateSeuratObject(counts=rna,min.features = 200, min.cells = 3,project=folder)
  seurat.obj[["percent.mt"]]<-PercentageFeatureSet(seurat.obj,pattern="^MT-")
  seurat.filtered<-subset(seurat.obj,subset=nFeature_RNA>200 & nCount_RNA < 7000 & percent.mt < 8 )
  
  sct <- SCTransform(seurat.filtered, verbose = FALSE)
  return(sct)
}

#' \code{preproc.cellranger.with.adt}
#' @param folder character. Folder where the 10XGenomics data to read lies.
#' @param sct logical (default TRUE). Whether or not to apply SCT transformation.
#' @details Reads the gene data and creates a seurat object, then applies some pre-processing on the cells, genes, and row/column names. Includes pre-processing of ADT data.
#' @return Returns the pre-processed seurat object.
#' @export

preproc.cellranger.with.adt <- function(folder, sct=T){
  data<-Read10X(data.dir = paste(DATAPATH, folder, sep=''))
  seurat.obj<-CreateSeuratObject(counts=data$`Gene Expression`,min.features = 200, min.cells = 3,project=folder)
  seurat.adt<-CreateAssayObject(counts=data$`Antibody Capture`[, colnames(seurat.obj)])

  seurat.obj[["ADT"]]<-seurat.adt

  DefaultAssay(seurat.obj) <- "ADT"
  seurat.obj <- NormalizeData(seurat.obj, normalization.method = "CLR",margin=2)
  
  DefaultAssay(seurat.obj) <- "RNA"
  
  # Filters out cells that express human genes, only keeps cells that express murine genes
  rna <- seurat.obj@assays$RNA
  tmp.genes<-gsub("-.*","",rownames(rna))
  hum <- names(table(tmp.genes))[1]
  mur <- names(table(tmp.genes))[2]
  hcts <- colSums(rna[which(tmp.genes==hum),])
  mcts <- colSums(rna[which(tmp.genes==mur),])
  ratio <- mcts/hcts
  plot(sort(ratio),yaxp  = c(0, 2500, 50))
  grid()
  cutoff <- 100 # can change this cutoff based on what the above plot shows
  rna@counts <- rna@counts[,which(ratio>cutoff)] # remove cells expressing human genes
  rna@data <- rna@data[,which(ratio>cutoff)]
  rna@counts <- rna@counts[-which(tmp.genes==hum),] # filter human genes out
  rna@data <- rna@data[-which(tmp.genes==hum),]
  rna@data@Dimnames[[1]] <- gsub("mm10---","",rna@data@Dimnames[[1]]) # remove gene prefix
  rna@counts@Dimnames[[1]] <- gsub("mm10---","",rna@counts@Dimnames[[1]])
  # rna@Dimnames[[2]] <- sapply(rna@Dimnames[[2]], function(x){substr(x, 1, nchar(x)-2)}, USE.NAMES = F) # remove -1 at end of cell name
  
  seurat.obj@assays$RNA <- rna
  seurat.obj[["percent.mt"]]<-PercentageFeatureSet(seurat.obj,pattern="^MT-")
  seurat.filtered<-subset(seurat.obj,subset=nFeature_RNA>500 & nCount_RNA < 32000 & percent.mt < 10)
  if(sct) seurat.filtered <- SCTransform(seurat.filtered, verbose = FALSE)
  return(seurat.filtered)
}

#' \code{vae.seurat}
#' @param seurat.obj Seurat object. Output of one of the previous functions.
#' @param scale logical (default FALSE). Whether or not to scale data before applying dimensionality reduction.
#' @param center logical (default FALSE). Whether or not to center data before applying dimensionality reduction.
#' @details Applies vaevictis dimensionality reduction to the seurat data.
#' @return Returns the seurat object enriched with the embeddings.
#' @export

vae.seurat <- function(seurat.obj, scale = F, center = F){
  dt <- unname(t(as.matrix(seurat.obj@assays[[DefaultAssay(seurat.obj)]]@data)))
  dt <- scale(dt, scale = scale, center = center)
  red <- vaevictis(dt,ww=c(10., 0., 10., 5., 2., 1.),epochs=50L,patience=8.)
  rownames(red$layout) <- seurat.obj@assays[[DefaultAssay(seurat.obj)]]@data@Dimnames[[2]]
  seurat.obj[["vae"]] <- CreateDimReducObject(embeddings = red$layout, key = "vae_", assay = DefaultAssay(seurat.obj))
  return(list(srt = seurat.obj, model = red))
}

#' \code{annot.singleR}
#' @param seurat.obj Seurat object. Output of one of the previous functions.
#' @param ref.data SummarizedExperiment object. Reference data for automatic annotation of cell populations.
#' @details Uses SingleR to automatically annotate the cells of a seurat object, based on a provided reference dataset.
#' @return Returns a named list with the seurat object enriched with the annotations, and the annotations themselves.
#' @export
annot.singleR <- function(seurat.obj, ref.data){
  sce.obj <- as.SingleCellExperiment(seurat.obj, assay = DefaultAssay(seurat.obj))
  
  # using log transformation to match with the reference data
  counts <- assay(sce.obj, "counts")
  libsizes <- colSums(counts)
  size.factors <- libsizes/mean(libsizes)
  logcounts(sce.obj) <- log2(t(t(counts)/size.factors) + 1)
  
  pred.hesc <- SingleR(test = assay(sce.obj, 'counts'), ref=ref.data,
                       labels = ref.data$label.fine)
  
  seurat.obj@meta.data["singleR.annot"] <- pred.hesc$pruned.labels
  return(list(data=seurat.obj, annot=pred.hesc))
}

#' \code{export.indices}
#' @param seurat.obj Seurat object. Output of one of the previous functions.
#' @param fileprefix character. Path prefix to specify where to save the indices.
#' @details Exports the cell and gene indices of a Seurat object to .csv, to be used with another tool.
#' @export
export.indices <- function(seurat.obj, fileprefix){
  gene.indices <- seurat.obj@assays[[DefaultAssay(seurat.obj)]]@counts@Dimnames[[1]]
  cell.indices <- seurat.obj@assays[[DefaultAssay(seurat.obj)]]@counts@Dimnames[[2]]
  
  # cell.indices <- sapply(cell.indices, function(x){substr(x, 1, nchar(x)-2)}, USE.NAMES = F)
  
  write.table(gene.indices, 
              file=paste(fileprefix, 'index_genes.csv', sep=''), 
              sep=',', row.names=F, col.names=F)
  write.table(cell.indices, 
              file=paste(fileprefix, 'index_cells.csv', sep=''), 
              sep=',', row.names=F, col.names=F)
}

#' \code{remove.celltypes}
#' @param seurat.obj Seurat object. Output of one of the previous functions.
#' @param types character vector. Cell populations to filter out.
#' @param meta character. Key to find the cell annotations in the Seurat object.
#' @details Removes specified cell populations from a Seurat object. Needs existing cell annotations.
#' @return Returns the filtered seurat object.
#' @export
remove.celltypes <- function(seurat.obj, types, meta = "singleR.annot"){
  labels <- seurat.obj[[meta]]
  to.rm <- grepl(paste(types,collapse="|"), x=seurat.obj$singleR.annot)
  to.rm <- !to.rm
  new <- subset(seurat.obj, cells = which(to.rm))
  return(new)
}

#' \code{filter.cells.adt}
#' @param seurat.obj Seurat object. Output of one of the previous functions.
#' @param marker character. ADT marker to use to filter cells.
#' @param threshold numeric. Cells expressing in the specified ADT marker above this threshold will be filtered out.
#' @details Filters cells in a Seurat object based on their expression level in a specified ADT marker.
#' @return Returns the filtered seurat object.
#' @export
filter.cells.adt <- function(seurat.obj, marker, threshold = 0.){
  levels <- seurat.obj@assays$ADT[marker,]
  brcd <- seurat.obj@assays[[DefaultAssay(seurat.obj)]]@data@Dimnames[[2]]
  seurat.obj <- subset(seurat.obj, cells = brcd[which(levels<=threshold)])
  return(seurat.obj)
}

# Reference data for annotation
immgen <- ImmGenData()

## Most recent iteration of the pre-processing pipeline. Earlier iterations are commented below.

# Read and preprocess each sample
dnwt <- preproc.cellranger.with.adt('08-23/DNWT/outs/per_sample_outs/DNWT/count/sample_feature_bc_matrix/')
dnmut <- preproc.cellranger.with.adt('08-23/DNMUT/outs/per_sample_outs/DNMUT/count/sample_feature_bc_matrix/')
totwt <- preproc.cellranger.with.adt('08-23/TOTWT/outs/per_sample_outs/TOTWT/count/sample_feature_bc_matrix/')
totmut <- preproc.cellranger.with.adt('08-23/TOTMUT/outs/per_sample_outs/TOTMUT/count/sample_feature_bc_matrix/')

dnwt <- filter.cells.adt(dnwt, marker='ADT-CD161', threshold = 0.)
dnmut <- filter.cells.adt(dnmut, marker='ADT-CD161', threshold = 0.)
totwt <- filter.cells.adt(totwt, marker='ADT-CD161', threshold = 0.)
totmut <- filter.cells.adt(totmut, marker='ADT-CD161', threshold = 0.)

# Merge the samples together
combined <- merge(x = dnwt, y = dnmut,)
combined <- merge(x = combined, y = totwt,)
combined <- merge(x = combined, y = totmut,)

# Use SCT transformation again on the merged dataset
# Not doing so results in an error because of the mismatched SCT models
combined <- SCTransform(combined, verbose = FALSE)

# Find the most variable features and filter the Seurat object using them
tmp.obj<-FindVariableFeatures(combined,assay="SCT",nfeatures=3000)
variable.genes<-VariableFeatures(tmp.obj,assay="SCT")

combined@assays$SCT@counts <- combined@assays$SCT@counts[variable.genes,]
combined@assays$SCT@data <- combined@assays$SCT@data[variable.genes,]

# Apply MAGIC on the data and replace the data slot in the Seurat object with the enhanced data
cbdt <- t(as.matrix(combined@assays[[DefaultAssay(combined)]]@data))
cb.magic <- Rmagic::magic(cbdt)

newd <- as(t(as.matrix(cb.magic$result)), "sparseMatrix")
combined <- SetAssayData(combined, slot = "data", as(newd, "sparseMatrix"), assay = 'SCT')
combined <- SetAssayData(combined, slot = "counts", as(newd, "sparseMatrix"), assay = 'SCT')

# Annotate cells using Rmagic
rescombined <- annot.singleR(combined, immgen)
combined.a <- rescombined$data

# Reduce dimensionality using vaevictis
# WARNING : a conflict between Rmagic and vaevictis can cause the session to crash
# Rmagic and vaevictis should be used on separate containers until a fix is found
# Save the seurat object enhanced with Rmagic as a .RData file and resume the analysis in another container 
res.va <- vae.seurat(combined.a)
combined.va <- res.va$srt
res.va$model$vae$save()

# Plot cells with their annotations and save to pdf
p = DimPlot(combined.va, reduction = "vae", group.by = "singleR.annot", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('All')
pdf(file = "all4.pdf",width = 20,height = 20)
p
dev.off()

# Add SingleR cell annotations to the original samples
n1 <- dim(dnwt)[2]
n2 <- dim(dnmut)[2] + n1
n3 <- dim(totwt)[2] + n2
n4 <- dim(totmut)[2] + n3

rr <- combined.va@reductions$vae@cell.embeddings[1:n1,]
rownames(rr) <- dnwt@assays[[DefaultAssay(dnwt)]]@data@Dimnames[[2]]
dnwt[["vae"]] <- CreateDimReducObject(embeddings = rr, key = "vae_", assay = DefaultAssay(dnwt))
dnwt@meta.data[["singleR.annot.fromc"]] <- combined.va@meta.data$singleR.annot[1:n1]

rr <- combined.va@reductions$vae@cell.embeddings[(n1+1):n2,]
rownames(rr) <- dnmut@assays[[DefaultAssay(dnmut)]]@data@Dimnames[[2]]
dnmut[["vae"]] <- CreateDimReducObject(embeddings = rr, key = "vae_", assay = DefaultAssay(dnmut))
dnmut@meta.data[["singleR.annot.fromc"]] <- combined.va@meta.data$singleR.annot[(n1+1):n2]

rr <- combined.va@reductions$vae@cell.embeddings[(n2+1):n3,]
rownames(rr) <- totwt@assays[[DefaultAssay(totwt)]]@data@Dimnames[[2]]
totwt[["vae"]] <- CreateDimReducObject(embeddings = rr, key = "vae_", assay = DefaultAssay(totwt))
totwt@meta.data[["singleR.annot.fromc"]] <- combined.va@meta.data$singleR.annot[(n2+1):n3]

rr <- combined.va@reductions$vae@cell.embeddings[(n3+1):n4,]
rownames(rr) <- totmut@assays[[DefaultAssay(totmut)]]@data@Dimnames[[2]]
totmut[["vae"]] <- CreateDimReducObject(embeddings = rr, key = "vae_", assay = DefaultAssay(totmut))
totmut@meta.data[["singleR.annot.fromc"]] <- combined.va@meta.data$singleR.annot[(n3+1):n4]

# Plot samples separately with their annotations and save to pdf
p1 = DimPlot(dnwt, reduction = "vae", group.by = "singleR.annot.fromc", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('dnwt')+xlim(-5, 5)+ylim(-4,5)
p2 = DimPlot(dnmut, reduction = "vae", group.by = "singleR.annot.fromc", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('dnmut')+xlim(-5, 5)+ylim(-4,5)
p3 = DimPlot(totwt, reduction = "vae", group.by = "singleR.annot.fromc", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('totwt')+xlim(-5, 5)+ylim(-4,5)
p4 = DimPlot(totmut, reduction = "vae", group.by = "singleR.annot.fromc", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('totmut')+xlim(-5, 5)+ylim(-4,5)

pdf(file = "4samples.pdf",width = 20,height = 20)
p1+p2+p3+p4
dev.off()

# Group small populations to make plot easier to read
most.abundant <- names(table(combined.va@meta.data$singleR.annot))[unname(which(table(combined.va@meta.data$singleR.annot)>100))]
combined.va@meta.data[["singleR.annot.broad"]] <- sapply(combined.va@meta.data$singleR.annot, function(x){ifelse(x%in%most.abundant, x, "other")})

# add labels for samples to the merged object
cc <- c("DNWT", "DNMUT", "TOTWT", "TOTMUT")
sl <- rep(cc[1], n1)
sl <- c(sl, rep(cc[2], n2-n1))
sl <- c(sl, rep(cc[3], n3-n2))
sl <- c(sl, rep(cc[4], n4-n3))

combined.va@meta.data[["sample.labels"]] <- sl

# remove non-Tcells and Tgd cells from the object
Tcellids <- combined.va@assays$SCT@data@Dimnames[[2]][grep('T cells', combined.va@meta.data$singleR.annot)]
combined.vaT <- subset(combined.va, cells = Tcellids)
p = DimPlot(combined.vaT, reduction = "vae", group.by = "singleR.annot", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('All')
pdf(file = "allT.pdf",width = 20,height = 20)
p
dev.off()

# save cell barcodes to use with the scvelo/velocyto tools
# velocyto only accepts .tsv.gz files, so these files still need to be gzipped via command line
# 4 original samples
write.table(colnames(dnwt@assays[[DefaultAssay(dnwt)]]), 
            file='velocyto_barcodes/DNWT_filt.tsv', 
            quote=FALSE, sep='\t', col.names = NA)
write.table(colnames(dnmut@assays[[DefaultAssay(dnmut)]]), 
            file='velocyto_barcodes/DNMUT_filt.tsv', 
            quote=FALSE, sep='\t', col.names = NA)
write.table(colnames(totwt@assays[[DefaultAssay(totwt)]]), 
            file='velocyto_barcodes/TOTWT_filt.tsv', 
            quote=FALSE, sep='\t', col.names = NA)
write.table(colnames(totmut@assays[[DefaultAssay(totmut)]]), 
            file='velocyto_barcodes/TOTMUT_filt.tsv', 
            quote=FALSE, sep='\t', col.names = NA)
# Full merged sample
write.table(colnames(combined.va@assays[[DefaultAssay(combined.va)]]),
            file='velocyto_barcodes/ALL_filt.tsv',
            quote=FALSE, sep='\t', col.names = NA)
# Merged sample with no non-Tcells and no Tgd
write.table(Tcellids, file='velocyto_barcodes/ALL_filt_onlyT.tsv', 
            quote=FALSE, sep='\t', col.names=NA)



# ### first 
# 
# dnwt <- preproc.cellranger('DNWT/outs/filtered_feature_bc_matrix')
# resdnwt <- annot.singleR(dnwt, immgen)
# dnwt.a <- resdnwt$data
# export.indices(dnwt.a, 'seurat_mice_out/DNWT/')
# dnwt.va <- vae.seurat(dnwt.a)
# write.table(dnwt.va@reductions$vae@cell.embeddings, file='seurat_mice_out/DNWT/vae_red.csv', sep=',',col.names=F)
# 
# p1 = DimPlot(dnwt.va, reduction = "vae", group.by = "singleR.annot", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('DNWT')
# p1
# 
# dnmut <- preproc.cellranger('DNMUT/outs/filtered_feature_bc_matrix')
# resdnmut <- annot.singleR(dnmut, immgen)
# dnmut.a <- resdnmut$data
# export.indices(dnmut.a, 'seurat_mice_out/DNMUT/')
# dnmut.va <- vae.seurat(dnmut.a)
# write.table(dnmut.va@reductions$vae@cell.embeddings, file='seurat_mice_out/DNMUT/vae_red.csv', sep=',',col.names=F)
# 
# p2 = DimPlot(dnmut.va, reduction = "vae", group.by = "singleR.annot", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('DNMUT')
# p2
# 
# 
# totwt <- preproc.cellranger('TOTWT/outs/filtered_feature_bc_matrix')
# restotwt <- annot.singleR(totwt, immgen)
# totwt.a <- restotwt$data
# export.indices(totwt.a, 'seurat_mice_out/TOTWT/')
# totwt.va <- vae.seurat(totwt.a)
# write.table(totwt.va@reductions$vae@cell.embeddings, file='seurat_mice_out/TOTWT/vae_red.csv', sep=',',col.names=F)
# 
# p3 = DimPlot(totwt.va, reduction = "vae", group.by = "singleR.annot", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('TOTWT')
# p3
# 
# 
# totmut <- preproc.cellranger('TOTMUT/outs/filtered_feature_bc_matrix')
# restotmut <- annot.singleR(totmut, immgen)
# totmut.a <- restotmut$data
# export.indices(totmut.a, 'seurat_mice_out/TOTMUT/')
# totmut.va <- vae.seurat(totmut.a)
# write.table(totmut.va@reductions$vae@cell.embeddings, file='seurat_mice_out/TOTMUT/vae_red.csv', sep=',',col.names=F)
# 
# p4 = DimPlot(totmut.va, reduction = "vae", group.by = "singleR.annot", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('TOTMUT')
# p4
# 
# p1 + p2
# p3 + p4
# 
# write.table(colnames(dnwt.va@assays[[DefaultAssay(dnwt.va)]]), 
#             file='velocyto_barcodes/DNWT_filt.tsv', 
#             quote=FALSE, sep='\t', col.names = NA)
# write.table(colnames(dnmut.va@assays[[DefaultAssay(dnmut.va)]]), 
#             file='velocyto_barcodes/DNMUT_filt.tsv', 
#             quote=FALSE, sep='\t', col.names = NA)
# write.table(colnames(totwt.va@assays[[DefaultAssay(totwt.va)]]), 
#             file='velocyto_barcodes/TOTWT_filt.tsv', 
#             quote=FALSE, sep='\t', col.names = NA)
# write.table(colnames(totmut.va@assays[[DefaultAssay(totmut.va)]]), 
#             file='velocyto_barcodes/TOTMUT_filt.tsv', 
#             quote=FALSE, sep='\t', col.names = NA)
# 
# 
# # merge all 4 (seurat merge)
# # preprocess seurat
# # filter using ADT
# # Rmagic
# # dimred vae
# 
# combined <- merge(x = dnwt, y = dnmut,)
# combined <- merge(x = combined, y = totwt,)
# combined <- merge(x = combined, y = totmut,)
# 
# cbdt <- t(as.matrix(combined@assays[[DefaultAssay(combined)]]@data))
# cb.magic <- Rmagic::magic(cbdt)
# 
# newd <- as(t(as.matrix(cb.magic$result)), "sparseMatrix")
# combined <- SetAssayData(combined, slot = "data", as(newd, "sparseMatrix"), assay = 'SCT')
# 
# rescombined <- annot.singleR(combined, immgen)
# combined.a <- rescombined$data
# 
# drop <- c('Epithelial', 'Endothelial', 'Stromal', 'Fibroblast', 'DC ')
# combined.a.trim <- remove.celltypes(combined.va, drop)
# 
# combined.va <- vae.seurat(combined.a.trim, scale=T, center=T)
# 
# p = DimPlot(combined.va, reduction = "vae", group.by = "singleR.annot", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('All')
# p
# 
# 
# ### second
# dnwt <- preproc.cellranger.with.adt('08-23/DNWT/outs/per_sample_outs/DNWT/count/sample_feature_bc_matrix/')
# dnmut <- preproc.cellranger.with.adt('08-23/DNMUT/outs/per_sample_outs/DNMUT/count/sample_feature_bc_matrix/')
# totwt <- preproc.cellranger.with.adt('08-23/TOTWT/outs/per_sample_outs/TOTWT/count/sample_feature_bc_matrix/')
# totmut <- preproc.cellranger.with.adt('08-23/TOTMUT/outs/per_sample_outs/TOTMUT/count/sample_feature_bc_matrix/')
# 
# write.table(colnames(dnwt@assays[[DefaultAssay(dnwt)]]), 
#             file='velocyto_barcodes/DNWT_filt.tsv', 
#             quote=FALSE, sep='\t', col.names = NA)
# write.table(colnames(dnmut@assays[[DefaultAssay(dnmut)]]), 
#             file='velocyto_barcodes/DNMUT_filt.tsv', 
#             quote=FALSE, sep='\t', col.names = NA)
# write.table(colnames(totwt@assays[[DefaultAssay(totwt)]]), 
#             file='velocyto_barcodes/TOTWT_filt.tsv', 
#             quote=FALSE, sep='\t', col.names = NA)
# write.table(colnames(totmut@assays[[DefaultAssay(totmut)]]), 
#             file='velocyto_barcodes/TOTMUT_filt.tsv', 
#             quote=FALSE, sep='\t', col.names = NA)
# 
# combined <- merge(x = dnwt, y = dnmut,)
# combined <- merge(x = combined, y = totwt,)
# combined <- merge(x = combined, y = totmut,)
# 
# tmp.obj<-FindVariableFeatures(combined,assay="RNA",nfeatures=3000)
# variable.genes<-VariableFeatures(tmp.obj,assay="RNA")
# 
# combined@assays$SCT@counts <- combined@assays$SCT@counts[intersect(variable.genes, combined@assays$SCT@counts@Dimnames[[1]]),]
# combined@assays$SCT@data <- combined@assays$SCT@data[intersect(variable.genes, combined@assays$SCT@data@Dimnames[[1]]),]
# 
# cbdt <- t(as.matrix(combined@assays[[DefaultAssay(combined)]]@data))
# cb.magic <- Rmagic::magic(cbdt)
# 
# newd <- as(t(as.matrix(cb.magic$result)), "sparseMatrix")
# combined <- SetAssayData(combined, slot = "data", as(newd, "sparseMatrix"), assay = 'SCT')
# combined <- SetAssayData(combined, slot = "counts", as(newd, "sparseMatrix"), assay = 'SCT')
# 
# rescombined <- annot.singleR(combined, immgen)
# combined.a <- rescombined$data
# 
# drop <- c('Epithelial', 'Endothelial', 'Stromal', 'Fibroblast', 'DC ')
# combined.a.trim <- remove.celltypes(combined.a, drop)
# combined.a.trim <- subset(combined.a.trim, cells=which(!is.na(combined.a.trim$singleR.annot)))
# 
# combined.va <- vae.seurat(combined.a, scale=T)
# 
# p = DimPlot(combined.va, reduction = "vae", group.by = "singleR.annot.light", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('All')
# p
# 
# n1 <- dim(dnwt)[2]
# n2 <- dim(dnmut)[2] + n1
# n3 <- dim(totwt)[2] + n2
# n4 <- dim(totmut)[2] + n3
# 
# rr <- combined.va@reductions$vae@cell.embeddings[1:n1,]
# rownames(rr) <- dnwt@assays[[DefaultAssay(dnwt)]]@data@Dimnames[[2]]
# dnwt[["vae"]] <- CreateDimReducObject(embeddings = rr, key = "vae_", assay = DefaultAssay(dnwt))
# dnwt@meta.data[["singleR.annot"]] <- combined.va@meta.data$singleR.annot[1:n1]
# dnwt@meta.data[["singleR.annot.light"]] <- combined.va@meta.data$singleR.annot.light[1:n1]
# 
# rr <- combined.va@reductions$vae@cell.embeddings[(n1+1):n2,]
# rownames(rr) <- dnmut@assays[[DefaultAssay(dnmut)]]@data@Dimnames[[2]]
# dnmut[["vae"]] <- CreateDimReducObject(embeddings = rr, key = "vae_", assay = DefaultAssay(dnmut))
# dnmut@meta.data[["singleR.annot"]] <- combined.va@meta.data$singleR.annot[(n1+1):n2]
# dnmut@meta.data[["singleR.annot.light"]] <- combined.va@meta.data$singleR.annot.light[(n1+1):n2]
# 
# rr <- combined.va@reductions$vae@cell.embeddings[(n2+1):n3,]
# rownames(rr) <- totwt@assays[[DefaultAssay(totwt)]]@data@Dimnames[[2]]
# totwt[["vae"]] <- CreateDimReducObject(embeddings = rr, key = "vae_", assay = DefaultAssay(totwt))
# totwt@meta.data[["singleR.annot"]] <- combined.va@meta.data$singleR.annot[(n2+1):n3]
# totwt@meta.data[["singleR.annot.light"]] <- combined.va@meta.data$singleR.annot.light[(n2+1):n3]
# 
# rr <- combined.va@reductions$vae@cell.embeddings[(n3+1):n4,]
# rownames(rr) <- totmut@assays[[DefaultAssay(totmut)]]@data@Dimnames[[2]]
# totmut[["vae"]] <- CreateDimReducObject(embeddings = rr, key = "vae_", assay = DefaultAssay(totmut))
# totmut@meta.data[["singleR.annot"]] <- combined.va@meta.data$singleR.annot[(n3+1):n4]
# totmut@meta.data[["singleR.annot.light"]] <- combined.va@meta.data$singleR.annot.light[(n3+1):n4]
# 
# 
# 
# p1 = DimPlot(dnwt, reduction = "vae", group.by = "singleR.annot.light", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('DNWT')
# p2 = DimPlot(dnmut, reduction = "vae", group.by = "singleR.annot.light", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('DNMUT')
# p3 = DimPlot(totwt, reduction = "vae", group.by = "singleR.annot.light", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('TOTWT')
# p4 = DimPlot(totmut, reduction = "vae", group.by = "singleR.annot.light", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('TOTMUT')
# 
# p1+p2+p3+p4
# 
# # third
# 
# dnwt <- preproc.cellranger.with.adt('08-23/DNWT/outs/per_sample_outs/DNWT/count/sample_feature_bc_matrix/')
# dnmut <- preproc.cellranger.with.adt('08-23/DNMUT/outs/per_sample_outs/DNMUT/count/sample_feature_bc_matrix/')
# totwt <- preproc.cellranger.with.adt('08-23/TOTWT/outs/per_sample_outs/TOTWT/count/sample_feature_bc_matrix/')
# totmut <- preproc.cellranger.with.adt('08-23/TOTMUT/outs/per_sample_outs/TOTMUT/count/sample_feature_bc_matrix/')
# 
# # dnwt
# tmp.obj<-FindVariableFeatures(dnwt,assay="SCT",nfeatures=3000)
# variable.genes<-VariableFeatures(tmp.obj,assay="SCT")
# 
# dnwt@assays$SCT@counts <- dnwt@assays$SCT@counts[variable.genes,]
# dnwt@assays$SCT@data <- dnwt@assays$SCT@data[variable.genes,]
# 
# dnwt <- filter.cells.adt(dnwt, marker='ADT-CD161', threshold = 0.)
# 
# dt <- t(as.matrix(dnwt@assays[[DefaultAssay(dnwt)]]@data))
# dt.magic <- Rmagic::magic(dt)
# 
# newd <- as(t(as.matrix(dt.magic$result)), "sparseMatrix")
# dnwt <- SetAssayData(dnwt, slot = "data", as(newd, "sparseMatrix"), assay = 'SCT')
# dnwt <- SetAssayData(dnwt, slot = "counts", as(newd, "sparseMatrix"), assay = 'SCT')
# 
# resdnwt <- annot.singleR(dnwt, immgen)
# dnwt.a <- resdnwt$data
# 
# # dnmut
# tmp.obj<-FindVariableFeatures(dnmut,assay="SCT",nfeatures=3000)
# variable.genes<-VariableFeatures(tmp.obj,assay="SCT")
# 
# dnmut@assays$SCT@counts <- dnmut@assays$SCT@counts[variable.genes,]
# dnmut@assays$SCT@data <- dnmut@assays$SCT@data[variable.genes,]
# 
# dnmut <- filter.cells.adt(dnmut, marker='ADT-CD161', threshold = 0.)
# 
# dt <- t(as.matrix(dnmut@assays[[DefaultAssay(dnmut)]]@data))
# dt.magic <- Rmagic::magic(dt)
# 
# newd <- as(t(as.matrix(dt.magic$result)), "sparseMatrix")
# dnmut <- SetAssayData(dnmut, slot = "data", as(newd, "sparseMatrix"), assay = 'SCT')
# dnmut <- SetAssayData(dnmut, slot = "counts", as(newd, "sparseMatrix"), assay = 'SCT')
# 
# resdnmut <- annot.singleR(dnmut, immgen)
# dnmut.a <- resdnmut$data
# 
# # totwt
# tmp.obj<-FindVariableFeatures(totwt,assay="SCT",nfeatures=3000)
# variable.genes<-VariableFeatures(tmp.obj,assay="SCT")
# 
# totwt@assays$SCT@counts <- totwt@assays$SCT@counts[variable.genes,]
# totwt@assays$SCT@data <- totwt@assays$SCT@data[variable.genes,]
# 
# totwt <- filter.cells.adt(totwt, marker='ADT-CD161', threshold = 0.)
# 
# dt <- t(as.matrix(totwt@assays[[DefaultAssay(totwt)]]@data))
# dt.magic <- Rmagic::magic(dt)
# 
# newd <- as(t(as.matrix(dt.magic$result)), "sparseMatrix")
# totwt <- SetAssayData(totwt, slot = "data", as(newd, "sparseMatrix"), assay = 'SCT')
# totwt <- SetAssayData(totwt, slot = "counts", as(newd, "sparseMatrix"), assay = 'SCT')
# 
# restotwt <- annot.singleR(totwt, immgen)
# totwt.a <- restotwt$data
# 
# # totmut
# tmp.obj<-FindVariableFeatures(totmut,assay="SCT",nfeatures=3000)
# variable.genes<-VariableFeatures(tmp.obj,assay="SCT")
# 
# totmut@assays$SCT@counts <- totmut@assays$SCT@counts[variable.genes,]
# totmut@assays$SCT@data <- totmut@assays$SCT@data[variable.genes,]
# 
# totmut <- filter.cells.adt(totmut, marker='ADT-CD161', threshold = 0.)
# 
# dt <- t(as.matrix(totmut@assays[[DefaultAssay(totmut)]]@data))
# dt.magic <- Rmagic::magic(dt)
# 
# newd <- as(t(as.matrix(dt.magic$result)), "sparseMatrix")
# totmut <- SetAssayData(totmut, slot = "data", as(newd, "sparseMatrix"), assay = 'SCT')
# totmut <- SetAssayData(totmut, slot = "counts", as(newd, "sparseMatrix"), assay = 'SCT')
# 
# restotmut <- annot.singleR(totmut, immgen)
# totmut.a <- restotmut$data
# 
# 
# # dnwt.va <- vae.seurat(dnwt.a)
# # 
# # p1 = DimPlot(dnwt.va, reduction = "vae", group.by = "singleR.annot", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('All')
# # p1
# # 
# # dnmut.va <- vae.seurat(dnmut.a)
# # 
# # p2 = DimPlot(dnmut.va, reduction = "vae", group.by = "singleR.annot", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('All')
# # p2
# # 
# # totwt.va <- vae.seurat(totwt.a)
# # 
# # p3 = DimPlot(totwt.va, reduction = "vae", group.by = "singleR.annot", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('All')
# # p3
# # 
# # totmut.va <- vae.seurat(totmut.a)
# # 
# # p4 = DimPlot(totmut.va, reduction = "vae", group.by = "singleR.annot", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('All')
# # p4
# 
# # fourth
# 
# dnwt <- preproc.cellranger.with.adt('08-23/DNWT/outs/per_sample_outs/DNWT/count/sample_feature_bc_matrix/')
# dnmut <- preproc.cellranger.with.adt('08-23/DNMUT/outs/per_sample_outs/DNMUT/count/sample_feature_bc_matrix/')
# totwt <- preproc.cellranger.with.adt('08-23/TOTWT/outs/per_sample_outs/TOTWT/count/sample_feature_bc_matrix/')
# totmut <- preproc.cellranger.with.adt('08-23/TOTMUT/outs/per_sample_outs/TOTMUT/count/sample_feature_bc_matrix/')
# 
# # dnwt
# tmp.obj<-FindVariableFeatures(dnwt,assay="SCT",nfeatures=3000)
# variable.genes<-VariableFeatures(tmp.obj,assay="SCT")
# 
# dnwt@assays$SCT@counts <- dnwt@assays$SCT@counts[variable.genes,]
# dnwt@assays$SCT@data <- dnwt@assays$SCT@data[variable.genes,]
# 
# dnwt <- filter.cells.adt(dnwt, marker='ADT-CD161', threshold = 0.)
# 
# # dnmut
# tmp.obj<-FindVariableFeatures(dnmut,assay="SCT",nfeatures=3000)
# variable.genes<-VariableFeatures(tmp.obj,assay="SCT")
# 
# dnmut@assays$SCT@counts <- dnmut@assays$SCT@counts[variable.genes,]
# dnmut@assays$SCT@data <- dnmut@assays$SCT@data[variable.genes,]
# 
# dnmut <- filter.cells.adt(dnmut, marker='ADT-CD161', threshold = 0.)
# 
# # totwt
# tmp.obj<-FindVariableFeatures(totwt,assay="SCT",nfeatures=3000)
# variable.genes<-VariableFeatures(tmp.obj,assay="SCT")
# 
# totwt@assays$SCT@counts <- totwt@assays$SCT@counts[variable.genes,]
# totwt@assays$SCT@data <- totwt@assays$SCT@data[variable.genes,]
# 
# totwt <- filter.cells.adt(totwt, marker='ADT-CD161', threshold = 0.)
# 
# # totmut
# tmp.obj<-FindVariableFeatures(totmut,assay="SCT",nfeatures=3000)
# variable.genes<-VariableFeatures(tmp.obj,assay="SCT")
# 
# totmut@assays$SCT@counts <- totmut@assays$SCT@counts[variable.genes,]
# totmut@assays$SCT@data <- totmut@assays$SCT@data[variable.genes,]
# 
# totmut <- filter.cells.adt(totmut, marker='ADT-CD161', threshold = 0.)
# 
# write.table(colnames(dnwt@assays[[DefaultAssay(dnwt)]]), 
#             file='velocyto_barcodes/DNWT_filt.tsv', 
#             quote=FALSE, sep='\t', col.names = NA)
# write.table(colnames(dnmut@assays[[DefaultAssay(dnmut)]]), 
#             file='velocyto_barcodes/DNMUT_filt.tsv', 
#             quote=FALSE, sep='\t', col.names = NA)
# write.table(colnames(totwt@assays[[DefaultAssay(totwt)]]), 
#             file='velocyto_barcodes/TOTWT_filt.tsv', 
#             quote=FALSE, sep='\t', col.names = NA)
# write.table(colnames(totmut@assays[[DefaultAssay(totmut)]]), 
#             file='velocyto_barcodes/TOTMUT_filt.tsv', 
#             quote=FALSE, sep='\t', col.names = NA)
# 
# combined <- merge(x = dnwt, y = dnmut,)
# combined <- merge(x = combined, y = totwt,)
# combined <- merge(x = combined, y = totmut,)
# 
# combined <- filter.cells.adt(combined, marker='ADT-CD161', threshold = 0.)
# 
# cbdt <- t(as.matrix(combined@assays[[DefaultAssay(combined)]]@data))
# cb.magic <- Rmagic::magic(cbdt)
# 
# newd <- as(t(as.matrix(cb.magic$result)), "sparseMatrix")
# combined <- SetAssayData(combined, slot = "data", as(newd, "sparseMatrix"), assay = 'SCT')
# combined <- SetAssayData(combined, slot = "counts", as(newd, "sparseMatrix"), assay = 'SCT')
# 
# rescombined <- annot.singleR(combined, immgen)
# combined.a <- rescombined$data
# 
# combined.va <- vae.seurat(combined.a)
# 
# p = DimPlot(combined.va, reduction = "vae", group.by = "singleR.annot", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('All')
# pdf(file = "all2.pdf",width = 20,height = 20)
# p
# dev.off()
# 
# 
# n1 <- dim(dnwt)[2]
# n2 <- dim(dnmut)[2] + n1
# n3 <- dim(totwt)[2] + n2
# n4 <- dim(totmut)[2] + n3
# 
# rr <- combined.va@reductions$vae@cell.embeddings[1:n1,]
# rownames(rr) <- dnwt@assays[[DefaultAssay(dnwt)]]@data@Dimnames[[2]]
# dnwt[["vae"]] <- CreateDimReducObject(embeddings = rr, key = "vae_", assay = DefaultAssay(dnwt))
# dnwt@meta.data[["singleR.annot.fromc"]] <- combined.va@meta.data$singleR.annot[1:n1]
# 
# rr <- combined.va@reductions$vae@cell.embeddings[(n1+1):n2,]
# rownames(rr) <- dnmut@assays[[DefaultAssay(dnmut)]]@data@Dimnames[[2]]
# dnmut[["vae"]] <- CreateDimReducObject(embeddings = rr, key = "vae_", assay = DefaultAssay(dnmut))
# dnmut@meta.data[["singleR.annot.fromc"]] <- combined.va@meta.data$singleR.annot[(n1+1):n2]
# 
# rr <- combined.va@reductions$vae@cell.embeddings[(n2+1):n3,]
# rownames(rr) <- totwt@assays[[DefaultAssay(totwt)]]@data@Dimnames[[2]]
# totwt[["vae"]] <- CreateDimReducObject(embeddings = rr, key = "vae_", assay = DefaultAssay(totwt))
# totwt@meta.data[["singleR.annot.fromc"]] <- combined.va@meta.data$singleR.annot[(n2+1):n3]
# 
# rr <- combined.va@reductions$vae@cell.embeddings[(n3+1):n4,]
# rownames(rr) <- totmut@assays[[DefaultAssay(totmut)]]@data@Dimnames[[2]]
# totmut[["vae"]] <- CreateDimReducObject(embeddings = rr, key = "vae_", assay = DefaultAssay(totmut))
# totmut@meta.data[["singleR.annot.fromc"]] <- combined.va@meta.data$singleR.annot[(n3+1):n4]
# 
# p1 = DimPlot(dnwt, reduction = "vae", group.by = "singleR.annot.fromc", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('dnwt')+xlim(-5, 5)+ylim(-4,5)
# p2 = DimPlot(dnmut, reduction = "vae", group.by = "singleR.annot.fromc", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('dnmut')+xlim(-5, 5)+ylim(-4,5)
# p3 = DimPlot(totwt, reduction = "vae", group.by = "singleR.annot.fromc", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('totwt')+xlim(-5, 5)+ylim(-4,5)
# p4 = DimPlot(totmut, reduction = "vae", group.by = "singleR.annot.fromc", label = TRUE, label.size = 3 ,repel = TRUE) + NoLegend() + ggtitle('totmut')+xlim(-5, 5)+ylim(-4,5)
# 
# pdf(file = "dnwt.pdf",width = 20,height = 20)
# p1+p2+p3+p4
# dev.off()