library(Biobase)
library(flowCore)
library(flowViz)

## Creates a flowframe .fcs from a seurat object data
## Replace pbmc3k by the name of the seurat object to analyze

# data
# can be any assay
adt.full <- pbmc3k@assays$ADT

# conf <- read.csv('seurat/seurat_data/Thymus/ADT_feature_ref.csv')
# filt.adt <- adt.full[gsub('_', '-', conf$name),]

dta <- t(adt.full[,])

# add dimred

dta <- cbind(dta, pbmc3k@reductions$vae@cell.embeddings)
colnames(dta)[c(ncol(dta)-1, ncol(dta))] <- c('vae1', 'vae2')

# you need to prepare some metadata
meta <- data.frame(name=dimnames(dta)[[2]],
                   desc=dimnames(dta)[[2]]
)
meta$range <- apply(apply(dta,2,range),2,diff)
meta$minRange <- apply(dta,2,min)
meta$maxRange <- apply(dta,2,max)

head(meta)


# meta <- rbind(meta, c('predicted.celltype', 'predicted.celltype', 0, 0, 0))

# a flowFrame is the internal representation of a FCS file
ff <- new("flowFrame",
          exprs=dta,
          parameters=AnnotatedDataFrame(meta)
)

# now you can save it back to the filesystem
write.FCS(ff,'thymus2_adt.fcs')
