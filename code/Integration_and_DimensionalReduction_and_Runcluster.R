# R code
# Yandong Zheng


# Harmony integration
library(Seurat)
library(hdf5r)
library(DoubletFinder)
library(RColorBrewer)
library(dplyr)
library(patchwork)
library(cowplot)
library(reshape2)
library(magrittr)

library(future)
future::plan(strategy = 'multicore', workers = 10)
options(future.globals.maxSize = 1000*1024^3)

save_wd <- "./downstream_analysis/all_sample_rds/"
file.list <- list.files(save_wd)


scRNAlist <- list()
for (i in file.list) {
  scRNAlist[[i]] <- readRDS(paste0(save_wd, file.list))
  DefaultAssay(scRNAlist[[i]]) <- "SCT"
  scRNAlist[[i]] <- subset(scRNAlist[[i]], subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & 
                             percent.mt < 5 & Double %in% "Singlet")
}
saveRDS(scRNAlist, paste0(save_wd, "Human_heart_scRNAlist_combined.rds"))


## RunPCA
var.features <- SelectIntegrationFeatures(object.list = scRNAlist)

scRNA_harmony <- merge(x = scRNAlist[[1]], y = scRNAlist[2:length(scRNAlist)], merge.data = TRUE)
saveRDS(scRNA_harmony,paste0(save_wd,"Human_heart_scRNAlist_merge.rds"))
gc()

scRNA_harmony <- RunPCA(object = scRNA_harmony, assay = "SCT", features = var.features, 
                        npcs = 100)

##integration
system.time({
  scRNA_harmony <- RunHarmony(scRNA_harmony, assay.use = "SCT", group.by.vars = "orig.ident", 
                              plot_convergence = TRUE, dims.use = 1:100)
})

saveRDS(scRNA_harmony, paste0(save_wd,"Human_heart_merge_object_RunHarmony.rds"))
gc()



# Dimensional Reduction and Runcluster
scRNA_harmony@meta.data$sample <- scRNA_harmony@meta.data$orig.ident
p1 <- DimPlot(object = scRNA_harmony, reduction = "harmony", pt.size = .1, group.by = "sample")
p2 <- VlnPlot(object = scRNA_harmony, features = "harmony_1", group.by = "sample", pt.size = 0)


png(paste(save_wd,"Human_heart_scRNA_harmony_DimHeatmap.png",sep = ""),height = 8000,width = 3000)
DimHeatmap(scRNA_harmony, dims = 1:100, cells = 500, balanced = TRUE)
dev.off()


scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:50)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:50)
scRNA_harmony <- FindClusters(scRNA_harmony, resolution = 2.5)
gc()
saveRDS(scRNA_harmony,paste0(save_wd, "Human_heart_Dimensional_Reduction_Runcluster.rds"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# Differential expression and marker selection 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
DefaultAssay(seurat.object) <- "RNA" 
seurat.object <- NormalizeData(seurat.object) 

Idents(seurat.object) <- "seurat_clusters"
cluster.marker <- FindAllMarkers(seurat.object, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
cluster.marker.final <- subset(cluster.marker, p_val_adj < 0.05 & avg_log2FC > 0.5)

