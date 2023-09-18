# R code
# Yandong Zheng
# Remove doublecell


# Load packages and data
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
plan("multiprocess", workers = 20)

args = commandArgs(trailingOnly = T)
type <- args[1]
sample <- args[2]

work_wd <- paste0("./h5file_after_cellranger/", type, "/")
save_wd <- paste0("./downstream_analysis/", type, "/")
dir.create(save_wd)

QC_wd <- paste(save_wd,"QC_plot/",sep = "")
dir.create(QC_wd)

# Load the dataset
print("Load the dataset")
mat.object <- Read10X_h5(filename = paste0(work_wd, sample, "/output_filtered.h5"))
seurat.object <- CreateSeuratObject(counts = mat.object, project = paste0(type, "_", sample), min.cells = 3, min.features = 200)
seurat.object <- RenameCells(seurat.object, add.cell.id = paste0(type, "_", sample))
seurat.object[["percent.mt"]] <- PercentageFeatureSet(seurat.object, pattern = "^MT-")

doublecell_save_wd <- paste0(save_wd, "DoublrFider_output/")
dir.create(doublecell_save_wd)

print("Run SCTransform")
seurat.object <- SCTransform(seurat.object, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = T)
seurat.object <- RunPCA(seurat.object, verbose = T)
pc.num = 1:10
seurat.object <- FindNeighbors(seurat.object, dims = pc.num, verbose = T)
seurat.object <- FindClusters(seurat.object, resolution = 1)
seurat.object <- RunUMAP(seurat.object, dims = pc.num)

print("Run doublefinder")

sweep.res.list <- paramSweep_v3(seurat.object, PCs = pc.num, sct = TRUE) 
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
print("The step1 is over")

DoubletRate = ncol(seurat.object) * 7.5 * 1e-6 
annotations <- seurat.object@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(DoubletRate*ncol(seurat.object))
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
print("The step2 is over")

pN_value <- 0.25
pK_value <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
seurat.object <- doubletFinder_v3(seurat.object, PCs = pc.num, pN = pN_value, pK = pK_value, 
                                  nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
seurat.object <- doubletFinder_v3(seurat.object, PCs = pc.num, pN = pN_value, pK = pK_value, 
                                  nExp = nExp_poi.adj, reuse.pANN = pANN_value, sct = T)
print("The step3 is over")

double_adi_colname <- paste0("DF.classifications_", pN_value, "_", pK_value, '_', nExp_poi.adj)
seurat.object@meta.data$Double <- seurat.object[[]][,double_adi_colname]
print(table(seurat.object@meta.data$Double))
saveRDS(seurat.object, paste0(doublecell_save_wd, sample, "_anno_doublecell.rds"))


seurat.object <- subset(seurat.object, Double %in% "Singlet")
delete_col <- c(grep(pattern = "^pANN", x = colnames(seurat.object@meta.data), value = F),
                grep(pattern = "^DF.classifications", x = colnames(seurat.object@meta.data), value = F))
seurat.object@meta.data %<>% .[,-delete_col]
saveRDS(seurat.object, paste0(doublecell_save_wd, sample, "_remove_doublecell.rds"))

sub_QC_wd <- paste0(QC_wd, sample, "/")
dir.create(sub_QC_wd)

pdf(paste0(sub_QC_wd, sample, ".pdf"), height = 5, width = 12)
print(VlnPlot(seurat.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
plot1 <- FeatureScatter(seurat.object, feature1 = "nFeature_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(plot_grid(plot1, plot2))
dev.off()




