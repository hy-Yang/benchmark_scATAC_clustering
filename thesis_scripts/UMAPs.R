#plot clustering results
library(CellBench)
library(Signac)
library(Seurat)

library(ggplot2)
library(patchwork)
set.seed(1234)
set1_res <- strip_timing(clustering_result[1:8, ])


sce <- set1_res$result[[1]]
sce
counts <- counts(sce)
metadata <- as.data.frame(colData(sce))

metadata$cicero_log_LSI = set1_res$result[[1]]$clustering_res
metadata$cicero_size_LSI = set1_res$result[[2]]$clustering_res
metadata$cicero_log_PCA  = set1_res$result[[3]]$clustering_res
metadata$cicero_size_PCA = set1_res$result[[4]]$clustering_res
metadata$destin  = set1_res$result[[5]]$clustering_res
metadata$scABC = set1_res$result[[6]]$clustering_res
metadata$signac  = set1_res$result[[7]]$clustering_res
metadata$epiconv = set1_res$result[[8]]$clustering_res
pbmc$sc3= sce_1$sc3_7_clusters



pbmc <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'cellmix_669',
  min.cells = 1,
  meta.data = metadata
)
pbmc$sc3= sce_1$sc3_7_clusters
pbmc$seurat=  Idents(seurat_f1_669)
pbmc$RaceID =  sc@cluster$kpart

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(object = pbmc, reduction = 'svd', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'svd', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
pbmc$celltype <- as.character(pbmc$celltype)
pbmc$celltype[pbmc$celltype == "H8383"] <- "H838"
pbmc$celltype <- factor(pbmc$celltype)
par(mfrow=c(3,3))
DimPlot(object = pbmc, group.by = "celltype", label = F)
DimPlot(object = pbmc, group.by = "cicero_log_LSI", label = F)
DimPlot(object = pbmc, group.by = "cicero_size_LSI", label = F)
DimPlot(object = pbmc, group.by = "cicero_log_PCA", label = F)
DimPlot(object = pbmc, group.by = "cicero_size_PCA", label = F)
DimPlot(object = pbmc, group.by = "destin", label = F)
DimPlot(object = pbmc, group.by = "signac", label = F)
DimPlot(object = pbmc, group.by = "scABC", label = F)
DimPlot(object = pbmc, group.by = "epiconv", label = F)
DimPlot(object = pbmc, group.by = "sc3", label = F)
DimPlot(object = pbmc, group.by = "seurat", label = F)
DimPlot(object = pbmc, group.by = "RaceID", label = F)
res_ARI <- set1_res %>%
  apply_methods(clustering_evaluation)


set3_res <- strip_timing(res_c_v2[17:24, ])

sce <- set3_res$result[[1]]
sce
counts <- counts(sce)
metadata <- as.data.frame(colData(sce))
metadata$cicero_log_LSI = set3_res$result[[1]]$clustering_res
metadata$cicero_size_LSI = set3_res$result[[2]]$clustering_res
metadata$cicero_log_PCA  = set3_res$result[[3]]$clustering_res
metadata$cicero_size_PCA = set3_res$result[[4]]$clustering_res
metadata$destin  = set3_res$result[[5]]$clustering_res
metadata$scABC = set3_res$result[[6]]$clustering_res
metadata$signac  = set3_res$result[[7]]$clustering_res
metadata$epiconv = set3_res$result[[8]]$clustering_res

pbmc <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'cellmix_669',
  min.cells = 1,
  meta.data = metadata
)
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

pbmc <- RunUMAP(object = pbmc, reduction = 'svd', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'svd', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
pbmc$celltype <- as.character(pbmc$celltype)
pbmc$celltype[pbmc$celltype == "H8383"] <- "H838"
pbmc$celltype <- factor(pbmc$celltype)
DimPlot(object = pbmc, group.by = "celltype", label = T)
DimPlot(object = pbmc, group.by = "cicero_log_LSI", label = T)
DimPlot(object = pbmc, group.by = "cicero_size_LSI", label = T)
DimPlot(object = pbmc, group.by = "cicero_log_PCA", label = T)
DimPlot(object = pbmc, group.by = "cicero_size_PCA", label = T)
DimPlot(object = pbmc, group.by = "destin", label = T)
DimPlot(object = pbmc, group.by = "scABC", label = T)
DimPlot(object = pbmc, group.by = "signac", label = T)
DimPlot(object = pbmc, group.by = "epiconv", label = T)

set4_res <- strip_timing(res_c_v2[25:32, ])

sce <- set4_res$result[[1]]
sce
counts <- counts(sce)
metadata <- as.data.frame(colData(sce))
metadata$cicero_log_LSI = set4_res$result[[1]]$clustering_res
metadata$cicero_size_LSI = set4_res$result[[2]]$clustering_res
metadata$cicero_log_PCA  = set4_res$result[[3]]$clustering_res
metadata$cicero_size_PCA = set4_res$result[[4]]$clustering_res
metadata$destin  = set4_res$result[[5]]$clustering_res
metadata$scABC = set4_res$result[[6]]$clustering_res
metadata$signac  = set4_res$result[[7]]$clustering_res
metadata$epiconv = set4_res$result[[8]]$clustering_res

pbmc <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'cellmix_669',
  min.cells = 1,
  meta.data = metadata
)
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

pbmc <- RunUMAP(object = pbmc, reduction = 'svd', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'svd', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
pbmc$celltype <- as.character(pbmc$celltype)
pbmc$celltype[pbmc$celltype == "H8383"] <- "H838"
pbmc$celltype <- factor(pbmc$celltype)
DimPlot(object = pbmc, group.by = "celltype", label = F)
DimPlot(object = pbmc, group.by = "cicero_log_LSI", label = F)
DimPlot(object = pbmc, group.by = "cicero_size_LSI", label = F)
DimPlot(object = pbmc, group.by = "cicero_log_PCA", label = F)
DimPlot(object = pbmc, group.by = "cicero_size_PCA", label = F)
DimPlot(object = pbmc, group.by = "destin", label = F)
DimPlot(object = pbmc, group.by = "scABC", label = F)
DimPlot(object = pbmc, group.by = "signac", label = F)
DimPlot(object = pbmc, group.by = "epiconv", label = F)
pbmc$scABC
set4_res$result[[6]]$clustering_res

res_ARI <- set4_res %>%
  apply_methods(clustering_evaluation)
round(res_ARI$result,2)

