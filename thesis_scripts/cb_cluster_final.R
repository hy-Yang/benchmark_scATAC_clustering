# clustering using cellbench
# input:  peak by cell matrix in sce


library(CellBench)

# ----cicero----
library(cicero)
library(monocle3)
cicero_log_LSI_UMAP <- function(sce) {

  indata <- counts(sce)
  # format cell info
  cellinfo <- data.frame(sce$barcode, stringsAsFactors = FALSE)
  row.names(cellinfo) <- cellinfo$sce.barcode
  names(cellinfo) <- "cells"
  # format peak info
  peakinfo <- rowData(sce)
  # make CDS
  input_cds <- suppressWarnings(new_cell_data_set(indata,
                                                  cell_metadata = cellinfo,
                                                  gene_metadata = peakinfo))
  input_cds$celltype <- colData(sce)$celltype
  # binarize the matrix
  counts(input_cds)@x[counts(input_cds)@x > 0] <- 1
  # Ensure there are no peaks included with zero reads
  input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]
  
  set.seed(2017)
  input_cds <- detect_genes(input_cds)
  input_cds <- estimate_size_factors(input_cds)
  input_cds <- preprocess_cds(input_cds, method = "LSI")
  input_cds <- reduce_dimension(input_cds, reduction_method = "UMAP",
                                preprocess_method = "LSI")
  input_cds <- cluster_cells(input_cds)
  
  colData(sce)$clustering_res <- factor(input_cds@clusters$UMAP$cluster_result$optim_res$membership)
  print('cicero_log_LSI_UMAP done!')

  return(sce)
}
cicero_size_LSI_UMAP <- function(sce) {
  indata <- counts(sce)
  # format cell info
  cellinfo <- data.frame(sce$barcode, stringsAsFactors = FALSE)
  row.names(cellinfo) <- cellinfo$sce.barcode
  names(cellinfo) <- "cells"
  # format peak info
  peakinfo <- rowData(sce)
  # make CDS
  input_cds <- suppressWarnings(new_cell_data_set(indata,
                                                  cell_metadata = cellinfo,
                                                  gene_metadata = peakinfo))
  input_cds$celltype <- colData(sce)$celltype
  # binarize the matrix
  counts(input_cds)@x[counts(input_cds)@x > 0] <- 1
  #Ensure there are no peaks included with zero reads
  input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]
  
  set.seed(2017)
  input_cds <- detect_genes(input_cds)
  input_cds <- estimate_size_factors(input_cds)
  input_cds <- preprocess_cds(input_cds, norm_method = 'size_only', method = "LSI")
  input_cds <- reduce_dimension(input_cds, reduction_method = "UMAP",
                                preprocess_method = "LSI")
  input_cds <- cluster_cells(input_cds)
  
  colData(sce)$clustering_res <- factor(input_cds@clusters$UMAP$cluster_result$optim_res$membership)
  
  print('cicero_size_LSI_UMAP done!')
  return(sce)
}
cicero_log_PCA_UMAP <- function(sce) {
  indata <- counts(sce)
  # format cell info
  cellinfo <- data.frame(sce$barcode, stringsAsFactors = FALSE)
  row.names(cellinfo) <- cellinfo$sce.barcode
  names(cellinfo) <- "cells"
  # format peak info
  peakinfo <- rowData(sce)
  # make CDS
  input_cds <- suppressWarnings(new_cell_data_set(indata,
                                                  cell_metadata = cellinfo,
                                                  gene_metadata = peakinfo))
  input_cds$celltype <- colData(sce)$celltype
  # binarize the matrix
  counts(input_cds)@x[counts(input_cds)@x > 0] <- 1
  #Ensure there are no peaks included with zero reads
  input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]
  
  set.seed(2017)
  input_cds <- detect_genes(input_cds)
  input_cds <- estimate_size_factors(input_cds)
  input_cds <- preprocess_cds(input_cds, norm_method = 'log', method = "PCA")
  input_cds <- reduce_dimension(input_cds, reduction_method = "UMAP",
                                preprocess_method = "PCA")
  input_cds <- cluster_cells(input_cds)
  
  colData(sce)$clustering_res <- factor(input_cds@clusters$UMAP$cluster_result$optim_res$membership)
  
  print('cicero_log_PCA_UMAP done!')
  return(sce)
}
cicero_size_PCA_UMAP <- function(sce) {
  indata <- counts(sce)
  # format cell info
  cellinfo <- data.frame(sce$barcode, stringsAsFactors = FALSE)
  row.names(cellinfo) <- cellinfo$sce.barcode
  names(cellinfo) <- "cells"
  # format peak info
  peakinfo <- rowData(sce)
  # make CDS
  input_cds <- suppressWarnings(new_cell_data_set(indata,
                                                  cell_metadata = cellinfo,
                                                  gene_metadata = peakinfo))
  input_cds$celltype <- colData(sce)$celltype
  # binarize the matrix
  counts(input_cds)@x[counts(input_cds)@x > 0] <- 1
  #Ensure there are no peaks included with zero reads
  input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]
  
  set.seed(2017)
  input_cds <- detect_genes(input_cds)
  input_cds <- estimate_size_factors(input_cds)
  input_cds <- preprocess_cds(input_cds, norm_method = 'size_only', method = "PCA")
  input_cds <- reduce_dimension(input_cds, reduction_method = "UMAP",
                                preprocess_method = "PCA")
  input_cds <- cluster_cells(input_cds)
  
  colData(sce)$clustering_res <- factor(input_cds@clusters$UMAP$cluster_result$optim_res$membership)
  
  print('cicero_size_PCA_UMAP done!')
  return(sce)
}
# ----destin----
library(destin)
destin_elbow_2to20 <- function(sce){
  countsMat = counts(sce)
  countsMat = as(countsMat, "dgCMatrix")
  countsMat[countsMat > 1] = 1
  bedData = with(rowData(sce), GRanges(chr, IRanges(bp1, bp2)))
  colData = colData(sce)
  rse = SummarizedExperiment::SummarizedExperiment(assays=list(counts=countsMat),rowRanges=bedData, colData=colData)
  # filter out Y chrom
  #yIndex = seqnames(rowRanges(rse)) == "chrY"
  #rse = rse[!yIndex]
  
  model = "hg38" # choose from hg19, hg38, mm10
  rse = annotateRSE(rse, model)
  rse = doQC(rse, regionSumCutoff = 5, cellSumCutoffSDs = 10) # cellSumCutoffSDs = 3, set to 3 to remove some cells
  clusterEst = estimateNClusters(rse, nClustersRange = 2:20)
  nClusters = clusterEst$nClustersList$logLikeElbow
  
  nCores = 8
  clusterResults = destinGrid(rse, nClusters = nClusters, nCores = nCores)
  
  colData(sce)$clustering_res <- clusterResults$cluster$cluster
  print('destin_elbow_2to20 done!')
  return(sce)
}
# ----scABC--------------
library(scABC)
scABC_5clusters <- function(sce) {
  InSilicoSCABCForeGroundMatrix = counts(sce)
  InSilicoSCABCPeaks = rowData(sce)
  # always filter out low-represented peaks # keep peaks with minimum 2 reads in at least 10 cells
  InSilicoSCABCForeGroundMatrixPeaksFiltered = filterPeaks(InSilicoSCABCForeGroundMatrix, InSilicoSCABCPeaks, nreads_thresh = 2, ncells_thresh = 10) 
  # always filter out low-represented samples ()
  # manually set "readsFGthresh" to 1000 (minimum number of reads in each cell across all peaks). BackGround with all elements equal to NAs since BackGround is not available
  # set tp 100 yo keep all barcodes
  InSilicoSCABCForeGroundMatrixSamplesFiltered = filterSamples(ForeGround = InSilicoSCABCForeGroundMatrixPeaksFiltered$ForeGroundMatrix, 
                                                               BackGround = matrix(nrow = dim(InSilicoSCABCForeGroundMatrixPeaksFiltered$ForeGroundMatrix)[1], 
                                                                                   ncol = dim(InSilicoSCABCForeGroundMatrixPeaksFiltered$ForeGroundMatrix)[2]), 
                                                               readsFGthresh = 100) 
  # filtered foreground matrix and peaks
  InSilicoSCABCForeGroundMatrix = InSilicoSCABCForeGroundMatrixSamplesFiltered$ForeGroundMatrix
  InSilicoSCABCPeaks = InSilicoSCABCForeGroundMatrixPeaksFiltered$peaks
  # define weights to compute landmarks
  weights = apply(InSilicoSCABCForeGroundMatrix, 2, mean)
  InSilicoSCABCLandmarksWithoutBam = computeLandmarks(InSilicoSCABCForeGroundMatrix, weights = weights, nCluster = 5, nTop = 5000)
  # true cluster assignments (not necessary)
  InSilicoClusterAssignments = colData(sce)[, c('barcode', 'celltype')]
  InSilicoClusterAssignments = InSilicoClusterAssignments[match(colnames(InSilicoSCABCForeGroundMatrix), InSilicoClusterAssignments[,1]),]
  # clustering results 
  InSilicoLandMarkAssignmentsWithoutBam = assign2landmarks(InSilicoSCABCForeGroundMatrix, InSilicoSCABCLandmarksWithoutBam)
  
  colData(sce)$clustering_res <- factor(InSilicoLandMarkAssignmentsWithoutBam)
  
  return(sce)
  print('scABC done')
}

# ----Signac--------------

library(Signac)
library(Seurat)   
library(EnsDb.Hsapiens.v75)

signac_norm_SVD_SLM_0.8 <- function(sce){
  counts <- counts(sce)
  metadata <- as.data.frame(colData(sce))
  
  pbmc <- CreateSeuratObject(
    counts = counts,
    assay = 'peaks',
    project = 'cellmix_669',
    min.cells = 1,
    meta.data = metadata
  )
  
  pbmc <- RunTFIDF(pbmc)
  pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
  pbmc <- RunSVD(
    object = pbmc,
    assay = 'peaks',
    reduction.key = 'LSI_',
    reduction.name = 'lsi'
  )
  #FeatureScatter(pbmc, 'LSI_1', 'nCount_peaks')
  # lsi 1 removed
  pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
  pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
  pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)
  
  colData(sce)$clustering_res <- pbmc@meta.data$seurat_clusters
  
  print('signac done!')
  return(sce)
}



# plotCluster(clusterResults, type = "t-SNE") 
# DimPlot(object = pbmc, label = TRUE, shape.by = 'celltype', pt.size = 2)
# ----epiConv---
epiconv_c <- function(sce){
  source("/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/sc_ATAC_seq/github_entries/epiConv/source/epiConv_functions.R")
  # extract count matrix and add colnames and rownames
  mat<-counts(sce)
  barcode <- colnames(sce)
  rownames(mat)<-paste(rowData(sce)$chr,":",
                       rowData(sce)$bp1,"-",
                       rowData(sce)$bp2,sep="")
  
  lib_size<-lib.estimate(mat) # binarised colsum
  freq<-freq.estimate(mat) # rowsum/number of barcodes
  
  # no filter used, create epiconv object
  res_epiConv<-create.epiconv(meta.features=data.frame(barcode=barcode,
                                                       lib_size=lib_size)) # add metadata
  res_epiConv<-add.mat(obj=res_epiConv,x=mat[freq!=0,],name="peak")
  
  # filter libsize > 1000
  #freq<-freq.estimate(mat[,lib_size>1000])
  #res_epiConv<-create.epiconv(meta.features=data.frame(barcode=barcode[lib_size>1000],
  #                                                     lib_size=lib_size[lib_size>1000]))
  #res_epiConv<-add.mat(obj=res_epiConv,x=mat[freq!=0,lib_size>1000],name="peak")
  
  # tfidf normlisation
  mat<-tfidf.norm(mat=res_epiConv@mat[["peak"]],lib_size=res_epiConv$lib_size)
  infv<-inf.estimate(mat[,sample(1:ncol(mat),size=500)],
                     sample_size=0.125,nsim=30) # learn a small value to replace infinite value in the analysis below
  
  # similarities between single cells
  sample_size<-floor(nrow(mat)/8)
  nsim<-30
  sample_matrix<-lapply(1:nsim,function(x) sample(1:nrow(mat),size=sample_size))
  
  Smat<-matrix(0,res_epiConv@ncell,res_epiConv@ncell)
  for(i in sample_matrix){
    Smat<-Smat+epiConv.matrix(mat=mat[i,],inf_replace=infv)
  }
  Smat<-Smat/nsim
  res_epiConv<-add.similarity(obj=res_epiConv,x=Smat,name="sampl_simp")
  #denoise the similarities between single cells
  Smat<-batch.blur(Smat=res_epiConv[["sampl_simp"]],
                   batch=NULL,
                   knn=50)
  res_epiConv<-add.similarity(res_epiConv,x=Smat,name="sampl_simp_denoise")
  #umap to learn the low-dimensional embedding of the data
  umap_settings<-umap::umap.defaults
  umap_settings$input<-"dist"
  umap_settings$n_components<-2
  umap_res<-umap::umap(max(Smat)-Smat,config=umap_settings)$layout
  res_epiConv<-add.embedding(obj=res_epiConv,x=umap_res,name="sampl_simp_denoise")
  
  library(densityClust)
  ncluster<-5
  dclust_obj<-densityClust(res_epiConv@embedding[["sampl_simp_denoise"]],gaussian=T)
  rho_cut<-quantile(dclust_obj$rho,0.5)
  delta_cut<-sort(dclust_obj$delta[dclust_obj$rho>=rho_cut],decreasing=T)[ncluster+1]
  clust<-findClusters(dclust_obj,rho=rho_cut,delta=delta_cut)$clusters
  
  colData(sce)$clustering_res <- clust
  
  print('epiconv done!')
  return(sce)
}

#######################


lib90_sce_669_dgT <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/sc_ATAC_seq/data/lib90_sce_669_dgT.rds")
lib90_sce_libsize1000_698_dgT <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/sc_ATAC_seq/data/lib90_sce_libsize1000_698_dgT.rds")
lib90_sce_libsize100_893_dgT <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/sc_ATAC_seq/data/lib90_sce_libsize100_893_dgT.rds")
lib90_sce_libsize10_1259_dgT <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/sc_ATAC_seq/data/lib90_sce_libsize10_1259_dgT.rds")
lib90_sce_demuxlet_1283_dgT <- readRDS("/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/sc_ATAC_seq/data/lib90_sce_demuxlet_1283_dgT.rds")


# ----dataset---------

dataset <- list(
  cm_clean_669 = lib90_sce_669_dgT,
  cm_1000_698 = lib90_sce_libsize1000_698_dgT,
  cm_100_893 = lib90_sce_libsize100_893_dgT,
  cm_10_1259 = lib90_sce_libsize10_1259_dgT,
  cm_all = lib90_sce_demuxlet_1283_dgT
)

# ----clustering method-------

clustering_method <- list(
  cicero_log_LSI = cicero_log_LSI_UMAP,
  cicero_size_LSI = cicero_size_LSI_UMAP,
  cicero_log_PCA = cicero_log_PCA_UMAP,
  cicero_size_PCA = cicero_size_PCA_UMAP,
  destin = destin_elbow_2to20,
  scABC = scABC_5clusters,
  signac = signac_norm_SVD_SLM_0.8,
  epiconv = epiconv_c
)

res_c <- dataset %>%
  time_methods(clustering_method)

itself <- function(sce){
  return(sce)
  }
n <- list(n=itself)
res_n <- dataset %>% apply_methods(n)



saveRDS(res_c, "/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/sc_ATAC_seq/data/clustering_result.rds")
res_c %>%
  unpack_timing()
# ----apply metrics for evaluation------

library(mclust)
ARI_matric = function(sce){
  if(!("clustering_res" %in% colnames(colData(sce)))){
    return(NA)
  }
  if ("group" %in% colnames(colData(sce))){
    ari_val = adjustedRandIndex(sce$group, sce$clustering_res)
  }else{
    ari_val = adjustedRandIndex(sce$celltype, sce$clustering_res)
  }
  return(ari_val)
}

clustering_evaluation <- list(
  ARI=ARI_matric
)

res_ARI <- res_c %>%
  strip_timing() %>%
  apply_methods(clustering_evaluation)

# ----plot----
#ggplot(data=res_ARI,aes(x=clustering_method,y=result,col=impute_method))+geom_boxplot()+theme_bw()
ggplot(data=res_ARI,aes(x=clustering_method,y=result)) +
  geom_bar(stat = 'identity', aes(fill = clustering_method)) +
  geom_text(aes(label = round(result, 3)), vjust = 0) +
  theme_bw()

ggplot(data=res_ARI,aes(x=clustering_method,y=result,col=dataset))+geom_boxplot()+theme_bw()
#ggplot(data=res_ARI,aes(x=clustering_method,y=result,col=norm_method))+geom_boxplot()+theme_bw()





#cat(paste(format(Sys.time(), "%a %b %d %X %Y"), "SC_all Done! saved the result to file.\n"), file = log_file, append = TRUE)

