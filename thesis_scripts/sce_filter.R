# create sce and keep barcodes in (cellranger, demuxlet) and get stats plot
# create sce using raw_feature matrix

library(DropletUtils)
library(ggplot2)

# Inputs
# lib 90
dir <- "/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/cellmix_ATACseq/outs/raw_peak_bc_matrix"
hdf5 <- "/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/cellmix_ATACseq/outs/filtered_peak_bc_matrix.h5"
barcodes_cr <- "/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/cellmix_ATACseq/outs/filtered_peak_bc_matrix/barcodes.tsv"
celltype <- "/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/sc_ATAC_seq/data/short_read_demuxlet_90/short_90_barcodes_celline.tab"
demux_single <- "/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/sc_ATAC_seq/data/short_read_demuxlet_90/lib90_single.txt"
cr_cluster <- "/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/cellmix_ATACseq/outs/analysis/clustering/graphclust/clusters.csv"
barcode_stats <- "/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/cellmix_ATACseq/outs/singlecell.csv"

createSceATAC <- function(hdf5, dir, barcodes_cr, celltype, barcodes_stats, cr_cluster) {
  
  # read raw matrix in hdf5
  sce <- read10xCounts(hdf5, sample.names = names(hdf5), col.names = TRUE)
  
  # add peak info to rowData and rowRanges
  bedData <-  rtracklayer::import(file.path(dir, "peaks.bed"), format = "BED")
  peakRegion <- read.table(file.path(dir, "peaks.bed"), stringsAsFactors = F)
  peakRegion$V4 <- paste0(peakRegion$V1,'_', peakRegion$V2, '_', peakRegion$V3)
  rowRanges(sce) <- bedData
  rowData(sce) <- data.frame(peaks = peakRegion$V4) 
  
  # union of barcodes in CellrangerATAC and demuxlet, keep these rows/barcodes
  barcodes_cr <- read.delim(barcodes_cr,header = F, stringsAsFactors = F)[ , 1]
  celltype <- read.table(celltype, stringsAsFactors = F, col.names = c("Barcode", "celltype"))[c(-1, -2), ]
  celltype$celltype <- as.factor(celltype$celltype)
  barcodes <- sort(union(barcodes_cr, celltype[,1]),)
  sce <- sce[,colnames(sce) %in% barcodes]
  
  # read barcode statistics into dataframe
  barcode_stats <- read.csv(barcode_stats, header = T, stringsAsFactors = F)
  barcode_stats <- barcode_stats[barcode_stats$barcode %in% barcodes, -9]
  
  # add barcode statistics, celltype, cellranger cluster, fragment ends count to colData
  # barcode stats
  barcode_stats <- merge(colData(sce), barcode_stats, by.x = "Barcode", by.y = "barcode", sort = F)
  barcode_stats$id  <- 1:nrow(barcode_stats)
  # celltype
  barcode_stats <- merge(barcode_stats, celltype, by="Barcode", all = T, sort = F)
  # cellranger cluster
  cr_cluster <- read.csv(cr_cluster, header = T, col.names = c("Barcode", "cr_cluster"), stringsAsFactors = F)
  cr_cluster$cr_cluster <- as.factor(cr_cluster$cr_cluster)
  barcode_stats <- merge(barcode_stats, cr_cluster, by="Barcode", all = T, sort = F)
  # fragment_ends_count
  fragment_ends_count <- colSums(counts(sce))
  fragment_ends_count <- as.data.frame(fragment_ends_count, row.names = names(fragment_ends_count))
  fragment_ends_count$Barcode <- row.names(fragment_ends_count)
  barcode_stats <- merge(barcode_stats, fragment_ends_count, by="Barcode", all = T, sort = F)
  # fragment_ends_count_binary
  fragment_ends_count_binary <- colSums(counts(sce) > 0)
  fragment_ends_count_binary <- as.data.frame(fragment_ends_count_binary, row.names = names(fragment_ends_count_binary))
  fragment_ends_count_binary$Barcode <- row.names(fragment_ends_count_binary)
  barcode_stats <- merge(barcode_stats, fragment_ends_count_binary, by="Barcode", all = T, sort = F)
  # demuxlet_single
  demuxlet_single <- read.table(demux_single, stringsAsFactors = F, header = F, col.names = c("Barcode", "demux_single"))
  demuxlet_single$demux_single <- as.factor(demuxlet_single$demux_single)
  barcode_stats <- merge(barcode_stats, demuxlet_single, by="Barcode", all = T, sort = F)
  # Final
  barcode_stats <- barcode_stats[order(barcode_stats$id), ]
  rownames(barcode_stats) <- barcode_stats$Barcode
  colData(sce) <- barcode_stats
  
  return(sce)
}

sce <- createSceATAC(hdf5, dir, barcodes_cr, celltype, barcodes_stats, cr_cluster)
saveRDS(sce, "/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/sc_ATAC_seq/data/lib90_sce_demux_cr.rds")
# counts(sce)
# colData(sce)
# rowData(sce)

# plot
stats <- data.frame(colData(sce))

# 1. colsum in matrix, i.e. fragment ends in peaks count

ggplot(stats, aes(celltype, log10(fragment_ends_count))) +
  geom_violin(aes(fill = celltype)) +
  geom_jitter(size = 0.5) +
  geom_hline(yintercept = c(log10(1000), log10(5000))) +
  theme(axis.title.y = element_blank()) + 
  ggtitle("log10(No. of fragment-ends in peaks)") # y > 2.1 or log>5

ggplot(stats, aes(cr_cluster, log10(fragment_ends_count))) +
  geom_violin(aes(fill = cr_cluster)) +
  geom_jitter(size = .5, ) +
  geom_hline(yintercept = 2.24)
  theme(axis.title.y = element_blank()) +
  ggtitle("log10(No. of fragment-ends in peaks) - cr_cluster")


ggplot(stats, aes(cr_cluster, log10(fragment_ends_count))) +
  geom_violin(aes(fill = cr_cluster)) +
  geom_jitter(size = 1, aes(color = celltype)) +
  theme(axis.title.y = element_blank()) +
  ggtitle("log10(No. of fragment-ends in peaks) - cr_cluste") # y > 2.1 or log>5

ggplot(stats, aes(celltype, log10(fragment_ends_count))) +
  geom_violin(aes(fill = celltype)) +
  geom_jitter(size = 1, aes(color = factor(cr_cluster, levels = c(5,1,6,4,3,2,NA)))) +
  theme(axis.title.y = element_blank()) +
  labs(title='Log10(Library Size) i.e. colSum of matrix',x='Cell Type',color = "CellRanger", fill='Cell Type') # y > 2.1 or log>5



# 2. fragments overlapping peaks - similar to previous
ggplot(stats, aes(celltype, log10(peak_region_fragments))) +
  geom_violin(aes(fill = celltype)) +
  geom_jitter(size = 0.2) +
  theme(axis.title.y = element_blank()) +
  ggtitle("log10(number of fragments overlapping peaks)") # y > 3 or > 1.7/1.8

ggplot(stats, aes(cr_cluster, log10(peak_region_fragments))) +
  geom_violin(aes(fill = cr_cluster)) +
  geom_jitter(size = 0.2) +
  geom_hline(yintercept = 1.94) +
  theme(axis.title.y = element_blank()) +
  ggtitle("log10(number of fragments overlapping peaks - cr_cluster)")

ggplot(stats, aes(demux_single, log10(peak_region_fragments))) +
  geom_violin(aes(fill = demux_single)) +
  geom_jitter(size = 1.0) +
  theme(axis.title.y = element_blank()) +
  ggtitle("log10(number of fragments overlapping peaks - demux_single)")


# 3. Total non-dup read-pairs passed filters
ggplot(stats, aes(celltype, log10(passed_filters)) ) +
  geom_violin(aes(fill = celltype)) +
  geom_jitter(size = 0.5) +
  theme(axis.title.y = element_blank()) +
  ggtitle("log10(Total Non-dup Read-pairs Passed Filters)")

ggplot(stats, aes(cr_cluster, log10(passed_filters)) ) +
  geom_violin(aes(fill = cr_cluster)) +
  geom_jitter(size = 0.5) +
  theme(axis.title.y = element_blank()) +
  ggtitle("log10(Total Non-dup Read-pairs Passed Filters - cr_cluster)")

ggplot(stats, aes(celltype, log10(passed_filters)) ) +
  geom_violin(aes(fill = celltype)) +
  geom_jitter(size = 1, aes(color=factor(cr_cluster, levels = c(5,1,6,4,3,2,NA)))) +
  theme(axis.title.y = element_blank()) +
  labs(title="log10(Total Non-dup Read-pairs Passed Filters)",x='Cell Type',color = "CellRanger", fill='Cell Type')

# 4. total read-pairs
ggplot(stats, aes(celltype, log10(total)) ) +
  geom_violin(aes(fill = celltype)) +
  geom_jitter(size = 0.5) +
  theme(axis.title.y = element_blank()) +
  ggtitle("log10(Total Read-pairs)")

ggplot(stats, aes(celltype, log10(total)) ) +
  geom_violin(aes(fill = celltype)) +
  geom_jitter(size =1.5, aes(color = factor(cr_cluster, levels = c(5,1,6,4,3,2,NA)))) +
  theme(axis.title.y = element_blank()) +
  labs(title="log10(Total Read-pairs)",x='Cell Type',color = "CellRanger", fill='Cell Type')


ggplot(stats, aes(cr_cluster, log10(total)) ) +
  geom_violin(aes(fill = cr_cluster)) +
  geom_jitter(size = 0.5) +
  theme(axis.title.y = element_blank()) +
  ggtitle("log10(Total Read-pairs - cr_cluster)")

ggplot(stats, aes(demux_single, log10(total)) ) +
  geom_violin(aes(fill = demux_single)) +
  geom_jitter(size = 0.5) +
  theme(axis.title.y = element_blank()) +
  ggtitle("log10(Total Read-pairs - demux_single)")

# 5. % mito

ggplot(stats, aes(celltype, log10(mitochondrial/total))) +
  geom_violin(aes(fill = celltype)) +
  geom_jitter(size = 0.5) + 
  theme(axis.title.y = element_blank()) +
  ggtitle("log10(% Mitochondria and Non-nuclear Read-pairs)")






# filtering

ggplot(stats, aes(celltype, log10(passed_filters)) ) +
  geom_violin(aes(fill = celltype)) +
  geom_jitter(size = 0.5, aes(color = cr_cluster)) +
  geom_hline(yintercept = log10(5000), color = 'red') +
  theme(axis.title.y = element_blank()) +
  ggtitle("log10(Total Non-dup Read-pairs Passed Filters)")

ggplot(stats, aes(celltype, log10(fragment_ends_count_binary)) ) +
  geom_violin(aes(fill = celltype)) +
  geom_jitter(size = 1.0, aes(color = factor(cr_cluster, levels = c(5,1,6,4,3,2,NA)))) +
  geom_hline(yintercept = c(log10(100), log10(1000), log10(10)), color = 'red') +
  theme(axis.title.y = element_blank()) +
  labs(title="log10(Binary Library Size)",x='Cell type',color = "CellRanger", fill='Cell Type')

stats1 <- stats[(log10(stats$passed_filters) > log10(5000) & !is.na(stats$celltype)), ]
sum(log10(stats$passed_filters) > log10(5000) & !is.na(stats$celltype)) # > 669
table(stats1$celltype)
# A549  H1975  H2228  H8383 HCC827 
# 130    212     91     94    142 

# filter by binary lib size > 1000
stats2 <- stats[(stats$fragment_ends_count_binary > 1000 & !is.na(stats$celltype)), ]
sum(stats$fragment_ends_count_binary > 1000 & !is.na(stats$celltype)) # > 698
table(stats2$celltype)
#A549  H1975  H2228  H8383 HCC827 
#132    222    101     92    151 

ggplot(stats1, aes(celltype, log10(fragment_ends_count))) +
  geom_violin(aes(fill = celltype)) +
  geom_jitter(size = 1, aes(color = cr_cluster)) +
  theme(axis.title.y = element_blank()) +
  ggtitle("log10(No. of fragment-ends in peaks)") 

ggplot(stats, aes(celltype, log10(total))) +
  geom_violin(aes(fill = celltype)) +
  geom_jitter(size = 1) +
  theme(axis.title.y = element_blank()) +
  ggtitle("log10(non-dup total read-pairs)") 


sce_clean <- sce[, (colData(sce)$passed_filters>5000 & !is.na(colData(sce)$celltype))]
barcodes_clean <- sce_clean$Barcode
saveRDS(sce_clean, "/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/sc_ATAC_seq/data/lib90_sce_clean_669.rds" )
saveRDS(barcodes_clean,  "/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/sc_ATAC_seq/data/lib90_barcodes_clean_vec_669.rds")
#rm(list = c('hdf5', 'dir', 'celltype'), envir = .GlobalEnv)
# A549  H1975  H2228  H8383 HCC827 
# 130    212     91     94    142 

ggplot(stats, aes(celltype, log10(fragment_ends_count_binary)) ) +
  geom_violin(aes(fill = celltype)) +
  geom_jitter(size = 0.5, aes(color = factor(cr_cluster, levels = c(5,1,6,4,3,2,NA)))) +
  geom_hline(yintercept = c(log10(100), log10(1000), log10(10)), color = 'red') +
  theme(axis.title.y = element_blank()) +
  labs(title="log10(Binary Library Size)",x='Cell type',color = "CellRanger", fill='Cell Type')

sce_libsize1000 <- sce[, (colData(sce)$fragment_ends_count_binary > 1000 & !is.na(colData(sce)$celltype))]
barcodes_libsize1000 <- sce_libsize1000$Barcode
saveRDS(sce_libsize1000, "/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/sc_ATAC_seq/data/lib90_sce_libsize1000_698.rds" )
saveRDS(barcodes_libsize1000,  "/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/sc_ATAC_seq/data/lib90_barcodes_libsize1000_698.rds")
#celltype <- "/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/sc_ATAC_seq/data/short_read_demuxlet_90/lib90_single.txt"
#A549  H1975  H2228  H8383 HCC827 
#132    222    101     92    151 

# filter by binary lib size > 100
stats2 <- stats[(stats$fragment_ends_count_binary > 100 & !is.na(stats$celltype)), ]
sum(stats$fragment_ends_count_binary > 100 & !is.na(stats$celltype)) # > 893
table(stats2$celltype)
#A549  H1975  H2228  H8383 HCC827 
#164    256    128    155    190
sce_libsize100 <- sce[, (colData(sce)$fragment_ends_count_binary > 100 & !is.na(colData(sce)$celltype))]
barcodes_libsize100 <- sce_libsize100$Barcode
saveRDS(sce_libsize100, "/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/sc_ATAC_seq/data/lib90_sce_libsize100_893.rds" )
saveRDS(barcodes_libsize100,  "/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/sc_ATAC_seq/data/lib90_barcodes_libsize100_893.rds")

# filter by binary lib size > 10
stats2 <- stats[(stats$fragment_ends_count_binary > 10 & !is.na(stats$celltype)), ]
sum(stats$fragment_ends_count_binary > 10 & !is.na(stats$celltype)) # > 893
table(stats2$celltype)
#A549  H1975  H2228  H8383 HCC827 
#204    288    162    384    221
sce_libsize10 <- sce[, (colData(sce)$fragment_ends_count_binary > 10 & !is.na(colData(sce)$celltype))]
barcodes_libsize10 <- sce_libsize10$Barcode
saveRDS(sce_libsize10, "/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/sc_ATAC_seq/data/lib90_sce_libsize10_1259.rds" )
saveRDS(barcodes_libsize10,  "/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/sc_ATAC_seq/data/lib90_barcodes_libsize10_1259.rds")

# all with demuxlet results
stats2 <- stats[!is.na(stats$celltype), ]
sum( !is.na(stats$celltype)) # > 893
table(stats2$celltype)
#A549  H1975  H2228  H8383 HCC827 
#207    290    164    397    225    
sce_demux <- sce[, !is.na(colData(sce)$celltype)]
barcodes_demux <- sce_demux$Barcode
saveRDS(sce_demux, "/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/sc_ATAC_seq/data/lib90_sce_demux_1283.rds" )
saveRDS(barcodes_demux,  "/stornext/General/data/user_managed/grpu_mritchie_1/HaoyuYang/sc_ATAC_seq/data/lib90_barcodes_demux_1283.rds")
