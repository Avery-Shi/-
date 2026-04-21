# ==============================================================================
# 酵母菌S288C转录组 - 表达矩阵构建
# 实验组：acid | 对照组：blank
# ==============================================================================

rm(list = ls())

# ===================== 处理 实验组（acid） =====================
setwd('C:/Users/shi60/Desktop')

if(!dir.exists('./raw_data')){
  dir.create('./raw_data')
}
setwd('./raw_data/')

# 导入count矩阵
train_data <- read.table('./acid.count', sep = "\t", header=TRUE, col.names = c("ID", "Count"))

# 导入GTF注释
library(rtracklayer)
library(dplyr)
library(tibble)

gene_annotation <- import('./genomic.gtf')
gene_annotation <- as.data.frame(gene_annotation)
gene_annotation <- gene_annotation[gene_annotation$type == 'gene', ]
gene_annotation <- gene_annotation[gene_annotation$gene_biotype == 'protein_coding', ]

# 提取基因ID与symbol
gene_symbol <- data.frame(gene_annotation$gene_id, gene_annotation$gene)
colnames(gene_symbol) <- c("ID", "symbol")

# 合并count与symbol
data_count <- inner_join(train_data, gene_symbol, by = "ID")
data_count <- dplyr::select(data_count, -ID)
data_count <- dplyr::select(data_count, symbol, everything())
data_count <- arrange(data_count, desc(Count))
data_count <- distinct(data_count, symbol, .keep_all = TRUE)
data_count <- as.data.frame(data_count)
data_count <- na.omit(data_count)
rownames(data_count) <- NULL
data_count <- column_to_rownames(data_count, var = "symbol")

# 输出 count 文件（保留你原来的文件名）
write.csv(data_count, file = 'dat.count.csv', row.names = TRUE)

# 提取基因长度
gene_length <- dplyr::select(gene_annotation, gene, width)
colnames(gene_length) <- c('symbol', 'length')
gene_length <- gene_length[gene_length$symbol %in% rownames(data_count), ]
gene_length <- gene_length[!duplicated(gene_length$symbol), ]

# FPKM计算
countToFpkm <- function(counts, effLen){
  N <- sum(counts)
  exp(log(counts) + log(1e9) - log(effLen) - log(N))
}

dat <- data_count %>% filter(rownames(.) %in% gene_length$symbol)
rownames(dat) <- rownames(data_count)[rownames(data_count) %in% gene_length$symbol]
dat <- as.data.frame(lapply(dat, as.numeric))
rownames(dat) <- rownames(data_count)

expr_fpkm <- apply(dat, 2, countToFpkm, effLen = gene_length$length)

# 自动Log2判断
qx <- as.numeric(quantile(expr_fpkm, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=TRUE))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { 
  expr_fpkm[which(expr_fpkm <= 0)] <- NaN
  expr_fpkm <- log2(expr_fpkm + 1)
  print("acid log2 transform finished")
} else {
  print("acid log2 not needed")
}

# 输出 FPKM（保留原名）
expr_fpkm_acid <- data.frame(expr_fpkm)
write.csv(expr_fpkm_acid, file = 'acid.fpkm.csv')

# 输出 NCBI_count.csv（保留原名）
gene_symbol <- data.frame(gene_id = gene_annotation$gene_id, symbol = gene_annotation$gene)
data_count <- train_data
data_count <- inner_join(data_count, gene_symbol, by = c("ID" = "gene_id"))
result_data <- dplyr::select(data_count, ID, symbol, Count)
colnames(result_data) <- c("NCBI_id", "gene", "Count")
write.csv(result_data, file = 'NCBI_count.csv', row.names = FALSE)

# ===================== 处理 对照组（blank） =====================
setwd('C:/Users/shi60/Desktop')

if(!dir.exists('./raw_data_2')){
  dir.create('./raw_data_2')
}
setwd('./raw_data_2/')

# 导入count
train_data <- read.table('./blank.count', sep = "\t", header=TRUE, col.names = c("ID", "Count"))

# GTF注释
gene_annotation <- import('./genomic.gtf')
gene_annotation <- as.data.frame(gene_annotation)
gene_annotation <- gene_annotation[gene_annotation$type == 'gene', ]
gene_annotation <- gene_annotation[gene_annotation$gene_biotype == 'protein_coding', ]

gene_symbol <- data.frame(gene_annotation$gene_id, gene_annotation$gene)
colnames(gene_symbol) <- c("ID", "symbol")

# 合并与整理
data_count <- inner_join(train_data, gene_symbol, by = "ID")
data_count <- dplyr::select(data_count, -ID)
data_count <- dplyr::select(data_count, symbol, everything())
data_count <- arrange(data_count, desc(Count))
data_count <- distinct(data_count, symbol, .keep_all = TRUE)
data_count <- as.data.frame(data_count)
data_count <- na.omit(data_count)
rownames(data_count) <- NULL
data_count <- column_to_rownames(data_count, var = "symbol")

# 输出 count（保留原名）
write.csv(data_count, file = 'dat.count.csv', row.names = TRUE)

# 基因长度
gene_length <- dplyr::select(gene_annotation, gene, width)
colnames(gene_length) <- c('symbol', 'length')
gene_length <- gene_length[gene_length$symbol %in% rownames(data_count), ]
gene_length <- gene_length[!duplicated(gene_length$symbol), ]

# FPKM
countToFpkm <- function(counts, effLen){
  N <- sum(counts)
  exp(log(counts) + log(1e9) - log(effLen) - log(N))
}

dat <- data_count %>% filter(rownames(.) %in% gene_length$symbol)
rownames(dat) <- rownames(data_count)[rownames(data_count) %in% gene_length$symbol]
dat <- as.data.frame(lapply(dat, as.numeric))
rownames(dat) <- rownames(data_count)

expr_fpkm <- apply(dat, 2, countToFpkm, effLen = gene_length$length)

# Log2
qx <- as.numeric(quantile(expr_fpkm, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=TRUE))
LogC <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { 
  expr_fpkm[which(expr_fpkm <= 0)] <- NaN
  expr_fpkm <- log2(expr_fpkm + 1)
  print("blank log2 transform finished")
} else {
  print("blank log2 not needed")
}

# 输出 FPKM（保留原名）
expr_fpkm_blank <- data.frame(expr_fpkm)
write.csv(expr_fpkm_blank, file = 'blank.fpkm.csv')

# 输出 NCBI_2_count.csv（保留你最终文件名）
gene_symbol <- data.frame(gene_id = gene_annotation$gene_id, symbol = gene_annotation$gene)
data_count <- train_data
data_count <- inner_join(data_count, gene_symbol, by = c("ID" = "gene_id"))
result_data <- dplyr::select(data_count, ID, symbol, Count)
colnames(result_data) <- c("NCBI_id", "gene", "Count")
write.csv(result_data, file = 'NCBI_2_count.csv', row.names = FALSE)

# ===================== 完成 =====================
print("所有表达矩阵已生成，全部使用你原来的文件名！")