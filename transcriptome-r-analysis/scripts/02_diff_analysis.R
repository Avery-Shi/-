rm(list = ls())
setwd('C:/Users/shi60/Desktop/01_DEGs/')

# 安装包（仅首次运行）
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR", "EnhancedVolcano", "pheatmap", "ggplot2"), force=TRUE)

# 加载包
library(edgeR)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)

# 读入 count 矩阵
data <- read.csv("NCBI_2_count.csv", row.names = 1)
data <- na.omit(data)
count_data <- data

# 过滤全0基因
keep <- rowSums(count_data) > 0
count_data <- count_data[keep,]

# 分组信息（control + experiment）
group <- factor(c("Control", "Experiment"))

# 创建 DGEList 对象
dge <- DGEList(counts = count_data, group = group)

# 过滤低表达基因
keep <- rowSums(cpm(dge) > 1) >= 1
dge <- dge[keep, , keep.lib.sizes=FALSE]

# TMM 标准化
dge <- calcNormFactors(dge)

# 设置离散度
dge$common.dispersion <- 0.01

# 差异表达分析
et <- exactTest(dge)
results <- topTags(et, n = Inf)
results_df <- as.data.frame(results$table)

# 标记显著 DEGs (FDR < 0.05 & |logFC| > 1)
results_df$significant <- ifelse(results_df$FDR < 0.05 & abs(results_df$logFC) > 1, "YES", "NO")
degs <- subset(results_df, significant == "YES")

# 统计
up <- subset(degs, logFC > 1)
down <- subset(degs, logFC < -1)
cat("上调基因：", nrow(up), "\n")
cat("下调基因：", nrow(down), "\n")
cat("总差异基因：", nrow(degs), "\n")

# 保存结果
write.csv(results_df, "all_diff_results.csv", row.names = TRUE)
write.csv(degs, "DEGs_significant.csv", row.names = TRUE)
write.csv(up, "DEGs_up.csv", row.names = TRUE)
write.csv(down, "DEGs_down.csv", row.names = TRUE)

# ---------------------- 可视化 ----------------------
# MA 图
plotMD(results_df, main = "MA Plot", ylim = c(-5, 5))

# 火山图
EnhancedVolcano(results_df,
                lab = rownames(results_df),
                x = "logFC",
                y = "FDR",
                pCutoff = 0.05,
                FCcutoff = 1,
                title = "Volcano Plot",
                pointSize = 2,
                labSize = 3)

# 热图（差异基因表达量）
sig_genes <- rownames(degs)
exp_matrix <- count_data[sig_genes,]
pheatmap(exp_matrix, scale = "row", show_rownames = F, main = "DEGs Expression Heatmap")