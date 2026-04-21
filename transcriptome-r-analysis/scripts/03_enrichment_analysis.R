
rm(list = ls())
setwd('C:/Users/shi60/Desktop/富集分析')

# 安装包（仅首次运行）
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Sc.sgd.db", "enrichplot", "pathview", "ggplot2"), force=TRUE)

# 加载包
library(clusterProfiler)
library(org.Sc.sgd.db)
library(enrichplot)
library(ggplot2)
library(pathview)

# ==============================================================================
# 1. 读取差异基因
# ==============================================================================
up_genes <- read.csv("edgeR_upregulated_genes.csv")
down_genes <- read.csv("edgeR_downregulated_genes.csv")

gene_up <- up_genes$gene_id
gene_down <- down_genes$gene_id

# ==============================================================================
# 2. GO 富集分析（BP / CC / MF）
# ==============================================================================
go_up <- enrichGO(
  gene = gene_up,
  OrgDb = org.Sc.sgd.db,
  keyType = "ENSEMBL",
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

go_down <- enrichGO(
  gene = gene_down,
  OrgDb = org.Sc.sgd.db,
  keyType = "ENSEMBL",
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

# 保存 GO 结果
write.csv(as.data.frame(go_up), "GO_upregulated.csv", row.names = FALSE)
write.csv(as.data.frame(go_down), "GO_downregulated.csv", row.names = FALSE)

# 可视化 GO
dotplot(go_up, showCategory=10, title="GO Upregulated")
dotplot(go_down, showCategory=10, title="GO Downregulated")

# ==============================================================================
# 3. KEGG 富集分析
# ==============================================================================
kegg_up <- enrichKEGG(
  gene = gene_up,
  organism = "sce",
  pvalueCutoff = 0.05
)

kegg_down <- enrichKEGG(
  gene = gene_down,
  organism = "sce",
  pvalueCutoff = 0.05
)

# 保存 KEGG 结果
write.csv(as.data.frame(kegg_up), "KEGG_upregulated.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_down), "KEGG_downregulated.csv", row.names = FALSE)

# 可视化 KEGG
dotplot(kegg_up, showCategory=10, title="KEGG Upregulated")
dotplot(kegg_down, showCategory=10, title="KEGG Downregulated")
barplot(kegg_up, showCategory=15)
barplot(kegg_down, showCategory=15)

# ==============================================================================
# 4. 合并所有差异基因做整体富集
# ==============================================================================
merged_data <- read.csv("merged_genes.csv")
gene_all <- merged_data$gene_id

kegg_all <- enrichKEGG(gene = gene_all, organism = "sce", pvalueCutoff = 0.05)
write.csv(as.data.frame(kegg_all), "KEGG_all_genes.csv", row.names = FALSE)

dotplot(kegg_all, showCategory=15, title="KEGG All DEGs")
barplot(kegg_all, showCategory=20)

# ==============================================================================
# 5. 通路图绘制（pathview）
# ==============================================================================
if (!is.null(kegg_all) && nrow(as.data.frame(kegg_all)) > 0) {
  pathway_id <- kegg_all$ID[1]
  pathview(gene.data = gene_all, pathway.id = pathway_id, species = "sce")
}

# ============================= 完成 =============================