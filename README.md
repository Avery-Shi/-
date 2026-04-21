# transcriptome-r-analysis
本项目为转录组数据的完整生物信息学分析流程，包含从数据预处理、差异基因筛选到功能富集的全流程代码与结果，实验组使用0.5%乙酸胁迫的酵母菌转录组。

## 研究内容
- 原始测序数据质控与过滤（FastQC、trim_galore）
- 序列比对与表达定量（HISAT2、featureCounts）
- 差异表达基因分析（edgeR，无重复设计）
- GO/KEGG功能富集分析（clusterProfiler）
- 可视化：火山图、热图、MA图、通路图

## 软件环境
- Linux：Ubuntu 20.04
- R：4.5.0
- 核心工具：FastQC、trim_galore、HISAT2、featureCounts、edgeR、clusterProfiler

## 项目结构
transcriptome-r-analysis/
├── data/                 # 数据文件（只放小文件）
│   ├── rawdata/          # 原始数据（只放处理后的矩阵，不放fq）
│   ├── filtered_data/    # 数据过滤结果
│   └── gene_id_map/      # Name_ID 映射文件
├── scripts/              # 你的 R 代码（关键！）
│   ├── 01_expression_quant.R    # 表达定量
│   ├── 02_diff_analysis.R       # 差异基因分析
│   ├── 03_enrichment_analysis.R # KEGG/GO富集
│   └── utils.R                  # 通用函数
├── results/              # 分析结果（图片、表格）
    ├── degs/             # 差异基因列表、火山图
    ├── kegg/             # KEGG富集结果
    └── plots/            # 热图、PCA图等可视化
