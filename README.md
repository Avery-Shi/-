## transcriptome-r-analysis
本项目为转录组数据的完整生物信息学分析流程，包含从数据预处理、差异基因筛选到功能富集的全流程代码与结果，实验组使用0.5%乙酸胁迫的酵母菌转录组。

##研究对象与胁迫条件
- Organism: Saccharomyces cerevisiae S288C (budding yeast, reference strain)
- Stress: Acetic acid stress
  
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
├── data/                 # 数据文件
│   ├── rawdata/          # 处理好的表达矩阵与计算好的FPKM值
│   ├── filtered_data/    # 合并后的count矩阵和用于差异分析的表格
│   └── gene_id_map/      # Name_ID 映射文件
│
├── scripts/              # R 代码
│   ├── 01_expression_quant.R    # 表达定量
│   ├── 02_diff_analysis.R       # 差异基因分析
│   └── 03_enrichment_analysis.R # KEGG/GO富集
│   
└── results/              # 分析结果
    ├── filter/                # 数据过滤与质控
    ├── expression/            # 表达定量
    ├── different/             # 差异基因分析
    ├── go_kegg/               # GO与KEGG富集分析
    ├── metabolic_pathway/     # 重要代谢通路
    └── final_results/         # 酵母菌相应胁迫的机制解析
