# R_reproduce — 论文方法复现代码


- 运行环境：**R 4.4.2**（Windows，ucrt）



## 1. 依赖包（实测版本）

```
Seurat 5.2.1        SingleR 2.8.0      celldex 1.16.0
scDblFinder 1.20.2  DoubletFinder 2.0.6 SingleCellExperiment
monocle3 1.4.26     Scissor 2.0.0      CellChat 1.6.1
lme4 1.1.36         lmerTest 3.2.1     survival 3.8.6   survminer 0.5.0
data.table 1.16.4   dplyr 1.1.4        ggplot2 4.0.2    tidyr 1.3.1
patchwork 1.3.0     readxl 1.4.3       Matrix 1.7.5
```

- `02` 里的 ARI 是**自带函数**，不依赖 `mclust`（本机未安装该包）。
- `10` 需要 `geo_cache/CIBERSORT.R` 与 `geo_cache/LM22.txt`（CIBERSORT 官方脚本与签名矩阵）。

---

## 2. 数据前置

| 数据 | 路径 | 说明 |
|---|---|---|
| 10x 原始矩阵 | `data/{1N,2N,3N,1T,2T,3T}` | GSE196756，3 癌旁 + 3 肿瘤 |
| TCGA-ESCC | `TCGA_ESCC_expr.rds`、`TCGA_ESCC_survival.csv` | n=95，事件 32 |
| GEO 验证集 | `GSE53622_*`、`GSE53624_*` | n=60 / 119 |
| 中间对象 | `output/rdata/`、`output/myeloid_subtypes/tables/myeloid_scored.rds`、`output/caf_v2/...` | 若已存在可跳过 01 |

---

## 3. 复现深度说明（重要）

这些脚本**已通过语法检查**，但**没有全部端到端重跑**：
01 需要 10x 原始数据 + SingleR 参考集（耗时长），含随机过程。


- 脚本的作用是**让方法与代码一一对应、可追溯**，不是声称"重跑必然一字不差"
- 涉及随机性的步骤（UMAP、聚类、Scissor、Monocle3）重跑会得到数值相近但非逐位相同的输出

---

