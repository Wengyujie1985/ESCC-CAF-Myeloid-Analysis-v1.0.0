###############################################################################
# 07  HLA-high ECM-CAF: threshold, identity and apCAF decomposition
#                                                            [Methods 2.8]
#
###############################################################################
suppressPackageStartupMessages({library(Seurat); library(dplyr)})

OUT <- "output"
dir.create(file.path(OUT, "tables"), showWarnings = FALSE, recursive = TRUE)

fib <- readRDS("CAF_4subtype_object.rds")
ecm <- subset(fib, cells = colnames(fib)[fib$CAF_subtype == "ECM_CAF"])
DefaultAssay(ecm) <- "RNA"
dat <- GetAssayData(ecm, assay = "RNA", layer = "data")
cat("ECM_CAF:", ncol(ecm), "\n")                         # 383

# ---- 1. threshold -----------------------------------------------------------
hla  <- as.numeric(dat["HLA-DRA", ])
thr  <- 2 * mean(hla)
ecm$HLA_class <- ifelse(hla >= thr, "HLA_high", "HLA_low")
ecm$HLA_class <- factor(ecm$HLA_class, levels = c("HLA_high", "HLA_low"))
cat(sprintf("threshold (log1p CP10k) = %.4f | back-transformed = %.4f CP10k\n",
            thr, expm1(thr)))
print(table(ecm$HLA_class))
cat(sprintf("zero-expressing cells: %.1f%%\n", 100 * mean(hla == 0)))

# ---- 2. identity: detection rates -------------------------------------------
genes <- c("COL1A1","DCN","LUM","PTPRC","CD68","LYZ","CD1C","CLEC10A")
det <- sapply(intersect(genes, rownames(dat)), function(g) {
  as.numeric(GetAssayData(ecm, assay = "RNA", layer = "counts")[g, ] > 0)
})
id_tab <- data.frame(gene = colnames(det),
                     pct_high = round(100 * colMeans(det[ecm$HLA_class == "HLA_high", , drop = FALSE]), 1),
                     pct_low  = round(100 * colMeans(det[ecm$HLA_class == "HLA_low",  , drop = FALSE]), 1))
print(id_tab)
write.csv(id_tab, file.path(OUT, "tables", "HLA_identity_detection.csv"), row.names = FALSE)

n_immune <- rowSums(det[, intersect(c("CD68","LYZ","CLEC10A","CD1C"), colnames(det)), drop = FALSE])
cat(sprintf(">=2 immune markers: HLA-high %.1f%% vs HLA-low %.1f%%\n",
            100 * mean(n_immune[ecm$HLA_class == "HLA_high"] >= 2),
            100 * mean(n_immune[ecm$HLA_class == "HLA_low"]  >= 2)))

# ---- 3. depth and doublet context -------------------------------------------
cat(sprintf("median nCount: HLA-high %.0f vs HLA-low %.0f\n",
            median(ecm$nCount_RNA[ecm$HLA_class == "HLA_high"]),
            median(ecm$nCount_RNA[ecm$HLA_class == "HLA_low"])))

# ---- 4. apCAF score decomposition -------------------------------------------
mhc2   <- c("HLA-DPA1","HLA-DPB1","HLA-DRB1","CD74")
costim <- c("CD40","CD80","CD86")
ecm$apCAF_full   <- colMeans(dat[intersect(c(mhc2, costim), rownames(dat)), , drop = FALSE])
ecm$apCAF_mhc2   <- colMeans(dat[intersect(mhc2, rownames(dat)), , drop = FALSE])
ecm$apCAF_costim <- colMeans(dat[intersect(costim, rownames(dat)), , drop = FALSE])

res <- lapply(c("apCAF_full","apCAF_mhc2","apCAF_costim"), function(v) {
  a <- ecm[[v]][ecm$HLA_class == "HLA_high"]; b <- ecm[[v]][ecm$HLA_class == "HLA_low"]
  data.frame(component = v, mean_high = mean(a), mean_low = mean(b),
             p = wilcox.test(a, b)$p.value)
}) %>% bind_rows()
print(res, digits = 4)
write.csv(res, file.path(OUT, "tables", "apCAF_decomposition.csv"), row.names = FALSE)

# ---- 5. sensitivity: >=2 UMI floor ------------------------------------------
umi <- as.numeric(GetAssayData(ecm, assay = "RNA", layer = "counts")["HLA-DRA", ])
cat(sprintf("HLA-high cells with exactly 1 HLA-DRA UMI: %d/%d\n",
            sum(umi[ecm$HLA_class == "HLA_high"] == 1),
            sum(ecm$HLA_class == "HLA_high")))
write.csv(data.frame(cell = colnames(ecm), HLA_class = ecm$HLA_class,
                     HLA_DRA_norm = hla, HLA_DRA_umi = umi, nCount = ecm$nCount_RNA),
          file.path(OUT, "tables", "HLA_percell.csv"), row.names = FALSE)

saveRDS(ecm, file.path(OUT, "tables", "ecm_hla_object.rds"))
