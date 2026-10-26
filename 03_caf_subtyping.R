###############################################################################
# 03  CAF subtypes by cluster-level argmax over four programs   [Methods 2.4]
#
###############################################################################
suppressPackageStartupMessages({library(Seurat); library(dplyr)})
set.seed(42)

OUT <- "output"
dir.create(file.path(OUT, "tables"), showWarnings = FALSE, recursive = TRUE)

fib <- readRDS("CAF_subtyped_object.rds")
DefaultAssay(fib) <- "RNA"
fib <- JoinLayers(fib, assay = "RNA")
fib <- NormalizeData(fib, verbose = FALSE)
fib <- FindVariableFeatures(fib, nfeatures = 2000, verbose = FALSE)
fib <- ScaleData(fib, verbose = FALSE)
fib <- RunPCA(fib, npcs = 20, verbose = FALSE)
fib <- RunUMAP(fib, dims = 1:15, verbose = FALSE)
fib <- FindNeighbors(fib, dims = 1:15, verbose = FALSE)
fib <- FindClusters(fib, resolution = 0.5, verbose = FALSE)
cat("CAF cells:", ncol(fib), "| clusters:", nlevels(fib), "\n")     # 2,749

# ---- four literature-based programs (exactly the lists in Supplementary S1) --
caf_programs <- list(
  ECM_CAF     = c("MMP1","MMP11","MMP13","MMP3","COL1A1","COL3A1","COL5A1",
                  "COL10A1","COL11A1","POSTN","TNC","FN1","LOX","LOXL1",
                  "FAP","CTSK","SPP1","COL6A1","COL6A2","COL6A3"),
  myCAF       = c("ACTA2","TAGLN","CNN1","MYH11","MYL9","MYLK","CALD1",
                  "TPM1","TPM2","LMOD1","HOPX","DES","PDGFRB"),
  apCAF       = c("HLA-DRA","CD74","HLA-DPA1","HLA-DPB1","HLA-DQA1",
                  "HLA-DQB1","HLA-DRB1","CIITA","CD40"),
  resting_CAF = c("DCN","LUM","PI16","C7","CFH","CLU","CST3",
                  "FBLN1","FBLN2","COL14A1","CXCL12","PTGDS","DPT"))

expr <- GetAssayData(fib, assay = "RNA", layer = "data")
score_list <- lapply(names(caf_programs), function(p) {
  g <- intersect(caf_programs[[p]], rownames(expr))
  cat(sprintf("  %-12s %2d/%2d genes\n", p, length(g), length(caf_programs[[p]])))
  colMeans(expr[g, , drop = FALSE])
})
score_mat <- do.call(cbind, score_list); colnames(score_mat) <- names(caf_programs)

# ---- argmax at cluster level -------------------------------------------------
clusters <- as.character(Idents(fib))
clust_program <- do.call(rbind, lapply(sort(unique(clusters)), function(cl) {
  colMeans(score_mat[clusters == cl, , drop = FALSE])
}))
rownames(clust_program) <- sort(unique(clusters))
print(round(clust_program, 3))

cluster_subtype <- setNames(apply(clust_program, 1, function(x) names(which.max(x))),
                            rownames(clust_program))
print(cluster_subtype)
fib$CAF_v2 <- cluster_subtype[clusters]

cat("\nsubtype counts:\n"); print(table(fib$CAF_v2))
cat("\nsubtype x group:\n"); print(table(fib$CAF_v2, fib$group))
saveRDS(fib, file.path(OUT, "tables", "CAF_4subtype_object.rds"))

comp <- as.data.frame(table(fib$CAF_v2, fib$group)) %>%
  rename(subtype = Var1, group = Var2, n = Freq) %>%
  group_by(group) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup()
write.csv(comp, file.path(OUT, "tables", "CAF_subtype_composition.csv"), row.names = FALSE)
print(comp)
