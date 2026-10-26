###############################################################################
# 04  CAF trajectory (Monocle3)   [Methods 2.5 / Results 3.4]
#
###############################################################################
suppressPackageStartupMessages({library(Seurat); library(monocle3); library(dplyr)})
set.seed(42)
OUT <- "output"
dir.create(file.path(OUT, "tables"), showWarnings = FALSE, recursive = TRUE)

fib <- readRDS("CAF_subtyped_object.rds")
DefaultAssay(fib) <- "RNA"
cat("CAF cells:", ncol(fib), "| old subtypes:",
    paste(unique(fib$CAF_subtype), collapse = ", "), "\n")

expr_mat  <- GetAssayData(fib, assay = "RNA", layer = "counts")
cell_meta <- data.frame(cell = colnames(fib), fib@meta.data, row.names = 1)
gene_meta <- data.frame(gene = rownames(expr_mat), row.names = rownames(expr_mat))

cds <- new_cell_data_set(expr_mat, cell_metadata = cell_meta, gene_metadata = gene_meta)
cds <- preprocess_cds(cds, num_dim = 30, method = "PCA")   # 30 PCs, not the default 50
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")
cds <- cluster_cells(cds, resolution = 1e-5, cluster_method = "leiden")
cds <- learn_graph(cds, use_partition = FALSE,
                   learn_graph_control = list(ncenter = 200))

# ---- root: cluster with the highest resting_CAF proportion ------------------
cd <- as.data.frame(colData(cds))
rest_prop <- tapply(cd$CAF_subtype == "resting_CAF", cd$cluster, mean)
root_cluster <- names(which.max(rest_prop))
cat(sprintf("root cluster = %s (resting_CAF %.1f%%)\n", root_cluster, 100 * max(rest_prop)))

clust_cells <- colnames(cds)[cd$cluster == root_cluster]
clust_cells <- clust_cells[order(match(clust_cells, colnames(cds)))]
cds <- order_cells(cds, reduction_method = "UMAP",
                   root_cells = clust_cells[1:min(50, length(clust_cells))])

pst <- pseudotime(cds)
res <- data.frame(cell = names(pst), pseudotime = pst,
                  subtype = cd$CAF_subtype, partition = cd$partition,
                  cluster = cd$cluster, group = cd$group)
write.csv(res, file.path(OUT, "tables", "pseudotime.csv"), row.names = FALSE)

# ---- terminal branch composition -------------------------------------------
g   <- principal_graph(cds)[["UMAP"]]
deg <- igraph::degree(g)
ends <- names(deg)[deg == 1]
cat("terminal nodes:", paste(ends, collapse = ", "), "\n")

node_pos <- t(principal_graph(cds)[["UMAP"]][ends, , drop = FALSE])
nodes <- data.frame(node = ends, UMAP1 = node_pos[, 1], UMAP2 = node_pos[, 2])
emb <- reducedDims(cds)[["UMAP"]]
nearest <- apply(emb, 1, function(p) {
  d <- sqrt((nodes$UMAP1 - p[1])^2 + (nodes$UMAP2 - p[2])^2)
  nodes$node[which.min(d)]
})
tb <- as.data.frame.matrix(table(res$subtype, nearest))
tb$subtype <- rownames(tb)
write.csv(tb, file.path(OUT, "tables", "terminal_branch_composition.csv"), row.names = FALSE)
print(tb)

# ---- descriptive surrogates ------------------------------------------------
fib <- AddModuleScore(fib,
  features = list(intersect(c(cc.genes.updated.2019$s.genes,
                              cc.genes.updated.2019$g2m.genes), rownames(fib))),
  name = "Proliferation")
sur <- fib@meta.data %>% group_by(CAF_subtype) %>%
  summarise(n = n(), proliferation = mean(Proliferation1), nFeature = mean(nFeature_RNA))
write.csv(sur, file.path(OUT, "tables", "surrogate_metrics.csv"), row.names = FALSE)
print(as.data.frame(sur))
