###############################################################################
# 06  CellChat: ECM_CAF -> Myeloid communication   [Methods 2.7]
#
###############################################################################
suppressPackageStartupMessages({
  library(Seurat); library(CellChat); library(dplyr)
})
set.seed(42)
OUT <- "output"
dir.create(file.path(OUT, "tables"), showWarnings = FALSE, recursive = TRUE)

# ---- 1. load and define the 8 groups ---------------------------------------
obj <- readRDS("ESCC_Seurat_final.rds")
if ("RNA" %in% Assays(obj)) DefaultAssay(obj) <- "RNA"

caf_obj   <- readRDS("CAF_4subtype_object.rds")
ecm_cells <- colnames(caf_obj)[caf_obj$CAF_v2 == "ECM_CAF"]
cat("ECM_CAF barcodes:", length(ecm_cells), "\n")
stopifnot(mean(ecm_cells %in% colnames(obj)) >= 0.90)

group_cells <- list(
  "ECM_CAF"     = ecm_cells,
  "Fibroblast"  = setdiff(colnames(obj)[obj$celltype == "Fibroblasts"], ecm_cells),
  "Myeloid"     = colnames(obj)[obj$celltype == "Myeloid"],
  "Epithelial"  = colnames(obj)[obj$celltype == "Epithelial"],
  "T_cells"     = colnames(obj)[obj$celltype == "T_cells"],
  "B_cells"     = colnames(obj)[obj$celltype == "B_cells"],
  "Endothelial" = colnames(obj)[obj$celltype == "Endothelial"],
  "NK_cells"    = colnames(obj)[obj$celltype == "NK_cells"])
print(sapply(group_cells, length))

# no barcode may belong to two groups
for (g in setdiff(names(group_cells), "ECM_CAF")) {
  stopifnot(length(intersect(group_cells[["ECM_CAF"]], group_cells[[g]])) == 0)
}
cells_use <- unlist(group_cells, use.names = FALSE)
stopifnot(!anyDuplicated(cells_use))

# ---- 2. matrix + meta, re-indexed by name ----------------------------------
sub_obj <- subset(obj, cells = cells_use); rm(obj); gc()

data_input <- as.matrix(GetAssayData(sub_obj, assay = "RNA", layer = "data")[, cells_use])
stopifnot(identical(colnames(data_input), cells_use))          # gate 1

meta_df <- data.frame(cell  = cells_use,
                      group = rep(names(group_cells), times = lengths(group_cells)),
                      stringsAsFactors = FALSE)
rownames(meta_df) <- meta_df$cell
stopifnot(identical(rownames(meta_df), colnames(data_input)))  # gate 2

cellchat <- createCellChat(object = data_input, meta = meta_df, group.by = "group")
cellchat@DB <- CellChatDB.human

# gate 3: label -> cell mapping must be the intended one
obj_cols <- colnames(cellchat@data)
for (g in names(group_cells)) {
  idx <- which(cellchat@idents == g)
  if (!identical(obj_cols[idx], group_cells[[g]])) stop("GROUP MISMATCH: ", g)
}
cat("alignment gates: PASS\n")

# ---- 3. marker QC: ECM_CAF must look like fibroblasts ----------------------
tri <- function(x) mean(quantile(x, c(.25, .5, .5, .75)))
qc_genes <- intersect(c("COL1A1","COL1A2","DCN","LUM","EPCAM","PTPRC","CD3D","CD68","PECAM1"),
                      rownames(cellchat@data))
qc <- sapply(levels(cellchat@idents), function(g)
  sapply(qc_genes, function(gene) tri(as.numeric(cellchat@data[gene, cellchat@idents == g]))))
write.csv(round(qc, 3), file.path(OUT, "tables", "QC_marker_trimean_by_group.csv"))
print(round(qc, 2))
stopifnot(qc["COL1A1","ECM_CAF"] > 1.5, qc["DCN","ECM_CAF"] > 1.5,
          qc["EPCAM","ECM_CAF"] < 0.5, qc["PTPRC","ECM_CAF"] < 0.8)

# ---- 4. CellChat ------------------------------------------------------------
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = FALSE)   # nboot 100, trim 0.1
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file.path(OUT, "tables", "cellchat_all_groups.rds"))

# ---- 5. ECM_CAF -> Myeloid pairs -------------------------------------------
src <- which(levels(cellchat@idents) == "ECM_CAF")
tgt <- which(levels(cellchat@idents) == "Myeloid")
p <- cellchat@net$prob[src, tgt, ]; q <- cellchat@net$pval[src, tgt, ]
pairs <- data.frame(LR = names(p), prob = as.numeric(p), pval = as.numeric(q))
write.csv(pairs, file.path(OUT, "tables", "ECM_CAF_to_Myeloid_ALL_pairs.csv"), row.names = FALSE)
sig <- subset(pairs, pval < 0.05)
cat("significant LR pairs:", nrow(sig), "\n")           # 312
cat("top pair:", sig$LR[which.max(sig$prob)], max(sig$prob), "\n")

# ---- 6. outgoing / incoming strength ---------------------------------------
strength <- data.frame(celltype = levels(cellchat@idents),
                       outgoing = rowSums(cellchat@net$weight, dims = 2),
                       incoming = colSums(cellchat@net$weight, dims = 2))
write.csv(strength, file.path(OUT, "tables", "signal_role_strength.csv"), row.names = FALSE)
print(strength[order(-strength$outgoing), ])
