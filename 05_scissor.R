###############################################################################
# 05  Scissor: phenotype-associated cells from bulk survival   [Methods 2.6]
#
###############################################################################
suppressPackageStartupMessages({
  library(Seurat); library(survival); library(dplyr); library(ggplot2)
})
source("scissor_patched.R")
set.seed(42)

BULK <- "批量数据"
OUT  <- "output"
dir.create(file.path(OUT, "tables"), showWarnings = FALSE, recursive = TRUE)

# ---- data -------------------------------------------------------------------
obj <- readRDS("ESCC_Seurat_light.rds")
DefaultAssay(obj) <- Assays(obj)[1]

caf_obj  <- readRDS("CAF_4subtype_object.rds")
caf_meta <- caf_obj@meta.data[, c("CAF_v2", "group", "sample"), drop = FALSE]
obj_meta <- obj@meta.data
obj_meta$CAF_v2 <- caf_meta$CAF_v2[match(rownames(obj_meta), rownames(caf_meta))]
obj_meta$CAF_v2[is.na(obj_meta$CAF_v2)] <- "Non-CAF"
obj <- AddMetaData(obj, metadata = obj_meta$CAF_v2, col.name = "CAF_v2")

bulk_expr <- as.matrix(readRDS(file.path(BULK, "TCGA_ESCC_expr.rds")))
surv      <- read.csv(file.path(BULK, "TCGA_ESCC_survival.csv"), stringsAsFactors = FALSE)

phenotype <- as.matrix(data.frame(time = surv$OS.time, status = surv$OS.event))
rownames(phenotype) <- surv$sample
bulk_expr <- bulk_expr[, intersect(colnames(bulk_expr), surv$sample), drop = FALSE]
phenotype <- phenotype[colnames(bulk_expr), , drop = FALSE]
cat("bulk:", ncol(bulk_expr), "samples,", sum(phenotype[, "status"]), "events\n")

# ---- cell selection ---------------------------------------------------------
caf_cells_all <- colnames(obj)[obj$CAF_v2 != "Non-CAF"]
noncaf_cells  <- setdiff(colnames(obj), caf_cells_all)
noncaf_sample <- sample(noncaf_cells, min(length(noncaf_cells), 8000))
cells_use     <- c(caf_cells_all, noncaf_sample)
obj_sub       <- subset(obj, cells = cells_use)
cat("cells:", length(cells_use), "(", length(caf_cells_all), "CAF +",
    length(noncaf_sample), "non-CAF )\n")

common   <- intersect(rownames(bulk_expr), rownames(obj))
vargenes <- intersect(VariableFeatures(obj_sub), common)
if (length(vargenes) < 2000) vargenes <- common[1:min(5000, length(common))]
cat("variable genes used:", length(vargenes), "\n")

sc_mat <- as.matrix(GetAssayData(obj_sub, assay = Assays(obj_sub)[1],
                                 layer = "data")[vargenes, , drop = FALSE])

# ---- Scissor ----------------------------------------------------------------
res <- scissor_patched(
  bulk_dataset = as.matrix(bulk_expr),
  sc_dataset   = sc_mat,
  phenotype    = phenotype,
  tag          = NULL,
  alpha        = c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
  cutoff       = 0.2,
  family       = "cox",
  Save_file    = file.path(OUT, "tables", "Scissor_inputs_95.RData"))

saveRDS(res, file.path(OUT, "tables", "Scissor_result_95.rds"))

# Scissor returns the selected cells as two ID vectors - there is NO per-cell
# `res$Scissor` column. Assigning the labels yourself is the only way.
scissor_pos <- res$Scissor_pos
scissor_neg <- res$Scissor_neg
cat(sprintf("Scissor+ %d cells | Scissor- %d cells\n",
            length(scissor_pos), length(scissor_neg)))     # 291 / 61

obj$Scissor <- "Background"
obj$Scissor[colnames(obj) %in% scissor_pos] <- "Scissor+"
obj$Scissor[colnames(obj) %in% scissor_neg] <- "Scissor-"

write.csv(data.frame(cell = colnames(obj), sample = obj$sample, celltype = obj$celltype,
                     caf_v2 = obj$CAF_v2, Scissor = obj$Scissor),
          file.path(OUT, "tables", "scissor_cell_labels.csv"), row.names = FALSE)

# ---- (a) Scissor+ fraction within each CAF subtype --------------------------
caf_cells <- colnames(obj)[obj$CAF_v2 != "Non-CAF"]
frac <- as.data.frame.matrix(table(obj@meta.data[caf_cells, "CAF_v2"],
                                   obj@meta.data[caf_cells, "Scissor"]))
if (!"Scissor+" %in% colnames(frac)) frac[["Scissor+"]] <- 0
frac$Scissor_plus_pct <- frac[["Scissor+"]] / rowSums(frac) * 100
frac$n <- rowSums(frac)
write.csv(frac, file.path(OUT, "tables", "Scissor_fraction_CAF_subtypes_95.csv"))
print(frac)                       # every CAF subtype: 0 Scissor+ (the honest result)

# ---- (b) cell-type composition of the Scissor+ cells ------------------------
comp <- obj@meta.data[obj$Scissor == "Scissor+", ] %>%
  count(celltype, name = "n") %>% mutate(pct = 100 * n / sum(n), run = "95p")
write.csv(comp, file.path(OUT, "tables", "scissor_pos_composition.csv"), row.names = FALSE)
print(comp)                       # Myeloid 237/291 = 81.4%

# ---- (c) signature genes: top-100 enriched in Scissor+ cells ---------------
# enrichment = mean single-cell expression in Scissor+ cells minus the mean over
# ALL modelled cells, then restricted to genes present in the bulk matrix.
pos_avail <- intersect(scissor_pos, colnames(sc_mat))
gene_sig  <- rowMeans(sc_mat[, pos_avail, drop = FALSE]) - rowMeans(sc_mat)
sig_genes <- intersect(names(sort(gene_sig, decreasing = TRUE))[1:100], rownames(bulk_expr))
cat("signature genes:", length(sig_genes), "\n")

# ---- (d) patient score from BULK expression, then KM and Cox ----------------
# NOTE: the patient-level score is the bulk mean of the signature genes, NOT the
# mean of the single-cell values - that is what makes it a bulk phenotype.
sig_score <- colMeans(bulk_expr[sig_genes, , drop = FALSE])
sv <- data.frame(sample  = colnames(bulk_expr),
                 OS.time = phenotype[, "time"],
                 OS.event = phenotype[, "status"],
                 Scissor_score = as.numeric(sig_score))
sv$group <- ifelse(sv$Scissor_score >= median(sv$Scissor_score), "High", "Low")
write.csv(sv, file.path(OUT, "tables", "survival_by_Scissor_95.csv"), row.names = FALSE)

lr <- survdiff(Surv(OS.time, OS.event) ~ group, data = sv)
p_lr <- 1 - pchisq(lr$chisq, 1)
cx <- coxph(Surv(OS.time, OS.event) ~ Scissor_score, data = sv)
s  <- summary(cx)
sd_score <- sd(sv$Scissor_score)
cat(sprintf("log-rank p = %.4f | HR/unit %.3f (%.3f-%.3f), p = %.3g | HR/SD %.3f\n",
            p_lr, exp(coef(cx)), exp(confint(cx))[1], exp(confint(cx))[2],
            s$coefficients[1, "Pr(>|z|)"], exp(coef(cx) * sd_score)))

write.csv(data.frame(n = nrow(sv), events = sum(sv$OS.event),
                     n_scissor_pos = length(scissor_pos), n_sig_genes = length(sig_genes),
                     logrank_p = p_lr, HR_per_unit = exp(coef(cx)),
                     HR_lo = exp(confint(cx))[1], HR_hi = exp(confint(cx))[2],
                     cox_p = s$coefficients[1, "Pr(>|z|)"],
                     HR_per_SD = exp(coef(cx) * sd_score),
                     HR_per_SD_lo = exp(confint(cx)[1] * sd_score),
                     HR_per_SD_hi = exp(confint(cx)[2] * sd_score),
                     score_SD = sd_score),
          file.path(OUT, "tables", "Scissor_KM_stats_95.csv"), row.names = FALSE)
