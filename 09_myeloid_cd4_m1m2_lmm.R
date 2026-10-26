###############################################################################
# 09  CD4 and M1/M2 in tumour myeloid cells: LMM, descriptives, LOO
#                                                          [Methods 2.10]
###############################################################################
suppressPackageStartupMessages({
  library(Seurat); library(lme4); library(lmerTest); library(dplyr)
})
OUT1 <- "output1"
OUT2 <- "output2"
dir.create(file.path(OUT1, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT2, "tables"), showWarnings = FALSE, recursive = TRUE)

d   <- readRDS("myeloid_scored.rds")
dat <- GetAssayData(d, assay = "RNA", layer = "data")

m1_genes <- intersect(c("NOS2","TNF","IL1B","IL6","CXCL9","CXCL10","CXCL11","CD80",
                        "CD86","IDO1","IRF5","CCL3","CCL4"), rownames(d))
m2_genes <- intersect(c("CD163","MRC1","MS4A4A","C1QA","C1QB","C1QC","TGM2","IL10",
                        "ARG1","MERTK","FOLR2","CCL18","TGFB1"), rownames(d))

my <- data.frame(cell    = colnames(d),
                 sample  = d$sample,
                 group   = d$group,
                 subtype = as.character(d$Myeloid_subtype),
                 CD4     = as.numeric(dat["CD4", ]),
                 HLA_DRA = as.numeric(dat["HLA-DRA", ]),
                 M1      = colMeans(as.matrix(dat[m1_genes, ])),
                 M2      = colMeans(as.matrix(dat[m2_genes, ])),
                 nCount  = d$nCount_RNA,
                 row.names = colnames(d))
write.csv(my, file.path(OUT1, "tables", "myeloid_CD4_M1M2.csv"), row.names = FALSE)
cat(sprintf("myeloid %d | M1 genes %d | M2 genes %d\n", nrow(my), length(m1_genes), length(m2_genes)))

tum <- my[my$group == "Tumor", ]
tum$CD4_class <- factor(ifelse(tum$CD4 > 0, "CD4+", "CD4-"), levels = c("CD4-", "CD4+"))
cat(sprintf("tumour myeloid n = %d | CD4+ %d (%.1f%%)\n",
            nrow(tum), sum(tum$CD4_class == "CD4+"), 100 * mean(tum$CD4_class == "CD4+")))

# =============================================================================
# (A) ORIGINAL model: all tumour myeloid, continuous CD4, random intercept
# =============================================================================
cat("\n=== (A) all tumour myeloid ===\n")
m_m2 <- lmer(M2 ~ CD4 + (1 | sample), data = tum)
m_m1 <- lmer(M1 ~ CD4 + (1 | sample), data = tum)
s_m2 <- summary(m_m2); s_m1 <- summary(m_m1)

stats_df <- data.frame(
  response     = c("M2", "M1"),
  fixed_effect = c(fixef(m_m2)["CD4"], fixef(m_m1)["CD4"]),
  KR_p         = c(s_m2$coefficients["CD4", "Pr(>|t|)"], s_m1$coefficients["CD4", "Pr(>|t|)"]),
  singular     = c(isSingular(m_m2), isSingular(m_m1)),
  population   = "Tumor Myeloid (receiver)")
print(stats_df, digits = 5)
write.csv(stats_df, file.path(OUT1, "tables", "mixed_model_stats.csv"), row.names = FALSE)

desc_df <- data.frame(
  comparison   = c("M2_CD4pos_vs_neg", "M1_CD4pos_vs_neg"),
  CD4pos_mean  = c(mean(tum$M2[tum$CD4_class == "CD4+"]), mean(tum$M1[tum$CD4_class == "CD4+"])),
  CD4neg_mean  = c(mean(tum$M2[tum$CD4_class == "CD4-"]), mean(tum$M1[tum$CD4_class == "CD4-"])),
  Wilcoxon_p   = c(wilcox.test(M2 ~ CD4_class, tum)$p.value,
                   wilcox.test(M1 ~ CD4_class, tum)$p.value),
  note         = "DESCRIPTIVE only (cells nested in samples)")
print(desc_df, digits = 4)
write.csv(desc_df, file.path(OUT1, "tables", "descriptive_stats.csv"), row.names = FALSE)

cor_df <- data.frame(
  comparison   = c("CD4_vs_M2", "CD4_vs_M1", "CD4_vs_HLA_DRA"),
  spearman_rho = c(cor(tum$CD4, tum$M2, method = "spearman"),
                   cor(tum$CD4, tum$M1, method = "spearman"),
                   cor(tum$CD4, tum$HLA_DRA, method = "spearman")),
  note         = "DESCRIPTIVE only")
print(cor_df, digits = 4)
write.csv(cor_df, file.path(OUT1, "tables", "correlations_descriptive.csv"), row.names = FALSE)

# HLA-DRA+ receiver subset (step 4 of the original script)
hla_pos <- tum[tum$HLA_DRA > 0, ]
hla_pos$CD4_class <- factor(ifelse(hla_pos$CD4 > 0, "CD4+", "CD4-"), levels = c("CD4-", "CD4+"))
hla_df <- data.frame(
  population = "Tumor HLA-DRA+ Myeloid", n = nrow(hla_pos),
  n_CD4pos   = sum(hla_pos$CD4_class == "CD4+"),
  M2_CD4pos  = mean(hla_pos$M2[hla_pos$CD4_class == "CD4+"]),
  M2_CD4neg  = mean(hla_pos$M2[hla_pos$CD4_class == "CD4-"]),
  M2_p       = wilcox.test(M2 ~ CD4_class, hla_pos)$p.value,
  M1_CD4pos  = mean(hla_pos$M1[hla_pos$CD4_class == "CD4+"]),
  M1_CD4neg  = mean(hla_pos$M1[hla_pos$CD4_class == "CD4-"]),
  M1_p       = wilcox.test(M1 ~ CD4_class, hla_pos)$p.value)
print(hla_df, digits = 4)
write.csv(hla_df, file.path(OUT1, "tables", "HLA_DR_pos_receiver_subset.csv"), row.names = FALSE)

# =============================================================================
# (B) MANUSCRIPT primary: TAM only
# =============================================================================
cat("\n=== (B) TAM only ===\n")
tam <- tum[tum$subtype == "TAM", ]
cat(sprintf("TAM n = %d | patients %d | CD4+ %.1f%%\n",
            nrow(tam), length(unique(tam$sample)), 100 * mean(tam$CD4_class == "CD4+")))

extract <- function(m, label, term) {
  s  <- summary(m)$coefficients
  ci <- tryCatch(confint(m, parm = term, method = "profile", oldNames = FALSE),
                 error = function(e) matrix(NA, 2, 1))
  vc <- as.data.frame(VarCorr(m))
  data.frame(model = label, n = nrow(m@frame), estimate = s[term, "Estimate"],
             ci_lo = ci[1], ci_hi = ci[2], p = s[term, "Pr(>|t|)"],
             singular = isSingular(m), ICC = vc$vcov[1] / sum(vc$vcov))
}
tab <- bind_rows(
  extract(lmer(M2 ~ CD4 + (1|sample), data = tam), "LMM M2 ~ CD4", "CD4"),
  extract(lmer(M1 ~ CD4 + (1|sample), data = tam), "LMM M1 ~ CD4", "CD4"),
  extract(lmer(M2 ~ CD4_class + (1|sample), data = tam), "LMM M2 ~ CD4_class", "CD4_classCD4+"),
  extract(lmer(M1 ~ CD4_class + (1|sample), data = tam), "LMM M1 ~ CD4_class", "CD4_classCD4+"))
print(tab, digits = 4)
write.csv(tab, file.path(OUT2, "tables", "TAM_LMM_table.csv"), row.names = FALSE)

w <- data.frame(
  score = c("M2", "M1"),
  high  = c(mean(tam$M2[tam$CD4_class == "CD4+"]), mean(tam$M1[tam$CD4_class == "CD4+"])),
  low   = c(mean(tam$M2[tam$CD4_class == "CD4-"]), mean(tam$M1[tam$CD4_class == "CD4-"])),
  p     = c(wilcox.test(M2 ~ CD4_class, tam)$p.value,
            wilcox.test(M1 ~ CD4_class, tam)$p.value))
print(w, digits = 4)
write.csv(w, file.path(OUT2, "tables", "TAM_wilcoxon.csv"), row.names = FALSE)

# ---- (B3) per-cell Spearman within TAM (descriptive, Methods 2.10) ----------
# Reported in Results 3.9 as rho = 0.45 (M2) and 0.08 (M1). Cells are nested in
# 3 patients, so these are descriptive only - no patient-level inference.
sp <- function(x, y) suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
partial_sp <- function(x, y, z) {          # rank-residual partial correlation
  suppressWarnings(cor.test(residuals(lm(rank(x) ~ rank(z))),
                            residuals(lm(rank(y) ~ rank(z))), method = "spearman"))
}

s_m2 <- sp(tam$CD4, tam$M2); s_m1 <- sp(tam$CD4, tam$M1)
cor_tab <- data.frame(
  comparison   = c("CD4_vs_M2", "CD4_vs_M1"),
  n            = nrow(tam),
  spearman_rho = c(unname(s_m2$estimate), unname(s_m1$estimate)),
  p            = c(s_m2$p.value, s_m1$p.value),
  rho_partial_nCount = c(partial_sp(tam$CD4, tam$M2, log10(tam$nCount))$estimate,
                         partial_sp(tam$CD4, tam$M1, log10(tam$nCount))$estimate),
  note         = "cell-level descriptive (294 TAM cells, 3 patients)")
print(cor_tab, digits = 4)
write.csv(cor_tab, file.path(OUT2, "tables", "TAM_CD4_correlation.csv"), row.names = FALSE)

# per-patient direction (the honest version: 3 independent units, not 294)
per <- tam %>% group_by(sample) %>%
  summarise(n = n(),
            rho_M2 = cor(CD4, M2, method = "spearman"),
            rho_M1 = cor(CD4, M1, method = "spearman"), .groups = "drop")
print(as.data.frame(per), digits = 4)
write.csv(per, file.path(OUT2, "tables", "TAM_CD4_correlation_per_patient.csv"), row.names = FALSE)
cat(sprintf("patients with positive CD4~M2 rho: %d/%d\n", sum(per$rho_M2 > 0), nrow(per)))

loo <- bind_rows(lapply(unique(tam$sample), function(s) {
  sub <- tam[tam$sample != s, ]
  m <- lmer(M2 ~ CD4 + (1|sample), data = sub)
  ci <- confint(m, parm = "CD4", method = "profile", oldNames = FALSE)
  data.frame(omitted = s, n = nrow(sub), estimate = fixef(m)["CD4"],
             lo = ci[1], hi = ci[2], p = summary(m)$coefficients["CD4", "Pr(>|t|)"])
}))
print(loo, digits = 4)
write.csv(loo, file.path(OUT2, "tables", "TAM_LMM_leave_one_out.csv"), row.names = FALSE)

# fixed-effect sensitivity (sample as a fixed factor)
fe <- bind_rows(
  data.frame(model = "FE M2 ~ CD4 + sample",
             estimate = coef(lm(M2 ~ CD4 + factor(sample), data = tam))["CD4"],
             p = summary(lm(M2 ~ CD4 + factor(sample), data = tam))$coefficients["CD4", 4]),
  data.frame(model = "FE M1 ~ CD4 + sample",
             estimate = coef(lm(M1 ~ CD4 + factor(sample), data = tam))["CD4"],
             p = summary(lm(M1 ~ CD4 + factor(sample), data = tam))$coefficients["CD4", 4]))
print(fe, digits = 4)

# residual diagnostics for the TAM M2 model
m_tam <- lmer(M2 ~ CD4 + (1|sample), data = tam)
sw <- shapiro.test(residuals(m_tam) / sigma(m_tam))
cat(sprintf("Shapiro-Wilk on TAM M2 residuals: W = %.4f, p = %.3f\n", sw$statistic, sw$p.value))
