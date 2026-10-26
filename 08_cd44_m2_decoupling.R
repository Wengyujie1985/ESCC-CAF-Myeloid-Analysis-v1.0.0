###############################################################################
# 08  CD44 vs the M2 program: Simpson's paradox   [Methods 2.9]
#
###############################################################################
suppressPackageStartupMessages({library(Seurat); library(dplyr)})

OUT <- "output"
dir.create(file.path(OUT, "tables"), showWarnings = FALSE, recursive = TRUE)

d <- readRDS("myeloid_scored.rds")
dat <- d@meta.data
dat$CD44    <- as.numeric(GetAssayData(d, assay = "RNA", layer = "data")["CD44", ])
dat$subtype  <- as.character(dat$Myeloid_subtype)
dat$log10_nCount <- log10(dat$nCount_RNA)
cat("myeloid cells:", nrow(dat), "| CD44 detected:",
    sprintf("%.1f%%\n", 100 * mean(dat$CD44 > 0)))

cliff <- function(a, b) {
  r <- rank(c(a, b)); n1 <- length(a)
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * length(b)) * 2 - 1
}
partial_sp <- function(x, y, z) {
  cor(residuals(lm(rank(x) ~ rank(z))), residuals(lm(rank(y) ~ rank(z))))
}

# ---- 1. whole-myeloid median split ------------------------------------------
cut <- median(dat$CD44)
dat$CD44_group <- factor(ifelse(dat$CD44 > cut, "CD44-high", "CD44-low"),
                         levels = c("CD44-low", "CD44-high"))
w <- wilcox.test(M2_score ~ CD44_group, data = dat)
cat(sprintf("median split at %.3f | Wilcoxon p = %.2e | Cliff's d = %+.3f\n",
            cut, w$p.value,
            cliff(dat$M2_score[dat$CD44_group == "CD44-high"],
                  dat$M2_score[dat$CD44_group == "CD44-low"])))
cat("subtype composition differs (this is the confounder):\n")
print(round(100 * prop.table(table(dat$CD44_group, dat$subtype), 1), 1))

# ---- 2. within-subtype ------------------------------------------------------
res <- lapply(c("Monocyte", "TAM", "DC"), function(s) {
  x <- filter(dat, subtype == s)
  ct <- cor.test(x$CD44, x$M2_score, method = "spearman")
  data.frame(subtype = s, n = nrow(x), rho = unname(ct$estimate),
             p = ct$p.value,
             partial_rho_nCount = partial_sp(x$CD44, x$M2_score, x$log10_nCount))
}) %>% bind_rows()
res$p_adj <- p.adjust(res$p, "BH")
print(res, digits = 4)
write.csv(res, file.path(OUT, "tables", "CD44_M2_correlation_by_subtype.csv"), row.names = FALSE)

# ---- 3. within-patient direction (descriptive) ------------------------------
per <- dat %>% group_by(orig.ident) %>% filter(n() >= 30) %>%
  summarise(n = n(), rho = cor(CD44, M2_score, method = "spearman"), .groups = "drop")
print(as.data.frame(per))
cat(sprintf("patients with positive rho: %d/%d\n", sum(per$rho > 0), nrow(per)))
