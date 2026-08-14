#!/usr/bin/env Rscript
# MMSage M2: Noback根选择 (复用已有Spearman相关矩阵)

MC <- "E:/Onedrive/mon345/02_TRAjMM/heart/metacard"
OUT <- file.path(MC, "mmsage_out")
dir.create(file.path(OUT, "roots"), showWarnings=FALSE, recursive=TRUE)

cat("========== MMSage M2: 根选择 (Noback) ==========\n")

cache <- readRDS(file.path(OUT, "M1_cache.rds"))
clr_met <- cache$clr_met
all_target_cags <- cache$all_target_cags

cor_mat <- as.matrix(read.table(file.path(MC, "results/method02_spearman_matrix.tsv"),
                                header=TRUE, sep="\t", row.names=1, check.names=FALSE))
pval_mat <- as.matrix(read.table(file.path(MC, "results/method02_spearman_pvalues.tsv"),
                                 header=TRUE, sep="\t", row.names=1, check.names=FALSE))
cat("Spearman矩阵:", nrow(cor_mat), "x", ncol(cor_mat), "\n")

clr_met_mean <- rowMeans(clr_met, na.rm=TRUE)

NArank_grid <- c(1, 2, 3, 5)
p_thresh <- 0.05

all_roots <- list()

for (cag in all_target_cags) {
  if (!cag %in% rownames(cor_mat)) {
    cat("  ", cag, "不在Spearman矩阵中，跳过\n")
    next
  }
  cors <- cor_mat[cag, ]
  pvals <- pval_mat[cag, ]

  valid <- which(!is.na(cors) & !is.na(pvals) & cors > 0 & pvals < p_thresh)
  if (length(valid) == 0) {
    cat("  ", cag, "无正相关代谢物(P<0.05)，跳过\n")
    next
  }

  valid_mets <- names(cors)[valid]
  valid_means <- clr_met_mean[valid_mets]
  sorted_mets <- valid_mets[order(valid_means, decreasing=TRUE)]

  for (nar in NArank_grid) {
    roots <- head(sorted_mets, nar)
    for (r in roots) {
      all_roots[[length(all_roots)+1]] <- data.frame(
        CAG=cag, Root=r, NArank=nar,
        Cor=cors[r], Pval=pvals[r], CLR_mean=clr_met_mean[r],
        stringsAsFactors=FALSE)
    }
  }
  cat("  ", cag, ": 正相关代谢物", length(valid), "个, 根候选top5:",
      paste(head(sorted_mets,5), collapse=","), "\n")
}

roots_df <- do.call(rbind, all_roots)
write.csv(roots_df, file.path(OUT, "roots/all_roots_noback.csv"), row.names=FALSE)
cat("\n根选择完成:", nrow(roots_df), "条 (", length(unique(roots_df$CAG)), "个CAG )\n")
cat("NArank分布:\n")
print(table(roots_df$NArank))
cat("已保存: mmsage_out/roots/all_roots_noback.csv\n")
cat("========== M2 完成 ==========\n")
