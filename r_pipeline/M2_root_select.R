#!/usr/bin/env Rscript
# MMSage M2: Noback根选择 (复用已有Spearman相关矩阵)

MC <- Sys.getenv("MMSAGE_JOB_DIR")
if (!nzchar(MC)) stop("MMSAGE_JOB_DIR must identify the job input directory")
OUT <- Sys.getenv("MMSAGE_OUTPUT_DIR", unset=file.path(MC, "mmsage_out"))
dir.create(file.path(OUT, "roots"), showWarnings=FALSE, recursive=TRUE)

cat("========== MMSage M2: 根选择 (Noback) ==========\n")

cache <- readRDS(file.path(OUT, "M1_cache.rds"))
clr_met <- cache$clr_met
all_target_cags <- cache$all_target_cags

cor_mat <- as.matrix(read.table(file.path(MC, "results/method02_spearman_matrix.tsv"),
                                header=TRUE, sep="\t", row.names=1, check.names=FALSE, quote='"', comment.char=""))
pval_mat <- as.matrix(read.table(file.path(MC, "results/method02_spearman_pvalues.tsv"),
                                 header=TRUE, sep="\t", row.names=1, check.names=FALSE, quote='"', comment.char=""))
cat("Spearman矩阵:", nrow(cor_mat), "x", ncol(cor_mat), "\n")

# Root selection follows TreeMM's root_metabolite_analysis, which uses
# SpiecEasi CLR on the abundance table (separately from the RCLR tensor).
met_root <- read.table(file.path(MC, "data/metabolites_wide.tsv"), header=TRUE,
                       sep="\t", check.names=FALSE, row.names=1, quote='"', comment.char="")
met_root <- as.matrix(met_root)
positive_root <- met_root[is.finite(met_root) & met_root > 0]
if (!length(positive_root)) stop("Metabolite table has no positive values")
met_root[met_root == 0] <- min(positive_root)
# Equivalent CLR implementation avoids requiring SpiecEasi solely for this
# transform and works on the server's current R/Bioconductor versions.
met_root <- t(met_root)
met_root <- met_root / pmax(rowSums(met_root), 1e-12)
log_root <- log(met_root)
clr_met_root <- sweep(log_root, 1, rowMeans(log_root), "-")
clr_met_mean <- colMeans(clr_met_root, na.rm=TRUE)

NArank_grid <- c(1, 2, 3, 5)
# TreeMM's root analysis defaults to P_threshold=1; expose it for configuration
# while retaining the same all-positive-edge behavior for arbitrary inputs.
p_thresh <- suppressWarnings(as.numeric(Sys.getenv("MMSAGE_ROOT_P_THRESHOLD", unset="1")))
if (!is.finite(p_thresh) || p_thresh < 0 || p_thresh > 1) p_thresh <- 1

all_roots <- list()

for (cag in all_target_cags) {
  if (!cag %in% rownames(cor_mat)) {
    cat("  ", cag, "不在Spearman矩阵中，跳过\n")
    next
  }
  cors <- cor_mat[cag, ]
  pvals <- pval_mat[cag, ]

  # Match TreeMM root_metabolite_analysis: first restrict to the original
  # metabolite universe, then retain positive edges under P_threshold.
  valid <- which(names(cors) %in% colnames(met_root) & !is.na(cors) &
                 !is.na(pvals) & cors > 0 & pvals < p_thresh)
  if (length(valid) == 0) {
    # Keep the pipeline total for arbitrary uploaded tables. When a feature has
    # no nominally significant positive edge, use its strongest finite positive
    # correlation as a deterministic fallback root.
    valid <- which(names(cors) %in% colnames(met_root) & !is.na(cors) & is.finite(cors) & cors > 0)
    if (length(valid) == 0) {
      valid <- which(names(cors) %in% colnames(met_root) & !is.na(cors) & is.finite(cors))
    }
    if (length(valid) == 0) {
      cat("  ", cag, "没有可用相关性，跳过\n")
      next
    }
    cat("  ", cag, "无显著正相关，使用最强相关代谢物作为回退根\n")
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
