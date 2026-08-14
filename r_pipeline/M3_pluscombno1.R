#!/usr/bin/env Rscript
# MMSage M3: per-microbe 交互矩阵 (comb_XY)
# sCor[j,s] = clr_metabolite[j,s] * clr_microbe[m,s]

MC <- "E:/Onedrive/mon345/02_TRAjMM/heart/metacard"
OUT <- file.path(MC, "mmsage_out")
dir.create(file.path(OUT, "pluscomb"), showWarnings=FALSE, recursive=TRUE)

cat("========== MMSage M3: pluscombno1 ==========\n")

cache <- readRDS(file.path(OUT, "M1_cache.rds"))
# M1保存的是 特征×样本 方向，这里转成 样本×特征
clr_mic <- t(cache$clr_mic)   # samples x microbes
clr_met <- t(cache$clr_met)   # samples x metabolites
target_cags <- cache$all_target_cags  # CAG字符向量

commrow <- rownames(clr_mic)
stopifnot(identical(commrow, rownames(clr_met)))

cat("Samples:", length(commrow), "| Metabolites:", ncol(clr_met), "| Target CAGs:", length(target_cags), "\n")

for (m in target_cags) {
  vecY <- clr_mic[, m]  # 微生物向量
  # 一次性算所有代谢物: sCor[s,j] = clr_met[s,j] * vecY[s]
  scor_mat <- clr_met * vecY  # 样本 x 代谢物
  # 长格式
  long_df <- data.frame(
    samples = rep(commrow, times=ncol(clr_met)),
    rowcolumn = rep(paste(colnames(clr_met), m, sep="-"), each=length(commrow)),
    sCor = as.vector(scor_mat),
    stringsAsFactors = FALSE
  )
  out_file <- file.path(OUT, "pluscomb", paste0("pluscombno1_", m, ".csv"))
  write.csv(long_df, out_file, row.names=FALSE)
  cat("  ", m, "-> rows:", nrow(long_df), "\n")
}

cat("\nM3完成: ", file.path(OUT, "pluscomb/"), "\n")
