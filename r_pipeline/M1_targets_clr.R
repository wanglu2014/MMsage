#!/usr/bin/env Rscript
# MMSage M1: 目标菌CAG解析 + CLR变换
# 所有输入/输出均在 heart/metacard 下

suppressMessages(library(SpiecEasi))

MC <- "E:/Onedrive/mon345/02_TRAjMM/heart/metacard"
OUT <- file.path(MC, "mmsage_out")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

cat("========== MMSage M1: 目标菌 + CLR ==========\n")

# 1. 读taxonomy映射(CAG→species)
taxon <- read.table(file.path(MC,"data/taxon_names.tsv"), header=TRUE, sep="\t",
                    stringsAsFactors=FALSE, comment.char="", quote="", fill=TRUE,
                    col.names=c("CAG","species"))

# 2. 目标菌
targets <- list(
  list(microbe="Dorea",                        metabolite="vanillactate",     direction=-1),
  list(microbe="Faecalibacterium prausnitzii", metabolite="PAGln",            direction=-1),
  list(microbe="Faecalibacterium prausnitzii", metabolite="4-cresyl sulfate", direction=-1),
  list(microbe="Roseburia faecis",             metabolite="PAGln",            direction=-1),
  list(microbe="Roseburia faecis",             metabolite="4-cresyl sulfate", direction=-1)
)

# 3. 读abundance
mic <- read.table(file.path(MC,"data/microbes_wide.tsv"), header=TRUE, sep="\t",
                 check.names=FALSE, row.names=1)
met <- read.table(file.path(MC,"data/metabolites_wide.tsv"), header=TRUE, sep="\t",
                 check.names=FALSE, row.names=1)
cat("微生物:", nrow(mic),"x",ncol(mic)," 代谢物:",nrow(met),"x",ncol(met),"\n")

# 4. 解析目标CAG(只保留在矩阵中的)
target_cags <- list()
for (t in unique(sapply(targets, function(x) x$microbe))) {
  idx <- grep(t, taxon$species, ignore.case=TRUE)
  cags <- intersect(taxon$CAG[idx], rownames(mic))
  target_cags[[t]] <- cags
  cat("  ", t, ":", length(cags), "CAG:", paste(cags, collapse=","), "\n")
}
all_target_cags <- unique(unlist(target_cags))
cat("目标CAG总数:", length(all_target_cags), "\n")

# 5. 目标代谢物列匹配
met_syn <- list(
  vanillactate = c("vanillactate"),
  PAGln = c("phenylacetylglutamine"),
  "4-cresyl sulfate" = c("p-cresol sulfate","4-cresyl sulfate")
)
target_mets <- list()
for (nm in names(met_syn)) {
  hit <- unique(unlist(lapply(met_syn[[nm]], function(s)
    grep(s, rownames(met), ignore.case=TRUE, value=TRUE))))
  target_mets[[nm]] <- hit
  cat("  代谢物", nm, "->", paste(hit, collapse=" | "), "\n")
}

# 6. CLR transformation within each sample; output remains features x samples.
mic_z <- as.matrix(mic)
mic_z[mic_z == 0] <- min(mic_z[mic_z != 0])
clr_mic <- clr(mic_z)

met_z <- as.matrix(met)
# CLR requires positive values, while the integrated metabolite matrix contains negatives.
if (min(met_z, na.rm=TRUE) <= 0) met_z <- met_z - min(met_z, na.rm=TRUE) + 1e-6
met_z[met_z == 0] <- min(met_z[met_z != 0])
clr_met <- clr(met_z)

stopifnot(identical(dim(clr_mic), dim(mic_z)),
          identical(dim(clr_met), dim(met_z)),
          max(abs(colMeans(clr_mic))) < 1e-10,
          max(abs(colMeans(clr_met))) < 1e-10)
cat("CLR complete clr_mic:", dim(clr_mic), "clr_met:", dim(clr_met), "\n")

# 7. 保存
saveRDS(list(targets=targets, target_cags=target_cags, all_target_cags=all_target_cags,
             target_mets=target_mets, clr_mic=clr_mic, clr_met=clr_met, taxon=taxon),
        file.path(OUT,"M1_cache.rds"))
cat("已保存: mmsage_out/M1_cache.rds\n")
cat("========== M1 完成 ==========\n")
