#!/usr/bin/env Rscript
# MMSage M4: monocle3 trajectories over a comprehensive parameter grid.
# Each metabolite is a monocle3 cell and each study sample is a feature.

suppressMessages({
  library(monocle3)
  library(tidyr)
  library(Matrix)
  library(data.table)
  library(parallel)
})

MC <- "E:/Onedrive/mon345/02_TRAjMM/heart/metacard"
OUT <- file.path(MC, "mmsage_out")
COORD_DIR <- file.path(OUT, "coords")
LOG_DIR <- file.path(OUT, "logs")
dir.create(COORD_DIR, showWarnings=FALSE, recursive=TRUE)
dir.create(LOG_DIR, showWarnings=FALSE, recursive=TRUE)

args <- commandArgs(trailingOnly=TRUE)
SMOKE <- "smoke" %in% args
RESUME <- !SMOKE && !("no-resume" %in% args)
SEED <- 1L

read_positive_int <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset=as.character(default))))
  if (length(value) != 1L || is.na(value) || value < 1L) default else value
}

UMAP_CORES <- if (SMOKE) 1L else read_positive_int("MMSAGE_UMAP_CORES", 1L)

if (SMOKE) {
  grid <- expand.grid(
    num_dim=c(2, 3),
    neighbors=c(5, 10),
    min_dist=0.1,
    metric="cosine",
    cluster_res=0.1,
    KEEP.OUT.ATTRS=FALSE,
    stringsAsFactors=FALSE
  )
} else {
  grid <- expand.grid(
    num_dim=c(2, 3, 5, 8, 10, 15, 20),
    neighbors=c(2, 3, 5, 10, 15, 20),
    min_dist=c(0.1, 0.5, 1),
    metric=c("cosine", "euclidean"),
    cluster_res=c(0.1, 0.5),
    KEEP.OUT.ATTRS=FALSE,
    stringsAsFactors=FALSE
  )
}

roots_file <- file.path(OUT, "roots", "all_roots_noback.csv")
cache_file <- file.path(OUT, "M1_cache.rds")
stopifnot(file.exists(roots_file), file.exists(cache_file))

roots_df <- read.csv(roots_file, stringsAsFactors=FALSE, check.names=FALSE)
required_root_columns <- c("CAG", "Root", "NArank")
if (!all(required_root_columns %in% names(roots_df))) {
  stop("Root manifest is missing required columns: ",
       paste(setdiff(required_root_columns, names(roots_df)), collapse=", "))
}

cache <- readRDS(cache_file)
cags <- cache$all_target_cags
if (SMOKE) cags <- cags[1]

format_number <- function(x) format(x, trim=TRUE, scientific=FALSE)

config_tag <- function(num_dim, neighbors, min_dist, metric, cluster_res) {
  sprintf(
    "dim%d_nbr%d_dist%s_met%s_clu%s",
    as.integer(num_dim),
    as.integer(neighbors),
    format_number(min_dist),
    metric,
    format_number(cluster_res)
  )
}

coord_path <- function(cag, tag, narank) {
  file.path(
    COORD_DIR,
    sprintf("coords_%s_seed%d_%s_NArank%d.csv", cag, SEED, tag, as.integer(narank))
  )
}

write_coord_file <- function(x, path) {
  temporary <- sprintf("%s.tmp.%d", path, Sys.getpid())
  on.exit(unlink(temporary), add=TRUE)
  write.csv(x, temporary, row.names=FALSE)
  if (file.exists(path) && !file.remove(path)) {
    stop("Unable to replace coordinate file: ", path)
  }
  if (!file.rename(temporary, path)) {
    stop("Unable to finalize coordinate file: ", path)
  }
}

dimension_groups <- list(
  low=c(2, 3, 5, 8),
  high=c(10, 15, 20)
)

if (SMOKE) {
  tasks <- list(list(
    cag=cags[[1]],
    task_id="smoke",
    result_label=cags[[1]],
    grid=grid,
    progress_file=file.path(LOG_DIR, paste0("M4_progress_", cags[[1]], ".tsv"))
  ))
} else {
  task_spec <- expand.grid(
    metric=c("cosine", "euclidean"),
    dimension_group=names(dimension_groups),
    cag=cags,
    KEEP.OUT.ATTRS=FALSE,
    stringsAsFactors=FALSE
  )
  tasks <- lapply(seq_len(nrow(task_spec)), function(i) {
    spec <- task_spec[i, , drop=FALSE]
    task_id <- sprintf("met%s_dim%s", spec$metric, spec$dimension_group)
    task_grid <- grid[
      grid$metric == spec$metric &
        grid$num_dim %in% dimension_groups[[spec$dimension_group]],
      ,
      drop=FALSE
    ]
    list(
      cag=spec$cag,
      task_id=task_id,
      result_label=sprintf("%s [%s]", spec$cag, task_id),
      grid=task_grid,
      progress_file=file.path(
        LOG_DIR,
        sprintf("M4_progress_%s_%s.tsv", spec$cag, task_id)
      )
    )
  })

  expected_task_count <- length(cags) * length(dimension_groups) *
    length(unique(grid$metric))
  if (length(tasks) != expected_task_count ||
      length(unique(vapply(tasks, `[[`, character(1), "progress_file"))) !=
        expected_task_count ||
      sum(vapply(tasks, function(task) nrow(task$grid), integer(1))) !=
        length(cags) * nrow(grid)) {
    stop("Full task grid is not a complete, non-overlapping partition")
  }
}

cat("========== MMSage M4 ==========" , "\n")
cat("Mode:", if (SMOKE) "SMOKE" else "FULL", "\n")
cat("Configurations per CAG:", nrow(grid), "\n")
cat("Target CAGs:", length(cags), "\n")
cat("Scheduled tasks:", length(tasks), "\n")
cat("UMAP cores per worker:", UMAP_CORES, "\n")
cat("Resume enabled:", RESUME, "\n\n")

progress_files_for_cag <- function(cag) {
  legacy_name <- paste0("M4_progress_", cag, ".tsv")
  task_prefix <- paste0("M4_progress_", cag, "_")
  candidates <- list.files(
    LOG_DIR,
    pattern="^M4_progress_.*[.]tsv$",
    full.names=TRUE
  )
  candidate_names <- basename(candidates)
  candidates[candidate_names == legacy_name | startsWith(candidate_names, task_prefix)]
}

read_completed_configs <- function(cag) {
  if (!RESUME) return(character(0))

  progress_files <- progress_files_for_cag(cag)
  progress_lines <- unlist(lapply(progress_files, function(path) {
    tryCatch(
      readLines(path, warn=FALSE),
      error=function(e) character(0)
    )
  }), use.names=FALSE)
  unique(sub(
    "^OK\t([^\t]+).*$",
    "\\1",
    grep("^OK\t", progress_lines, value=TRUE)
  ))
}

run_cag <- function(task) {
  cag <- task$cag
  task_grid <- task$grid
  progress_file <- task$progress_file
  result_label <- task$result_label
  pluscomb_file <- file.path(OUT, "pluscomb", paste0("pluscombno1_", cag, ".csv"))
  if (!file.exists(pluscomb_file)) {
    return(sprintf("SKIP %s: pluscomb file is missing", result_label))
  }

  if (SMOKE || !file.exists(progress_file)) {
    writeLines("status\tconfig\ttime\tmessage", progress_file)
  }
  completed <- read_completed_configs(cag)

  log_status <- function(status, tag, message="") {
    line <- paste(status, tag, format(Sys.time(), "%Y-%m-%dT%H:%M:%S"), message, sep="\t")
    cat(line, "\n", file=progress_file, append=TRUE, sep="")
  }

  raw <- data.table::fread(pluscomb_file, data.table=FALSE, showProgress=FALSE)
  required_input_columns <- c("samples", "rowcolumn", "sCor")
  if (!all(required_input_columns %in% names(raw))) {
    return(sprintf("FAIL %s: pluscomb columns are invalid", result_label))
  }

  wide <- tidyr::pivot_wider(raw, names_from=samples, values_from=sCor)
  mat <- as.matrix(wide[, -1, drop=FALSE])
  rownames(mat) <- wide$rowcolumn
  mat[is.na(mat)] <- 0
  rm(raw, wide)
  gc(verbose=FALSE)

  metabolite_names <- sub(paste0("-", cag, "$"), "", rownames(mat))
  cag_roots <- roots_df[roots_df$CAG == cag, required_root_columns, drop=FALSE]
  naranks <- sort(unique(as.integer(cag_roots$NArank)))
  if (nrow(cag_roots) == 0L || length(naranks) == 0L) {
    return(sprintf("FAIL %s: no roots in manifest", result_label))
  }

  cell_meta <- data.frame(
    samples=rownames(mat),
    gene_short_name=rownames(mat),
    row.names=rownames(mat),
    stringsAsFactors=FALSE
  )
  gene_meta <- data.frame(
    gene_short_name=colnames(mat),
    row.names=colnames(mat),
    stringsAsFactors=FALSE
  )
  base_cds <- new_cell_data_set(
    t(mat),
    cell_metadata=cell_meta,
    gene_metadata=gene_meta
  )
  base_cds <- base_cds[, Matrix::colSums(exprs(base_cds)) != 0]
  suppressWarnings(base_cds <- estimate_size_factors(base_cds))

  root_cells_by_rank <- lapply(naranks, function(narank) {
    candidates <- paste0(cag_roots$Root[cag_roots$NArank == narank], "-", cag)
    intersect(candidates, colnames(base_cds))
  })
  names(root_cells_by_rank) <- as.character(naranks)
  missing_root_ranks <- naranks[lengths(root_cells_by_rank) == 0L]
  if (length(missing_root_ranks) > 0L) {
    return(sprintf(
      "FAIL %s: no usable roots for NArank %s",
      result_label,
      paste(missing_root_ranks, collapse=",")
    ))
  }

  is_complete <- function(tag) {
    tag %in% completed && all(file.exists(vapply(
      naranks,
      function(narank) coord_path(cag, tag, narank),
      character(1)
    )))
  }

  n_ok <- 0L
  n_skip <- 0L
  n_fail <- 0L
  error_messages <- character(0)

  for (dim_value in unique(task_grid$num_dim)) {
    dim_grid <- task_grid[task_grid$num_dim == dim_value, , drop=FALSE]
    dim_tags <- mapply(
      config_tag,
      dim_grid$num_dim,
      dim_grid$neighbors,
      dim_grid$min_dist,
      dim_grid$metric,
      dim_grid$cluster_res,
      USE.NAMES=FALSE
    )
    if (RESUME && all(vapply(dim_tags, is_complete, logical(1)))) {
      n_skip <- n_skip + nrow(dim_grid)
      next
    }

    preprocessed <- tryCatch(
      suppressWarnings(preprocess_cds(base_cds, num_dim=dim_value, norm_method="none")),
      error=function(e) e
    )
    if (inherits(preprocessed, "error")) {
      message_text <- conditionMessage(preprocessed)
      n_fail <- n_fail + nrow(dim_grid)
      error_messages <- c(error_messages, message_text)
      for (tag in dim_tags) log_status("FAIL", tag, message_text)
      next
    }

    umap_grid <- unique(dim_grid[, c("neighbors", "min_dist", "metric"), drop=FALSE])
    for (u in seq_len(nrow(umap_grid))) {
      neighbors_value <- umap_grid$neighbors[u]
      min_dist_value <- umap_grid$min_dist[u]
      metric_value <- umap_grid$metric[u]
      cluster_values <- sort(unique(dim_grid$cluster_res[
        dim_grid$neighbors == neighbors_value &
          dim_grid$min_dist == min_dist_value &
          dim_grid$metric == metric_value
      ]))

      tags <- vapply(cluster_values, function(cluster_value) {
        config_tag(dim_value, neighbors_value, min_dist_value, metric_value, cluster_value)
      }, character(1))
      complete_flags <- if (RESUME) {
        vapply(tags, is_complete, logical(1))
      } else {
        rep(FALSE, length(tags))
      }
      if (RESUME && all(complete_flags)) {
        n_skip <- n_skip + length(tags)
        next
      }

      set.seed(SEED)
      embedded <- tryCatch(
        reduce_dimension(
          preprocessed,
          umap.n_neighbors=neighbors_value,
          umap.min_dist=min_dist_value,
          umap.metric=metric_value,
          preprocess_method="PCA",
          cores=UMAP_CORES,
          verbose=FALSE
        ),
        error=function(e) e
      )
      if (inherits(embedded, "error")) {
        message_text <- conditionMessage(embedded)
        failed_tags <- tags[!complete_flags]
        n_fail <- n_fail + length(failed_tags)
        error_messages <- c(error_messages, message_text)
        for (tag in failed_tags) log_status("FAIL", tag, message_text)
        next
      }

      for (j in seq_along(cluster_values)) {
        tag <- tags[j]
        if (RESUME && complete_flags[j]) {
          n_skip <- n_skip + 1L
          next
        }

        cluster_value <- cluster_values[j]
        result <- tryCatch({
          clustered <- cluster_cells(
            embedded,
            resolution=cluster_value,
            random_seed=SEED,
            verbose=FALSE
          )
          graph_cds <- learn_graph(
            clustered,
            use_partition=FALSE,
            verbose=FALSE
          )
          umap_coordinates <- reducedDims(graph_cds)$UMAP
          keep_index <- match(colnames(graph_cds), rownames(mat))
          if (anyNA(keep_index)) stop("Cell names no longer match the interaction matrix")

          for (narank in naranks) {
            ordered <- order_cells(
              graph_cds,
              root_cells=root_cells_by_rank[[as.character(narank)]]
            )
            pt <- pseudotime(ordered)
            if (length(pt) != ncol(graph_cds) || any(!is.finite(pt))) {
              stop("Non-finite or incomplete pseudotime values")
            }
            output <- data.frame(
              metabolite=metabolite_names[keep_index],
              rowcolumn=colnames(graph_cds),
              UMAP1=umap_coordinates[, 1],
              UMAP2=umap_coordinates[, 2],
              Pseudotime=as.numeric(pt[colnames(graph_cds)]),
              stringsAsFactors=FALSE
            )
            write_coord_file(output, coord_path(cag, tag, narank))
          }
          TRUE
        }, error=function(e) e)

        if (identical(result, TRUE)) {
          n_ok <- n_ok + 1L
          completed <- unique(c(completed, tag))
          log_status("OK", tag)
        } else {
          message_text <- conditionMessage(result)
          n_fail <- n_fail + 1L
          error_messages <- c(error_messages, message_text)
          log_status("FAIL", tag, message_text)
        }
      }
    }
    rm(preprocessed)
    gc(verbose=FALSE)
  }

  errors <- if (length(error_messages) > 0L) {
    paste0(" | errors: ", paste(unique(error_messages), collapse=" ; "))
  } else {
    ""
  }
  sprintf("%s: ok=%d skip=%d fail=%d%s", result_label, n_ok, n_skip, n_fail, errors)
}

run_task_safely <- function(task) {
  tryCatch(
    run_cag(task),
    error=function(e) sprintf(
      "FAIL %s: unhandled task error: %s",
      task$result_label,
      conditionMessage(e)
    )
  )
}

physical_cores <- suppressWarnings(as.integer(parallel::detectCores(logical=FALSE)))
if (length(physical_cores) != 1L || is.na(physical_cores) || physical_cores < 1L) {
  physical_cores <- 1L
}
default_workers <- min(length(tasks), physical_cores)
n_workers <- if (SMOKE) 1L else min(
  length(tasks),
  read_positive_int("MMSAGE_WORKERS", default_workers)
)
cat("Task workers:", n_workers, "\n")

start_time <- Sys.time()
if (n_workers > 1L) {
  cluster <- makeCluster(n_workers)
  clusterExport(
    cluster,
    c(
      "OUT", "COORD_DIR", "LOG_DIR", "roots_df", "required_root_columns",
      "SEED", "SMOKE",
      "RESUME", "UMAP_CORES", "format_number", "config_tag", "coord_path",
      "write_coord_file", "progress_files_for_cag", "read_completed_configs",
      "run_cag", "run_task_safely"
    )
  )
  clusterEvalQ(cluster, {
    suppressMessages({
      library(monocle3)
      library(tidyr)
      library(Matrix)
      library(data.table)
    })
  })
  results <- parLapply(cluster, tasks, run_task_safely)
  stopCluster(cluster)
} else {
  results <- lapply(tasks, run_task_safely)
}

elapsed_minutes <- as.numeric(difftime(Sys.time(), start_time, units="mins"))
result_lines <- unlist(results, use.names=FALSE)
summary_lines <- c(
  sprintf("mode=%s", if (SMOKE) "SMOKE" else "FULL"),
  sprintf("configurations_per_cag=%d", nrow(grid)),
  sprintf("scheduled_tasks=%d", length(tasks)),
  sprintf("task_workers=%d", n_workers),
  sprintf("elapsed_minutes=%.2f", elapsed_minutes),
  result_lines
)
writeLines(summary_lines, file.path(LOG_DIR, "M4_summary.log"))

cat("\nElapsed minutes:", round(elapsed_minutes, 2), "\n")
cat("M4 summary:\n", paste(result_lines, collapse="\n"), "\n", sep="")

has_failures <- length(result_lines) != length(tasks) ||
  anyNA(result_lines) ||
  any(grepl("fail=[1-9]|^FAIL|^SKIP", result_lines), na.rm=TRUE)
if (has_failures) quit(status=1L)
