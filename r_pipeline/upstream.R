# Extracted and adapted from monocle3 1.4.27 (MIT; Cole Trapnell).
# See LICENSE and README.md in this directory.
`new_cell_data_set` <- function (expression_data, cell_metadata = NULL, gene_metadata = NULL, verbose = FALSE) 
{
    assertthat::assert_that(methods::is(expression_data, "matrix") || is_sparse_matrix(expression_data) || is(expression_data, 
        "IterableMatrix"), msg = paste("Argument expression_data must be a", "matrix - either sparse from the", "Matrix package, dense,", 
        "or a BPCells matrix"))
    if (!is.null(cell_metadata)) {
        assertthat::assert_that(nrow(cell_metadata) == ncol(expression_data), msg = paste("cell_metadata must be NULL or have", 
            "the same number of rows as columns", "in expression_data"))
        assertthat::assert_that(!is.null(row.names(cell_metadata)) & all(row.names(cell_metadata) == colnames(expression_data)), 
            msg = paste("row.names of cell_metadata must be equal to colnames of", "expression_data"))
    }
    if (!is.null(gene_metadata)) {
        assertthat::assert_that(nrow(gene_metadata) == nrow(expression_data), msg = paste("gene_metadata must be NULL or have", 
            "the same number of rows as rows", "in expression_data"))
        assertthat::assert_that(!is.null(row.names(gene_metadata)) & all(row.names(gene_metadata) == row.names(expression_data)), 
            msg = paste("row.names of gene_metadata must be equal to row.names of", "expression_data"))
    }
    if (is.null(cell_metadata)) {
        cell_metadata <- data.frame(cell = colnames(expression_data), row.names = colnames(expression_data))
    }
    if (!("gene_short_name" %in% colnames(gene_metadata))) {
        warning("gene_metadata must contain a column verbatim ", "named 'gene_short_name' for certain functions.")
    }
    matrix_info <- get_matrix_info(mat = expression_data)
    if (!(matrix_info[["matrix_class"]] %in% c("r_dense_matrix", "dgCMatrix", "dgTMatrix", "BPCells"))) {
        stop("new_cell_data_set: invalid expression_data matrix class")
    }
    if (!(matrix_info[["matrix_class"]] %in% c("dgCMatrix", "BPCells"))) {
        expression_data <- methods::as(expression_data, "CsparseMatrix")
    }
    sce <- SingleCellExperiment(list(counts = expression_data), rowData = gene_metadata, colData = cell_metadata)
    cds <- methods::new("mmsage_cell_data_set", assays = SummarizedExperiment::Assays(list(counts = expression_data)), colData = colData(sce), 
        int_elementMetadata = SingleCellExperiment::int_elementMetadata(sce), int_colData = SingleCellExperiment::int_colData(sce), 
        int_metadata = SingleCellExperiment::int_metadata(sce), metadata = S4Vectors::metadata(sce), NAMES = NULL, elementMetadata = elementMetadata(sce)[, 
            0], rowRanges = rowRanges(sce))
    if (is(counts(cds), "IterableMatrix")) {
        cds <- set_cds_row_order_matrix(cds = cds)
    }
    S4Vectors::metadata(cds)$cds_version <- "1.4.27-extracted"
    clusters <- stats::setNames(S4Vectors::SimpleList(), character(0))
    cds <- estimate_size_factors(cds)
    cds
}

`estimate_size_factors` <- function (cds, round_exprs = TRUE, method = c("mean-geometric-mean-total", "mean-geometric-mean-log-total")) 
{
    method <- match.arg(method)
    if (methods::is(SingleCellExperiment::counts(cds), "IterableMatrix")) {
        if (any(BPCells::colSums(SingleCellExperiment::counts(cds)) == 0)) {
            warning("Your CDS object contains cells with zero reads. ", "This causes size factor calculation to fail. Please remove ", 
                "the zero read cells using ", "cds <- cds[,Matrix::colSums(counts(cds)) != 0] and then ", "run cds <- estimate_size_factors(cds)")
            return(cds)
        }
    }
    else {
        if (any(Matrix::colSums(SingleCellExperiment::counts(cds)) == 0)) {
            warning("Your CDS object contains cells with zero reads. ", "This causes size factor calculation to fail. Please remove ", 
                "the zero read cells using ", "cds <- cds[,Matrix::colSums(counts(cds)) != 0] and then ", "run cds <- estimate_size_factors(cds)")
            return(cds)
        }
    }
    if (is_sparse_matrix(SingleCellExperiment::counts(cds))) {
        size_factors(cds) <- estimate_sf_sparse(SingleCellExperiment::counts(cds), round_exprs = round_exprs, method = method)
    }
    else if (methods::is(SingleCellExperiment::counts(cds), "IterableMatrix")) {
        size_factors(cds) <- estimate_sf_bpcells(SingleCellExperiment::counts(cds), round_exprs = round_exprs, method = method)
    }
    else {
        size_factors(cds) <- estimate_sf_dense(SingleCellExperiment::counts(cds), round_exprs = round_exprs, method = method)
    }
    return(cds)
}

`preprocess_cds` <- function (cds, method = c("PCA", "LSI"), num_dim = 50, norm_method = c("log", "size_only", "none"), use_genes = NULL, pseudo_count = NULL, 
    scaling = TRUE, verbose = FALSE, build_nn_index = FALSE, nn_control = list()) 
{
    assertthat::assert_that(tryCatch(expr = ifelse(match.arg(method) == "", TRUE, TRUE), error = function(e) FALSE), msg = "method must be one of 'PCA' or 'LSI'")
    method <- match.arg(method)
    assertthat::assert_that(tryCatch(expr = ifelse(match.arg(norm_method) == "", TRUE, TRUE), error = function(e) FALSE), 
        msg = "norm_method must be one of 'log', 'size_only' or 'none'")
    norm_method <- match.arg(norm_method)
    assertthat::assert_that(assertthat::is.count(num_dim))
    if (!is.null(use_genes)) {
        assertthat::assert_that(is.character(use_genes))
        assertthat::assert_that(all(use_genes %in% row.names(rowData(cds))), msg = paste("use_genes must be NULL, or all must", 
            "be present in the row.names of rowData(cds)"))
    }
    assertthat::assert_that(!is.null(size_factors(cds)), msg = paste("You must call estimate_size_factors before calling", 
        "preprocess_cds."))
    assertthat::assert_that(sum(is.na(size_factors(cds))) == 0, msg = paste("One or more cells has a size factor of", "NA."))
    if (build_nn_index) {
        nn_control <- set_nn_control(mode = 1, nn_control = nn_control, nn_control_default = get_global_variable("nn_control_annoy_cosine"), 
            nn_index = NULL, k = NULL, verbose = verbose)
    }
    set.seed(2016)
    FM <- SingleCellExperiment::counts(cds)
    iterable_matrix_flag <- is(FM, "IterableMatrix")
    FM <- normalize_expr_data(FM = FM, size_factors = size_factors(cds), norm_method = norm_method, pseudo_count = pseudo_count)
    if (nrow(FM) == 0) {
        stop("all rows have standard deviation zero")
    }
    if (!is.null(use_genes)) {
        FM <- FM[use_genes, ]
    }
    if (method == "PCA") {
        cds <- initialize_reduce_dim_metadata(cds, "PCA")
        cds <- initialize_reduce_dim_model_identity(cds, "PCA")
        if (verbose) 
            message("Remove noise by PCA ...")
        if (!iterable_matrix_flag) {
            if (verbose) {
                message("preprocess_cds: FM matrix class: ", class(FM))
                message()
                message("preprocess_cds: str(FM):")
                message(str(FM))
                message()
            }
            fm_rowsums = Matrix::rowSums(FM)
            FM <- FM[is.finite(fm_rowsums) & fm_rowsums != 0, ]
            irlba_res <- sparse_prcomp_irlba(Matrix::t(FM), n = min(num_dim, min(dim(FM)) - 1), center = scaling, scale. = scaling, 
                verbose = verbose)
        }
        else {
            if (verbose) {
                message("preprocess_cds: FM matrix info:")
                message(show_matrix_info(matrix_info = get_matrix_info(mat = FM), "  "), appendLF = FALSE)
                message()
            }
            fm_rowsums = BPCells::rowSums(FM)
            FM <- FM[is.finite(fm_rowsums) & fm_rowsums != 0, ]
            irlba_res <- bpcells_prcomp_irlba(BPCells::t(FM), n = min(num_dim, min(dim(FM)) - 1), center = scaling, scale. = scaling, 
                verbose = verbose)
        }
        preproc_res <- irlba_res$x
        row.names(preproc_res) <- colnames(cds)
        SingleCellExperiment::reducedDims(cds)[[method]] <- as.matrix(preproc_res)
        irlba_rotation <- irlba_res$rotation
        row.names(irlba_rotation) <- rownames(FM)
        cds@reduce_dim_aux[["PCA"]][["model"]][["num_dim"]] <- num_dim
        cds@reduce_dim_aux[["PCA"]][["model"]][["norm_method"]] <- norm_method
        cds@reduce_dim_aux[["PCA"]][["model"]][["use_genes"]] <- use_genes
        cds@reduce_dim_aux[["PCA"]][["model"]][["pseudo_count"]] <- pseudo_count
        cds@reduce_dim_aux[["PCA"]][["model"]][["svd_v"]] <- irlba_rotation
        cds@reduce_dim_aux[["PCA"]][["model"]][["svd_sdev"]] <- irlba_res$sdev
        cds@reduce_dim_aux[["PCA"]][["model"]][["svd_center"]] <- irlba_res$center
        cds@reduce_dim_aux[["PCA"]][["model"]][["svd_scale"]] <- irlba_res$svd_scale
        cds@reduce_dim_aux[["PCA"]][["model"]][["prop_var_expl"]] <- irlba_res$sdev^2/sum(irlba_res$sdev^2)
        matrix_id <- get_unique_id(SingleCellExperiment::reducedDims(cds)[["PCA"]])
        counts_identity <- get_counts_identity(cds)
        cds <- set_reduce_dim_matrix_identity(cds, "PCA", "matrix:PCA", matrix_id, counts_identity[["matrix_type"]], counts_identity[["matrix_id"]], 
            "matrix:PCA", matrix_id)
        cds <- set_reduce_dim_model_identity(cds, "PCA", "matrix:PCA", matrix_id, "none", "none")
        if (build_nn_index) {
            nn_index <- make_nn_index(subject_matrix = SingleCellExperiment::reducedDims(cds)[[method]], nn_control = nn_control, 
                verbose = verbose)
            cds <- tryCatch(set_cds_nn_index(cds = cds, reduction_method = method, nn_index = nn_index, verbose = verbose), 
                error = function(c) {
                  stop(paste0(trimws(c), "\n* error in preprocess_cds"))
                })
        }
        else {
            cds <- tryCatch(clear_cds_nn_index(cds = cds, reduction_method = method, nn_method = "all"), error = function(c) {
                stop(paste0(trimws(c), "\n* error in preprocess_cds"))
            })
        }
    }
    else if (method == "LSI") {
        cds <- initialize_reduce_dim_metadata(cds, "LSI")
        cds <- initialize_reduce_dim_model_identity(cds, "LSI")
        if (!iterable_matrix_flag) {
            fm_rowsums <- Matrix::rowSums(FM)
        }
        else {
            fm_rowsums <- BPCells::rowSums(FM)
        }
        FM <- FM[is.finite(fm_rowsums) & fm_rowsums != 0, ]
        tfidf_res <- tfidf(count_matrix = FM, iterable_matrix_flag = iterable_matrix_flag)
        preproc_res <- tfidf_res[["tf_idf_counts"]]
        num_col <- ncol(preproc_res)
        if (!iterable_matrix_flag) {
            irlba_res <- irlba::irlba(A = Matrix::t(preproc_res), nv = min(num_dim, min(dim(FM)) - 1))
        }
        else {
            matrix_control_res <- set_matrix_control_pca(mat = FM, verbose = verbose)
            preproc_res_commit <- set_matrix_class(mat = BPCells::t(preproc_res), matrix_control = matrix_control_res)
            irlba_res <- irlba::irlba(A = BPCells:::linear_operator(preproc_res_commit), nv = min(num_dim, min(dim(FM)) - 
                1))
            rm_bpcells_dir(mat = preproc_res_commit)
            gc()
        }
        if (verbose) {
            message("singular values (head)")
            message(paste(head(irlba_res$d), collapse = " "))
            message("")
            message("umat: ", paste(dim(irlba_res$u), collapse = " "))
            message("vtmat: ", paste(dim(irlba_res$v), collapse = " "))
        }
        preproc_res <- irlba_res$u %*% diag(irlba_res$d)
        row.names(preproc_res) <- colnames(cds)
        SingleCellExperiment::reducedDims(cds)[[method]] <- as.matrix(preproc_res)
        irlba_rotation = irlba_res$v
        row.names(irlba_rotation) = rownames(FM)
        cds@reduce_dim_aux[["LSI"]][["model"]][["num_dim"]] <- num_dim
        cds@reduce_dim_aux[["LSI"]][["model"]][["norm_method"]] <- norm_method
        cds@reduce_dim_aux[["LSI"]][["model"]][["use_genes"]] <- use_genes
        cds@reduce_dim_aux[["LSI"]][["model"]][["pseudo_count"]] <- pseudo_count
        cds@reduce_dim_aux[["LSI"]][["model"]][["log_scale_tf"]] <- tfidf_res[["log_scale_tf"]]
        cds@reduce_dim_aux[["LSI"]][["model"]][["frequencies"]] <- tfidf_res[["frequencies"]]
        cds@reduce_dim_aux[["LSI"]][["model"]][["scale_factor"]] <- tfidf_res[["scale_factor"]]
        cds@reduce_dim_aux[["LSI"]][["model"]][["col_sums"]] <- tfidf_res[["col_sums"]]
        cds@reduce_dim_aux[["LSI"]][["model"]][["row_sums"]] <- tfidf_res[["row_sums"]]
        cds@reduce_dim_aux[["LSI"]][["model"]][["num_cols"]] <- tfidf_res[["num_cols"]]
        cds@reduce_dim_aux[["LSI"]][["model"]][["svd_v"]] <- irlba_rotation
        cds@reduce_dim_aux[["LSI"]][["model"]][["svd_sdev"]] <- irlba_res$d/sqrt(max(1, num_col - 1))
        matrix_id <- get_unique_id(SingleCellExperiment::reducedDims(cds)[["LSI"]])
        counts_identity <- get_counts_identity(cds)
        cds <- set_reduce_dim_matrix_identity(cds, "LSI", "matrix:LSI", matrix_id, counts_identity[["matrix_type"]], counts_identity[["matrix_id"]], 
            "matrix:LSI", matrix_id)
        cds <- set_reduce_dim_model_identity(cds, "LSI", "matrix:LSI", matrix_id, "none", "none")
        if (build_nn_index) {
            nn_index <- make_nn_index(subject_matrix = SingleCellExperiment::reducedDims(cds)[[method]], nn_control = nn_control, 
                verbose = verbose)
            cds <- tryCatch(set_cds_nn_index(cds = cds, reduction_method = method, nn_index = nn_index, verbose = verbose), 
                error = function(c) {
                  stop(paste0(trimws(c), "\n* error in preprocess_cds"))
                })
        }
        else {
            cds <- tryCatch(clear_cds_nn_index(cds = cds, reduction_method = method, nn_method = "all"), error = function(c) {
                stop(paste0(trimws(c), "\n* error in preprocess_cds"))
            })
        }
    }
    if (!is.null(cds@reduce_dim_aux[["Aligned"]]) && !is.null(cds@reduce_dim_aux[["Aligned"]][["model"]][["beta"]])) {
        cds@reduce_dim_aux[["Aligned"]][["model"]][["beta"]] <- NULL
    }
    cds
}

`reduce_dimension` <- function (cds, max_components = 2, reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"), preprocess_method = NULL, 
    umap.metric = "cosine", umap.min_dist = 0.1, umap.n_neighbors = 15L, umap.fast_sgd = FALSE, umap.nn_method = "annoy", 
    verbose = FALSE, cores = 1, build_nn_index = FALSE, nn_control = list(), ...) 
{
    extra_arguments <- list(...)
    reduce_dim_preprocess_method_check = get_global_variable("reduce_dim_preprocess_method_check")
    assertthat::assert_that(tryCatch(expr = ifelse(match.arg(reduction_method) == "", TRUE, TRUE), error = function(e) FALSE), 
        msg = "reduction_method must be one of 'UMAP', 'PCA', 'tSNE', 'LSI', 'Aligned'")
    reduction_method <- match.arg(reduction_method)
    assertthat::assert_that(is.logical(build_nn_index), msg = paste("build_nn_index must be either TRUE or FALSE"))
    if (reduce_dim_preprocess_method_check) {
        if (is.null(preprocess_method)) {
            if ("Aligned" %in% names(SingleCellExperiment::reducedDims(cds))) {
                preprocess_method = "Aligned"
                message("No preprocess_method specified, and aligned coordinates ", "have been computed previously. Using preprocess_method = 'Aligned'")
            }
            else {
                preprocess_method = "PCA"
                message("No preprocess_method specified, using preprocess_method = 'PCA'")
            }
        }
        else {
            assertthat::assert_that(preprocess_method %in% c("PCA", "LSI", "Aligned"), msg = "preprocess_method must be one of 'PCA' or 'LSI'")
        }
    }
    if (build_nn_index) {
        if (reduction_method == "tSNE" || reduction_method == "UMAP") 
            nn_control_default <- get_global_variable("nn_control_annoy_euclidean")
        else nn_control_default <- get_global_variable("nn_control_annoy_cosine")
        nn_control <- set_nn_control(mode = 1, nn_control = nn_control, nn_control_default = nn_control_default, nn_index = NULL, 
            k = NULL, verbose = verbose)
    }
    assertthat::assert_that(assertthat::is.count(max_components))
    assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[preprocess_method]]), msg = paste("Data has not been preprocessed with", 
        "chosen method:", preprocess_method, "Please run preprocess_cds with", "method =", preprocess_method, "before running reduce_dimension."))
    if (reduce_dim_preprocess_method_check) {
        if (reduction_method == "PCA") {
            assertthat::assert_that(preprocess_method == "PCA", msg = paste("preprocess_method must be 'PCA' when", "reduction_method = 'PCA'"))
            assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[["PCA"]]), msg = paste("When reduction_method = 'PCA', the", 
                "cds must have been preprocessed for", "PCA. Please run preprocess_cds with", "method = 'PCA' before running", 
                "reduce_dimension with", "reduction_method = 'PCA'."))
        }
        if (reduction_method == "LSI") {
            assertthat::assert_that(preprocess_method == "LSI", msg = paste("preprocess_method must be 'LSI' when", "reduction_method = 'LSI'"))
            assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[["LSI"]]), msg = paste("When reduction_method = 'LSI', the", 
                "cds must have been preprocessed for", "LSI. Please run preprocess_cds with", "method = 'LSI' before running", 
                "reduce_dimension with", "reduction_method = 'LSI'."))
        }
        if (reduction_method == "Aligned") {
            assertthat::assert_that(preprocess_method == "Aligned", msg = paste("preprocess_method must be 'Aligned' when", 
                "reduction_method = 'Aligned'"))
            assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[["Aligned"]]), msg = paste("When reduction_method = 'Aligned', the", 
                "cds must have been aligned.", "Please run align_cds before running", "reduce_dimension with", "reduction_method = 'Aligned'."))
        }
    }
    set.seed(2016)
    if (reduction_method == "UMAP" && (umap.fast_sgd == TRUE || cores > 1)) {
        message("Note: reduce_dimension will produce slightly different ", "output each time you run it unless you set ", 
            "'umap.fast_sgd = FALSE' and 'cores = 1'")
    }
    preprocess_mat <- SingleCellExperiment::reducedDims(cds)[[preprocess_method]]
    if (reduction_method == "PCA") {
        if (build_nn_index && is.null(cds@reduce_dim_aux[[reduction_method]][["nn_index"]])) {
            nn_index <- make_nn_index(subject_matrix = SingleCellExperiment::reducedDims(cds)[[reduction_method]], nn_control = nn_control, 
                verbose = verbose)
            cds <- tryCatch(set_cds_nn_index(cds = cds, reduction_method = reduction_method, nn_index = nn_index, verbose = verbose), 
                error = function(c) {
                  stop(paste0(trimws(c), "\n* error in reduce_dimension"))
                })
        }
        if (verbose) 
            message("Returning preprocessed PCA matrix")
    }
    else if (reduction_method == "LSI") {
        if (build_nn_index && is.null(cds@reduce_dim_aux[[reduction_method]][["nn_index"]])) {
            nn_index <- make_nn_index(subject_matrix = SingleCellExperiment::reducedDims(cds)[[reduction_method]], nn_control = nn_control, 
                verbose = verbose)
            cds <- tryCatch(set_cds_nn_index(cds = cds, reduction_method = reduction_method, nn_index = nn_index, verbose = verbose), 
                error = function(c) {
                  stop(paste0(trimws(c), "\n* error in reduce_dimension"))
                })
        }
        if (verbose) 
            message("Returning preprocessed LSI matrix")
    }
    else if (reduction_method == "Aligned") {
        if (build_nn_index && is.null(cds@reduce_dim_aux[[reduction_method]][["nn_index"]])) {
            nn_index <- make_nn_index(subject_matrix = SingleCellExperiment::reducedDims(cds)[[reduction_method]], nn_control = nn_control, 
                verbose = verbose)
            cds <- tryCatch(set_cds_nn_index(cds = cds, reduction_method = reduction_method, nn_index = nn_index, verbose = verbose), 
                error = function(c) {
                  stop(paste0(trimws(c), "\n* error in reduce_dimension"))
                })
        }
        if (verbose) 
            message("Returning preprocessed Aligned matrix")
    }
    else if (reduction_method == "tSNE") {
        if (verbose) 
            message("Reduce dimension by tSNE ...")
        cds <- initialize_reduce_dim_metadata(cds, "tSNE")
        cds <- initialize_reduce_dim_model_identity(cds, "tSNE")
        tsne_res <- Rtsne::Rtsne(as.matrix(preprocess_mat), dims = max_components, pca = FALSE, check_duplicates = FALSE, 
            ...)
        if (max_components < 1) 
            warning("bad loop: max_components < 1")
        tsne_data <- tsne_res$Y[, 1:max_components]
        row.names(tsne_data) <- colnames(tsne_data)
        SingleCellExperiment::reducedDims(cds)$tSNE <- tsne_data
        matrix_id <- get_unique_id(SingleCellExperiment::reducedDims(cds)[["tSNE"]])
        reduce_dim_matrix_identity <- get_reduce_dim_matrix_identity(cds, preprocess_method)
        set_reduce_dim_matrix_identity(cds, "tSNE", "matrix:tSNE", matrix_id, reduce_dim_matrix_identity[["matrix_type"]], 
            reduce_dim_matrix_identity[["matrix_id"]], "matrix:tSNE", matrix_id)
        reduce_dim_model_identity <- get_reduce_dim_model_identity(cds, preprocess_method)
        set_reduce_dim_model_identity(cds, "tSNE", "matrix:tSNE", matrix_id, reduce_dim_model_identity[["model_type"]], reduce_dim_model_identity[["model_id"]])
        if (build_nn_index) {
            nn_index <- make_nn_index(subject_matrix = SingleCellExperiment::reducedDims(cds)[[reduction_method]], nn_control = nn_control, 
                verbose = verbose)
            cds <- tryCatch(set_cds_nn_index(cds = cds, reduction_method = reduction_method, nn_index = nn_index, verbose = verbose), 
                error = function(c) {
                  stop(paste0(trimws(c), "\n* error in reduce_dimension"))
                })
        }
        else {
            cds <- tryCatch(clear_cds_nn_index(cds = cds, reduction_method = reduction_method, nn_method = "all"), error = function(c) {
                stop(paste0(trimws(c), "\n* error in reduce_dimension"))
            })
        }
    }
    else if (reduction_method == c("UMAP")) {
        cds <- add_citation(cds, "UMAP")
        if (verbose) 
            message("Running Uniform Manifold Approximation and Projection")
        cds <- initialize_reduce_dim_metadata(cds, "UMAP")
        cds <- initialize_reduce_dim_model_identity(cds, "UMAP")
        umap_model <- uwot::umap(as.matrix(preprocess_mat), n_components = max_components, metric = umap.metric, min_dist = umap.min_dist, 
            n_neighbors = umap.n_neighbors, fast_sgd = umap.fast_sgd, n_threads = cores, verbose = verbose, nn_method = umap.nn_method, 
            ret_model = TRUE, ...)
        set.seed(2016)
        umap_res <- uwot::umap_transform(X = as.matrix(preprocess_mat), model = umap_model, n_threads = 1)
        row.names(umap_res) <- colnames(cds)
        SingleCellExperiment::reducedDims(cds)[["UMAP"]] <- umap_res
        cds@reduce_dim_aux[["UMAP"]][["model"]][["umap_preprocess_method"]] <- preprocess_method
        cds@reduce_dim_aux[["UMAP"]][["model"]][["max_components"]] <- max_components
        cds@reduce_dim_aux[["UMAP"]][["model"]][["umap_metric"]] <- umap.metric
        cds@reduce_dim_aux[["UMAP"]][["model"]][["umap_min_dist"]] <- umap.min_dist
        cds@reduce_dim_aux[["UMAP"]][["model"]][["umap_n_neighbors"]] <- umap.n_neighbors
        cds@reduce_dim_aux[["UMAP"]][["model"]][["umap_fast_sgd"]] <- umap.fast_sgd
        cds@reduce_dim_aux[["UMAP"]][["model"]][["umap_model"]] <- umap_model
        matrix_id <- get_unique_id(SingleCellExperiment::reducedDims(cds)[["UMAP"]])
        reduce_dim_matrix_identity <- get_reduce_dim_matrix_identity(cds, preprocess_method)
        cds <- set_reduce_dim_matrix_identity(cds, "UMAP", "matrix:UMAP", matrix_id, reduce_dim_matrix_identity[["matrix_type"]], 
            reduce_dim_matrix_identity[["matrix_id"]], "matrix:UMAP", matrix_id)
        reduce_dim_model_identity <- get_reduce_dim_model_identity(cds, preprocess_method)
        cds <- set_reduce_dim_model_identity(cds, "UMAP", "matrix:UMAP", matrix_id, reduce_dim_model_identity[["model_type"]], 
            reduce_dim_model_identity[["model_id"]])
        if (build_nn_index) {
            nn_index <- make_nn_index(subject_matrix = SingleCellExperiment::reducedDims(cds)[[reduction_method]], nn_control = nn_control, 
                verbose = verbose)
            cds <- tryCatch(set_cds_nn_index(cds = cds, reduction_method = reduction_method, nn_index = nn_index, verbose = verbose), 
                error = function(c) {
                  stop(paste0(trimws(c), "\n* error in reduce_dimension"))
                })
        }
        else {
            cds <- tryCatch(clear_cds_nn_index(cds = cds, reduction_method = reduction_method, nn_method = "all"), error = function(c) {
                stop(paste0(trimws(c), "\n* error in reduce_dimension"))
            })
        }
    }
    cds@principal_graph_aux[[reduction_method]] <- NULL
    cds@principal_graph[[reduction_method]] <- NULL
    cds@clusters[[reduction_method]] <- NULL
    cds
}

`cluster_cells` <- function (cds, reduction_method = c("UMAP", "tSNE", "PCA", "LSI", "Aligned"), k = 20, cluster_method = c("leiden", "louvain"), 
    num_iter = 2, partition_qval = 0.05, weight = FALSE, resolution = NULL, random_seed = 42, verbose = FALSE, nn_control = list(), 
    ...) 
{
    assertthat::assert_that(tryCatch(expr = ifelse(match.arg(reduction_method) == "", TRUE, TRUE), error = function(e) FALSE), 
        msg = "reduction_method must be one of 'UMAP', 'tSNE', 'PCA', 'LSI', 'Aligned'")
    reduction_method <- match.arg(reduction_method)
    assertthat::assert_that(tryCatch(expr = ifelse(match.arg(cluster_method) == "", TRUE, TRUE), error = function(e) FALSE), 
        msg = "cluster_method must be one of 'leiden', 'louvain'")
    cluster_method <- match.arg(cluster_method)
    assertthat::assert_that(methods::is(cds, "mmsage_cell_data_set"))
    assertthat::assert_that(is.character(reduction_method))
    assertthat::assert_that(is.logical(weight))
    assertthat::assert_that(assertthat::is.count(num_iter))
    assertthat::assert_that(assertthat::is.count(k))
    assertthat::assert_that(!is.null(colnames(cds)), msg = message("cluster_cells: the cds is missing cell names, which are required by cluster_cells."))
    if (!is.null(resolution) & cluster_method == "louvain") {
        message("Resolution can only be used when cluster_method is ", "'leiden'. Switching to leiden clustering.")
        cluster_method <- "leiden"
    }
    if (!is.null(resolution)) {
        assertthat::assert_that(is.numeric(resolution))
    }
    assertthat::assert_that(is.numeric(partition_qval))
    assertthat::assert_that(is.logical(verbose))
    assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]), msg = paste("No dimensionality reduction for", 
        reduction_method, "calculated.", "Please run reduce_dimension with", "reduction_method =", reduction_method, "before running cluster_cells"))
    if (reduction_method == "tSNE" || reduction_method == "UMAP") 
        nn_control_default <- get_global_variable("nn_control_annoy_euclidean")
    else nn_control_default <- get_global_variable("nn_control_annoy_cosine")
    nn_control <- set_nn_control(mode = 3, nn_control = nn_control, nn_control_default = nn_control_default, nn_index = NULL, 
        k = k, verbose = verbose)
    nn_method <- nn_control[["method"]]
    nn_index <- NULL
    reduced_dim_res <- SingleCellExperiment::reducedDims(cds)[[reduction_method]]
    if (is.null(random_seed)) {
        random_seed <- sample.int(.Machine$integer.max, 1)
    }
    if (verbose) 
        message("Running ", cluster_method, " clustering algorithm ...")
    if (cluster_method == "louvain") {
        cluster_result <- tryCatch(louvain_clustering(data = reduced_dim_res, pd = colData(cds), weight = weight, nn_index = nn_index, 
            k = k, nn_control = nn_control, louvain_iter = num_iter, random_seed = random_seed, verbose = verbose), error = function(c) {
            stop(paste0(trimws(c), "\n* error in cluster_cells"))
        })
        if (length(unique(cluster_result$optim_res$membership)) > 1) {
            cluster_graph_res <- compute_partitions(cluster_result$g, cluster_result$optim_res, partition_qval, verbose)
            partitions <- igraph::components(cluster_graph_res$cluster_g)$membership[cluster_result$optim_res$membership]
            partitions <- as.factor(partitions)
        }
        else {
            partitions <- rep(1, nrow(colData(cds)))
        }
        names(partitions) <- row.names(reduced_dim_res)
        clusters <- factor(igraph::membership(cluster_result$optim_res))
        cds@clusters[[reduction_method]] <- list(cluster_result = cluster_result, partitions = partitions, clusters = clusters)
    }
    else if (cluster_method == "leiden") {
        cds <- add_citation(cds, "leiden")
        cluster_result <- tryCatch(leiden_clustering(data = reduced_dim_res, pd = colData(cds), weight = weight, nn_index = nn_index, 
            k = k, nn_control = nn_control, num_iter = num_iter, resolution_parameter = resolution, random_seed = random_seed, 
            verbose = verbose, ...), error = function(c) {
            stop(paste0(trimws(c), "\n* error in cluster_cells"))
        })
        if (length(unique(cluster_result$optim_res$membership)) > 1) {
            cluster_graph_res <- compute_partitions(cluster_result$g, cluster_result$optim_res, partition_qval, verbose)
            partitions <- igraph::components(cluster_graph_res$cluster_g)$membership[cluster_result$optim_res$membership]
            partitions <- as.factor(partitions)
        }
        else {
            partitions <- rep(1, nrow(colData(cds)))
        }
        names(partitions) <- row.names(reduced_dim_res)
        clusters <- factor(igraph::membership(cluster_result$optim_res))
        cds@clusters[[reduction_method]] <- list(cluster_result = cluster_result, partitions = partitions, clusters = clusters)
    }
    cds <- add_citation(cds, "clusters")
    cds <- add_citation(cds, "partitions")
    return(cds)
}

`learn_graph` <- function (cds, use_partition = TRUE, close_loop = TRUE, learn_graph_control = NULL, verbose = FALSE) 
{
    reduction_method <- "UMAP"
    if (!is.null(learn_graph_control)) {
        assertthat::assert_that(methods::is(learn_graph_control, "list"))
        assertthat::assert_that(all(names(learn_graph_control) %in% c("euclidean_distance_ratio", "geodesic_distance_ratio", 
            "minimal_branch_len", "orthogonal_proj_tip", "prune_graph", "scale", "ncenter", "nn.k", "rann.k", "maxiter", 
            "eps", "L1.gamma", "L1.sigma", "nn.method", "nn.metric", "nn.n_trees", "nn.search_k", "nn.M", "nn.ef_construction", 
            "nn.ef", "nn.grain_size", "nn.cores")), msg = "Unknown variable in learn_graph_control")
    }
    if (!is.null(learn_graph_control[["rann.k"]]) && !is.null(learn_graph_control[["nn.k"]])) {
        assertthat::assert_that(learn_graph_control[["rann.k"]] == learn_graph_control[["nn.k"]], msg = paste0("both learn_graph_control$nn.k and learn_graph_control$rann.k are", 
            " defined and are unequal. See help(learn_graph) for more", " information."))
    }
    if (is.null(learn_graph_control[["nn.k"]]) && !is.null(learn_graph_control[["rann.k"]])) 
        learn_graph_control[["nn.k"]] <- learn_graph_control[["rann.k"]]
    euclidean_distance_ratio <- ifelse(is.null(learn_graph_control$euclidean_distance_ratio), 1, learn_graph_control$euclidean_distance_ratio)
    geodesic_distance_ratio <- ifelse(is.null(learn_graph_control$geodesic_distance_ratio), 1/3, learn_graph_control$geodesic_distance_ratio)
    minimal_branch_len <- ifelse(is.null(learn_graph_control$minimal_branch_len), 10, learn_graph_control$minimal_branch_len)
    orthogonal_proj_tip <- ifelse(is.null(learn_graph_control$orthogonal_proj_tip), FALSE, learn_graph_control$orthogonal_proj_tip)
    prune_graph <- ifelse(is.null(learn_graph_control$prune_graph), TRUE, learn_graph_control$prune_graph)
    ncenter <- learn_graph_control$ncenter
    scale <- ifelse(is.null(learn_graph_control$scale), FALSE, learn_graph_control$scale)
    nn.k <- ifelse(is.null(learn_graph_control[["nn.k"]]), 25, learn_graph_control[["nn.k"]])
    maxiter <- ifelse(is.null(learn_graph_control$maxiter), 10, learn_graph_control$maxiter)
    eps <- ifelse(is.null(learn_graph_control$eps), 1e-05, learn_graph_control$eps)
    L1.gamma <- ifelse(is.null(learn_graph_control$L1.gamma), 0.5, learn_graph_control$L1.gamma)
    L1.sigma <- ifelse(is.null(learn_graph_control$L1.sigma), 0.01, learn_graph_control$L1.sigma)
    assertthat::assert_that(methods::is(cds, "mmsage_cell_data_set"))
    assertthat::assert_that(reduction_method %in% c("UMAP"), msg = paste0("unsupported or invalid reduction method '", reduction_method, 
        "'"))
    assertthat::assert_that(is.logical(use_partition))
    assertthat::assert_that(is.logical(close_loop))
    assertthat::assert_that(is.logical(verbose))
    assertthat::assert_that(is.logical(orthogonal_proj_tip))
    assertthat::assert_that(is.logical(prune_graph))
    assertthat::assert_that(is.logical(scale))
    assertthat::assert_that(is.numeric(euclidean_distance_ratio))
    assertthat::assert_that(is.numeric(geodesic_distance_ratio))
    assertthat::assert_that(is.numeric(minimal_branch_len))
    if (!is.null(ncenter)) {
        assertthat::assert_that(assertthat::is.count(ncenter))
    }
    assertthat::assert_that(assertthat::is.count(maxiter))
    assertthat::assert_that(assertthat::is.count(nn.k))
    assertthat::assert_that(is.numeric(eps))
    assertthat::assert_that(is.numeric(L1.sigma))
    assertthat::assert_that(is.numeric(L1.sigma))
    assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]), msg = paste("No dimensionality reduction for", 
        reduction_method, "calculated.", "Please run reduce_dimension with", "reduction_method =", reduction_method, "and cluster_cells before running", 
        "learn_graph."))
    assertthat::assert_that(!is.null(cds@clusters[[reduction_method]]), msg = paste("No cell clusters for", reduction_method, 
        "calculated.", "Please run cluster_cells with", "reduction_method =", reduction_method, "before running learn_graph."))
    nn_control <- list()
    if (!is.null(learn_graph_control[["nn.method"]])) 
        nn_control[["method"]] <- learn_graph_control[["nn.method"]]
    if (!is.null(learn_graph_control[["nn.metric"]])) 
        nn_control[["metric"]] <- learn_graph_control[["nn.metric"]]
    if (!is.null(learn_graph_control[["nn.n_trees"]])) 
        nn_control[["n_trees"]] <- learn_graph_control[["nn.n_trees"]]
    if (!is.null(learn_graph_control[["nn.search_k"]])) 
        nn_control[["search_k"]] <- learn_graph_control[["nn.search_k"]]
    if (!is.null(learn_graph_control[["nn.M"]])) 
        nn_control[["M"]] <- learn_graph_control[["nn.M"]]
    if (!is.null(learn_graph_control[["nn.ef_construction"]])) 
        nn_control[["ef_construction"]] <- learn_graph_control[["nn.ef_construction"]]
    if (!is.null(learn_graph_control[["nn.ef"]])) 
        nn_control[["ef"]] <- learn_graph_control[["nn.ef"]]
    if (!is.null(learn_graph_control[["nn.grain_size"]])) 
        nn_control[["grain_size"]] <- learn_graph_control[["nn.grain_size"]]
    if (!is.null(learn_graph_control[["nn.cores"]])) 
        nn_control[["cores"]] <- learn_graph_control[["nn.cores"]]
    if (verbose) 
        report_nn_control("nn_control: ", nn_control)
    nn_control_default <- get_global_variable("nn_control_annoy_euclidean")
    nn_control <- set_nn_control(mode = 3, nn_control = nn_control, nn_control_default = nn_control_default, nn_index = NULL, 
        k = nn.k, verbose = verbose)
    if (use_partition) {
        partition_list <- cds@clusters[[reduction_method]]$partitions
    }
    else {
        partition_list <- rep(1, nrow(colData(cds)))
    }
    multi_tree_DDRTree_res <- tryCatch(multi_component_RGE(cds, scale = scale, reduction_method = reduction_method, partition_list = partition_list, 
        irlba_pca_res = SingleCellExperiment::reducedDims(cds)[[reduction_method]], max_components = max_components, ncenter = ncenter, 
        nn.k = nn.k, nn_control = nn_control, maxiter = maxiter, eps = eps, L1.gamma = L1.gamma, L1.sigma = L1.sigma, close_loop = close_loop, 
        euclidean_distance_ratio = euclidean_distance_ratio, geodesic_distance_ratio = geodesic_distance_ratio, prune_graph = prune_graph, 
        minimal_branch_len = minimal_branch_len, verbose = verbose), error = function(c) {
        stop(paste0(trimws(c), "\n* error in learn_graph"))
    })
    rge_res_W <- multi_tree_DDRTree_res$ddrtree_res_W
    rge_res_Z <- multi_tree_DDRTree_res$ddrtree_res_Z
    rge_res_Y <- multi_tree_DDRTree_res$ddrtree_res_Y
    cds <- multi_tree_DDRTree_res$cds
    dp_mst <- multi_tree_DDRTree_res$dp_mst
    principal_graph(cds)[[reduction_method]] <- dp_mst
    cds@principal_graph_aux[[reduction_method]]$dp_mst <- rge_res_Y
    cds <- project2MST(cds, project_point_to_line_segment, orthogonal_proj_tip, verbose, reduction_method, rge_res_Y)
    cds
}

`order_cells` <- function (cds, reduction_method = "UMAP", root_pr_nodes = NULL, root_cells = NULL, verbose = FALSE) 
{
    assertthat::assert_that(methods::is(cds, "mmsage_cell_data_set"))
    assertthat::assert_that(assertthat::are_equal("UMAP", reduction_method), msg = paste("Currently only 'UMAP' is accepted as a", 
        "reduction_method."))
    assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]), msg = paste0("No dimensionality reduction for ", 
        reduction_method, " calculated. ", "Please run reduce_dimension with ", "reduction_method = ", reduction_method, 
        ", cluster_cells, and learn_graph ", "before running order_cells."))
    assertthat::assert_that(!is.null(cds@clusters[[reduction_method]]), msg = paste("No cell clusters for", reduction_method, 
        "calculated.", "Please run cluster_cells with", "reduction_method =", reduction_method, "and run learn_graph before running", 
        "order_cells."))
    assertthat::assert_that(!is.null(principal_graph(cds)[[reduction_method]]), msg = paste("No principal graph for", reduction_method, 
        "calculated.", "Please run learn_graph with", "reduction_method =", reduction_method, "before running order_cells."))
    assertthat::assert_that(igraph::vcount(principal_graph(cds)[[reduction_method]]) < 10000, msg = paste("principal graph is too large. order_cells doesn't support", 
        "more than 10 thousand centroids."))
    if (!is.null(root_pr_nodes)) {
        assertthat::assert_that(all(root_pr_nodes %in% igraph::V(principal_graph(cds)[[reduction_method]])$name), msg = paste("All provided root_pr_nodes must be present in the", 
            "principal graph."))
    }
    if (!is.null(root_cells)) {
        assertthat::assert_that(all(root_cells %in% row.names(colData(cds))), msg = paste("All provided root_cells must be", 
            "present in the cell data set."))
    }
    if (is.null(root_cells) & is.null(root_pr_nodes)) {
        assertthat::assert_that(interactive(), msg = paste("When not in interactive mode, either", "root_pr_nodes or root_cells", 
            "must be provided."))
    }
    assertthat::assert_that(!all(c(!is.null(root_cells), !is.null(root_pr_nodes))), msg = paste("Please specify either root_pr_nodes", 
        "or root_cells, not both."))
    if (is.null(root_pr_nodes) & is.null(root_cells)) {
        if (interactive()) {
            root_pr_nodes <- select_trajectory_roots(cds, reduction_method = reduction_method)
            if (length(root_pr_nodes) == 0) {
                stop("No root node was chosen!")
            }
        }
    }
    else if (!is.null(root_cells)) {
        closest_vertex <- cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_closest_vertex
        root_pr_nodes <- unique(paste("Y_", closest_vertex[root_cells, ], sep = ""))
    }
    cds@principal_graph_aux[[reduction_method]]$root_pr_nodes <- root_pr_nodes
    cc_ordering <- extract_general_graph_ordering(cds, root_pr_nodes, verbose, reduction_method)
    cds@principal_graph_aux[[reduction_method]]$pseudotime <- cc_ordering[row.names(colData(cds)), ]$pseudo_time
    names(cds@principal_graph_aux[[reduction_method]]$pseudotime) <- row.names(colData(cds))
    cds
}

`size_factors<-` <- function (cds, value) 
{
    stopifnot(methods::is(cds, "mmsage_cell_data_set"))
    colData(cds)$Size_Factor <- value
    methods::validObject(cds)
    cds
}

`is_sparse_matrix` <- function (x) 
{
    any(class(x) %in% c("dgCMatrix", "dgTMatrix", "lgCMatrix", "CsparseMatrix"))
}

`get_matrix_info` <- function (mat) 
{
    matrix_info <- tryCatch(get_matrix_class(mat = mat), error = function(c) {
        stop(paste0(trimws(c), "\n* error in get_matrix_info"))
    })
    if (is.null(matrix_info[["matrix_class"]])) {
        stop("get_matrix_info: unable to infer matrix_class")
    }
    if (matrix_info[["matrix_class"]] != "BPCells") {
        return(matrix_info)
    }
    bmat <- bpcells_find_base_matrix(mat = mat)
    if (!(class(bmat) %in% c("UnpackedMatrixMem_uint32_t", "UnpackedMatrixMem_float", "UnpackedMatrixMem_double", "PackedMatrixMem_uint32_t", 
        "PackedMatrixMem_float", "PackedMatrixMem_double", "MatrixDir", "Iterable_dgCMatrix_wrapper"))) {
        stop("get_matrix_info: unrecognized BPCells matrix class \"", class(mat), "\"")
        return(NULL)
    }
    matrix_info[["matrix_class"]] <- "BPCells"
    if (class(bmat) == "Iterable_dgCMatrix_wrapper") {
        matrix_info[["matrix_mode"]] <- "dgCMatrix"
    }
    else if (class(bmat) == "UnpackedMatrixMem_uint32_t") {
        matrix_info[["matrix_mode"]] <- "mem"
        matrix_info[["matrix_type"]] <- "uint32_t"
        matrix_info[["matrix_compress"]] <- FALSE
    }
    else if (class(bmat) == "UnpackedMatrixMem_float") {
        matrix_info[["matrix_mode"]] <- "mem"
        matrix_info[["matrix_type"]] <- "float"
        matrix_info[["matrix_compress"]] <- FALSE
    }
    else if (class(bmat) == "UnpackedMatrixMem_double") {
        matrix_info[["matrix_mode"]] <- "mem"
        matrix_info[["matrix_type"]] <- "double"
        matrix_info[["matrix_compress"]] <- FALSE
    }
    else if (class(bmat) == "PackedMatrixMem_uint32_t") {
        matrix_info[["matrix_mode"]] <- "mem"
        matrix_info[["matrix_type"]] <- "uint32_t"
        matrix_info[["matrix_compress"]] <- TRUE
    }
    else if (class(bmat) == "PackedMatrixMem_float") {
        matrix_info[["matrix_mode"]] <- "mem"
        matrix_info[["matrix_type"]] <- "float"
        matrix_info[["matrix_compress"]] <- TRUE
    }
    else if (class(bmat) == "PackedMatrixMem_double") {
        matrix_info[["matrix_mode"]] <- "mem"
        matrix_info[["matrix_type"]] <- "double"
        matrix_info[["matrix_compress"]] <- TRUE
    }
    else if (class(bmat) == "MatrixDir") {
        matrix_info[["matrix_mode"]] <- "dir"
        matrix_info[["matrix_type"]] <- bmat@type
        matrix_info[["matrix_compress"]] <- bmat@compressed
        matrix_info[["matrix_buffer_size"]] <- bmat@buffer_size
        matrix_info[["matrix_path"]] <- bmat@dir
    }
    return(matrix_info)
}

`size_factors` <- function (cds) 
{
    stopifnot(methods::is(cds, "mmsage_cell_data_set"))
    sf <- colData(cds)$Size_Factor
    names(sf) <- colnames(SingleCellExperiment::counts(cds))
    sf
}

`estimate_sf_sparse` <- function (counts, round_exprs = TRUE, method = "mean-geometric-mean-total") 
{
    if (round_exprs) 
        counts <- round(counts)
    if (method == "mean-geometric-mean-total") {
        cell_total <- Matrix::colSums(counts)
        sfs <- cell_total/exp(mean(log(cell_total)))
    }
    else if (method == "mean-geometric-mean-log-total") {
        cell_total <- Matrix::colSums(counts)
        sfs <- log(cell_total)/exp(mean(log(log(cell_total))))
    }
    sfs[is.na(sfs)] <- 1
    sfs
}

`estimate_sf_dense` <- function (counts, round_exprs = TRUE, method = "mean-geometric-mean-total") 
{
    CM <- counts
    if (round_exprs) 
        CM <- round(CM)
    if (method == "mean-geometric-mean-log-total") {
        cell_total <- apply(CM, 2, sum)
        sfs <- log(cell_total)/exp(mean(log(log(cell_total))))
    }
    else if (method == "mean-geometric-mean-total") {
        cell_total <- apply(CM, 2, sum)
        sfs <- cell_total/exp(mean(log(cell_total)))
    }
    sfs[is.na(sfs)] <- 1
    sfs
}

`set_nn_control` <- function (mode, nn_control = list(), nn_control_default = list(), nn_index = NULL, k = NULL, verbose = FALSE) 
{
    default_method <- "annoy"
    default_metric <- "euclidean"
    default_k <- 25
    default_n_trees <- 50
    default_M <- 48
    default_ef_construction <- 200
    default_ef <- 150
    default_grain_size <- 1
    default_cores <- 1
    default_annoy_random_seed <- 42
    assertthat::assert_that(methods::is(nn_control, "list"))
    assertthat::assert_that(methods::is(nn_control_default, "list"))
    allowed_control_parameters <- c("method", "metric", "n_trees", "search_k", "M", "ef_construction", "ef", "grain_size", 
        "cores", "show_values", "annoy_random_seed")
    assertthat::assert_that(all(names(nn_control) %in% allowed_control_parameters), msg = "set_nn_control: unknown variable in nn_control")
    assertthat::assert_that(all(names(nn_control_default) %in% allowed_control_parameters), msg = "set_nn_control: unknown variable in nn_control_default")
    assertthat::assert_that(assertthat::is.count(mode) && mode >= 1 && mode <= 3, msg = paste0("set_nn_control: invalid mode value. Mode must be an integer with the value 1, 2, or 3."))
    nn_control_out <- list()
    nn_control_out[["method"]] <- select_nn_parameter_value("method", nn_control, nn_control_default, default_method)
    if (nn_control_out[["method"]] == "nn2") {
    }
    else if (nn_control_out[["method"]] == "annoy") {
        nn_control_out[["metric"]] <- select_nn_parameter_value("metric", nn_control, nn_control_default, default_metric)
        assertthat::assert_that(nn_control_out[["metric"]] %in% c("euclidean", "cosine", "manhattan", "hamming"), msg = paste0("set_nn_control: nearest neighbor metric for annoy must be one of 'euclidean', 'cosine', 'manhattan', or 'hamming'"))
        if (bitwAnd(mode, 1)) {
            nn_control_out[["n_trees"]] <- select_nn_parameter_value("n_trees", nn_control, nn_control_default, default_n_trees)
            nn_control_out[["annoy_random_seed"]] <- select_nn_parameter_value("annoy_random_seed", nn_control, nn_control_default, 
                default_annoy_random_seed)
            assertthat::assert_that(assertthat::is.count(nn_control_out[["n_trees"]]))
        }
        if (bitwAnd(mode, 2)) {
            nn_control_out[["search_k"]] <- tryCatch(select_annoy_search_k(mode, nn_control, nn_control_default, nn_index, 
                k, default_n_trees, default_k), error = function(c) {
                stop(paste0(trimws(c), "\n* error in set_nn_control"))
            })
            nn_control_out[["grain_size"]] <- select_nn_parameter_value("grain_size", nn_control, nn_control_default, default_grain_size)
            nn_control_out[["cores"]] <- select_nn_parameter_value("cores", nn_control, nn_control_default, default_cores)
            assertthat::assert_that(assertthat::is.count(nn_control_out[["search_k"]]))
            assertthat::assert_that(assertthat::is.count(nn_control_out[["grain_size"]]))
            assertthat::assert_that(assertthat::is.count(nn_control_out[["cores"]]))
        }
    }
    else if (nn_control_out[["method"]] == "hnsw") {
        nn_control_out[["metric"]] <- select_nn_parameter_value("metric", nn_control, nn_control_default, default_metric)
        assertthat::assert_that(nn_control_out[["metric"]] %in% c("euclidean", "l2", "cosine", "ip"), msg = paste0("set_nn_control: nearest neighbor metric for HNSW must be one of 'euclidean', 'l2', 'cosine', or 'ip'"))
        nn_control_out[["grain_size"]] <- select_nn_parameter_value("grain_size", nn_control, nn_control_default, default_grain_size)
        nn_control_out[["cores"]] <- select_nn_parameter_value("cores", nn_control, nn_control_default, default_cores)
        assertthat::assert_that(assertthat::is.count(nn_control_out[["grain_size"]]))
        assertthat::assert_that(assertthat::is.count(nn_control_out[["cores"]]))
        if (bitwAnd(mode, 1)) {
            nn_control_out[["M"]] <- select_nn_parameter_value("M", nn_control, nn_control_default, default_M)
            nn_control_out[["ef_construction"]] <- select_nn_parameter_value("ef_construction", nn_control, nn_control_default, 
                default_ef_construction)
            assertthat::assert_that(assertthat::is.count(nn_control_out[["M"]]))
            assertthat::assert_that(assertthat::is.count(nn_control_out[["ef_construction"]]))
            assertthat::assert_that(nn_control_out[["M"]] >= 2, msg = paste0("set_nn_control: HNSW nearest neighbor M index build parameter must be >= 2"))
        }
        if (bitwAnd(mode, 2)) {
            assertthat::assert_that(assertthat::is.count(k), msg = paste0("set_nn_control: parameter k must be set for method='hnsw' and modes 2 and 3."))
            nn_control_out[["ef"]] <- select_nn_parameter_value("ef", nn_control, nn_control_default, default_ef)
            assertthat::assert_that(assertthat::is.count(nn_control_out[["ef"]]))
            assertthat::assert_that(nn_control_out[["ef"]] >= k, msg = paste0("set_nn_control: HNSW nearest neighbor ef index search parameter must be >= k (", 
                k, ")"))
        }
    }
    else stop("set_nn_control: unsupported nearest neighbor method '", nn_control_out[["method"]], "'")
    if (verbose) {
        cs <- get_call_stack_as_string()
        message("set_nn_control: call stack: ", cs)
        report_nn_control("  nn_control: ", nn_control_out)
    }
    if (!is.null(nn_control[["show_values"]]) && nn_control[["show_values"]] == TRUE) {
        report_nn_control("  nn_control: ", nn_control = nn_control_out)
        stop_no_noise()
    }
    return(nn_control_out)
}

`get_global_variable` <- function (variable_name = NULL) 
{
    value <- tryCatch({
        v <- get("guard_element", envir = ._._global_variable_env_._.)
        if (v != "sanity_check") 
            stop()
        v
    }, error = function(msg) {
        message("Global variable storage is compromised.")
        return(NA)
    })
    if (is.na(value)) {
        return(NA)
    }
    if (!is.null(variable_name)) {
        value <- tryCatch({
            get(variable_name, envir = ._._global_variable_env_._.)
        }, error = function(msg) {
            message("'", variable_name, "'", " is not a global variable.")
            return(NA)
        })
    }
    else {
        value = list()
        variable_names <- ls(envir = ._._global_variable_env_._.)
        for (variable_name in variable_names) {
            value[[variable_name]] <- get(variable_name, envir = ._._global_variable_env_._.)
        }
    }
    return(value)
}

`normalize_expr_data` <- function (FM, size_factors = NULL, norm_method = c("log", "size_only", "none"), pseudo_count = NULL) 
{
    assertthat::assert_that(!is.null(size_factors))
    assertthat::assert_that(length(size_factors) == ncol(FM))
    norm_method <- match.arg(norm_method)
    if (is.null(pseudo_count)) {
        if (norm_method == "log") 
            pseudo_count <- 1
        else pseudo_count <- 0
    }
    if (!is(FM, "IterableMatrix")) {
        if (norm_method == "log") {
            FM <- Matrix::t(Matrix::t(FM)/size_factors)
            if (pseudo_count != 1 || is_sparse_matrix(FM) == FALSE) {
                FM <- FM + pseudo_count
                FM <- log2(FM)
            }
            else {
                FM@x = log2(FM@x + 1)
            }
        }
        else if (norm_method == "size_only") {
            FM <- Matrix::t(Matrix::t(FM)/size_factors)
            FM <- FM + pseudo_count
        }
    }
    else {
        if (norm_method == "log") {
            FM <- BPCells::t(BPCells::t(FM)/size_factors)
            if (pseudo_count == 1) {
                FM <- log1p(FM)/log(2)
            }
            else {
                FM <- log1p(FM + pseudo_count - 1)/log(2)
            }
        }
        else if (norm_method == "size_only") {
            FM <- BPCells::t(BPCells::t(FM)/size_factors)
            FM <- FM + pseudo_count
        }
    }
    return(FM)
}

`initialize_reduce_dim_metadata` <- function (cds, reduction_method = c("PCA", "LSI", "Aligned", "tSNE", "UMAP")) 
{
    assertthat::assert_that(methods::is(cds, "mmsage_cell_data_set"), msg = paste("cds parameter is not a cell_data_set"))
    assertthat::assert_that(tryCatch(expr = ifelse(match.arg(reduction_method) == "", TRUE, TRUE), error = function(e) FALSE), 
        msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
    reduction_method <- match.arg(reduction_method)
    if (is.null(SingleCellExperiment::int_metadata(cds))) {
        SingleCellExperiment::int_metadata(cds) <- list()
    }
    if (is.null(SingleCellExperiment::int_metadata(cds)[["reduce_dim_metadata"]])) {
        SingleCellExperiment::int_metadata(cds)[["reduce_dim_metadata"]] <- list()
    }
    SingleCellExperiment::int_metadata(cds)[["reduce_dim_metadata"]][[reduction_method]] <- list()
    SingleCellExperiment::int_metadata(cds)[["reduce_dim_metadata"]][[reduction_method]][["identity"]] <- list()
    return(cds)
}

`initialize_reduce_dim_model_identity` <- function (cds, reduction_method = c("PCA", "LSI", "Aligned", "tSNE", "UMAP")) 
{
    assertthat::assert_that(methods::is(cds, "mmsage_cell_data_set"), msg = paste("cds parameter is not a cell_data_set"))
    assertthat::assert_that(tryCatch(expr = ifelse(match.arg(reduction_method) == "", TRUE, TRUE), error = function(e) FALSE), 
        msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
    reduction_method <- match.arg(reduction_method)
    if (is.null(cds@reduce_dim_aux[[reduction_method]])) {
        cds@reduce_dim_aux[[reduction_method]] <- S4Vectors::SimpleList()
    }
    cds@reduce_dim_aux[[reduction_method]][["model"]] <- S4Vectors::SimpleList()
    cds@reduce_dim_aux[[reduction_method]][["model"]][["identity"]] <- S4Vectors::SimpleList()
    return(cds)
}

`sparse_prcomp_irlba` <- function (x, n = 3, retx = TRUE, center = TRUE, scale. = FALSE, verbose = FALSE, ...) 
{
    if (verbose) {
        message("pca: sparse_prcomp_irlba: matrix class: ", class(x))
        message(paste0(show_matrix_info(matrix_info = get_matrix_info(mat = x), indent = "  ")), appendLF = FALSE)
        message()
    }
    a <- names(as.list(match.call()))
    ans <- list(scale = scale.)
    if ("tol" %in% a) 
        warning("The `tol` truncation argument from `prcomp` is not supported by\n            `prcomp_irlba`. If specified, `tol` is passed to the `irlba`\n            function to control that algorithm's convergence tolerance. See\n            `?prcomp_irlba` for help.")
    orig_x <- x
    if (!methods::is(x, "DelayedMatrix")) {
        x = DelayedArray::DelayedArray(x)
    }
    args <- list(A = orig_x, nv = n)
    if (is.logical(center)) {
        if (center) 
            args$center <- DelayedMatrixStats::colMeans2(x)
    }
    else args$center <- center
    if (is.logical(scale.)) {
        if (is.numeric(args$center)) {
            scale. <- sqrt(DelayedMatrixStats::colVars(x))
            if (ans$scale) 
                ans$totalvar <- ncol(x)
            else ans$totalvar <- sum(scale.^2)
        }
        else {
            if (ans$scale) {
                scale. <- sqrt(DelayedMatrixStats::colSums2(x^2)/(max(1, nrow(x) - 1L)))
                ans$totalvar <- sum(sqrt(DelayedMatrixStats::colSums2(t(t(x)/scale.)^2)/(nrow(x) - 1L)))
            }
            else {
                ans$totalvar <- sum(DelayedMatrixStats::colSums2(x^2)/(nrow(x) - 1L))
            }
        }
        if (ans$scale) 
            args$scale <- scale.
    }
    else {
        args$scale <- scale.
        ans$totalvar <- sum(sqrt(DelayedMatrixStats::colSums2(t(t(x)/scale.)^2)/(nrow(x) - 1L)))
    }
    if (!missing(...)) 
        args <- c(args, list(...))
    if (verbose) {
        message("start irlba: ", Sys.time())
    }
    s <- do.call(irlba::irlba, args = args)
    if (verbose) {
        message("end irlba: ", Sys.time())
        message()
    }
    if (verbose) {
        message("singular values (head)")
        message(paste(head(s$d), collapse = " "))
        message()
    }
    ans$sdev <- s$d/sqrt(max(1, nrow(x) - 1))
    ans$rotation <- s$v
    colnames(ans$rotation) <- paste("PC", seq(1, ncol(ans$rotation)), sep = "")
    ans$center <- args$center
    ans$svd_scale <- args$scale
    if (retx) {
        ans <- c(ans, list(x = sweep(s$u, 2, s$d, FUN = `*`)))
        colnames(ans$x) <- paste("PC", seq(1, ncol(ans$rotation)), sep = "")
    }
    if (verbose) {
        message("umat: ", paste(dim(s$u), collapse = " "))
        message("vtmat: ", paste(dim(s$v), collapse = " "))
    }
    class(ans) <- c("irlba_prcomp", "prcomp")
    ans
}

`show_matrix_info` <- function (matrix_info, indent = "") 
{
    message("matrix_info:")
    if (!is.null(matrix_info[["matrix_class"]])) {
        message(paste0(indent, "class:       ", matrix_info[["matrix_class"]]))
    }
    if (!is.null(matrix_info[["matrix_mode"]])) {
        message(paste0(indent, "mode:        ", matrix_info[["matrix_mode"]]))
    }
    if (!is.null(matrix_info[["matrix_type"]])) {
        message(paste0(indent, "type:        ", matrix_info[["matrix_type"]]))
    }
    if (!is.null(matrix_info[["matrix_path"]])) {
        message(paste0(indent, "path:        ", matrix_info[["matrix_path"]]))
    }
    if (!is.null(matrix_info[["matrix_compress"]])) {
        message(paste0(indent, "compress:    ", matrix_info[["matrix_compress"]]))
    }
    if (!is.null(matrix_info[["matrix_buffer_size"]])) {
        message(paste0(indent, "buffer_size: ", matrix_info[["matrix_buffer_size"]]))
    }
}

`get_unique_id` <- function (object = NULL) 
{
    if (!is.null(object)) {
        object_dim <- dim(object)
        if (!methods::is(object, "IterableMatrix")) {
            object_checksum <- digest::digest(object)
        }
        else {
            object_checksum <- BPCells::checksum(object)
        }
        if (!is.null(object_dim)) 
            object_id <- list(checksum = object_checksum, dim = object_dim)
        else object_id <- list(checksum = object_checksum, dim = length(object))
    }
    else {
        id_count <- get_global_variable("id_count")
        rtime <- as.numeric(Sys.time()) * 1e+05 + id_count
        object_id <- openssl::md5(as.character(rtime))
        id_count <- id_count + 1
        set_global_variable("id_count", id_count)
    }
    return(object_id)
}

`get_counts_identity` <- function (cds) 
{
    assertthat::assert_that(methods::is(cds, "mmsage_cell_data_set"), msg = paste("cds parameter is not a cell_data_set"))
    if (is.null(SingleCellExperiment::int_metadata(cds)[["counts_metadata"]])) {
        initialize_counts_metadata(cds)
    }
    if (!is.null(SingleCellExperiment::int_metadata(cds)[["counts_metadata"]][["identity"]][["matrix_id"]])) {
        matrix_id <- SingleCellExperiment::int_metadata(cds)[["counts_metadata"]][["identity"]][["matrix_id"]]
    }
    else {
        matrix_id <- "none"
    }
    if (!is.null(SingleCellExperiment::int_metadata(cds)[["counts_metadata"]][["identity"]][["matrix_type"]])) {
        matrix_type <- SingleCellExperiment::int_metadata(cds)[["counts_metadata"]][["identity"]][["matrix_type"]]
    }
    else {
        matrix_type <- "matrix:counts"
    }
    return(list(matrix_id = matrix_id, matrix_type = matrix_type))
}

`set_reduce_dim_matrix_identity` <- function (cds, reduction_method = c("PCA", "LSI", "Aligned", "tSNE", "UMAP"), matrix_type, matrix_id, prev_matrix_type, prev_matrix_id, 
    model_type, model_id) 
{
    assertthat::assert_that(methods::is(cds, "mmsage_cell_data_set"), msg = paste("cds parameter is not a cell_data_set"))
    assertthat::assert_that(tryCatch(expr = ifelse(match.arg(reduction_method) == "", TRUE, TRUE), error = function(e) FALSE), 
        msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
    reduction_method <- match.arg(reduction_method)
    if (is.null(SingleCellExperiment::int_metadata(cds)[["reduce_dim_metadata"]][[reduction_method]])) {
        cds <- initialize_reduce_dim_metadata(cds = cds, reduction_method = reduction_method)
    }
    SingleCellExperiment::int_metadata(cds)[["reduce_dim_metadata"]][[reduction_method]][["identity"]][["matrix_type"]] <- matrix_type
    SingleCellExperiment::int_metadata(cds)[["reduce_dim_metadata"]][[reduction_method]][["identity"]][["matrix_id"]] <- matrix_id
    SingleCellExperiment::int_metadata(cds)[["reduce_dim_metadata"]][[reduction_method]][["identity"]][["prev_matrix_type"]] <- prev_matrix_type
    SingleCellExperiment::int_metadata(cds)[["reduce_dim_metadata"]][[reduction_method]][["identity"]][["prev_matrix_id"]] <- prev_matrix_id
    SingleCellExperiment::int_metadata(cds)[["reduce_dim_metadata"]][[reduction_method]][["identity"]][["model_type"]] <- model_type
    SingleCellExperiment::int_metadata(cds)[["reduce_dim_metadata"]][[reduction_method]][["identity"]][["model_id"]] <- model_id
    return(cds)
}

`set_reduce_dim_model_identity` <- function (cds, reduction_method = c("PCA", "LSI", "Aligned", "tSNE", "UMAP"), model_type, model_id, prev_model_type, prev_model_id, 
    model_path = "none") 
{
    assertthat::assert_that(methods::is(cds, "mmsage_cell_data_set"), msg = paste("cds parameter is not a cell_data_set"))
    assertthat::assert_that(tryCatch(expr = ifelse(match.arg(reduction_method) == "", TRUE, TRUE), error = function(e) FALSE), 
        msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
    reduction_method <- match.arg(reduction_method)
    if (is.null(cds@reduce_dim_aux[[reduction_method]][["model"]][["identity"]])) {
        cds <- initialize_reduce_dim_model_identity(cds = cds, reduction_method = reduction_method)
    }
    cds@reduce_dim_aux[[reduction_method]][["model"]][["identity"]][["model_type"]] <- model_type
    cds@reduce_dim_aux[[reduction_method]][["model"]][["identity"]][["model_id"]] <- model_id
    cds@reduce_dim_aux[[reduction_method]][["model"]][["identity"]][["prev_model_type"]] <- prev_model_type
    cds@reduce_dim_aux[[reduction_method]][["model"]][["identity"]][["prev_model_id"]] <- prev_model_id
    cds@reduce_dim_aux[[reduction_method]][["model"]][["identity"]][["model_path"]] <- model_path
    global_variable_name <- list(PCA = "reduce_dim_pca_model_version", LSI = "reduce_dim_lsi_model_version", Aligned = "reduce_dim_aligned_model_version", 
        tSNE = "reduce_dim_tsne_model_version", UMAP = "reduce_dim_umap_model_version")
    cds@reduce_dim_aux[[reduction_method]][["model"]][["identity"]][["model_version"]] <- get_global_variable(global_variable_name[[reduction_method]])
    return(cds)
}

`make_nn_index` <- function (subject_matrix, nn_control = list(), verbose = FALSE) 
{
    assertthat::assert_that(methods::is(subject_matrix, "matrix") || is_sparse_matrix(subject_matrix), msg = paste0("make_nn_matrix: the subject_matrix object must be of type matrix"))
    nn_control_default <- get_global_variable("nn_control_annoy_euclidean")
    nn_control <- set_nn_control(mode = 1, nn_control = nn_control, nn_control_default = nn_control_default, nn_index = NULL, 
        k = NULL, verbose = verbose)
    if (verbose) {
        message("make_nn_index:")
        report_nn_control("  nn_control: ", nn_control)
        tick("make_nn_index: build time")
    }
    nn_method <- nn_control[["method"]]
    metric = nn_control[["metric"]]
    num_row <- nrow(subject_matrix)
    num_col <- ncol(subject_matrix)
    if (!is.null(rownames(subject_matrix))) 
        checksum_rownames <- digest::digest(sort(rownames(subject_matrix)))
    else checksum_rownames <- NA_character_
    if (nn_method == "nn2") {
        stop("make_nn_index is not valid for method nn2")
    }
    else if (nn_method == "annoy") {
        monocle3_annoy_index_version <- get_global_variable("monocle3_annoy_index_version")
        annoy_index <- tryCatch(new_annoy_index(metric, num_col), error = function(c) {
            stop(paste0(trimws(c), "\n* error in make_nn_index"))
        })
        annoy_random_seed <- nn_control[["annoy_random_seed"]]
        annoy_index$setSeed(annoy_random_seed)
        n_trees <- nn_control[["n_trees"]]
        if (num_row > 0) {
            for (i in 1:num_row) annoy_index$addItem(i - 1, subject_matrix[i, ])
            annoy_index$build(n_trees)
        }
        annoy_index_version <- packageVersion("RcppAnnoy")
        nn_index <- list(method = "annoy", annoy_index = annoy_index, version = monocle3_annoy_index_version, annoy_index_version = annoy_index_version, 
            metric = metric, n_trees = n_trees, nrow = num_row, ncol = num_col, checksum_rownames = checksum_rownames, annoy_random_seed = annoy_random_seed)
    }
    else if (nn_method == "hnsw") {
        monocle3_hnsw_index_version <- get_global_variable("monocle3_hnsw_index_version")
        M <- nn_control[["M"]]
        ef_construction <- nn_control[["ef_construction"]]
        hnsw_index <- RcppHNSW::hnsw_build(X = subject_matrix, distance = metric, M = M, ef = ef_construction, verbose = verbose, 
            n_threads = nn_control[["cores"]], grain_size = nn_control[["grain_size"]])
        hnsw_index_version <- packageVersion("RcppHNSW")
        nn_index <- list(method = "hnsw", hnsw_index = hnsw_index, version = monocle3_hnsw_index_version, hnsw_index_version = hnsw_index_version, 
            metric = metric, M = M, ef_construction = ef_construction, nrow = num_row, ncol = num_col, checksum_rownames = checksum_rownames)
    }
    else stop("make_nn_index: unsupported nearest neighbor index type '", nn_method, "'")
    if (verbose) {
        tock()
    }
    return(nn_index)
}

`set_cds_nn_index` <- function (cds, reduction_method = c("UMAP", "PCA", "LSI", "Aligned", "tSNE"), nn_index, verbose = FALSE) 
{
    assertthat::assert_that(methods::is(cds, "mmsage_cell_data_set"), msg = paste("cds parameter is not a cell_data_set"))
    assertthat::assert_that(tryCatch(expr = ifelse(match.arg(reduction_method) == "", TRUE, TRUE), error = function(e) FALSE), 
        msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
    reduction_method <- match.arg(reduction_method)
    assertthat::assert_that(!is.null(SingleCellExperiment::reducedDims(cds)[[reduction_method]]), msg = paste0("When reduction_method = '", 
        reduction_method, "' the cds must have been processed for it.", " Please run the required processing function", " before this one."))
    nn_method <- nn_index[["method"]]
    if (nn_method == "nn2") {
        stop("set_cds_nn_index is not valid for method nn2")
    }
    else if (nn_method == "annoy") {
        cds@reduce_dim_aux[[reduction_method]][["nn_index"]][[nn_method]] <- S4Vectors::SimpleList()
        cds@reduce_dim_aux[[reduction_method]][["nn_index"]][[nn_method]][["nn_index"]] <- nn_index
        cds@reduce_dim_aux[[reduction_method]][["nn_index"]][[nn_method]][["matrix_id"]] <- get_reduce_dim_matrix_identity(cds, 
            reduction_method)[["matrix_id"]]
    }
    else if (nn_method == "hnsw") {
        cds@reduce_dim_aux[[reduction_method]][["nn_index"]][[nn_method]] <- S4Vectors::SimpleList()
        cds@reduce_dim_aux[[reduction_method]][["nn_index"]][[nn_method]][["nn_index"]] <- nn_index
        cds@reduce_dim_aux[[reduction_method]][["nn_index"]][[nn_method]][["matrix_id"]] <- get_reduce_dim_matrix_identity(cds, 
            reduction_method)[["matrix_id"]]
    }
    else stop("set_cds_nn_index: unsupported nearest neighbor index type '", nn_method, "'")
    return(cds)
}

`clear_cds_nn_index` <- function (cds, reduction_method = c("PCA", "LSI", "Aligned", "tSNE", "UMAP"), nn_method = c("annoy", "hnsw", "all")) 
{
    assertthat::assert_that(methods::is(cds, "mmsage_cell_data_set"), msg = paste("cds parameter is not a cell_data_set"))
    assertthat::assert_that(tryCatch(expr = ifelse(match.arg(reduction_method) == "", TRUE, TRUE), error = function(e) FALSE), 
        msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', 'UMAP'")
    reduction_method <- match.arg(reduction_method)
    assertthat::assert_that(tryCatch(expr = ifelse(match.arg(nn_method) == "", TRUE, TRUE), error = function(e) FALSE), msg = "nn_method must be one of 'annoy', 'hnsw', or 'all'")
    nn_method <- match.arg(nn_method)
    if (nn_method == "annoy") {
        if (!is.null(cds@reduce_dim_aux[[reduction_method]][["nn_index"]]) && !is.null(cds@reduce_dim_aux[[reduction_method]][["nn_index"]][[nn_method]])) {
            cds@reduce_dim_aux[[reduction_method]][["nn_index"]][[nn_method]] <- NULL
        }
    }
    else if (nn_method == "hnsw") {
        if (!is.null(cds@reduce_dim_aux[[reduction_method]][["nn_index"]]) && !is.null(cds@reduce_dim_aux[[reduction_method]][["nn_index"]][[nn_method]])) {
            cds@reduce_dim_aux[[reduction_method]][["nn_index"]][[nn_method]] <- NULL
        }
    }
    else if (nn_method == "all") {
        if (!is.null(cds@reduce_dim_aux[[reduction_method]][["nn_index"]])) {
            cds@reduce_dim_aux[[reduction_method]][["nn_index"]] <- S4Vectors::SimpleList()
        }
    }
    else stop("clear_cds_nn_index: unsupported nearest neighbor index type '", nn_method, "'")
    return(cds)
}

`tfidf` <- function (count_matrix, frequencies = TRUE, log_scale_tf = TRUE, scale_factor = 1e+05, block_size = 2e+09, iterable_matrix_flag = FALSE) 
{
    if (frequencies) {
        if (!iterable_matrix_flag) {
            col_sums <- Matrix::colSums(count_matrix)
            tf <- Matrix::t(Matrix::t(count_matrix)/col_sums)
        }
        else {
            col_sums <- BPCells::colSums(count_matrix)
            tf <- BPCells::t(BPCells::t(count_matrix)/col_sums)
        }
    }
    else {
        col_sums <- NA
        tf <- count_matrix
    }
    if (log_scale_tf) {
        if (!iterable_matrix_flag) {
            if (frequencies) {
                tf@x <- log1p(tf@x * scale_factor)
            }
            else {
                tf@x <- log1p(tf@x * 1)
            }
        }
        else {
            if (frequencies) {
                tf <- log1p(tf * scale_factor)
            }
            else {
                tf <- log1p(tf * 1)
            }
        }
    }
    num_cols <- ncol(count_matrix)
    if (!iterable_matrix_flag) {
        row_sums <- Matrix::rowSums(count_matrix > 0)
    }
    else {
        row_sums <- BPCells::rowSums(BPCells::binarize(count_matrix, threshold = 0))
    }
    idf <- log(1 + num_cols/row_sums)
    if (!iterable_matrix_flag) {
        tf_idf_counts = tryCatch({
            tf_idf_counts <- tf * idf
            tf_idf_counts
        }, error = function(e) {
            print(paste("TF*IDF multiplication too large for in-memory, falling back", "on DelayedArray."))
            options(DelayedArray.block.size = block_size)
            DelayedArray:::set_verbose_block_processing(TRUE)
            tf <- DelayedArray::DelayedArray(tf)
            idf <- as.matrix(idf)
            tf_idf_counts <- tf * idf
            tf_idf_counts
        })
    }
    else {
        tf_idf_counts <- tf * idf
    }
    rownames(tf_idf_counts) <- rownames(count_matrix)
    colnames(tf_idf_counts) <- colnames(count_matrix)
    if (!iterable_matrix_flag) {
        tf_idf_counts <- methods::as(tf_idf_counts, "sparseMatrix")
    }
    return(list(tf_idf_counts = tf_idf_counts, frequencies = frequencies, log_scale_tf = log_scale_tf, scale_factor = scale_factor, 
        col_sums = col_sums, row_sums = row_sums, num_cols = num_cols))
}

`get_reduce_dim_matrix_identity` <- function (cds, reduction_method = c("PCA", "LSI", "Aligned", "tSNE", "UMAP")) 
{
    assertthat::assert_that(methods::is(cds, "mmsage_cell_data_set"), msg = paste("cds parameter is not a cell_data_set"))
    assertthat::assert_that(tryCatch(expr = ifelse(match.arg(reduction_method) == "", TRUE, TRUE), error = function(e) FALSE), 
        msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
    reduction_method <- match.arg(reduction_method)
    if (is.null(SingleCellExperiment::int_metadata(cds)[["reduce_dim_metadata"]][[reduction_method]])) {
        cds <- initialize_reduce_dim_metadata(cds = cds, reduction_method = reduction_method)
        return(list(identity_exists = FALSE))
    }
    return(list(identity_exists = TRUE, matrix_type = SingleCellExperiment::int_metadata(cds)[["reduce_dim_metadata"]][[reduction_method]][["identity"]][["matrix_type"]], 
        matrix_id = SingleCellExperiment::int_metadata(cds)[["reduce_dim_metadata"]][[reduction_method]][["identity"]][["matrix_id"]], 
        prev_matrix_type = SingleCellExperiment::int_metadata(cds)[["reduce_dim_metadata"]][[reduction_method]][["identity"]][["prev_matrix_type"]], 
        prev_matrix_id = SingleCellExperiment::int_metadata(cds)[["reduce_dim_metadata"]][[reduction_method]][["identity"]][["prev_matrix_id"]], 
        model_type = SingleCellExperiment::int_metadata(cds)[["reduce_dim_metadata"]][[reduction_method]][["identity"]][["model_type"]], 
        model_id = SingleCellExperiment::int_metadata(cds)[["reduce_dim_metadata"]][[reduction_method]][["identity"]][["model_id"]]))
}

`get_reduce_dim_model_identity` <- function (cds, reduction_method = c("PCA", "LSI", "Aligned", "tSNE", "UMAP")) 
{
    assertthat::assert_that(methods::is(cds, "mmsage_cell_data_set"), msg = paste("cds parameter is not a cell_data_set"))
    assertthat::assert_that(tryCatch(expr = ifelse(match.arg(reduction_method) == "", TRUE, TRUE), error = function(e) FALSE), 
        msg = "reduction_method must be one of 'PCA', 'LSI', 'Aligned', 'tSNE', or 'UMAP'")
    reduction_method <- match.arg(reduction_method)
    if (is.null(cds@reduce_dim_aux[[reduction_method]][["model"]][["identity"]])) {
        cds <- initialize_reduce_dim_model_identity(cds = cds, reduction_method = reduction_method)
        return(list(identity_exists = FALSE))
    }
    return(list(identity_exists = TRUE, model_type = cds@reduce_dim_aux[[reduction_method]][["model"]][["identity"]][["model_type"]], 
        model_id = cds@reduce_dim_aux[[reduction_method]][["model"]][["identity"]][["model_id"]], prev_model_type = cds@reduce_dim_aux[[reduction_method]][["model"]][["identity"]][["prev_model_type"]], 
        prev_model_id = cds@reduce_dim_aux[[reduction_method]][["model"]][["identity"]][["prev_model_id"]], model_path = cds@reduce_dim_aux[[reduction_method]][["model"]][["identity"]][["model_path"]], 
        version = cds@reduce_dim_aux[[reduction_method]][["model"]][["identity"]][["model_version"]]))
}

`add_citation` <- function (cds, citation_key) 
{
    citation_map <- list(UMAP = c("UMAP", "McInnes, L., Healy, J. & Melville, J. UMAP: Uniform Manifold Approximation and Projection for dimension reduction. Preprint at https://arxiv.org/abs/1802.03426 (2018)."), 
        MNN_correct = c("MNN Correct", "Haghverdi, L. et. al. Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors. Nat. Biotechnol. 36, 421-427 (2018). https://doi.org/10.1038/nbt.4091"), 
        partitions = c("partitioning", c("Levine, J. H., et. al. Data-driven phenotypic dissection of AML reveals progenitor-like cells that correlate with prognosis. Cell 162, 184-197 (2015). https://doi.org/10.1016/j.cell.2015.05.047", 
            "Wolf, F. A. et. al. PAGA: graph abstraction reconciles clustering with trajectory inference through a topology preserving map of single cells. Genome Biol. 20, 59 (2019). https://doi.org/10.1186/s13059-019-1663-x")), 
        clusters = c("clustering", "Levine, J. H. et. al. Data-driven phenotypic dissection of AML reveals progenitor-like cells that correlate with prognosis. Cell 162, 184-197 (2015). https://doi.org/10.1016/j.cell.2015.05.047"), 
        leiden = c("leiden", "Traag, V.A., Waltman, L. & van Eck, N.J. From Louvain to Leiden: guaranteeing well-connected communities. Scientific Reportsvolume 9, Article number: 5233 (2019). https://doi.org/10.1038/s41598-019-41695-z"), 
        bpcells = c("BPCells", "Parks, B & Abdi, I. BPCells: Single Cell Counts Matrices to PCA. https://bnprks.github.io/BPCells"))
    if (is.null(S4Vectors::metadata(cds)$citations) | citation_key == "Monocle") {
        S4Vectors::metadata(cds)$citations <- data.frame(method = c("Monocle", "Monocle", "Monocle"), citations = c("Trapnell C. et. al. The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. Nat. Biotechnol. 32, 381-386 (2014). https://doi.org/10.1038/nbt.2859", 
            "Qiu, X. et. al. Reversed graph embedding resolves complex single-cell trajectories. Nat. Methods 14, 979-982 (2017). https://doi.org/10.1038/nmeth.4402", 
            "Cao, J. et. al. The single-cell transcriptional landscape of mammalian organogenesis. Nature 566, 496-502 (2019). https://doi.org/10.1038/s41586-019-0969-x"))
    }
    S4Vectors::metadata(cds)$citations <- rbind(S4Vectors::metadata(cds)$citations, data.frame(method = citation_map[[citation_key]][1], 
        citations = citation_map[[citation_key]][2]))
    cds
}

`louvain_clustering` <- function (data, pd, weight = FALSE, nn_index = NULL, k = 20, nn_control = list(), louvain_iter = 1, random_seed = 0L, verbose = FALSE) 
{
    assertthat::assert_that(assertthat::is.count(k))
    cell_names <- row.names(pd)
    if (!identical(cell_names, row.names(pd))) 
        stop("Phenotype and row name from the data doesn't match")
    graph_result <- tryCatch(cluster_cells_make_graph(data = data, weight = weight, cell_names = cell_names, nn_index, k = k, 
        nn_control = nn_control, verbose = verbose), error = function(c) {
        stop(paste0(trimws(c), "\n* error in louvain_clustering"))
    })
    if (verbose) 
        message("  Run louvain clustering ...")
    t_start <- Sys.time()
    Qp <- -1
    optim_res <- NULL
    best_max_resolution <- "No resolution"
    if (louvain_iter >= 2) {
        random_seed <- NULL
    }
    if (louvain_iter < 1) 
        warning("bad loop: louvain_iter is < 1")
    for (iter in 1:louvain_iter) {
        if (verbose) {
            cat("Running louvain iteration ", iter, "...\n")
        }
        Q <- igraph::cluster_louvain(graph_result[["g"]])
        if (is.null(optim_res)) {
            Qp <- max(Q$modularity)
            optim_res <- Q
        }
        else {
            Qt <- max(Q$modularity)
            if (Qt > Qp) {
                optim_res <- Q
                Qp <- Qt
            }
        }
    }
    if (verbose) 
        message("Maximal modularity is ", Qp, "; corresponding resolution is ", best_max_resolution)
    t_end <- Sys.time()
    if (verbose) {
        message("\nRun kNN based graph clustering DONE, totally takes ", t_end - t_start, " s.")
        cat("  -Number of clusters:", length(unique(igraph::membership(optim_res))), "\n")
    }
    if (igraph::vcount(graph_result[["g"]]) < 3000) {
        coord <- NULL
        edge_links <- NULL
    }
    else {
        coord <- NULL
        edge_links <- NULL
    }
    igraph::V(graph_result[["g"]])$names <- as.character(igraph::V(graph_result[["g"]]))
    return(list(g = graph_result[["g"]], relations = graph_result[["relations"]], distMatrix = graph_result[["distMatrix"]], 
        coord = coord, edge_links = edge_links, optim_res = optim_res))
}

`compute_partitions` <- function (g, optim_res, qval_thresh = 0.05, verbose = FALSE) 
{
    cell_membership <- as.factor(igraph::membership(optim_res))
    membership_matrix <- Matrix::sparse.model.matrix(~cell_membership + 0)
    num_links <- Matrix::t(membership_matrix) %*% igraph::as_adjacency_matrix(g) %*% membership_matrix
    diag(num_links) <- 0
    louvain_modules <- levels(cell_membership)
    edges_per_module <- Matrix::rowSums(num_links)
    total_edges <- sum(num_links)
    theta <- (as.matrix(edges_per_module)/total_edges) %*% Matrix::t(edges_per_module/total_edges)
    var_null_num_links <- theta * (1 - theta)/total_edges
    num_links_ij <- num_links/total_edges - theta
    cluster_mat <- pnorm_over_mat(as.matrix(num_links_ij), var_null_num_links)
    num_links <- num_links_ij/total_edges
    num_links[is.nan(num_links)] <- 0
    cluster_mat[is.nan(cluster_mat)] <- 0
    cluster_mat <- matrix(stats::p.adjust(cluster_mat), nrow = length(louvain_modules), ncol = length(louvain_modules))
    sig_links <- as.matrix(num_links)
    row.names(sig_links) = colnames(sig_links) = louvain_modules
    sig_links[cluster_mat > qval_thresh] = 0
    diag(sig_links) <- 0
    cluster_g <- igraph::graph_from_adjacency_matrix(sig_links, weighted = TRUE, mode = "undirected")
    list(cluster_g = cluster_g, num_links = num_links, cluster_mat = cluster_mat)
}

`leiden_clustering` <- function (data, pd, weight = NULL, nn_index = NULL, k = 20, nn_control = list(), num_iter = 2, resolution_parameter = 1e-04, 
    random_seed = NULL, verbose = FALSE, ...) 
{
    extra_arguments <- list(...)
    if ("partition_type" %in% names(extra_arguments)) 
        partition_type <- extra_arguments[["partition_type"]]
    else partition_type <- "CPMVertexPartition"
    if ("initial_membership" %in% names(extra_arguments)) 
        initial_membership <- extra_arguments[["initial_membership"]]
    else initial_membership <- NULL
    if ("weights" %in% names(extra_arguments)) 
        edge_weights <- extra_arguments[["weights"]]
    else edge_weights <- NULL
    if ("node_sizes" %in% names(extra_arguments)) 
        node_sizes <- extra_arguments[["node_sizes"]]
    else node_sizes <- NULL
    assertthat::assert_that(assertthat::is.count(k))
    if (partition_type %in% c("ModularityVertexPartition", "SignificanceVertexPartition", "SurpriseVertexPartition")) {
        resolution_parameter = NA
    }
    else if (is.null(resolution_parameter)) {
        resolution_parameter = 1e-04
    }
    if (is.null(num_iter)) 
        num_iter = 2
    if (random_seed == 0L) 
        random_seed = NULL
    cell_names <- row.names(pd)
    if (!identical(cell_names, row.names(pd))) 
        stop("Phenotype and row name from the data don't match")
    graph_result <- tryCatch(cluster_cells_make_graph(data = data, weight = weight, cell_names = cell_names, nn_index, k = k, 
        nn_control = nn_control, verbose = verbose), error = function(c) {
        stop(paste0(trimws(c), "\n* error in leiden_clustering"))
    })
    if (verbose) 
        message("  Run leiden clustering ...")
    t_start <- Sys.time()
    if (verbose) {
        table_results <- data.frame(resolution_parameter = double(), quality = double(), modularity = double(), significance = double(), 
            number_clusters = integer())
    }
    best_modularity <- -1
    best_result <- NULL
    best_resolution_parameter <- "No resolution"
    if (length(resolution_parameter) < 1) 
        warning("bad loop: length(resolution_parameter) < 1")
    for (i in 1:length(resolution_parameter)) {
        cur_resolution_parameter <- resolution_parameter[i]
        cluster_result <- leidenbase::leiden_find_partition(graph_result[["g"]], partition_type = partition_type, initial_membership = initial_membership, 
            edge_weights = edge_weights, node_sizes = node_sizes, seed = random_seed, resolution_parameter = cur_resolution_parameter, 
            num_iter = num_iter, verbose = verbose)
        quality <- cluster_result[["quality"]]
        modularity <- cluster_result[["modularity"]]
        significance <- cluster_result[["significance"]]
        if (verbose) 
            table_results <- rbind(table_results, data.frame(resolution_parameter = cur_resolution_parameter, quality = quality, 
                modularity = modularity, significance = significance, cluster_count = max(cluster_result[["membership"]])))
        if (verbose) 
            message("    Current resolution is ", cur_resolution_parameter, "; Modularity is ", modularity, "; Quality is ", 
                quality, "; Significance is ", significance, "; Number of clusters is ", max(cluster_result[["membership"]]))
        if (modularity > best_modularity) {
            best_result <- cluster_result
            best_resolution_parameter <- cur_resolution_parameter
            best_modularity <- modularity
        }
        if (is.null(best_result)) {
            best_result <- cluster_result
            best_resolution_parameter <- NULL
            best_modularity <- cluster_result[["modularity"]]
        }
    }
    t_end <- Sys.time()
    if (verbose) {
        message("    Done. Run time: ", t_end - t_start, "s\n")
        message("  Clustering statistics")
        selected <- vector(mode = "character", length = length(resolution_parameter))
        if (length(resolution_parameter) < 1) 
            warning("bad loop: length(resolution_parameter) < 1")
        for (irespar in 1:length(resolution_parameter)) {
            if (identical(table_results[["resolution_parameter"]][irespar], best_resolution_parameter)) 
                selected[irespar] <- "*"
            else selected[irespar] <- " "
        }
        print(cbind(` ` = " ", table_results, selected), row.names = FALSE)
        message()
        message("  Cell counts by cluster")
        membership <- best_result[["membership"]]
        membership_frequency <- stats::aggregate(data.frame(cell_count = membership), list(cluster = membership), length)
        membership_frequency <- cbind(` ` = " ", membership_frequency, cell_fraction = sprintf("%.3f", membership_frequency[["cell_count"]]/sum(membership_frequency[["cell_count"]])))
        print(membership_frequency, row.names = FALSE)
        message()
        message("  Maximal modularity is ", best_modularity, " for resolution parameter ", best_resolution_parameter)
        message("\n  Run kNN based graph clustering DONE.\n  -Number of clusters: ", max(best_result[["membership"]]))
    }
    if (igraph::vcount(graph_result[["g"]]) < 3000) {
        coord <- NULL
        edge_links <- NULL
    }
    else {
        coord <- NULL
        edge_links <- NULL
    }
    igraph::V(graph_result[["g"]])$names <- as.character(igraph::V(graph_result[["g"]]))
    out_result <- list(membership = best_result[["membership"]], modularity = best_result[["modularity"]])
    names(out_result$membership) = cell_names
    return(list(g = graph_result[["g"]], relations = graph_result[["relations"]], distMatrix = graph_result[["distMatrix"]], 
        coord = coord, edge_links = edge_links, optim_res = out_result))
}

`report_nn_control` <- function (label = NULL, nn_control) 
{
    indent <- ""
    if (!is.null(label)) {
        indent <- "  "
    }
    message(ifelse(!is.null(label), label, ""))
    message(indent, "  method: ", ifelse(!is.null(nn_control[["method"]]), nn_control[["method"]], as.character(NA)))
    message(indent, "  metric: ", ifelse(!is.null(nn_control[["metric"]]), nn_control[["metric"]], as.character(NA)))
    if (is.null(nn_control[["method"]])) {
        return()
    }
    if (nn_control[["method"]] == "nn2") {
        message(indent, "  nn2 has no parameters")
    }
    else if (nn_control[["method"]] == "annoy") {
        message(indent, "  n_trees: ", ifelse(!is.null(nn_control[["n_trees"]]), nn_control[["n_trees"]], as.character(NA)))
        message(indent, "  search_k: ", ifelse(!is.null(nn_control[["search_k"]]), nn_control[["search_k"]], as.character(NA)))
        message(indent, "  cores: ", ifelse(!is.null(nn_control[["cores"]]), nn_control[["cores"]], as.character(NA)))
        message(indent, "  grain_size: ", ifelse(!is.null(nn_control[["grain_size"]]), nn_control[["grain_size"]], as.character(NA)))
    }
    else if (nn_control[["method"]] == "hnsw") {
        message(indent, "  M: ", ifelse(!is.null(nn_control[["M"]]), nn_control[["M"]], as.character(NA)))
        message(indent, "  ef_construction: ", ifelse(!is.null(nn_control[["ef_construction"]]), nn_control[["ef_construction"]], 
            as.character(NA)))
        message(indent, "  ef: ", ifelse(!is.null(nn_control[["ef"]]), nn_control[["ef"]], as.character(NA)))
        message(indent, "  cores: ", ifelse(!is.null(nn_control[["cores"]]), nn_control[["cores"]], as.character(NA)))
        message(indent, "  grain_size: ", ifelse(!is.null(nn_control[["grain_size"]]), nn_control[["grain_size"]], as.character(NA)))
    }
    else stop("report_nn_control: unsupported nearest neighbor method '", nn_control[["method"]], "'")
}

`multi_component_RGE` <- function (cds, scale = FALSE, reduction_method, partition_list, max_components, ncenter, irlba_pca_res, nn.k = 25, nn_control = list(), 
    maxiter, eps, L1.gamma, L1.sigma, close_loop = FALSE, euclidean_distance_ratio = 1, geodesic_distance_ratio = 1/3, prune_graph = TRUE, 
    minimal_branch_len = minimal_branch_len, verbose = FALSE) 
{
    cluster <- NULL
    X <- t(irlba_pca_res)
    dp_mst <- NULL
    pr_graph_cell_proj_closest_vertex <- NULL
    cell_name_vec <- NULL
    reducedDimK_coord <- NULL
    merge_rge_res <- NULL
    max_ncenter <- 0
    for (cur_comp in sort(unique(partition_list))) {
        if (verbose) {
            message("Processing partition component ", cur_comp)
        }
        X_subset <- X[, partition_list == cur_comp]
        if (verbose) 
            message("Current partition is ", cur_comp)
        if (scale) {
            X_subset <- t(as.matrix(scale(t(X_subset))))
        }
        if (is.null(ncenter)) {
            num_clusters_in_partition <- length(unique(clusters(cds, reduction_method)[colnames(X_subset)]))
            num_cells_in_partition = ncol(X_subset)
            curr_ncenter <- cal_ncenter(num_clusters_in_partition, num_cells_in_partition)
            if (is.null(curr_ncenter) || curr_ncenter >= ncol(X_subset)) {
                curr_ncenter <- ncol(X_subset) - 1
            }
        }
        else {
            curr_ncenter <- min(ncol(X_subset) - 1, ncenter)
        }
        if (verbose) 
            message("Using ", curr_ncenter, " nodes for principal graph")
        kmean_res <- NULL
        centers <- t(X_subset)[seq(1, ncol(X_subset), length.out = curr_ncenter), , drop = FALSE]
        centers <- centers + matrix(stats::rnorm(length(centers), sd = 1e-10), nrow = nrow(centers))
        kmean_res <- tryCatch({
            stats::kmeans(t(X_subset), centers = centers, iter.max = 100)
        }, error = function(err) {
            stats::kmeans(t(X_subset), centers = curr_ncenter, iter.max = 100)
        })
        if (kmean_res$ifault != 0) {
            message("kmeans returned ifault = ", kmean_res$ifault)
        }
        nearest_center <- find_nearest_vertex(t(kmean_res$centers), X_subset, process_targets_in_blocks = TRUE)
        medioids <- X_subset[, unique(nearest_center)]
        reduced_dim_res <- t(medioids)
        mat <- t(X_subset)
        if (is.null(nn.k)) {
            k <- round(sqrt(nrow(mat))/2)
            k <- max(10, k)
        }
        else {
            k <- nn.k
        }
        if (verbose) 
            message("Finding kNN with ", k, " neighbors")
        nn_method <- nn_control[["method"]]
        dx <- search_nn_matrix(subject_matrix = mat, query_matrix = mat, k = min(k, nrow(mat) - 1), nn_control = nn_control, 
            verbose = verbose)
        if (nn_method == "annoy" || nn_method == "hnsw") 
            dx <- swap_nn_row_index_point(nn_res = dx, verbose = verbose)
        nn.index <- dx$nn.idx[, -1]
        nn.dist <- dx$nn.dists[, -1]
        if (verbose) 
            message("Calculating the local density for each sample based on kNNs ...")
        rho <- exp(-rowMeans(nn.dist))
        mat_df <- as.data.frame(mat)
        tmp <- mat_df %>% tibble::rownames_to_column() %>% dplyr::mutate(cluster = kmean_res$cluster, density = rho) %>% 
            dplyr::group_by(cluster) %>% dplyr::top_n(n = 1, wt = density) %>% dplyr::arrange(-dplyr::desc(cluster))
        medioids <- X_subset[, tmp$rowname]
        reduced_dim_res <- t(medioids)
        graph_args <- list(X = X_subset, C0 = medioids, maxiter = maxiter, eps = eps, L1.gamma = L1.gamma, L1.sigma = L1.sigma, 
            verbose = verbose)
        rge_res <- do.call(calc_principal_graph, graph_args)
        names(rge_res)[c(2, 4, 5)] <- c("Y", "R", "objective_vals")
        stree <- rge_res$W
        stree_ori <- stree
        if (close_loop) {
            reduce_dims_old <- t(SingleCellExperiment::reducedDims(cds)[[reduction_method]])[, partition_list == cur_comp]
            connect_tips_res <- tryCatch(connect_tips(cds, pd = colData(cds)[partition_list == cur_comp, ], R = rge_res$R, 
                stree = stree, reducedDimK_old = rge_res$Y, reducedDimS_old = reduce_dims_old, k = 25, nn_control = nn_control, 
                kmean_res = kmean_res, euclidean_distance_ratio = euclidean_distance_ratio, geodesic_distance_ratio = geodesic_distance_ratio, 
                medioids = medioids, verbose = verbose), error = function(c) {
                stop(paste0(trimws(c), "\n* error in multi_component_RGE"))
            })
            stree <- connect_tips_res$stree
        }
        if (prune_graph) {
            if (verbose) {
                message("Running graph pruning ...")
            }
            stree <- prune_tree(stree_ori, as.matrix(stree), minimal_branch_len = minimal_branch_len)
            rge_res$Y <- rge_res$Y[, match(row.names(stree), row.names(stree_ori))]
            rge_res$R <- rge_res$R[, match(row.names(stree), row.names(stree_ori))]
            medioids <- medioids[, row.names(stree)]
        }
        if (is.null(merge_rge_res)) {
            if (ncol(rge_res$Y) < 1) 
                warning("bad loop: ncol(rge_res$Y) < 1")
            colnames(rge_res$Y) <- paste0("Y_", 1:ncol(rge_res$Y))
            merge_rge_res <- rge_res
            colnames(merge_rge_res$X) <- colnames(X_subset)
            row.names(merge_rge_res$R) <- colnames(X_subset)
            if (ncol(merge_rge_res$Y) < 1) 
                warning("bad loop: ncol(merge_rge_res$Y) < 1")
            colnames(merge_rge_res$R) <- paste0("Y_", 1:ncol(merge_rge_res$Y))
            merge_rge_res$R <- list(merge_rge_res$R)
            merge_rge_res$stree <- list(stree)
            merge_rge_res$objective_vals <- list(merge_rge_res$objective_vals)
        }
        else {
            colnames(rge_res$X) <- colnames(X_subset)
            row.names(rge_res$R) <- colnames(X_subset)
            colnames(rge_res$R) <- paste0("Y_", (ncol(merge_rge_res$Y) + 1):(ncol(merge_rge_res$Y) + ncol(rge_res$Y)), sep = "")
            colnames(rge_res$Y) <- paste("Y_", (ncol(merge_rge_res$Y) + 1):(ncol(merge_rge_res$Y) + ncol(rge_res$Y)), sep = "")
            merge_rge_res$Y <- cbind(merge_rge_res$Y, rge_res$Y)
            merge_rge_res$R <- c(merge_rge_res$R, list(rge_res$R))
            merge_rge_res$stree <- c(merge_rge_res$stree, list(stree))
            merge_rge_res$objective_vals <- c(merge_rge_res$objective_vals, list(rge_res$objective_vals))
        }
        if (is.null(reducedDimK_coord)) {
            if (ncol(rge_res$Y) < 1) 
                warning("bad loop: ncol(rge_res$Y) < 1")
            curr_cell_names <- paste("Y_", 1:ncol(rge_res$Y), sep = "")
            pr_graph_cell_proj_closest_vertex <- matrix(apply(rge_res$R, 1, which.max))
            cell_name_vec <- colnames(X_subset)
        }
        else {
            curr_cell_names <- paste("Y_", (ncol(reducedDimK_coord) + 1):(ncol(reducedDimK_coord) + ncol(rge_res$Y)), sep = "")
            pr_graph_cell_proj_closest_vertex <- rbind(pr_graph_cell_proj_closest_vertex, matrix(apply(rge_res$R, 1, which.max) + 
                ncol(reducedDimK_coord)))
            cell_name_vec <- c(cell_name_vec, colnames(X_subset))
        }
        curr_reducedDimK_coord <- rge_res$Y
        dimnames(stree) <- list(curr_cell_names, curr_cell_names)
        cur_dp_mst <- igraph::graph.adjacency(stree, mode = "undirected", weighted = TRUE)
        if (any(igraph::E(cur_dp_mst)$weight != 1)) {
            message("Warning: multi_component_RGE: not all weights are 1")
        }
        dp_mst <- igraph::graph.union(dp_mst, cur_dp_mst)
        reducedDimK_coord <- cbind(reducedDimK_coord, curr_reducedDimK_coord)
    }
    igraph::E(dp_mst)$weight <- rep(1, igraph::ecount(dp_mst))
    row.names(pr_graph_cell_proj_closest_vertex) <- cell_name_vec
    ddrtree_res_W <- as.matrix(rge_res$W)
    ddrtree_res_Z <- SingleCellExperiment::reducedDims(cds)[[reduction_method]]
    ddrtree_res_Y <- reducedDimK_coord
    R <- Matrix::sparseMatrix(i = 1, j = 1, x = 0, dims = c(ncol(cds), ncol(merge_rge_res$Y)))
    stree <- Matrix::sparseMatrix(i = 1, j = 1, x = 0, dims = c(ncol(merge_rge_res$Y), ncol(merge_rge_res$Y)))
    curr_row_id <- 1
    curr_col_id <- 1
    R_row_names <- NULL
    if (length(merge_rge_res$R) < 1) 
        warning("bad loop: length(merge_rge_res$R) < 1")
    for (i in 1:length(merge_rge_res$R)) {
        current_R <- merge_rge_res$R[[i]]
        stree[curr_col_id:(curr_col_id + ncol(current_R) - 1), curr_col_id:(curr_col_id + ncol(current_R) - 1)] <- merge_rge_res$stree[[i]]
        curr_row_id <- curr_row_id + nrow(current_R)
        curr_col_id <- curr_col_id + ncol(current_R)
        R_row_names <- c(R_row_names, row.names(current_R))
    }
    row.names(R) <- R_row_names
    R <- R[colnames(cds), ]
    cds@principal_graph_aux[[reduction_method]] <- list(stree = stree, Q = merge_rge_res$Q, R = R, objective_vals = merge_rge_res$objective_vals, 
        history = merge_rge_res$history)
    cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_closest_vertex <- as.data.frame(pr_graph_cell_proj_closest_vertex)[colnames(cds), 
        , drop = FALSE]
    if (ncol(ddrtree_res_Y) < 1) 
        warning("bad loop: ncol(ddrtree_res_Y) < 1")
    colnames(ddrtree_res_Y) <- paste0("Y_", 1:ncol(ddrtree_res_Y), sep = "")
    return(list(cds = cds, ddrtree_res_W = ddrtree_res_W, ddrtree_res_Z = ddrtree_res_Z, ddrtree_res_Y = ddrtree_res_Y, dp_mst = dp_mst))
}

`project2MST` <- function (cds, Projection_Method, orthogonal_proj_tip = FALSE, verbose, reduction_method, rge_res_Y) 
{
    target <- group <- distance_2_source <- rowname <- NULL
    dp_mst <- principal_graph(cds)[[reduction_method]]
    Z <- t(SingleCellExperiment::reducedDims(cds)[[reduction_method]])
    Y <- rge_res_Y
    cds <- findNearestPointOnMST(cds, reduction_method, rge_res_Y)
    closest_vertex <- cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_closest_vertex
    closest_vertex_names <- colnames(Y)[closest_vertex[, 1]]
    tip_leaves <- names(which(igraph::degree(dp_mst) == 1))
    if (!is.function(Projection_Method)) {
        P <- Y[, closest_vertex]
    }
    else {
        if (length(Z) < 1) 
            warning("bad loop: length(Z) < 1")
        P <- matrix(rep(0, length(Z)), nrow = nrow(Z))
        if (length(Z[1:2, ]) < 1) 
            warning("bad loop: length(Z[1:2, ]) < 1")
        nearest_edges <- matrix(rep(0, length(Z[1:2, ])), ncol = 2)
        row.names(nearest_edges) <- colnames(cds)
        if (length(closest_vertex) < 1) 
            warning("bad loop: length(closest_vertex) < 1")
        for (i in 1:length(closest_vertex)) {
            neighbors <- names(igraph::neighborhood(dp_mst, nodes = closest_vertex_names[i], mode = "all")[[1]])[-1]
            projection <- NULL
            distance <- NULL
            Z_i <- Z[, i]
            for (neighbor in neighbors) {
                if (closest_vertex_names[i] %in% tip_leaves) {
                  if (orthogonal_proj_tip) {
                    tmp <- projPointOnLine(Z_i, Y[, c(closest_vertex_names[i], neighbor)])
                  }
                  else {
                    tmp <- Projection_Method(Z_i, Y[, c(closest_vertex_names[i], neighbor)])
                  }
                }
                else {
                  tmp <- Projection_Method(Z_i, Y[, c(closest_vertex_names[i], neighbor)])
                }
                if (any(is.na(tmp))) {
                  tmp <- Y[, neighbor]
                }
                projection <- rbind(projection, tmp)
                distance <- c(distance, stats::dist(rbind(Z_i, tmp)))
            }
            if (class(projection)[1] != "matrix") {
                projection <- as.matrix(projection)
            }
            which_min <- which.min(distance)
            P[, i] <- projection[which_min, ]
            nearest_edges[i, ] <- c(closest_vertex_names[i], neighbors[which_min])
        }
    }
    colnames(P) <- colnames(Z)
    dp_mst_list <- igraph::decompose.graph(dp_mst)
    dp_mst_df <- NULL
    partitions <- cds@clusters[[reduction_method]]$partitions
    assertthat::assert_that(!is.null(names(cds@clusters[[reduction_method]]$partitions)), msg = "names(cds@clusters[[reduction_method]]$partitions) == NULL")
    assertthat::assert_that(!(length(colnames(cds)) != length(names(cds@clusters[[reduction_method]]$partitions))), msg = "length( colnames(cds) ) != length(names(cds@clusters[[reduction_method]]$partitions))")
    assertthat::assert_that(!any(colnames(cds) != names(cds@clusters[[reduction_method]]$partitions)), msg = "colnames(cds)!=names(cds@clusters[[reduction_method]]$partitions)")
    if (length(dp_mst_list) == 1 & length(unique(partitions)) > 1) {
        partitions[partitions != "1"] <- "1"
    }
    if (!is.null(partitions)) {
        for (cur_partition in sort(unique(partitions))) {
            data_df <- NULL
            if (verbose) {
                message("\nProjecting cells to principal points for partition: ", cur_partition)
            }
            subset_cds_col_names <- names(partitions[partitions == cur_partition])
            cur_z <- Z[, subset_cds_col_names]
            cur_p <- P[, subset_cds_col_names]
            if (ncol(cur_p) > 0 && nrow(cur_p) > 0) {
                cur_centroid_name <- igraph::V(dp_mst_list[[as.numeric(cur_partition)]])$name
                cur_nearest_edges <- nearest_edges[subset_cds_col_names, ]
                data_df <- cbind(as.data.frame(t(cur_p), stringsAsFactors = FALSE), apply(cur_nearest_edges, 1, sort) %>% 
                  t(), stringsAsFactors = FALSE)
                row.names(data_df) <- colnames(cur_p)
                if (nrow(cur_p) < 1) 
                  warning("bad loop: nrow(cur_p) < 1")
                colnames(data_df) <- c(paste0("P_", 1:nrow(cur_p)), "source", "target")
                data_df$distance_2_source <- sqrt(colSums((cur_p - rge_res_Y[, as.character(data_df[, "source"])])^2))
                data_df <- data_df %>% tibble::rownames_to_column() %>% dplyr::mutate(group = paste(source, target, sep = "_")) %>% 
                  dplyr::arrange(group, dplyr::desc(-distance_2_source))
                data_df <- data_df %>% dplyr::group_by(group) %>% dplyr::mutate(new_source = dplyr::lag(rowname), new_target = rowname)
                data_df[is.na(data_df$new_source), "new_source"] <- as.character(as.matrix(data_df[is.na(data_df$new_source), 
                  "source"]))
                added_rows <- which(is.na(data_df$new_source) & is.na(data_df$new_target))
                data_df <- as.data.frame(data_df, stringsAsFactors = FALSE)
                data_df <- as.data.frame(as.matrix(data_df), stringsAsFactors = FALSE)
                data_df[added_rows, c("new_source", "new_target")] <- data_df[added_rows - 1, c("rowname", "target")]
                aug_P = cbind(cur_p, rge_res_Y, stringsAsFactors = F)
                data_df$weight <- sqrt(colSums((aug_P[, data_df$new_source] - aug_P[, data_df$new_target]))^2)
                data_df$weight <- data_df$weight + min(data_df$weight[data_df$weight > 0])
                edge_list <- as.data.frame(igraph::get.edgelist(dp_mst_list[[as.numeric(cur_partition)]]), stringsAsFactors = FALSE)
                dp <- as.matrix(stats::dist(t(rge_res_Y)[cur_centroid_name, ]))
                edge_list$weight <- dp[cbind(edge_list[, 1], edge_list[, 2])]
                colnames(edge_list) <- c("new_source", "new_target", "weight")
                dp_mst_df <- Reduce(rbind, list(dp_mst_df, data_df[, c("new_source", "new_target", "weight")], edge_list))
            }
        }
    }
    dp_mst <- igraph::graph.data.frame(dp_mst_df, directed = FALSE)
    cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_tree <- dp_mst
    cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_dist <- P
    closest_vertex_df <- as.matrix(closest_vertex)
    row.names(closest_vertex_df) <- row.names(closest_vertex)
    cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df
    cds
}

`project_point_to_line_segment` <- function (p, df) 
{
    A <- df[, 1]
    B <- df[, 2]
    AB <- (B - A)
    AB_squared = sum(AB^2)
    if (AB_squared == 0) {
        q <- A
    }
    else {
        Ap <- (p - A)
        t <- sum(Ap * AB)/AB_squared
        if (t < 0) {
            q <- A
        }
        else if (t > 1) {
            q <- B
        }
        else {
            q <- A + t * AB
        }
    }
    return(q)
}

`extract_general_graph_ordering` <- function (cds, root_pr_nodes, verbose = TRUE, reduction_method) 
{
    Z <- t(SingleCellExperiment::reducedDims(cds)[[reduction_method]])
    Y <- cds@principal_graph_aux[[reduction_method]]$dp_mst
    pr_graph <- principal_graph(cds)[[reduction_method]]
    parents <- rep(NA, length(igraph::V(pr_graph)))
    states <- rep(NA, length(igraph::V(pr_graph)))
    if (any(is.na(igraph::E(pr_graph)$weight))) {
        igraph::E(pr_graph)$weight <- 1
    }
    closest_vertex <- find_nearest_vertex(Y[, root_pr_nodes, drop = FALSE], Z)
    closest_vertex_id <- colnames(cds)[closest_vertex]
    cell_wise_graph <- cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_tree
    cell_wise_distances <- igraph::distances(cell_wise_graph, v = closest_vertex_id)
    if (length(closest_vertex_id) > 1) {
        node_names <- colnames(cell_wise_distances)
        pseudotimes <- apply(cell_wise_distances, 2, min)
    }
    else {
        node_names <- names(cell_wise_distances)
        pseudotimes <- cell_wise_distances
    }
    names(pseudotimes) <- node_names
    ordering_df <- data.frame(sample_name = igraph::V(cell_wise_graph)$name, pseudo_time = as.vector(pseudotimes))
    row.names(ordering_df) <- ordering_df$sample_name
    return(ordering_df)
}

`get_matrix_class` <- function (mat) 
{
    nmatch <- 0
    matrix_info <- list()
    if (is(mat, "matrix")) {
        matrix_info[["matrix_class"]] <- "r_dense_matrix"
        nmatch <- nmatch + 1
    }
    if (any(class(mat) %in% c("dgCMatrix"))) {
        matrix_info[["matrix_class"]] <- "dgCMatrix"
        nmatch <- nmatch + 1
    }
    if (any(class(mat) %in% c("lgCMatrix"))) {
        matrix_info[["matrix_class"]] <- "lgCMatrix"
        nmatch <- nmatch + 1
    }
    if (is(mat, "dgTMatrix")) {
        matrix_info[["matrix_class"]] <- "dgTMatrix"
        nmatch <- nmatch + 1
    }
    if (is(mat, "dgeMatrix")) {
        matrix_info[["matrix_class"]] <- "dgeMatrix"
        nmatch <- nmatch + 1
    }
    if (is(mat, "IterableMatrix")) {
        matrix_info[["matrix_class"]] <- "BPCells"
        nmatch <- nmatch + 1
    }
    if (nmatch == 0) {
        stop("get_matrix_class: unrecognized matrix class")
    }
    if (nmatch > 1) {
        stop("get_matrix_class: ambiguous matrix class")
    }
    return(matrix_info)
}

`select_nn_parameter_value` <- function (parameter, nn_control, nn_control_default, default_value) 
{
    if (!is.null(nn_control[[parameter]])) {
        return(nn_control[[parameter]])
    }
    else if (!is.null(nn_control_default[[parameter]])) {
        return(nn_control_default[[parameter]])
    }
    return(default_value)
}

`select_annoy_search_k` <- function (mode, nn_control, nn_control_default, nn_index, k, default_n_trees, default_k) 
{
    if (!is.null(k)) {
        use_k <- k
        src_k <- "parameter"
    }
    else {
        use_k <- default_k
        src_k <- "default"
    }
    if (mode == 2 && !is.null(nn_index) && is.null(nn_index[["n_trees"]])) {
        stop("set_nn_control: unexpected condition: found old version of reduce_dim_aux")
    }
    if (!is.null(nn_control[["search_k"]])) {
        return(nn_control[["search_k"]])
    }
    else if (mode == 2 && !is.null(nn_index) && !is.null(nn_index[["n_trees"]]) && !is.null(k)) {
        n_trees <- nn_index[["n_trees"]]
        return(2 * n_trees * k)
    }
    else if (!is.null(nn_control[["n_trees"]])) {
        return(2 * nn_control[["n_trees"]] * use_k)
    }
    else if (!is.null(nn_control_default[["search_k"]])) {
        return(nn_control_default[["search_k"]])
    }
    else if (!is.null(nn_control_default[["n_trees"]])) {
        return(2 * nn_control_default[["n_trees"]] * use_k)
    }
    return(2 * default_n_trees * use_k)
}

`get_call_stack_as_string` <- function () 
{
    cs <- get_call_stack()
    scs <- ""
    for (i in seq(length(cs) - 1)) {
        csep <- ifelse(i == 1, "", " => ")
        scs <- sprintf("%s%s%s()", scs, csep, cs[[i]])
    }
    return(scs)
}

`stop_no_noise` <- function () 
{
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
}

`set_global_variable` <- function (variable_name, value) 
{
    assign(variable_name, value, envir = ._._global_variable_env_._.)
}

`initialize_counts_metadata` <- function (cds) 
{
    assertthat::assert_that(methods::is(cds, "mmsage_cell_data_set"), msg = paste("cds parameter is not a cell_data_set"))
    if (is.null(SingleCellExperiment::int_metadata(cds))) {
        SingleCellExperiment::int_metadata(cds) <- list()
    }
    SingleCellExperiment::int_metadata(cds)[["counts_metadata"]] <- list()
    return(cds)
}

`tick` <- function (msg = "") 
{
    set_global_variable("monocle3_timer_t0", Sys.time())
    set_global_variable("monocle3_timer_msg", msg)
}

`new_annoy_index` <- function (metric, ndim) 
{
    nn_class <- switch(metric, cosine = RcppAnnoy::AnnoyAngular, euclidean = RcppAnnoy::AnnoyEuclidean, hamming = RcppAnnoy::AnnoyHamming, 
        manhattan = RcppAnnoy::AnnoyManhattan, stop("unsupported annoy metric ", metric))
    nn_index <- methods::new(nn_class, ndim)
    return(nn_index)
}

`tock` <- function () 
{
    t1 <- Sys.time()
    t0 <- get_global_variable("monocle3_timer_t0")
    msg <- get_global_variable("monocle3_timer_msg")
    if (length(msg) > 0) {
        message(sprintf("%s %.2f seconds.", msg, difftime(t1, t0, units = "secs")))
    }
    else {
        return(t1 - t0)
    }
}

`cluster_cells_make_graph` <- function (data, weight, cell_names, nn_index = NULL, k = k, nn_control = list(), verbose) 
{
    if (is.data.frame(data)) 
        data <- as.matrix(data)
    if (!is.matrix(data)) 
        stop("Wrong input data, should be a data frame or matrix!")
    if (k < 1) {
        stop("k must be a positive integer!")
    }
    else if (k > nrow(data) - 2) {
        k <- nrow(data) - 2
        warning("The nearest neighbors includes the point itself, k must be smaller than\nthe ", "total number of points - 1 (all other points) - 1 ", 
            "(itself)! ", "Total number of points is ", nrow(data))
    }
    if (verbose) {
        message("Run kNN based graph clustering starts:", "\n", "  -Input data of ", nrow(data), " rows and ", ncol(data), 
            " columns", "\n", "  -k is set to ", k)
        message("  Finding nearest neighbors...")
    }
    nn_method <- nn_control[["method"]]
    if (nn_method == "nn2") {
        t1 <- system.time(tmp <- RANN::nn2(data, data, k + 1, searchtype = "standard"))
    }
    else {
        if (is.null(nn_index)) {
            nn_index <- make_nn_index(subject_matrix = data, nn_control = nn_control, verbose = verbose)
        }
        tmp <- tryCatch(search_nn_index(query_matrix = data, nn_index = nn_index, k = k + 1, nn_control = nn_control, verbose = verbose), 
            error = function(c) {
                stop(paste0(trimws(c), "\n* error in cluster_cells_make_graph"))
            })
        if (nn_method == "annoy" || nn_method == "hnsw") {
            tmp <- swap_nn_row_index_point(nn_res = tmp, verbose = verbose)
        }
    }
    neighborMatrix <- tmp[["nn.idx"]][, -1]
    distMatrix <- tmp[["nn.dists"]][, -1]
    if (verbose) {
        if (nn_method == "nn2") {
            message("DONE. Run time: ", t1[3], "s\n")
        }
        message("Compute jaccard coefficient between nearest-neighbor sets ...")
    }
    t2 <- system.time(links <- jaccard_coeff(neighborMatrix, weight))
    if (verbose) 
        message("DONE. Run time:", t2[3], "s\n", " Build undirected graph from the weighted links ...")
    links <- links[links[, 1] > 0, ]
    relations <- as.data.frame(links)
    colnames(relations) <- c("from", "to", "weight")
    relations$from <- cell_names[relations$from]
    relations$to <- cell_names[relations$to]
    t3 <- system.time(g <- igraph::graph.data.frame(relations, directed = FALSE))
    if (verbose) 
        message("DONE ~", t3[3], "s\n")
    return(list(g = g, distMatrix = distMatrix, relations = relations))
}

`cal_ncenter` <- function (num_cell_communities, ncells, nodes_per_log10_cells = 15) 
{
    round(num_cell_communities * nodes_per_log10_cells * log10(ncells))
}

`find_nearest_vertex` <- function (data_matrix, target_points, block_size = 50000, process_targets_in_blocks = FALSE) 
{
    closest_vertex = c()
    if (process_targets_in_blocks == FALSE) {
        num_blocks = ceiling(ncol(data_matrix)/block_size)
        if (num_blocks < 1) 
            warning("bad loop: num_blocks < 1")
        for (i in 1:num_blocks) {
            if (i < num_blocks) {
                block = data_matrix[, ((((i - 1) * block_size) + 1):(i * block_size))]
            }
            else {
                block = data_matrix[, ((((i - 1) * block_size) + 1):(ncol(data_matrix)))]
            }
            distances_Z_to_Y <- proxy::dist(t(block), t(target_points))
            closest_vertex_for_block <- apply(distances_Z_to_Y, 1, function(z) {
                which.min(z)
            })
            closest_vertex = append(closest_vertex, closest_vertex_for_block)
        }
    }
    else {
        num_blocks = ceiling(ncol(target_points)/block_size)
        dist_to_closest_vertex = rep(Inf, length(ncol(data_matrix)))
        closest_vertex = rep(NA, length(ncol(data_matrix)))
        if (num_blocks < 1) 
            warning("bad loop: num_blocks < 1")
        for (i in 1:num_blocks) {
            if (i < num_blocks) {
                block = target_points[, ((((i - 1) * block_size) + 1):(i * block_size))]
            }
            else {
                block = target_points[, ((((i - 1) * block_size) + 1):(ncol(target_points)))]
            }
            distances_Z_to_Y <- proxy::dist(t(data_matrix), t(block))
            closest_vertex_for_block <- apply(distances_Z_to_Y, 1, function(z) {
                which.min(z)
            })
            if (nrow(distances_Z_to_Y) < 1) 
                warning("bad loop: nrow(distances_Z_to_Y) < 1")
            new_block_distances <- distances_Z_to_Y[cbind(1:nrow(distances_Z_to_Y), closest_vertex_for_block)]
            updated_nearest_idx <- which(new_block_distances < dist_to_closest_vertex)
            closest_vertex[updated_nearest_idx] <- closest_vertex_for_block[updated_nearest_idx] + (i - 1) * block_size
            dist_to_closest_vertex[updated_nearest_idx] <- new_block_distances[updated_nearest_idx]
        }
    }
    stopifnot(length(closest_vertex) == ncol(data_matrix))
    return(closest_vertex)
}

`search_nn_matrix` <- function (subject_matrix, query_matrix, k = 25, nn_control = list(), verbose = FALSE) 
{
    assertthat::assert_that(methods::is(subject_matrix, "matrix") || is_sparse_matrix(subject_matrix), msg = paste0("search_nn_matrix: the subject_matrix object must be of type matrix"))
    assertthat::assert_that(methods::is(query_matrix, "matrix") || is_sparse_matrix(query_matrix), msg = paste0("search_nn_matrix: the query_matrix object must be of type matrix"))
    nn_control_default <- get_global_variable("nn_control_annoy_euclidean")
    nn_control <- set_nn_control(mode = 3, nn_control = nn_control, nn_control_default = nn_control_default, nn_index = NULL, 
        k = k, verbose = verbose)
    method <- nn_control[["method"]]
    k <- min(k, nrow(subject_matrix))
    if (verbose) {
        message("search_nn_matrix:")
        message("  k: ", k)
        report_nn_control("  nn_control: ", nn_control)
        tick("search_nn_matrix: search_time")
    }
    if (method == "nn2") {
        nn_res <- RANN::nn2(subject_matrix, query_matrix, k, searchtype = "standard")
    }
    else {
        nn_index <- tryCatch(make_nn_index(subject_matrix, nn_control = nn_control, verbose = verbose), error = function(c) {
            stop(paste0(trimws(c), "\n* error in search_nn_index"))
        })
        nn_res <- tryCatch(search_nn_index(query_matrix = query_matrix, nn_index = nn_index, k = k, nn_control = nn_control, 
            verbose = verbose), error = function(c) {
            stop(paste0(trimws(c), "\n* error in search_nn_matrix"))
        })
    }
    if (verbose) 
        tock()
    return(nn_res)
}

`swap_nn_row_index_point` <- function (nn_res, verbose = FALSE) 
{
    if (verbose) {
        count_nn_missing_self_index(nn_res, verbose)
    }
    if (check_nn_col1(nn_res$nn.idx)) 
        return(nn_res)
    idx <- nn_res[["nn.idx"]]
    dst <- nn_res[["nn.dists"]]
    diagnostics <- FALSE
    num_no_recall <- 0
    if (nrow(idx) == 0) 
        return(nn_res)
    for (irow in 1:nrow(idx)) {
        vidx <- idx[irow, ]
        vdst <- dst[irow, ]
        if (vidx[[1]] != irow) {
            if (diagnostics) {
                message("swap_nn_row_index_point: adjust nn matrix row: ", irow)
                message("swap_nn_row_index_point: idx row pre fix: ", paste(vidx, collapse = " "))
            }
            if (length(vidx) == 1) {
                num_no_recall <- num_no_recall + 1
                next
            }
            if (vidx[[2]] == irow) {
                vidx[[2]] <- vidx[[1]]
                vidx[[1]] <- irow
            }
            else {
                match <- FALSE
                for (i in seq(2, length(vidx), 1)) {
                  if (vidx[[i]] == irow) {
                    vidx[[i]] <- vidx[[1]]
                    vidx[[1]] <- irow
                    match <- TRUE
                    break
                  }
                }
                if (!match) {
                  if (diagnostics) {
                    message("swap_nn_row_index_point: dst row pre fix: ", paste(vdst, collapse = " "))
                  }
                  for (i in seq(1, length(vidx), 1)) {
                    if (vdst[[i]] > .Machine$double.xmin) {
                      num_no_recall <- num_no_recall + 1
                    }
                  }
                  vidx <- c(irow, vidx[1:(length(vidx) - 1)])
                  vdst <- c(0, vdst[1:(length(vdst) - 1)])
                  dst[irow, ] <- vdst
                  if (diagnostics) {
                    message("swap_nn_row_index_point: dst row post fix: ", paste(vdst, collapse = " "))
                  }
                }
            }
            idx[irow, ] <- vidx
            if (diagnostics) {
                message("swap_nn_row_index_point: idx row post fix: ", paste(vidx, collapse = " "))
            }
        }
    }
    if (num_no_recall > 0) {
        frac_recall <- (nrow(idx) - num_no_recall)/nrow(idx)
        format_recall <- sprintf("%3.1f", frac_recall * 100)
        message("The search result is expected to include the query row value (self)\n", "because the NN index includes the query objects; however, this search result\n", 
            "is missing ", num_no_recall, " self values (recall: ", format_recall, "%). Monocle3 has added the self\n", "values to the first column of the search result in order to allow further\n", 
            "analysis -- but it is missing important nearest neighbors so you need to\n", "increase the sensitivity for making and/or searching the index.")
    }
    nn_res_out <- list(nn.idx = idx, nn.dists = dst)
    return(nn_res_out)
}

`calc_principal_graph` <- function (X, C0, maxiter = 10, eps = 1e-05, L1.gamma = 0.5, L1.sigma = 0.01, verbose = TRUE) 
{
    C <- C0
    K <- ncol(C)
    objs <- c()
    if (maxiter < 1) 
        warning("bad loop: maxiter < 1")
    for (iter in 1:maxiter) {
        norm_sq <- repmat(t(colSums(C^2)), K, 1)
        Phi <- norm_sq + t(norm_sq) - 2 * t(C) %*% C
        g <- igraph::graph.adjacency(Phi, mode = "lower", diag = TRUE, weighted = TRUE)
        g_mst <- igraph::mst(g)
        stree <- igraph::get.adjacency(g_mst, attr = "weight", type = "lower")
        stree_ori <- stree
        stree <- as.matrix(stree)
        stree <- stree + t(stree)
        W <- stree != 0
        obj_W <- sum(sum(stree))
        res = soft_assignment(X, C, L1.sigma)
        P <- res$P
        obj_P <- res$obj
        obj <- obj_W + L1.gamma * obj_P
        objs = c(objs, obj)
        if (verbose) 
            message("iter = ", iter, " obj = ", obj)
        if (iter > 1) {
            relative_diff = abs(objs[iter - 1] - obj)/abs(objs[iter - 1])
            if (relative_diff < eps) {
                if (verbose) 
                  message("eps = ", relative_diff, ", converge.")
                break
            }
            if (iter >= maxiter) {
                if (verbose) 
                  message("eps = ", relative_diff, " reach maxiter.")
            }
        }
        C <- generate_centers(X, W, P, L1.gamma)
    }
    return(list(X = X, C = C, W = W, P = P, objs = objs))
}

`connect_tips` <- function (cds, pd, R, stree, reducedDimK_old, reducedDimS_old, k = 25, nn_control = nn_control, weight = FALSE, qval_thresh = 0.05, 
    kmean_res, euclidean_distance_ratio = 1, geodesic_distance_ratio = 1/3, medioids, verbose = FALSE) 
{
    random_seed <- 0L
    reduction_method <- "UMAP"
    if (is.null(row.names(stree)) & is.null(row.names(stree))) {
        if (ncol(stree) < 1) 
            warning("bad loop: ncol(stree) < 1")
        dimnames(stree) <- list(paste0("Y_", 1:ncol(stree)), paste0("Y_", 1:ncol(stree)))
    }
    stree <- as.matrix(stree)
    stree[stree != 0] <- 1
    mst_g_old <- igraph::graph_from_adjacency_matrix(stree, mode = "undirected")
    if (is.null(kmean_res)) {
        tmp <- matrix(apply(R, 1, which.max))
        row.names(tmp) <- colnames(reducedDimS_old)
        tip_pc_points <- which(igraph::degree(mst_g_old) == 1)
        data <- t(reducedDimS_old[, ])
        cluster_result <- tryCatch(louvain_clustering(data = data, pd = pd[, ], weight = weight, nn_index = NULL, k = k, 
            nn_control = nn_control, louvain_iter = 1, random_seed = 0L, verbose = verbose), error = function(c) {
            stop(paste0(trimws(c), "\n * error in connect_tips"))
        })
        cluster_result$optim_res$membership <- tmp[, 1]
    }
    else {
        tip_pc_points <- which(igraph::degree(mst_g_old) == 1)
        tip_pc_points_kmean_clusters <- sort(kmean_res$cluster[names(tip_pc_points)])
        data <- t(reducedDimS_old[, ])
        cluster_result <- tryCatch(louvain_clustering(data = data, pd = pd[row.names(data), ], weight = weight, nn_index = NULL, 
            k = k, nn_control = nn_control, louvain_iter = 1, random_seed = random_seed, verbose = verbose), error = function(c) {
            stop(paste0(trimws(c), "\n * error in connect_tips"))
        })
        cluster_result$optim_res$membership <- kmean_res$cluster
    }
    cluster_graph_res <- compute_partitions(cluster_result$g, cluster_result$optim_res, qval_thresh = qval_thresh, verbose = verbose)
    dimnames(cluster_graph_res$cluster_mat) <- dimnames(cluster_graph_res$num_links)
    valid_connection <- which(cluster_graph_res$cluster_mat < qval_thresh, arr.ind = TRUE)
    valid_connection <- valid_connection[apply(valid_connection, 1, function(x) {
        all(x %in% tip_pc_points)
    }), ]
    G <- cluster_graph_res$cluster_mat
    G[cluster_graph_res$cluster_mat < qval_thresh] <- -1
    G[cluster_graph_res$cluster_mat > 0] <- 0
    G <- -G
    if (all(G == 0, na.rm = TRUE)) {
        return(list(stree = igraph::get.adjacency(mst_g_old), Y = reducedDimK_old, G = G))
    }
    if (nrow(valid_connection) == 0) {
        return(list(stree = igraph::get.adjacency(mst_g_old), Y = reducedDimK_old, G = G))
    }
    mst_g <- mst_g_old
    diameter_dis <- igraph::diameter(mst_g_old)
    reducedDimK_df <- reducedDimK_old
    pb4 <- utils::txtProgressBar(max = length(nrow(valid_connection)), file = "", style = 3, min = 0)
    res <- stats::dist(t(reducedDimK_old))
    g <- igraph::graph_from_adjacency_matrix(as.matrix(res), weighted = TRUE, mode = "undirected")
    mst <- igraph::minimum.spanning.tree(g)
    max_node_dist <- max(igraph::E(mst)$weight)
    if (nrow(valid_connection) < 1) 
        warning("bad loop: nrow(valid_connection) < 1")
    for (i in 1:nrow(valid_connection)) {
        edge_vec <- sort(unique(cluster_result$optim_res$membership))[valid_connection[i, ]]
        edge_vec_in_tip_pc_point <- igraph::V(mst_g_old)$name[edge_vec]
        if (length(edge_vec_in_tip_pc_point) == 1) 
            next
        if (all(edge_vec %in% tip_pc_points) & (igraph::distances(mst_g_old, edge_vec_in_tip_pc_point[1], edge_vec_in_tip_pc_point[2]) >= 
            geodesic_distance_ratio * diameter_dis) & (euclidean_distance_ratio * max_node_dist > stats::dist(t(reducedDimK_old[, 
            edge_vec])))) {
            if (verbose) 
                message("edge_vec is ", edge_vec[1], "\t", edge_vec[2])
            if (verbose) 
                message("edge_vec_in_tip_pc_point is ", edge_vec_in_tip_pc_point[1], "\t", edge_vec_in_tip_pc_point[2])
            mst_g <- igraph::add_edges(mst_g, edge_vec_in_tip_pc_point)
        }
        utils::setTxtProgressBar(pb = pb4, value = pb4$getVal() + 1)
    }
    close(pb4)
    list(stree = igraph::get.adjacency(mst_g), Y = reducedDimK_df, G = G)
}

`prune_tree` <- function (stree_ori, stree_loop_closure, minimal_branch_len = 10) 
{
    if (ncol(stree_ori) < minimal_branch_len) 
        return(stree_loop_closure)
    dimnames(stree_loop_closure) <- dimnames(stree_ori)
    stree_ori[stree_ori != 0] <- 1
    stree_ori <- igraph::graph_from_adjacency_matrix(stree_ori, mode = "undirected", weighted = NULL)
    stree_loop_closure[stree_loop_closure != 0] <- 1
    stree_loop_closure <- igraph::graph_from_adjacency_matrix(stree_loop_closure, mode = "undirected", weighted = NULL)
    added_edges <- igraph::get.edgelist(stree_loop_closure - stree_ori)
    valid_edges <- matrix(ncol = 2, nrow = 0)
    edges_to_remove_df <- matrix(ncol = 2, nrow = 0)
    vertex_top_keep <- NULL
    if (nrow(added_edges) > 0) {
        edge_dists <- apply(added_edges, 1, function(x) igraph::distances(stree_ori, x[1], x[2]))
        valid_edges <- added_edges[which(edge_dists >= minimal_branch_len), , drop = FALSE]
        edges_to_remove_df <- added_edges[which(edge_dists < minimal_branch_len), , drop = FALSE]
    }
    if (nrow(valid_edges) > 0) {
        vertex_top_keep <- as.character(unlist(apply(valid_edges, 1, function(x) {
            igraph::shortest_paths(stree_ori, x[1], x[2])$vpath[[1]]$name
        })))
    }
    root_cell <- which(igraph::neighborhood.size(stree_ori) == 2)[1]
    if (is.na(root_cell)) {
        root_cell <- igraph::V(stree_ori)$name[1]
    }
    mst_traversal <- igraph::graph.dfs(stree_ori, root = root_cell, mode = "all", unreachable = FALSE, father = TRUE)
    mst_traversal$father <- as.numeric(mst_traversal$father)
    vertex_to_be_deleted <- c()
    if (length(mst_traversal$order) > 0) {
        for (i in 1:length(mst_traversal$order)) {
            curr_node <- tryCatch({
                mst_traversal$order[i]
            }, error = function(e) {
                NA
            })
            if (is.na(curr_node)) 
                next
            curr_node_name <- igraph::V(stree_ori)[curr_node]$name
            if (is.na(mst_traversal$father[curr_node]) == FALSE) {
                parent_node <- mst_traversal$father[curr_node]
                parent_node_name <- igraph::V(stree_ori)[parent_node]$name
                if (igraph::degree(stree_ori, v = parent_node_name) > 2) {
                  parent_neighbors <- igraph::neighbors(stree_ori, v = parent_node_name, mode = "all")
                  parent_neighbors_index <- sort(match(parent_neighbors$name, mst_traversal$order$name))
                  parent_neighbors <- mst_traversal$order$name[parent_neighbors_index]
                  tmp <- igraph::delete.edges(stree_ori, paste0(parent_node_name, "|", parent_neighbors))
                  tmp_decomposed <- igraph::decompose.graph(tmp)
                  comp_a <- tmp_decomposed[unlist(lapply(tmp_decomposed, function(x) {
                    parent_neighbors[2] %in% igraph::V(x)$name
                  }))][[1]]
                  comp_b <- tmp_decomposed[unlist(lapply(tmp_decomposed, function(x) {
                    parent_neighbors[3] %in% igraph::V(x)$name
                  }))][[1]]
                  diameter_len_a <- igraph::diameter(comp_a) + 1
                  diameter_len_b <- igraph::diameter(comp_b) + 1
                  if (diameter_len_a < minimal_branch_len) {
                    vertex_to_be_deleted <- c(vertex_to_be_deleted, igraph::V(comp_a)$name)
                  }
                  if (diameter_len_b < minimal_branch_len) {
                    vertex_to_be_deleted <- c(vertex_to_be_deleted, igraph::V(comp_b)$name)
                  }
                }
            }
        }
    }
    valid_vertex_to_be_deleted <- setdiff(vertex_to_be_deleted, vertex_top_keep)
    stree_loop_closure <- igraph::delete_vertices(stree_loop_closure, valid_vertex_to_be_deleted)
    tmp <- edges_to_remove_df[edges_to_remove_df[, 1] %in% igraph::V(stree_loop_closure)$name & edges_to_remove_df[, 2] %in% 
        igraph::V(stree_loop_closure)$name, , drop = FALSE]
    if (nrow(tmp) > 0) {
        edges_to_remove <- paste0(tmp[, 1], "|", tmp[, 2])
        stree_loop_closure <- igraph::delete.edges(stree_loop_closure, edges_to_remove)
    }
    return(igraph::get.adjacency(stree_loop_closure))
}

`findNearestPointOnMST` <- function (cds, reduction_method, rge_res_Y) 
{
    dp_mst <- principal_graph(cds)[[reduction_method]]
    dp_mst_list <- igraph::decompose.graph(dp_mst)
    if (length(unique(cds@clusters[[reduction_method]]$partitions)) != length(dp_mst_list)) {
        dp_mst_list <- list(dp_mst)
    }
    closest_vertex_df <- NULL
    cur_start_index <- 0
    if (length(dp_mst_list) < 1) 
        warning("bad loop: length(dp_mst_list) < 1")
    for (i in 1:length(dp_mst_list)) {
        cur_dp_mst <- dp_mst_list[[i]]
        if (length(dp_mst_list) == 1) {
            Z <- t(SingleCellExperiment::reducedDims(cds)[[reduction_method]])
        }
        else {
            Z <- t(SingleCellExperiment::reducedDims(cds)[[reduction_method]])[, cds@clusters[[reduction_method]]$partitions == 
                i]
        }
        Y <- rge_res_Y[, igraph::V(cur_dp_mst)$name]
        tip_leaves <- names(which(igraph::degree(cur_dp_mst) == 1))
        closest_vertex_ori <- find_nearest_vertex(Z, Y)
        closest_vertex <- closest_vertex_ori + cur_start_index
        closest_vertex_names <- colnames(Y)[closest_vertex_ori]
        cur_name <- names(closest_vertex)
        closest_vertex <- as.matrix(closest_vertex)
        row.names(closest_vertex) <- cur_name
        closest_vertex_df <- rbind(closest_vertex_df, closest_vertex)
        cur_start_index <- cur_start_index + igraph::vcount(cur_dp_mst)
    }
    closest_vertex_df <- closest_vertex_df[colnames(cds), , drop = FALSE]
    cds@principal_graph_aux[[reduction_method]]$pr_graph_cell_proj_closest_vertex <- closest_vertex_df
    cds
}

`projPointOnLine` <- function (point, line) 
{
    ap <- point - line[, 1]
    ab <- line[, 2] - line[, 1]
    res <- line[, 1] + c((ap %*% ab)/(ab %*% ab)) * ab
    return(res)
}

`get_call_stack` <- function () 
{
    cv <- as.vector(sys.calls())
    lcv <- length(cv)
    n <- lcv - 1
    ocv <- vector()
    for (i in seq(1, n, 1)) {
        elem <- stringr::str_split(as.character(cv[i]), "[(]", n = 2)[[1]][[1]]
        ocv <- c(ocv, elem)
    }
    return(ocv)
}

`search_nn_index` <- function (query_matrix, nn_index, k = 25, nn_control = list(), verbose = FALSE) 
{
    assertthat::assert_that(methods::is(query_matrix, "matrix") || is_sparse_matrix(query_matrix), msg = paste0("make_nn_matrix: the query_matrix object must be of type matrix"))
    assertthat::assert_that(assertthat::is.count(k))
    nn_control_default <- get_global_variable("nn_control_annoy_euclidean")
    nn_control <- set_nn_control(mode = 2, nn_control = nn_control, nn_control_default = nn_control_default, nn_index = nn_index, 
        k = k, verbose = verbose)
    nn_method <- nn_control[["method"]]
    assertthat::assert_that(nn_method %in% c("annoy", "hnsw"), msg = paste0("search_nn_index: unsupported nearest neighbor index type '", 
        nn_method, "'."))
    if (verbose) {
        message("search_nn_index:")
        message("  k: ", k)
        report_nn_control("  nn_control: ", nn_control)
        tick("search_nn_index: search time:")
    }
    k <- min(k, nrow(query_matrix))
    cores <- nn_control[["cores"]]
    if (nn_method == "nn2") {
        stop("search_nn_index is not valid for method nn2")
    }
    else if (nn_method == "annoy") {
        if (!test_annoy_index(nn_index = nn_index, verbose = verbose)) {
            stop("search_nn_index: the annoy nearest neighbor does not exist.")
        }
        num_row <- nrow(query_matrix)
        idx <- matrix(nrow = num_row, ncol = k)
        dists <- matrix(nrow = num_row, ncol = k)
        metric <- nn_control[["metric"]]
        search_k <- nn_control[["search_k"]]
        if (cores <= 1) {
            nn_res <- search_nn_annoy_index(query_matrix = query_matrix, nn_index = nn_index, metric = metric, k = k, search_k = search_k, 
                beg_row_index = 1, end_row_index = num_row)
            if (nn_res[["num_bad"]]) 
                stop("annoy was unable to find ", k, " nearest neighbors for ", nn_res[["num_bad"]], " rows. You may need to increase the n_trees and/or search_k parameter values.")
            nn_res <- list(nn.idx = nn_res[["idx"]], nn.dists = nn_res[["dists"]])
        }
        else {
            omp_num_threads <- get_global_variable("omp_num_threads")
            blas_num_threads <- get_global_variable("blas_num_threads")
            RhpcBLASctl::omp_set_num_threads(1L)
            RhpcBLASctl::blas_set_num_threads(1L)
            if (cores > num_row) 
                cores <- num_row
            tasks <- tasks_per_block(num_row, cores)
            beg_block <- c(0, cumsum(tasks))[1:cores] + 1
            end_block <- cumsum(tasks)
            nn_blocks <- list()
            inplan <- future::plan()
            future::plan(future::multicore, workers = cores)
            on.exit(future::plan(inplan), add = TRUE)
            for (iblock in seq(cores)) {
                nn_blocks[[iblock]] <- future::future({
                  search_nn_annoy_index(query_matrix = query_matrix, nn_index = nn_index, metric = metric, k = k, search_k = search_k, 
                    beg_row_index = beg_block[[iblock]], end_row_index = end_block[[iblock]])
                })
            }
            tot_bad <- 0
            for (iblock in seq(cores)) {
                nn_res <- future::value(nn_blocks[[iblock]])
                if (nrow(nn_res[["idx"]]) != end_block[[iblock]] - beg_block[[iblock]] + 1) {
                  stop("bad row count in nn_res")
                }
                idx[beg_block[[iblock]]:end_block[[iblock]], ] <- nn_res[["idx"]]
                dists[beg_block[[iblock]]:end_block[[iblock]], ] <- nn_res[["dists"]]
                tot_bad <- tot_bad + nn_res[["num_bad"]]
            }
            if (tot_bad) 
                stop("annoy was unable to find ", k, " nearest neighbors for ", tot_bad, " rows. You may need to increase the n_trees and/or search_k parameter values.")
            nn_res <- list(nn.idx = idx, nn.dists = dists)
            RhpcBLASctl::omp_set_num_threads(as.integer(omp_num_threads))
            RhpcBLASctl::blas_set_num_threads(as.integer(blas_num_threads))
        }
    }
    else if (nn_method == "hnsw") {
        if (!test_hnsw_index(nn_index = nn_index, verbose = verbose)) {
            stop("search_nn_index: the hnsw nearest neighbor does not exist.")
        }
        assertthat::assert_that(nn_control[["ef"]] >= k, msg = paste0("search_nn_index: ef must be >= k"))
        tmp <- RcppHNSW::hnsw_search(X = query_matrix, ann = nn_index[["hnsw_index"]], k = k, ef = nn_control[["ef"]], verbose = verbose, 
            n_threads = cores, grain_size = nn_control[["grain_size"]])
        nn_res <- list(nn.idx = tmp[["idx"]], nn.dists = tmp[["dist"]])
    }
    else stop("search_nn_index: unsupported nearest neighbor index type '", nn_method, "'")
    if (verbose) {
        tock()
    }
    return(nn_res)
}

`count_nn_missing_self_index` <- function (nn_res, verbose = FALSE) 
{
    idx <- nn_res[["nn.idx"]]
    dst <- nn_res[["nn.dists"]]
    len <- length(idx[1, ])
    num_missing <- 0
    if (nrow(idx) == 0) 
        return(0)
    for (irow in 1:nrow(idx)) {
        vidx <- idx[irow, ]
        vdst <- dst[irow, ]
        if (vidx[[1]] != irow) {
            dself <- FALSE
            dzero <- TRUE
            if (length(vidx) == 1) {
                num_missing <- num_missing + 1
                next
            }
            for (i in seq(1, len, 1)) {
                if (vidx[[i]] == irow) {
                  dself = TRUE
                  break
                }
                if (vdst[[i]] > .Machine$double.xmin) {
                  dzero <- FALSE
                }
            }
            if (dself == FALSE && dzero == FALSE) {
                num_missing <- num_missing + 1
            }
        }
    }
    if (verbose) {
        message("count_nn_missing_self_index:")
        message("  ", num_missing, " out of ", nrow(nn_res[["nn.idx"]]), " rows are missing the row index")
        message("  'recall': ", formatC((as.double(nrow(nn_res[["nn.idx"]]) - num_missing)/as.double(nrow(nn_res[["nn.idx"]]))) * 
            100), "%")
    }
    return(num_missing)
}

`check_nn_col1` <- function (mat) 
{
    irow <- seq(nrow(mat))
    if (any(mat[, 1] != irow)) {
        return(FALSE)
    }
    return(TRUE)
}

`repmat` <- function (X, m, n) 
{
    mx = dim(X)[1]
    nx = dim(X)[2]
    matrix(t(matrix(X, mx, nx * n)), mx * m, nx * n, byrow = TRUE)
}

`soft_assignment` <- function (X, C, sigma) 
{
    D <- nrow(X)
    N <- ncol(X)
    K <- ncol(C)
    norm_X_sq <- repmat(t(t(colSums(X^2))), 1, K)
    norm_C_sq <- repmat(t(colSums(C^2)), N, 1)
    dist_XC <- norm_X_sq + norm_C_sq - 2 * t(X) %*% C
    min_dist <- apply(dist_XC, 1, min)
    dist_XC <- dist_XC - repmat(t(t(min_dist)), 1, K)
    Phi_XC <- exp(-dist_XC/sigma)
    P <- Phi_XC/repmat(t(t(rowSums(Phi_XC))), 1, K)
    obj <- -sigma * sum(log(rowSums(exp(-dist_XC/sigma))) - min_dist/sigma)
    return(list(P = P, obj = obj))
}

`generate_centers` <- function (X, W, P, param.gamma) 
{
    D <- nrow(X)
    N <- nrow(X)
    K <- ncol(W)
    Q <- 2 * (diag(colSums(W)) - W) + param.gamma * diag(colSums(P))
    B <- param.gamma * X %*% P
    C <- B %*% solve(Q)
    return(C)
}

`test_annoy_index` <- function (nn_index, verbose = FALSE) 
{
    res <- TRUE
    index_obj <- NULL
    if (!is.null(nn_index[["annoy_index"]])) {
        index_obj <- nn_index[["annoy_index"]]
    }
    else if (!is.null(nn_index[["ann"]])) {
        index_obj <- nn_index[["ann"]]
    }
    else {
        if (verbose) {
            cs <- get_call_stack_as_string()
            message("test_annoy_index: the annoy nearest neighbor does not exist\ncall stack: ", cs)
        }
        return(FALSE)
    }
    dist_res <- tryCatch(index_obj$getDistance(0, 1), error = function(emsg) {
        if (verbose) {
            cs <- get_call_stack_as_string()
            message("test_annoy_index: the annoy nearest neighbor does not exist\ncall stack: ", cs)
        }
        res <<- FALSE
    })
    return(res)
}

`search_nn_annoy_index` <- function (query_matrix, nn_index, metric, k, search_k, beg_row_index, end_row_index) 
{
    assertthat::assert_that(beg_row_index <= end_row_index, msg = paste0("search_nn_annoy_index: beg_row_index must be <= end_row_index"))
    if (!is.null(nn_index[["version"]])) 
        annoy_index <- nn_index[["annoy_index"]]
    else if (!is.null(nn_index[["type"]])) 
        annoy_index <- nn_index[["ann"]]
    else annoy_index <- nn_index
    nrow <- end_row_index - beg_row_index + 1
    idx <- matrix(nrow = nrow, ncol = k)
    dists <- matrix(nrow = nrow, ncol = k)
    num_bad <- 0
    offset <- beg_row_index - 1
    for (i in seq(1, nrow)) {
        nn_list <- annoy_index$getNNsByVectorList(query_matrix[i + offset, ], k, search_k, TRUE)
        if (length(nn_list$item) != k) 
            num_bad <- num_bad + 1
        idx[i, ] <- nn_list$item
        dists[i, ] <- nn_list$distance
    }
    if (metric == "cosine") {
        dists <- 0.5 * dists * dists
    }
    return(list(idx = idx + 1, dists = dists, num_bad = num_bad))
}

`tasks_per_block` <- function (ntask = NULL, nblock = NULL) 
{
    tasks_block <- rep(trunc(ntask/nblock), nblock)
    remain <- ntask%%nblock
    if (remain) 
        for (i in seq(remain)) tasks_block[i] <- tasks_block[i] + 1
    return(tasks_block)
}

`test_hnsw_index` <- function (nn_index, verbose = FALSE) 
{
    res <- TRUE
    if (is.null(nn_index[["hnsw_index"]])) {
        if (verbose) {
            cs <- get_call_stack_as_string()
            message("test_hnsw_index: the hnsw nearest neighbor does not exist\ncall stack: ", cs)
        }
        return(FALSE)
    }
    size_res <- tryCatch(nn_index[["hnsw_index"]]$size(), error = function(emsg) {
        if (verbose) {
            cs <- get_call_stack_as_string()
            message("test_hnsw_index: the hnsw nearest neighbor does not exist\ncall stack: ", cs)
        }
        res <<- FALSE
    })
    return(res)
}

