# Runtime entry point. Only matrix-based PCA/UMAP trajectories are supported.
.source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))
.core_file <- normalizePath(tail(.source_files, 1)[[1]])
.core_dir <- dirname(.core_file)
trajectory_packages <- c('Matrix','tidyr','data.table','SingleCellExperiment',
 'assertthat','digest','dplyr','igraph','irlba','leidenbase','pbapply',
 'proxy','RANN','RcppAnnoy','RcppHNSW','RhpcBLASctl','uwot','openssl','plyr','future')
missing <- trajectory_packages[!vapply(trajectory_packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Missing trajectory packages: ',paste(missing,collapse=', '))
suppressPackageStartupMessages({
 library(SingleCellExperiment)
 library(Matrix)
})
`%>%` <- dplyr::`%>%`
`.` <- plyr::`.`
if(!methods::isClass('mmsage_cell_data_set')) methods::setClass(
 'mmsage_cell_data_set',contains='SingleCellExperiment',
 slots=c(reduce_dim_aux='SimpleList',principal_graph_aux='SimpleList',
 principal_graph='SimpleList',clusters='SimpleList'))
source(file.path(.core_dir,'settings.R'))
source(file.path(.core_dir,'upstream.R'))
exprs <- function(x) SingleCellExperiment::counts(x)
principal_graph <- function(x) x@principal_graph
`principal_graph<-` <- function(x,value) {x@principal_graph<-value;x}
principal_graph_aux <- function(x) x@principal_graph_aux
`principal_graph_aux<-` <- function(x,value) {x@principal_graph_aux<-value;x}
clusters <- function(x,reduction_method='UMAP') x@clusters[[reduction_method]]$clusters[colnames(x)]
partitions <- function(x,reduction_method='UMAP') x@clusters[[reduction_method]]$partitions[colnames(x)]
pseudotime <- function(x,reduction_method='UMAP') x@principal_graph_aux[[reduction_method]]$pseudotime[colnames(x)]
# Exact unweighted edge construction used by M4's default Leiden path.
jaccard_coeff <- function(R_idx,R_weight) {
 if(R_weight) stop('Weighted graph construction is outside the M4 runtime contract')
 cbind(rep(seq_len(nrow(R_idx)),each=ncol(R_idx)),as.vector(t(R_idx)),1)
}
pnorm_over_mat <- function(R_num_links_ij,R_var_null_num_links) {
 matrix(stats::pnorm(as.vector(R_num_links_ij),sd=sqrt(as.vector(R_var_null_num_links)),lower.tail=FALSE),nrow=nrow(R_num_links_ij))
}
# Unsupported branches fail explicitly; the production pipeline uses dgCMatrix,
# PCA, UMAP, unweighted Leiden and explicit roots.
.unsupported <- function(...) stop('This operation is outside the MMSage matrix trajectory runtime')
for(.name in c('select_trajectory_roots','bpcells_prcomp_irlba','estimate_sf_bpcells',
 'set_cds_row_order_matrix','set_matrix_control_pca','set_matrix_class','rm_bpcells_dir',
 'bpcells_find_base_matrix')) assign(.name,.unsupported)
