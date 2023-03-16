#' @import utils
#' @import stats
#' @import ggplot2
#' @import patchwork
#' @import dplyr 
#' @import tictoc

NULL

# Calculation ----
#' Generate a random variable with a defined correlation to an existing variables 
#' This function is adapted from the [post](https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables)
#' @param y The reference vector.
#' @param rho The target correlation
#' @param x A vector to start with. If specified, it must have the same length as `y`
#' @export
#' @importFrom magrittr `%>%`
#' 
complementCor = function(y, rho, x) {
  if (missing(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
  y.perp <- residuals(lm(x ~ y))                                                                                                                                                                                                                                              
  rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
}


calc_err = function(origin_mt, nmf_obj){
  if(is.matrix(origin_mt)){
    err_mt = origin_mt - prod(nmf_obj)
    err_value = sqrt(sum(err_mt ^ 2)) / length(origin_mt)
  } else { # for Matrix
    prod_nmf = as(prod(nmf_obj), 'dgCMatrix')
    err_mt = origin_mt - prod_nmf
    err_value = sqrt(sum(err_mt ^ 2)) / length(origin_mt)
  }
  return(err_value)
}

calc_jac = function(cor_mt, threshold = 0.9){
  n_intersect = sum(apply(cor_mt, 1, max) > threshold)
  return(
    n_intersect / (ncol(cor_mt) + nrow(cor_mt) - n_intersect)
  )
}


calc_jacind_from_cor = function(x_cor_mt, y_cor_mt){
  uni_genes = union(rownames(x_cor_mt), rownames(y_cor_mt))
  x_cor_exp_mt = rbind(x_cor_mt, matrix(NA, nrow = length(setdiff(
    uni_genes, rownames(x_cor_mt)
  )), ncol = ncol(x_cor_mt)))
  rownames(x_cor_exp_mt)[rownames(x_cor_exp_mt) == ""] = setdiff(uni_genes, rownames(x_cor_mt))
  colnames(x_cor_exp_mt) = paste0("x_", colnames(x_cor_exp_mt))
  
  y_cor_exp_mt = rbind(y_cor_mt, matrix(NA, nrow = length(setdiff(
    uni_genes, rownames(y_cor_mt)
  )), ncol = ncol(y_cor_mt)))
  rownames(y_cor_exp_mt)[rownames(y_cor_exp_mt) == ""] = setdiff(uni_genes, rownames(y_cor_exp_mt))
  colnames(y_cor_exp_mt) = paste0("y_", colnames(y_cor_exp_mt))
  
  # merge_w_mt = cbind(, y_cor_exp_mt[uni_genes, ])
  jac_mt = cor(
    x = x_cor_exp_mt[uni_genes,],
    y = y_cor_exp_mt[uni_genes,],
    use = 'pairwise.complete.obs',
    method = 'pearson'
  )
  
  return(calc_jac(jac_mt))
}

#' Quantile normalization
#'
#' Perform quantile normalization column-wisely on a matrix.
#' @param mt Matrix to norm by columns
#' @param verbose Logical,  print the log.
#' @param n_cores Number of cores to use. (to be added)
#' @return A quantile-normalized matrix.
#' @export
#'
quantile_normalize = function(mt, n_cores = 30, verbose = T){
  # can adapt the qnorm Rcpp from preprocessCore::normalize.quantiles,
  # however, noted that this function yield a slightly different results 
  # that in some columns the minimum values would differ from others.
  # Try S3 method in the future
  if(is.matrix(mt)) {
    if(verbose) message('Calculate ranks...')
    mt_rank = apply(mt, 2, rank, ties.method = "min")
    if(verbose) message('Sorting...')
    mt_sorted = apply(mt, 2, sort)
    mt_mean = rowMeans(mt_sorted)
    if(verbose) message('Normalizing...')
    mt_final = apply(mt_rank, 2, function(x) {
      mt_mean[x]
    })
    rownames(mt_final) = rownames(mt)
    return(mt_final)
    
  } else {
    if(verbose) message('Treat the input as a dgCMatrix...')
    mt = as(mt, 'dgCMatrix')
    # The 0s in the sparse matrix should be recognized as 1
    if(verbose) message('Calculate ranks...')
    n_nonzero = diff(mt@p)
    n_genes = nrow(mt)
    x_lst = split(mt@x, f = rep.int(1:ncol(mt), n_nonzero))
    rank_lst = parallel::mclapply(x_lst, mc.cores = n_cores, FUN = function(x_j){
      x_j_rank = rank(c(x_j, rep(0, n_genes - length(x_j))), ties.method = 'min')
      return(head(x_j_rank, length(x_j)))
    })
    # MM_rank = mt
    # MM_rank@x = as.double(unlist(rank_lst))
    # rm(rank_lst)
    
    if(verbose) message('Sorting...')
    sort_lst = parallel::mclapply(x_lst, mc.cores = n_cores, FUN = function(x_j){
      x_j_sorted = sort(c(x_j, rep(0, n_genes - length(x_j))))
      return(tail(x_j_sorted, length(x_j)))
    })
    idx_lst = parallel::mclapply(x_lst, mc.cores = n_cores, FUN = function(x_j){
      x_i = tail(0:(n_genes-1), length(x_j))
      return(x_i)
    })
    MM_sorted = mt
    rownames(MM_sorted) = NULL
    MM_sorted@x = unlist(sort_lst)
    MM_sorted@i = unlist(idx_lst)
    rm(x_lst, sort_lst, idx_lst)
    MM_mean = Matrix::rowMeans(MM_sorted)
    
    if(verbose) message('Normalizing...')
    # all.equal(mt_mean, MM_mean)
    reassign_lst = parallel::mclapply(rank_lst, mc.cores = n_cores, function(x_rank){
      return(MM_mean[x_rank])
    })
    MM_final = mt
    MM_final@x = unlist(reassign_lst)
    # all.equal(as.matrix(MM_final), mt_final)
    
    return(MM_final)
  }
}



sweep_MM <- function(x, margin, stats, fun = "*") {
# adapted from https://stackoverflow.com/questions/55407656/r-sweep-on-a-sparse-matrix
  x = as(x, 'dgTMatrix')
  f <- match.fun(fun)
  if (margin == 1) {
     idx <- x@i + 1
  } else {
     idx <- x@j + 1
  }
  x@x <- f(x@x, stats[idx])
  return(x)
}

#' Calculate normalized expresssion for UMI-based data
#'
#' 'calc_normExpr.Droplet' is a wrapper to calculate the normalized expression for droplet data.
#' It would calculate the normalized expresssion as: 'x/sum(x) * factor' for each sample.
#' @param x Matrix or sparse matrix. Each column specify a sample.
#' @param factor A number define the resulting library size for each sample after the normalization. Default 10000.
#' @param logarithmetic Logical. Should the results be log-transformed?
#' @return Matrix
#' @export
#'
calc_normExpr.Droplet = function(x,
                                 factor = 10000,
                                 logarithmetic = T) {
  if (is.matrix(x)) {
    y = apply(x, 2, function(i) {
      i / sum(i) * factor
    })
  } else {
    message('Try to deal as sparse matrices.')
    x = as(x, 'dgTMatrix')
    y = sweep_MM(x, 2, Matrix::colSums(x), fun = function(x,y){x/y * factor})
  }
  
  if (logarithmetic) {
    if(is.matrix(x)){
      return(log2(y + 1))
    } else {
      return(as(log2(y+1), 'dgTMatrix'))
    }
  } else {
    return(y)
  }
}

#' Calculate normalized expresssion for full-length RNA data
#'
#' 'calc_normExpr.TPM' is a wrapper to calculate the normalized expression for count data derived from full-length protocol.
#' It would calculate the normalized expresssion as TPM for each sample.
#' @param x Matrix or sparse matrix. Each column specify a sample.
#' @param logarithmetic Logical. Should the results be log-transformed?
#' @return Matrix
#' @export
#'
calc_normExpr.TPM = function(x, logarithmetic = T){
  y = apply(x, 2, FUN = calc_TPM)
  if(logarithmetic){
    return(log2(y+1))
  } else {
    return(y)
  }
}

#' Calculate TPM
#'
#' Calculate TPM, the length information is provided by R package `annotables`
#' @param x Named vector of gene counts
#' @return Named vector of TPM
#'
calc_TPM = function(x){
 # require(annotables)
  length_vec = setNames(nm = annotables::grch38$symbol, 
                        abs(annotables::grch38$start - annotables::grch38$end))
  length_vec = length_vec/1000 # per kilobass
  x_length = length_vec[names(x)]
  x_length = x_length[!is.na(x_length)]
  
  if(length(x_length) != length(x)){
    inter_genes = intersect(names(x_length), names(x))
    # warning(sprintf("Found length of %d / %d genes!", length(inter_genes), length(x))) 
    x = x[inter_genes]
    x_length = x_length[inter_genes]
  }
  
  RPK = x / x_length
  scaling_factor = sum(RPK)/1000000
  TPM = RPK/scaling_factor
  return(TPM)
}


#' Select features for single-cell data
#' 
#' Using differential expression to select feature genes for single-cell data. 
#' @param X Matrix, genes x cells, a library-size-normalized matrix or a UMI count matrix or a TPM matrix.
#' @param celltypes Cell-named vector to specify
#' @param batches Cell-named vector to specify the batch information, default = NULL, not run
#' @param method String. Can be 'ttest' or 'wilcox'. Default = 'wilcox'
#' @param n_cores Integer. Default using all cores.
#' @return A differential expression table.
#' @export
#' @importFrom magrittr `%>%`
#' 
select_features_sc = function(X, celltypes, batches = NULL, method = 'wilcox', n_cores = 0){
  # require(JasonToolBox)
  if(n_cores == 0)
    n_cores = parallel::detectCores()
  
  stopifnot(all(colnames(X) == names(celltypes)))
  
  if(!is.null(batches))
    stopifnot(all(colnames(X) == names(batches)))
  
  chosen_fun = switch (method,
    wilcox = dea.wilcox, 
    ttest = dea.ttest, 
    dea.wilcox
  )
  
  if(!is.null(batches)){
    message(sprintf('Run for %d batches.', length(unique(batches))))
    celltypes_de_lst = list()
    for(i_batch in unique(batches)){
      message(sprintf('Now running batch: %s', i_batch))
      tictoc::tic()
      # t_grp1_ind = celltypes == i_celltype
      # t_grp2_ind = celltypes != i_celltype
      t_celltypes = celltypes[batches == i_batch]
      t_batch_X = X[, batches == i_batch]
      gene_idx = Matrix::rowSums(t_batch_X > 0) > round(0.1* ncol(t_batch_X))
      message(sprintf('%d genes to test.', sum(gene_idx)))
      t_batch_X = t_batch_X[gene_idx, ]
      t_type_de_lst = list()
      for(i_celltype in unique(t_celltypes)){
        message(sprintf('Now running type: %s', i_celltype))
        t_de_tbl = chosen_fun(
          mt1_origin = t_batch_X[, t_celltypes == i_celltype, drop = F],
          mt2_origin = t_batch_X[, t_celltypes != i_celltype, drop = F],
          down_size = 5000, n_cores = n_cores, onlyPos = T)
        if(!is.null(t_de_tbl)){
          t_type_de_lst[[i_batch]] = t_de_tbl %>% mutate(batchID = i_batch, target_celltype = i_celltype) # %>% dplyr::filter(adj.P.Val < 0.05)
        } 
        invisible(gc())
      }
      t_type_de_tbl = bind_rows(t_type_de_lst)
      # t_remained_genes_tbl = t_type_de_tbl %>% group_by(id) %>% 
      #   summarise(Log2FC = max(Log2FC), min_signif = min(adj.P.Val), n_signif = length(unique(batchID))) 
      celltypes_de_lst[[i_celltype]] = t_type_de_tbl
      tictoc::toc()
    }
    celltypes_de_tbl = bind_rows(celltypes_de_lst)
    
  } else {
    message('Run for no batches.')
    celltypes_de_lst = list()
    for (i_celltype in unique(celltypes)) {
      message(sprintf('Now running type: %s', i_celltype))
      # target_cells = names(celltypes)[celltypes == i_celltype]
      t_de_tbl = chosen_fun(
        mt1_origin = X[, celltypes == i_celltype],
        mt2_origin = X[, celltypes != i_celltype],
        down_size = 10000,
        n_cores = n_cores,
        onlyPos = T
      )
      celltypes_de_lst[[i_celltype]] = t_de_tbl %>% mutate(target_celltype = i_celltype)
    }
    celltypes_de_tbl = bind_rows(celltypes_de_lst)
  }
  
  return(celltypes_de_tbl)
}

#' Differential expression analysis using wilcox test
#' 
#' @importFrom magrittr `%>%`
dea.wilcox = function (mt1_origin, mt2_origin, down_size = 5000, n_cores = 1, onlyPos = T,...){
  # mt1_origin = x[[1]]
  # mt2_origin = x[[2]]
  args = list(...)
  # progress = ifelse(is.null(args$progress), T, args$progress)
  # if mt is too large, downsample to 10000
  if(ncol(mt1_origin) < 3 || ncol(mt2_origin) < 3) {
    return(NULL)
  }
  
  message(sprintf('mt1: %d x %d; mt2: %d x %d', nrow(mt1_origin), ncol(mt1_origin), nrow(mt2_origin), ncol(mt2_origin)))
  
  mt1 = mt1_origin
  if(ncol(mt1) > down_size){
    mt1 = mt1_origin[, sample(1:ncol(mt1_origin), size = down_size)]
  }
  mt2 = mt2_origin
  if(ncol(mt2) > down_size){
    mt2 = mt2_origin[, sample(1:ncol(mt2_origin), size = down_size)]
  }
  
  stopifnot(nrow(mt1) == nrow(mt2))
  if(max(mt1) > 100){
    message('Not log-transformed? Transfering now.')
    mt1 = log1p(mt1)
    mt2 = log1p(mt2)
  }
  
  # message('Run')
  tbl_lst = pbmcapply::pbmclapply(mc.cores = n_cores, 1:nrow(mt1), function(i){
    X = mt1[i, ]
    Y = mt2[i, ]
    t_res = wilcox.test(X, Y)
    t_p_val = t_res$p.value
    t_mean1 = mean(expm1(X) + 1)
    t_mean2 = mean(expm1(Y) + 1)
    t_LFC = log2(t_mean1) - log2(t_mean2)
    return(tibble(id = rownames(mt1)[i], mean_grp1 = t_mean1, mean_grp2 = t_mean2, Log2FC = t_LFC, P.Val = t_p_val))
  })
  res_tbl = bind_rows(tbl_lst) %>% mutate(adj.P.Val = p.adjust(P.Val, method = 'BH'))
  if(onlyPos){
    res_tbl = res_tbl %>% dplyr::filter(Log2FC > 0)
  }
  return(res_tbl)
}


#' Select features for bulk data
#' 
#' Choose highly variable genes, using 'vst' algorithm adapted from Seurat.
#' @param X a library-size-normalized matrix or a UMI count matrix or a TPM matrix.
#' @param loess_span smooth factor for loess, default = 0.3
#' @param high_quantile quantile for selecting the variable genes, 
#' @return A list object with two named field: *fit* for loess fit object; *hvg_df* for a dataframe with highly variable genes.
#' @export
#' 
select_features = function(X, loess_span = 0.3, high_quantile = 0.5){
  # require(JasonToolBox)
  
  hvg_df = data.frame(featureID = rownames(X), means = rowMeans(X))
  if(is.matrix(X)) {
    hvg_df$vars = matrixStats::rowVars(X)
    hvg_df$notzeros = hvg_df$vars > 0
  } else {
    message('Calculating vars for the sparse matrix. This may take a while.')
    hvg_df$vars = JasonToolBox::apply_MM(X, 1, function(x){var(c(x, rep(0, ncol(X) - length(x))))})
    hvg_df$notzeros = hvg_df$vars > 0
  }

  fit = loess(formula = log10(vars) ~ log10(means),
              data = hvg_df[hvg_df$notzeros, ],
              span = loess_span)

  selected_genes = fit$residuals > quantile(fit$residuals[fit$residuals > 0], high_quantile)
  
  hvg_df$var_expected = 0
  hvg_df$var_expected[hvg_df$notzeros] = 10^fit$fitted
  hvg_df$selected = selected_genes[hvg_df$featureID]
  hvg_df$selected[is.na(hvg_df$selected)] = F
  
  return(list(fit = fit, hvg_df = hvg_df))
}



# #' Merge cells into metacell by clusters (deprecated)
# #' 
# #' Merge single cell data into meta cells. To be deprecated.
# #' @param expr_mt Expression matrix
# #' @param celltypes should be a cell-id-named vector of cell types
# #' @param target_number Number of cells to merge into one, default = 30.
# #' @param k_beta Numeric
# #' @param double_check Logical,
# #' @param n_cores Use when double_check, default using all cores available
# #' @param n_perm Use when double_check, default = 100
# #' @return A named list
# #' 
# mergeCells = function(expr_mt,
#                       celltypes,
#                       target_number = 30,
#                       K,
#                       k_beta = 2,
#                       double_check = F,
#                       ...) {
#   
#   # require(RcppParallel)
#   # require(proxyC)
#   # require(pbapply)
#   
#   args = list(...)
#   n_cores = ifelse(is.null(args$n_cores), RcppParallel::defaultNumThreads(), args$n_cores)
#   n_perm = ifelse(is.null(args$n_perm), 100, args$n_perm)
#   batch_num = ifelse(is.null(args$batch_num), 10000, args$batch_num)
#   
#   checkmate::assert(all(colnames(expr_mt) == names(celltypes)))
#   
#   RcppParallel::setThreadOptions(n_cores)
#   
#   message('Merging for metacells')
#   ttl_metacell_lst = list()
#   r = 1
#   redefine_celltypes = c()
#   tic.clear()
#   for (i_celltype in unique(celltypes)) {
#     n_cells = sum(celltypes == i_celltype)
#     target_type_expr_mt = expr_mt[, names(celltypes)[celltypes == i_celltype]]
#     
#     step_seq = unique(c(seq(1, ncol(target_type_expr_mt), by = batch_num), ncol(target_type_expr_mt)+1))
#     
#     tictoc::tic(msg = i_celltype)
#     for(i_step in 2:length(step_seq)){
#       i_idx1 = step_seq[i_step-1]
#       i_idx2 = step_seq[i_step]-1
#       
#       i_n_cells = i_idx2 - i_idx1 + 1
#       message('Batch processing: ', i_n_cells, ' cells of ', i_celltype)
#       i_target_mt = target_type_expr_mt[, i_idx1:i_idx2]
#       lognorm_mt = log1p(i_target_mt)
#     
#       cor_mat = proxyC::simil(lognorm_mt, method = 'cosine', margin = 2)
#       diag(cor_mat) = 0
#     
#       message('Calculating ranks')
#       rank_mat = pbapply::pbapply(cor_mat, 2, function(x) {
#         y = rank(-x, ties.method = 'max')
#         return(y)
#       })
#       S_mat = rank_mat * t(rank_mat)
#       S_mat = K ^ 2 - S_mat
#       S_mat[S_mat < 0] = 0
#       message('Regularizing ranks for in-node edges')
#       S_mat_i = t(pbapply::pbapply(S_mat, 1, function(x) {
#         y = k_beta * K - rank(-x, ties.method = 'max')
#         y[y < 0] = 0
#         return(y)
#       }))
#       message('Regularizing ranks for out-node edges')
#       adj_mt = pbapply::pbapply(S_mat_i, 2, function(x) {
#         y = K - rank(-x, ties.method = 'max')
#         y[y < 0] = 0
#         return(y)
#       })
#       leiden_obj = leiden::leiden(adj_mt, max_comm_size = target_number)
#       
#       lognorm_meta_mt = log1p(t(apply(i_target_mt, 1, function(x) {
#         tapply(x, INDEX = leiden_obj, mean)
#       })))
#       
#       tmp_expr_meta_mt = expm1(lognorm_meta_mt)
#       colnames(tmp_expr_meta_mt) = sprintf(
#         'metacell_%s_%d_%05.f',
#         stringr::str_replace(i_celltype, ' ', ''),
#         r,
#         as.numeric(colnames(tmp_expr_meta_mt))
#       )
#       
#       ttl_metacell_lst[[r]] = tmp_expr_meta_mt
#       r = r+1
#       redefine_celltypes = c(redefine_celltypes,
#                              setNames(nm = colnames(tmp_expr_meta_mt), rep(i_celltype, ncol(tmp_expr_meta_mt))))
#     }
#     
#     if (double_check) {
#       # Check meta distances
#       obs_DistinMeta = parallel::mclapply(split(t(lognorm_mt), leiden_obj), function(x) {
#         y_mt = matrix(x, ncol = nrow(lognorm_mt))
#         return(median(dist(y_mt), na.rm = T))
#       }, mc.cores = n_cores)
#       obs_DistinMeta_vec = unlist(obs_DistinMeta)
#       obs_DistinMeta_vec = obs_DistinMeta_vec[!is.na(obs_DistinMeta_vec)]
#       
#       message('Running permuation now, this step is computational expensive.')
#       # cl = parallel::makeCluster(n_cores)
#       perm_DistinMeta_lst = pbmcapply::pbmclapply(1:n_perm, mc.cores = n_cores, function(i_perm) {
#         sample_leiden_obj = sample(leiden_obj)
#         perm_DistinMeta = lapply(split(t(lognorm_mt), sample_leiden_obj), function(x) {
#           y_mt = matrix(x, ncol = nrow(lognorm_mt))
#           return(median(dist(y_mt), na.rm = T))
#         })
#         perm_DistinMeta_vec = unlist(perm_DistinMeta)
#         perm_DistinMeta_vec = perm_DistinMeta_vec[!is.na(perm_DistinMeta_vec)]
#         return(perm_DistinMeta_vec)
#       })
#       #   stopCluster(cl = cl)
#       plt_DistinMeta_tbl = obs_DistinMeta_vec %>% tibble::enframe
#       colnames(plt_DistinMeta_tbl) = c('meta_id', 'median_dist')
#       plt_DistinMeta_tbl = plt_DistinMeta_tbl %>% mutate(source = 'obs')
#       
#       perm_DistinMeta_tbl = lapply(1:length(perm_DistinMeta_lst), function(i_perm) {
#         x = perm_DistinMeta_lst[[i_perm]]
#         tmp_DistinMeta_tbl = enframe(x)
#         colnames(tmp_DistinMeta_tbl) = c('meta_id', 'median_dist')
#         tmp_DistinMeta_tbl = tmp_DistinMeta_tbl %>% mutate(source = sprintf('perm_%d', i_perm))
#         return(tmp_DistinMeta_tbl)
#       }) %>% bind_rows()
#       plt_DistinMeta_tbl = bind_rows(plt_DistinMeta_tbl, perm_DistinMeta_tbl) %>%
#         mutate(category = ifelse(source == 'obs', 'obs', 'perm'))
#       
#       p_distdistr = plt_DistinMeta_tbl %>% ggplot(aes(x = source, y = median_dist, color = category)) +
#         geom_violin() +
#         geom_boxplot(width = 0.5) +
#         labs(x = 'Source', y = 'Median distribution', title = 'Distance distribution between cells within each metacell') +
#         scale_color_manual(values = c('obs' = 'red', 'perm' = 'grey')) +
#         theme(axis.text.x = element_text(
#           angle = 90,
#           hjust = 1,
#           vjust = 0.5
#         ), legend.position = 'none')
#       # p_distdistr
#       # ggsave(
#       #   filename = paste0(getwd(), '/distance_within_metacell_', i_celltype,'_doublecheck.pdf'),
#       #   plot = p_distdistr,
#       #   width = 10,
#       #   height = 5,
#       #   device = 'pdf'
#       # )
#     }
#     
#     tictoc::toc()
#   }
#   
#   message('Merging matrices...')
#   expr_meta_mt = do.call(what = cbind, args = ttl_metacell_lst)
#   
#   return(list(mt = expr_meta_mt, celltypes = redefine_celltypes))
# }



# #' Merge single-cell and bulk data.ver1
# #' 
# #' Merge single cells and bulks/spots into one expression matrix.
# #' 
# #' @param sc_data Sparse matrix of normalized counts or TPM, genes x cells
# #' @param bulk_data Sparse matrix of TPM, genes x samples
# #' @param celltypes Named vector, should be a cell-id-named vector of cell types
# #' @param top_gene_ratio Numeric, top % of HVG genes to remain
# #' @param logarithmetic Logical, whether to log-transformed the data
# #' @param scaling Logical, whether to scale data by Frobenius norm
# #' @param quantile_normalization Logical, whether to perform quantile normalization
# #' @param downsample_size 
# #' @param draw_now Logical, default = T
# #' @param verbose Logical, default TRUE
# #' @return A named list
# #' @export
# #' 
# preprocessing = function(sc_data,
#                          bulk_data,
#                          celltypes, 
#                          top_gene_ratio = 0.5,
#                          logarithmetic = T, 
#                          scaling = T,
#                          quantile_normalization = F,
#                          downsample_size = NULL,
#                          draw_now = T,
#                          verbose = T, ...){
#   
#   stopifnot(all(colnames(sc_data) == names(celltypes)))
#   
#   if(!is.null(downsample_size)){
#     if (verbose) message('Downsampling')
#     celltypes2run = names((table(celltypes) > downsample_size) %>% .[.])
#     chosen_cells = names(celltypes)[!celltypes %in% celltypes2run]
#     for (i_celltype in celltypes2run) {
#       chosen_cells = c(chosen_cells, sample(names(celltypes)[celltypes == i_celltype], downsample_size))
#     }
#     celltypes = celltypes[chosen_cells]
#     sc_data = sc_data[, chosen_cells]
#   }
#   
#   
#   # rm 0 genes 
#   if(verbose) message('Now preprocessing')
#   genes2remain = Matrix::rowSums(sc_data > 1) > 0
#   sc_flt_data = sc_data[genes2remain,]
#   
#   
#   genes2remain = Matrix::rowSums(bulk_data > 1) > 0
#   bulk_norm_flt_mt = bulk_data[genes2remain,]
#   
#   # Extract intersect genes
#   inter_genes = intersect(rownames(sc_flt_data), rownames(bulk_norm_flt_mt))
#   
#   # sc_flt_data = calc_normExpr.TPM(sc_flt_data, logarithmetic = F)
#   # Takes about 1.5 min
#   # sc_flt_data = calc_normExpr.Droplet(sc_flt_data, logarithmetic = F)
#   
#   
#   # Identify HVG 
#   # if(verbose) message('Pass: Identifying HVGs for single cells')
#   # sc_hvg_obj = select_features(sc_flt_data, high_quantile = 1-top_gene_ratio)
#   # plt_obj = sc_hvg_obj
#   # 
#   # p_sc = plt_obj$hvg_df %>%
#   #   ggplot(aes(x = log10(means), y = log10(vars))) +
#   #   geom_point(shape = 'o', aes(color = selected)) +
#   #   geom_line(
#   #     data = data.frame(x = plt_obj$fit$x %>% as.numeric(),
#   #                       y = plt_obj$fit$fitted),
#   #     aes(x = x, y = y),
#   #     color = 'black'
#   #   ) +
#   #   theme_classic(base_size = 18) +
#   #   scale_color_manual(values = c("TRUE" = JasonToolBox::gg_color_hue(2)[1], "FALSE" = "grey50")) +
#   #   labs(title = sprintf(
#   #     "Single cells: CP10K, %d genes selected",
#   #     sum(plt_obj$hvg_df$selected)
#   #   ))
#   
#   if(verbose) message('Identifying HVGs for bulk data')
#   bulk_hvg_obj = select_features(bulk_norm_flt_mt, high_quantile = 1-top_gene_ratio)
#   plt_obj = bulk_hvg_obj
# 
#   p_bulk = plt_obj$hvg_df %>%
#     ggplot(aes(x = log10(means), y = log10(vars))) +
#     geom_point(shape = 'o', aes(color = selected)) +
#     geom_line(
#       data = data.frame(x = plt_obj$fit$x %>% as.numeric(), y = plt_obj$fit$fitted),
#       aes(x = x, y = y),
#       color = 'black'
#     ) +
#     theme_classic(base_size = 18) +
#     scale_color_manual(values = c("TRUE" = JasonToolBox::gg_color_hue(2)[1], "FALSE" = "grey50")) +
#     labs(title = sprintf("TCGA: TPM, %d genes selected", sum(plt_obj$hvg_df$selected)))
# 
#   if(draw_now){
#     # print(p_sc)
#     print(p_bulk)
#   }
#   
#   
#   # Filter for HVG for matrices 
#   if(verbose) message('Extracting HVGs')
#   inter_genes = intersect(rownames(bulk_norm_flt_mt), rownames(sc_flt_data))
#   all_variable_genes = union(bulk_hvg_obj$hvg_df$featureID[bulk_hvg_obj$hvg_df$selected],
#                              rownames(sc_flt_data)) %>% unique
#   # all_variable_genes = union(rownames(bulk_norm_flt_mt),
#                              # rownames(sc_flt_data)) %>% unique
#   
#   chosen_genes = unique(intersect(inter_genes, all_variable_genes))
#   if(verbose) message(sprintf("%d genes chosen.", length(chosen_genes)))
#   
#   if(verbose) message(sprintf(
#     "%d / %d TCGA bulk variable genes chosen",
#     sum(rownames(bulk_norm_flt_mt) %in% chosen_genes),
#     length(rownames(bulk_norm_flt_mt))
#   ))
#   
#   if(verbose) message(sprintf(
#     "%d / %d test single cells variable genes chosen",
#     sum(rownames(sc_flt_data) %in% chosen_genes),
#     length(rownames(sc_flt_data))
#   ))
#   
#   # normalize, select genes, scale by MSE
#   sc_norm_data = sc_flt_data[chosen_genes,]
#   
#   # 
#   # if(merge_cell){
#   #   if(verbose) message('Now downsampling')
#   #   sc_norm_data_obj = mergeCells(sc_norm_data, celltypes, target_number = 30, K = ncol(bulk_data), double_check = F)
#   #   sc_norm_data = sc_norm_data_obj$mt
#   #   celltypes = sc_norm_data_obj$celltypes
#   # } 
#   
#   sc_norm_data = as(sc_norm_data, 'dgCMatrix')
#   bulk_tpm_flt_MM = as(bulk_norm_flt_mt[chosen_genes,], 'dgCMatrix')
# 
#   # Normalize by library size
#   # scale_factor = colSums(sc_norm_data)
#   # sc_norm_data = sweep_MM(sc_norm_data, 2, scale_factor, function(x, y) {
#   #   x / y
#   # })
#   # 
#   # scale_factor = colSums(bulk_tpm_flt_MM)
#   # bulk_tpm_flt_MM = sweep_MM(bulk_tpm_flt_MM, 2, scale_factor, function(x, y) {
#   #   x / y
#   # })
#   if(logarithmetic){
#     if(verbose) message('Now log-transforming')
#     sc_norm_data@x = log1p(sc_norm_data@x)
#     bulk_tpm_flt_MM@x = log1p(bulk_tpm_flt_MM@x)
#   }
#   
#   # scale by source: Frobenius norm
#   if(scaling){
#     if(verbose) message('Now scaling')
#     scale_factor = sqrt(sum(sc_norm_data@x ^ 2)) # Frobenius norm for sc
#     sc_norm_data@x = sc_norm_data@x / scale_factor
#     
#     scale_factor = sqrt(sum(bulk_tpm_flt_MM@x ^ 2)) # Frobenius norm for bulk
#     bulk_tpm_flt_MM@x = bulk_tpm_flt_MM@x / scale_factor
#   }
#   
#   # bulk_tpm_flt_MM %>% apply(1, function(x){x/sqrt(sum(x^2)/length(x))}) %>% t() %>% as('dgTMatrix')
#   
#   merge_MM = Matrix::cbind2(sc_norm_data,
#                             bulk_tpm_flt_MM)
#   
#   meta.data = data.frame(
#     name = colnames(merge_MM),
#     type = ifelse(
#       colnames(merge_MM) %in% colnames(bulk_norm_flt_mt),
#       "bulk",
#       "singlecell"
#     )
#   ) %>% mutate(
#     celltype = ifelse(is.na(celltypes[name]), "bulk", celltypes[name])# ,
#     # major_celltype = ifelse(is.na(maj_celltype_vec[name]), "bulk", maj_celltype_vec[name])
#   )
#   
#   if(verbose) message(sprintf(
#     "Merged expression matrix is %d genes * %d samples.",
#     nrow(merge_MM),
#     ncol(merge_MM)
#   ))
#   
#   if (quantile_normalization == T) {
#     if(verbose) message('Now quantile normalizing')
#     merge_MM = quantile_normalize(merge_MM)
#   }
#   return(
#     list(
#       merge_MM = merge_MM,
#       meta.data = meta.data,
#       # plot_sc = p_sc,
#       plot_bulk = p_bulk,
#       logarithmetic = logarithmetic,
#       scale = scaling,
#       quantile_normalization = quantile_normalization
#     )
#   )
# }

#' Merge single-cell and bulk data.ver2
#' 
#' Merge single cells and bulks/spots into one expression matrix.
#' 
#' @param sc_data Sparse matrix of normalized counts or TPM, genes x cells
#' @param bulk_data Sparse matrix of TPM, genes x samples
#' @param celltypes Named vector, should be a cell-id-named vector of cell types
#' @param top_gene_ratio Numeric, top % of HVG genes to remain
#' @param logarithmetic Logical, whether to log-transformed the data
#' @param scaling Logical, whether to scale data by Frobenius norm
#' @param quantile_normalization Logical, whether to perform quantile normalization
#' @param downsample_size Number of cells to downsample for each cell type. Default = NULL, no sampling.
#' @param draw_now Logical, default = T
#' @param verbose Logical, default TRUE
#' @return A named list
#' @export
#' @importFrom magrittr `%>%`
#' 
preprocessing2 = function(sc_data,
                          bulk_data,
                          celltypes,
                          top_gene_ratio = 0.5,
                          logarithmetic = T,
                          scaling = F,
                          quantile_normalization = T,
                          downsample_size = NULL,
                          draw_now = T,
                          verbose = T) {
  
  
  
  stopifnot(all(colnames(sc_data) == names(celltypes)))
  
  if(!is.null(downsample_size)){
    if (verbose) message('Downsampling')
    celltypes2run = names((table(celltypes) > downsample_size) %>% .[.])
    chosen_cells = names(celltypes)[!celltypes %in% celltypes2run]
    for (i_celltype in celltypes2run) {
      chosen_cells = c(chosen_cells, sample(names(celltypes)[celltypes == i_celltype], downsample_size))
    }
    celltypes = celltypes[chosen_cells]
    sc_data = sc_data[, chosen_cells]
  }
  
  
  # rm 0 genes 
  if(verbose) message('Now preprocessing')
  genes2remain = Matrix::rowSums(sc_data > 1) > 0
  sc_flt_data = sc_data[genes2remain,]
  
  
  genes2remain = Matrix::rowSums(bulk_data > 1) > 0
  bulk_norm_flt_mt = bulk_data[genes2remain,]
  
  # Extract intersect genes
  inter_genes = intersect(rownames(sc_flt_data), rownames(bulk_norm_flt_mt))
  
  # sc_flt_data = calc_normExpr.TPM(sc_flt_data, logarithmetic = F)
  # Takes about 1.5 min
  # sc_flt_data = calc_normExpr.Droplet(sc_flt_data, logarithmetic = F)
  
  
  # Identify HVG 
  # if(verbose) message('Pass: Identifying HVGs for single cells')
  # sc_hvg_obj = select_features(sc_flt_data, high_quantile = 1-top_gene_ratio)
  # plt_obj = sc_hvg_obj
  # 
  # p_sc = plt_obj$hvg_df %>%
  #   ggplot(aes(x = log10(means), y = log10(vars))) +
  #   geom_point(shape = 'o', aes(color = selected)) +
  #   geom_line(
  #     data = data.frame(x = plt_obj$fit$x %>% as.numeric(),
  #                       y = plt_obj$fit$fitted),
  #     aes(x = x, y = y),
  #     color = 'black'
  #   ) +
  #   theme_classic(base_size = 18) +
  #   scale_color_manual(values = c("TRUE" = JasonToolBox::gg_color_hue(2)[1], "FALSE" = "grey50")) +
  #   labs(title = sprintf(
  #     "Single cells: CP10K, %d genes selected",
  #     sum(plt_obj$hvg_df$selected)
  #   ))
  
  if(verbose) message('Identifying HVGs for bulk data')
  bulk_hvg_obj = select_features(bulk_norm_flt_mt, high_quantile = 1 - top_gene_ratio)
  plt_obj = bulk_hvg_obj
  
  p_bulk = plt_obj$hvg_df %>%
    ggplot(aes(x = log10(means), y = log10(vars))) +
    geom_point(shape = 'o', aes(color = selected)) +
    geom_line(
      data = data.frame(x = plt_obj$fit$x %>% as.numeric(), y = plt_obj$fit$fitted),
      aes(x = x, y = y),
      color = 'black'
    ) +
    theme_classic(base_size = 18) +
    scale_color_manual(values = c("TRUE" = JasonToolBox::gg_color_hue(2)[1], "FALSE" = "grey50")) +
    labs(title = sprintf("Bulk: TPM, %d genes selected", sum(plt_obj$hvg_df$selected)))
  
  if(draw_now){
    # print(p_sc)
    print(p_bulk)
  }
  
  
  # Filter for HVG for matrices 
  if(verbose) message('Extracting HVGs')
  inter_genes = intersect(rownames(bulk_norm_flt_mt), rownames(sc_flt_data))
  all_variable_genes = union(bulk_hvg_obj$hvg_df$featureID[bulk_hvg_obj$hvg_df$selected],
                             rownames(sc_flt_data)) %>% unique
  # all_variable_genes = union(rownames(bulk_norm_flt_mt),
  # rownames(sc_flt_data)) %>% unique
  
  chosen_genes = unique(intersect(inter_genes, all_variable_genes))
  if(verbose) message(sprintf("%d genes chosen.", length(chosen_genes)))
  
  if(verbose) message(sprintf(
    "%d / %d bulk genes chosen",
    sum(rownames(bulk_norm_flt_mt) %in% chosen_genes),
    length(rownames(bulk_norm_flt_mt))
  ))
  
  if(verbose) message(sprintf(
    "%d / %d test single cells variable genes chosen",
    sum(rownames(sc_flt_data) %in% chosen_genes),
    length(rownames(sc_flt_data))
  ))
  
  # normalize, select genes, scale by MSE
  sc_norm_data = sc_flt_data[chosen_genes,]
  
  # 
  # if(merge_cell){
  #   if(verbose) message('Now downsampling')
  #   sc_norm_data_obj = mergeCells(sc_norm_data, celltypes, target_number = 30, K = ncol(bulk_data), double_check = F)
  #   sc_norm_data = sc_norm_data_obj$mt
  #   celltypes = sc_norm_data_obj$celltypes
  # } 
  
  sc_norm_data = as(sc_norm_data, 'dgCMatrix')
  bulk_tpm_flt_MM = as(bulk_norm_flt_mt[chosen_genes,], 'dgCMatrix')
  
  # Normalize by library size
  # scale_factor = colSums(sc_norm_data)
  # sc_norm_data = sweep_MM(sc_norm_data, 2, scale_factor, function(x, y) {
  #   x / y
  # })
  # 
  # scale_factor = colSums(bulk_tpm_flt_MM)
  # bulk_tpm_flt_MM = sweep_MM(bulk_tpm_flt_MM, 2, scale_factor, function(x, y) {
  #   x / y
  # })
  if(logarithmetic){
    if(verbose) message('Now log-transforming')
    sc_norm_data@x = log1p(sc_norm_data@x)
    bulk_tpm_flt_MM@x = log1p(bulk_tpm_flt_MM@x)
  }
  origin_merge_MM = Matrix::cbind2(sc_norm_data, bulk_tpm_flt_MM)
  
  # scale by source: Frobenius norm
  if(scaling){
    if(verbose) message('Now scaling')
    scale_factor = sqrt(sum(sc_norm_data@x ^ 2)) # Frobenius norm for sc
    sc_norm_data@x = sc_norm_data@x / scale_factor
    
    scale_factor = sqrt(sum(bulk_tpm_flt_MM@x ^ 2)) # Frobenius norm for bulk
    bulk_tpm_flt_MM@x = bulk_tpm_flt_MM@x / scale_factor
  }
  
  # bulk_tpm_flt_MM %>% apply(1, function(x){x/sqrt(sum(x^2)/length(x))}) %>% t() %>% as('dgTMatrix')
  
  merge_MM = Matrix::cbind2(sc_norm_data,
                            bulk_tpm_flt_MM)
  
  meta.data = data.frame(
    name = colnames(merge_MM),
    type = ifelse(
      colnames(merge_MM) %in% colnames(bulk_norm_flt_mt),
      "bulk",
      "singlecell"
    )
  ) %>% mutate(
    celltype = ifelse(is.na(celltypes[name]), "bulk", celltypes[name])# ,
    # major_celltype = ifelse(is.na(maj_celltype_vec[name]), "bulk", maj_celltype_vec[name])
  )
  
  if(verbose) message(sprintf(
    "Merged expression matrix is %d genes * %d samples.",
    nrow(merge_MM),
    ncol(merge_MM)
  ))
  
  if (quantile_normalization == T) {
    if(verbose) message('Now quantile normalizing')
    merge_MM = quantile_normalize(merge_MM)
  }
  return(
    list(
      origin_MM = origin_merge_MM, 
      merge_MM = merge_MM,
      meta.data = meta.data,
      # plot_sc = p_sc,
      plot_bulk = p_bulk,
      logarithmetic = logarithmetic,
      scale = scaling,
      quantile_normalization = quantile_normalization
    )
  )
}


#' Merge single-cell and bulk data.ver3
#' 
#' Merge single cells and bulks/spots into one expression matrix.
#' 
#' @param sc_data Sparse matrix of normalized counts or TPM, genes x cells
#' @param bulk_data Sparse matrix of TPM, genes x samples
#' @param celltypes Named vector, should be a cell-id-named vector of cell types
#' @param top_gene_ratio Numeric, top % of HVG genes to remain
#' @param logarithmetic Logical, whether to log-transformed the data
#' @param findHVG Logical, whether to find HVG for bulks
#' @param scaling Logical, whether to scale data by Frobenius norm
#' @param quantile_normalization Logical, whether to perform quantile normalization
#' @param downsample_size Number of cells to downsample for each cell type. Default = NULL, no sampling.
#' @param draw_now Logical, default = T
#' @param verbose Logical, default TRUE
#' @return A named list
#' @export
#' @importFrom magrittr `%>%`
#' 
preprocessing3 = function(sc_data,
                          bulk_data,
                          celltypes,
                          top_gene_ratio = 0.5,
                          logarithmetic = T,
                          findHVG = F, 
                          scaling = T,
                          quantile_normalization = F,
                          downsample_size = NULL,
                          draw_now = T,
                          verbose = T) {
  
  
  
  stopifnot(!is.null(colnames(sc_data)))
  stopifnot(!is.null(colnames(bulk_data)))
  stopifnot(!(length(sc_data) == 0))
  stopifnot(!(length(bulk_data) == 0))
  stopifnot(all(colnames(sc_data) == names(celltypes)))
  
  if(!is.null(downsample_size)){
    if (verbose) message('Downsampling')
    celltypes2run = names((table(celltypes) > downsample_size) %>% .[.])
    chosen_cells = names(celltypes)[!celltypes %in% celltypes2run]
    for (i_celltype in celltypes2run) {
      chosen_cells = c(chosen_cells, sample(names(celltypes)[celltypes == i_celltype], downsample_size))
    }
    celltypes = celltypes[chosen_cells]
    sc_data = sc_data[, chosen_cells]
  }
  
  
  # rm 0 genes 
  if(verbose) message('Now preprocessing')
  genes2remain = Matrix::rowSums(sc_data > 1) > 0
  sc_flt_data = sc_data[genes2remain,]
  
  
  genes2remain = Matrix::rowSums(bulk_data > 1) > 0
  bulk_norm_flt_mt = bulk_data[genes2remain,]
  
  # Extract intersect genes
  inter_genes = intersect(rownames(sc_flt_data), rownames(bulk_norm_flt_mt))
  
  # sc_flt_data = calc_normExpr.TPM(sc_flt_data, logarithmetic = F)
  # Takes about 1.5 min
  # sc_flt_data = calc_normExpr.Droplet(sc_flt_data, logarithmetic = F)
  
  
  # Identify HVG 
  if(findHVG) {
    if (verbose)
      message('Identifying HVGs for bulk data')
    bulk_hvg_obj = select_features(bulk_norm_flt_mt, high_quantile = 1 - top_gene_ratio)
    plt_obj = bulk_hvg_obj
    
    p_bulk = plt_obj$hvg_df %>%
      ggplot(aes(x = log10(means), y = log10(vars))) +
      geom_point(shape = 'o', aes(color = selected)) +
      geom_line(
        data = data.frame(
          x = plt_obj$fit$x %>% as.numeric(),
          y = plt_obj$fit$fitted
        ),
        aes(x = x, y = y),
        color = 'black'
      ) +
      theme_classic(base_size = 18) +
      scale_color_manual(values = c("TRUE" = JasonToolBox::gg_color_hue(2)[1], "FALSE" = "grey50")) +
      labs(title = sprintf("Bulk: TPM, %d genes highlighted", sum(plt_obj$hvg_df$selected)))
    
    if (draw_now) {
      # print(p_sc)
      print(p_bulk)
    }
    if(verbose) message('Extracting HVGs')
    
    # Filter for HVG for matrices
    all_variable_genes = union(bulk_hvg_obj$hvg_df$featureID[bulk_hvg_obj$hvg_df$selected],
                               rownames(sc_flt_data)) %>% unique
    
    chosen_genes = unique(intersect(inter_genes, all_variable_genes))
    if (verbose)
      message(sprintf("%d genes chosen.", length(chosen_genes)))
    
    if (verbose)
      message(sprintf(
        "%d / %d bulk genes chosen",
        sum(rownames(bulk_norm_flt_mt) %in% chosen_genes),
        length(rownames(bulk_norm_flt_mt))
      ))
    
    if (verbose)
      message(sprintf(
        "%d / %d single cells genes chosen",
        sum(rownames(sc_flt_data) %in% chosen_genes),
        length(rownames(sc_flt_data))
      ))
    
    # normalize, select genes, scale by MSE
    # sc_norm_data = sc_flt_data[chosen_genes, ]
    # bulk_norm_flt_mt = bulk_norm_flt_mt[chosen_genes,]
  } else {
    chosen_genes = unique(inter_genes)
  }
  
  sc_norm_data = sc_flt_data[chosen_genes, ]
  bulk_norm_flt_mt = bulk_norm_flt_mt[chosen_genes,]
  
  sc_norm_data = as(sc_norm_data, 'dgCMatrix')
  bulk_tpm_flt_MM = as(bulk_norm_flt_mt, 'dgCMatrix')
  
  # Normalize by library size
  # scale_factor = colSums(sc_norm_data)
  # sc_norm_data = sweep_MM(sc_norm_data, 2, scale_factor, function(x, y) {
  #   x / y
  # })
  # 
  # scale_factor = colSums(bulk_tpm_flt_MM)
  # bulk_tpm_flt_MM = sweep_MM(bulk_tpm_flt_MM, 2, scale_factor, function(x, y) {
  #   x / y
  # })
  if(logarithmetic){
    if(verbose) message('Now log-transforming')
    sc_norm_data@x = log1p(sc_norm_data@x)
    bulk_tpm_flt_MM@x = log1p(bulk_tpm_flt_MM@x)
  }
  origin_merge_MM = Matrix::cbind2(sc_norm_data, bulk_tpm_flt_MM)
  
  # scale by source: Frobenius norm
  if(scaling){
    if(verbose) message('Now scaling')
    scale_factor = sqrt(sum(sc_norm_data@x ^ 2)) # Frobenius norm for sc
    sc_norm_data@x = sc_norm_data@x / scale_factor
    
    scale_factor = sqrt(sum(bulk_tpm_flt_MM@x ^ 2)) # Frobenius norm for bulk
    bulk_tpm_flt_MM@x = bulk_tpm_flt_MM@x / scale_factor
  }
  
  # bulk_tpm_flt_MM %>% apply(1, function(x){x/sqrt(sum(x^2)/length(x))}) %>% t() %>% as('dgTMatrix')
  
  merge_MM = Matrix::cbind2(sc_norm_data,
                            bulk_tpm_flt_MM)
  
  meta.data = data.frame(
    name = colnames(merge_MM),
    type = ifelse(
      colnames(merge_MM) %in% colnames(bulk_norm_flt_mt),
      "bulk",
      "singlecell"
    )
  ) %>% mutate(
    celltype = ifelse(is.na(celltypes[name]), "bulk", celltypes[name])# ,
    # major_celltype = ifelse(is.na(maj_celltype_vec[name]), "bulk", maj_celltype_vec[name])
  )
  
  if(verbose) message(sprintf(
    "Merged expression matrix is %d genes * %d samples.",
    nrow(merge_MM),
    ncol(merge_MM)
  ))
  
  if (quantile_normalization == T) {
    if(verbose) message('Now quantile normalizing')
    merge_MM = quantile_normalize(merge_MM)
  }
  if(findHVG) {
    return(
      list(
        origin_MM = origin_merge_MM,
        merge_MM = merge_MM,
        meta.data = meta.data,
        # plot_sc = p_sc,
        plot_bulk = p_bulk,
        logarithmetic = logarithmetic,
        scale = scaling,
        quantile_normalization = quantile_normalization
      )
    )
  } else {
    return(
      list(
        origin_MM = origin_merge_MM,
        merge_MM = merge_MM,
        meta.data = meta.data,
        logarithmetic = logarithmetic,
        scale = scaling,
        quantile_normalization = quantile_normalization
      )
    )
    
  }
}


# Determine K ----

identifyLinetype = function(x, y) {
  line_df = dplyr::arrange(data.frame(x, y), x)
  start_point = head(line_df$y, 1)
  end_point = tail(line_df$y, 1)
  max_point = max(line_df$y)
  min_point = min(line_df$y)
  if (start_point == max_point && end_point == min_point) {
    return(list(
      type = 'decrease',
      crop_point = NULL
    ))
  } else if (start_point == min_point && end_point == max_point) {
    return(list(
      type = 'increase',
      crop_point = NULL
    ))
  } else if (max_point > start_point && max_point > end_point) {
    mid_point = line_df$x[line_df$y == max_point]
    return(list(
      type = 'cap',
      crop_point = mid_point
    ))
  } else if (min_point < start_point && min_point < end_point) {
    mid_point = line_df$x[line_df$y == min_point]
    return(list(
      type = 'cup',
      crop_point = mid_point
    ))
    
  } else {
    warning('This should not appear.')
    return(NULL)
  }
}

findKnee = function(x, y , ...) {
  kneer_obj = kneer::create_knee_locator(x, y, ...)
  kneer_obj = kneer::preprocess(kneer_obj)
  kneer_obj = kneer::find_knee(kneer_obj)
  return(kneer_obj)
}
 

#' Identify knee
#' 
#' @param x Numeric vector, x coordinates
#' @param y Numeric vector, y coordinates.
#' @param direction String, either 'up' or 'down'
#' @param draw_now Logical, if print the plot now
#' @param ... Other parameters for kneer::create_knee_locator
#' @return A named list
#' @export
#' @importFrom magrittr `%>%`
#'
identifyKnee = function(x, y, direction = 'down', draw_now = T, ...){
  
  # require(dplyr)
  # require(ggplot2)
  # require(kneer)
  # require(patchwork)
  
  n_cat = length(unique(x))
  plt_df = data.frame(x, y) %>% arrange(x)
  # Crop data according to direction
  loess_obj = loess(formula = y ~ x, plt_df)
  
  line_df = data.frame(x = unique(loess_obj$x[, 1]), y = unique(loess_obj$fitted))
  # fitted_y = setNames(line_df)
  
  linetype_obj = identifyLinetype(line_df$x, line_df$y)

    
  if (direction == 'up') {
    switch(
      linetype_obj$type,
      increase = {
      },
      decrease = {
        warning('The requested trend is inconsistent. Returning NULL')
        kneer_obj = NULL
      },
      cup = {
        crop_point = linetype_obj$crop_point
        plt_df = plt_df[plt_df$x > crop_point, ]
      },
      cap = {
      }
    )
  } else {
    switch(
      linetype_obj$type,
      increase = {
        warning('The requested trend is inconsistent. Returning NULL')
        kneer_obj = NULL
      },
      decrease = {
      },
      cup = {
      },
      cap = {
        crop_point = linetype_obj$crop_point
        plt_df = plt_df[plt_df$x > crop_point, ]
      }
    )
  }
  
  if(length(unique(plt_df$x)) < 0.1 * n_cat){
    warning('Fewer data points remain. The proper knee point cannot be found. Returning NULL')
    kneer_obj = NULL
  }
  
  tryCatch({kneer_obj = findKnee(plt_df$x, plt_df$y, ...)}, error = function(e){
    warning('Cannot find knee.')
    kneer_obj <<- NULL
  })
  
  if(!is.null(kneer_obj)){
    p_out = plotKnee(data.frame(x,y), fit_obj = kneer_obj)
  } else {
    p_out = plotLoess(data.frame(x,y), fit_obj = loess_obj)
  }
  if(draw_now){
    print(p_out)
  }
  if(is.null(kneer_obj)){
    return(list(kneer_point = NA, plot = p_out))
  } else {
    return(list(kneer_point = kneer_obj@knee, plot = p_out))
  }
}

#' Generate loess plot when knee cannot be found
#' 
#' @param plt_df A data frame with x and y to plot.
#' @param fit_obj A loess fit object generated by `loess`
#' @return A ggplot plot
#' 
plotLoess = function(plt_df, fit_obj){
  # require(ggplot2)
  p1 = ggplot(data = plt_df, aes(x = x, y = y)) +
    geom_point() +
    geom_line(data = data.frame(x = unique(fit_obj$x[, 1]), y = unique(fit_obj$fitted)), 
              linetype = 2)
  return(p1)
}

#' Generate knee plot
#' 
#' @param plt_df A data frame with x and y to plot.
#' @param fit_obj A knee object generated by `kneer`
#' @return A ggplot plot
#' 
plotKnee = function(plt_df, fit_obj) {
  # require(ggplot2)
  # require(patchwork)
  p1 = ggplot(data = plt_df, aes(x = x, y = y)) +
    geom_point() +
    geom_line(data = data.frame(x = fit_obj@x, fitted_y = fit_obj@y.smooth),
              aes(x = x, y = fitted_y)) +
    geom_vline(xintercept = fit_obj@knee, linetype = 2) +
    labs(x = 'x', y = 'y') +
    theme(axis.title.x = element_blank())
  # p1
  
  p2 = ggplot() +
    geom_line(data = data.frame(x = fit_obj@x, fitted_y = fit_obj@y.difference),
              aes(x = x, y = fitted_y)) +
    geom_vline(xintercept = fit_obj@knee, linetype = 2) +
    xlim(min(fit_obj@x), max(fit_obj@x)) +
    labs(x = 'x', y = 'y')
  
  p_merge = p1 + p2 + patchwork::plot_layout(ncol = 1) + patchwork::plot_annotation(title = sprintf(
    "Identify the heuristically optimal K:\nKnee point = %d",
    fit_obj@knee
  ))
  p_merge
  return(p_merge)
}


#' Run NMF for determination of K
#' 
#' A wrapper function to run NMF multiple times and calculate the index matrix helping determination of K
#' @param merge_mt Matrix or sparse matrix. Normalized expression matrix of single-cell and bulk. Genes x Samples.
#' @param Ks2run Integer vector. Ks to test.
#' @param n_reps Integer. Number of replication for each K.
#' @param n_cores Integer. Number of cores to use, recommended to use multicore. Default = 0, all cores possible.
#' @param draw_now Logical. If print out the plots.
#' @param verbose Logical. If print the logs.
#' @param return_res Logical. If return all results.
#' @return A named list with data to determine K
#' @export
#' 
determineK_runNMF = function(merge_mt,
                             Ks2run,
                             n_reps,
                             n_cores = 0,
                             draw_now = F,
                             verbose = T,
                             return_res = F) {
  
  # require(tictoc)
  # require(RcppML)
  # require(ggplot2)
  # require(dplyr)
  # require(checkmate)
  # args = list(...)
  # return_res = ifelse(is.null(args$return_res), F, T)
  
  # checkmate::assert_logical(verbose)
  # checkmate::assert_int(n_cores)
  
  options(RcppML.verbose = F)
  if (verbose) {
    message("Verbose mode.")
  } else {
    # message("Silence mode.")
  }
  
  if(n_cores == 0){
    options(RcppML.threads = n_cores)
    message('Using all cores possible: ', parallel::detectCores())
    n_cores = parallel::detectCores()
  } else {
    options(RcppML.threads = n_cores)
    message('Using ', n_cores, ' cores.')
  }
  
  # merge_mt = merge_mt$merge_MM
  # if (verbose) message("Now running k = ", k)
  if(return_res) res_lst = list()
  
  invisible(gc())
  res_monitor_df = data.frame()
  entropy_df = data.frame()
  jac_index_df = data.frame()
  
  for (k in Ks2run) {
    if (verbose) {
      message("Now running k = ", k)
    }
    ##### Run NMF #####
    if (verbose) {
      tictoc::tic(msg = "Fast NMF sequentially")
      JasonToolBox::initiatePB('.I1')
    }
    
    fnmf_model_lst = list()
    for (i_rep in 1:n_reps) {
      tictoc::tic()
      i_seed = 1000 * k + i_rep
      
      fnmf_model_lst[[i_rep]] = list(
        nmf_obj = RcppML::nmf(
          data = merge_mt,
          k = k,
          tol = 1e-5,
          seed = i_seed
        ))
      irep_time = toc(quiet = T)
      res_monitor_df = rbind(
        res_monitor_df,
        data.frame(
          k = k,
          seed = i_seed,
          time_sequentially = irep_time$toc - irep_time$tic
        )
      )
      if (verbose)
        JasonToolBox::processBar('.I1', i_rep, n_reps)
    }
    if (verbose)
      f_time_seq = tictoc::toc()
    if(return_res) res_lst[[as.character(k)]] = fnmf_model_lst
    invisible(gc())
    
    ########### Calculate Entropy #################
    # no need to store res_lst
    if (verbose)
      message("Calculating entropy...")
    tmp_df = bind_rows(parallel::mclapply(
      mc.cores = n_cores,
      X = 1:length(fnmf_model_lst),
      FUN = function(ith) {
        fnmf_obj = fnmf_model_lst[[ith]]
        fnmf_model = fnmf_obj$nmf_obj
        test_x = fnmf_model@h
        test_x[1:length(test_x)] = cut(fnmf_model@h, 50)
        w_entr_vec = apply(test_x, 2, function(x) {
          # y = entropy::entropy(table(x), unit = 'log2') / (-log2(1 / length(table(x)))) # normalized entropy
          y = entropy::entropy(x, unit = 'log2') / (-log2(1 / length(x))) # normalized entropy
          return(y)
        })
        
        # w_entr_vec = apply(fnmf_model@h, 2, function(x) {
        #   y = entropy::entropy(x, unit = 'log2') / (-log2(1 / length(x))) # normalized entropy
        #   return(y)
        # })
        
        return(data.frame(i_rep = ith,
                          median_entropy = median(w_entr_vec)))
      }
    ))
    tmp_df$k = k
    entropy_df = rbind(entropy_df, tmp_df)
    
    ######### Calculate Stability ########
    if (verbose)
      message('Calculating stability...')
    tmp_w_lst = lapply(
      fnmf_model_lst,
      FUN = function(fnmf_obj) {
        return(fnmf_obj$nmf_obj@w)
      }
    )
    
    combn_df = data.frame(t(combn(1:length(tmp_w_lst), 2)))
    combn_df$jac_ind = NA
    n_step = 5
    combn_df$index = (1:nrow(combn_df) %% n_step)
    JasonToolBox::initiatePB('.I1')
    
    for (grp_idx in unique(combn_df$index)) {
      t_combn_df = combn_df[combn_df$index == grp_idx,]
      tmp_jacs_lst = parallel::mclapply(
        mc.cores = n_cores,
        1:nrow(t_combn_df),
        FUN = function(i_row) {
          t_idx1 = t_combn_df[i_row, 1]
          t_idx2 = t_combn_df[i_row, 2]
          x_cor_mt = tmp_w_lst[[t_idx1]]
          y_cor_mt = tmp_w_lst[[t_idx2]]
          return(calc_jacind_from_cor(x_cor_mt, y_cor_mt))
        }
      )
      combn_df[combn_df$index == grp_idx,]$jac_ind = unlist(tmp_jacs_lst)
      
      JasonToolBox::processBar('.I1', which(grp_idx == unique(combn_df$index)), length(unique(combn_df$index)))
    }
    tmp_df = data.frame(k = k, sd_jac_ind = sd(combn_df$jac_ind))
    tmp_df$k = k
    jac_index_df = rbind(jac_index_df, tmp_df)
    
  }
  
  if(return_res){
    return(list(
      data_df = list(entropy_df = entropy_df, jacind_df = jac_index_df, time_df = res_monitor_df),
      fnmf_model_lst = fnmf_model_lst, 
      res_lst = res_lst
    ))
  }
  return(list(
    data_df = list(entropy_df = entropy_df, jacind_df = jac_index_df, time_df = res_monitor_df),
    fnmf_model_lst = fnmf_model_lst
  ))
  
}

#' Determine K to use
#' 
#' Determine the knee point using `kneer` to identify the optimum K to run cNMF.
#' @param data_obj Data object from `determineK_runNMF`
#' @param draw_now Logical. If print out the plots.
#' @param verbose Logical. If print the logs.
#' @param ... Other parameters for kneer::create_knee_locator
#' @return A named list with optimum K and data/objects to support it.
#' @export
#' 
determineK = function(data_obj, draw_now = F, verbose = T, ...) {
  # require(tictoc)
  # require(RcppML)
  # require(ggplot2)
  # require(dplyr)
  
  # return(list(
  #   data_df = list(entropy_df = entropy_df, jacind_df = jac_index_df, time_df = res_monitor_df),
  #   fnmf_model_lst = fnmf_model_lst
  # ))
  
  res_monitor_df = data_obj$data_df$time_df %>% as_tibble()
  res_monitor_df$k = as.factor(res_monitor_df$k)
  p_time = res_monitor_df %>%
    ggplot(aes(x = k, y = time_sequentially)) +
    geom_boxplot() +
    labs(y = 'Time elapsed(sec)')
  
  
  ith = 1
  fnmf_obj = data_obj$fnmf_model_lst[[ith]]
  fnmf_model = fnmf_obj$nmf_obj
  w_entr_vec = apply(fnmf_model@h, 2, function(x) {
    # freqs_x = freqs(x, method = 'ML')
    n_x = sum(x!=0)
    H = -log2(1/n_x)
    return(entropy::entropy(x, method = 'ML', unit = 'log2') / H)
  })
  
  p_entropy_dist = w_entr_vec %>% tibble::enframe() %>% ggplot(aes(x = value)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = median(w_entr_vec)) +
    labs(
      x = "Median of entropy for each sample",
      y = "Count",
      title = paste0(
        "The distribution of the median of entropy",
        "ith replication = ",
        ith
      )
    )
  
  entropy_df = data_obj$data_df$entropy_df %>% dplyr::filter(!is.na(median_entropy))
  knee_point_entropy = identifyKnee(
    x = entropy_df$k,
    y = entropy_df$median_entropy,
    direction = 'up',
    draw_now = draw_now, ...
  )
  
  jac_index_df = data_obj$data_df$jacind_df
  knee_point_jacind = identifyKnee(
    x = jac_index_df$k,
    y = jac_index_df$sd_jac_ind,
    direction = 'down',
    draw_now = draw_now, ...
  )
  
  if (draw_now) {
    print(p_entropy_dist)
    print(p_time)
    # print(knee_point_entropy$plot)
    # print(knee_point_jacind$plot)
  }
  
  return(list(
    K = max(knee_point_entropy$kneer_point, knee_point_jacind$kneer_point, na.rm = T), 
    data = list(entropy_df = entropy_df, jacind_df = jac_index_df),
    knee_objs = list(knee_entropy = knee_point_entropy, knee_jacind = knee_point_jacind),
    plots = list(
      p_entropy_dist, 
      p_time,
      knee_point_entropy$plot,
      knee_point_jacind$plot)
  ))
}

#' runcNMF 
#' 
#' run cNMF to extract programs
#' @param merge_mt Matrix or sparse matrix. Normalized expression matrix of single-cell and bulk. Genes x Samples.
#' @param k Integer. K programs to run.
#' @param n_reps Integer. Number of replication for each K.
#' @param n_cores Integer. Number of cores to use, recommended to use multicore. Default = 0, all cores possible.
#' @param row_subset Double, \(0-1\]. Proportion of rows to sample. Default = 1.
#' @param col_subset Double, \(0-1\]. Proportion of columns to sample. Default = 1.
#' @param verbose Logical. Whether to print the logs.
#' @return A list of all the NMF results.
#' @export
#' 
runcNMF = function(merge_mt,
                   k,
                   n_reps,
                   n_cores = 0,
                   row_subset = 1,
                   col_subset = 1,
                   verbose = T) {
  
  # require(tictoc)
  # require(RcppML)
  # require(ggplot2)
  # require(dplyr)
  
  
  options(RcppML.verbose = F)
  if (verbose) {
    message("Verbose mode.")
  } else {
    message("Silence mode.")
  }
  
  if (n_cores == 0) {
    options(RcppML.threads = n_cores)
    message('Using all cores possible: ', parallel::detectCores())
    n_cores = parallel::detectCores()
  } else {
    options(RcppML.threads = n_cores)
    message('Using ', n_cores, ' cores.')
  }
  
  message("Now running k = ", k)
  message("Using ", n_cores, " cores.")
  
  tictoc::tic(msg = "Fast NMF sequentially")
  fnmf_model_lst = list()
  JasonToolBox::initiatePB('.I1')
  for (i_rep in 1:n_reps) {
    tictoc::tic()
    i_seed = 1000 * k + i_rep
    
    test_mt = merge_mt[sample(1:nrow(merge_mt),
                              size = round(nrow(merge_mt) * row_subset),
                              replace = F),
                       sample(1:ncol(merge_mt),
                              size = round(ncol(merge_mt) * col_subset),
                              replace = F)]
    
    fnmf_model_lst[[i_rep]] = list(
      subset_cols = colnames(test_mt),
      subset_rows = rownames(test_mt),
      nmf_obj = RcppML::nmf(
        data = test_mt,
        k = k,
        tol = 1e-5,
        seed = i_seed
      )
    )
    irep_time = toc(quiet = T)
    
    
    
    JasonToolBox::processBar('.I1', i_rep, n_reps)
  }
  f_time_seq = tictoc::toc()
  
  
  
  invisible(gc())
  return(fnmf_model_lst)
}



# #' Identify consensus programs
# #' 
# #' A wrapper function to merge similar programs into one
# #' @param NMFs List of Results from repetitive NMF runs.
# #' @param n_cores Integer. Number of cores to use, recommended to use multicore. Default = 0, all cores possible.
# #' @param n_reps Integer. Number of replication for each K.
# #' @param reps_ratio Double between \(0, 1\). Keeping clusters that have at least `round(n_reps * reps_ratio)` programs.
# #' @param return_plot Logical. If return the plots for checking.
# #' @param verbose Logical. If print the logs.
# #' @return A named list with consensus programs
# #' @export
# #' 
# identifyConsensusProgram = function(NMFs,
#                                     n_reps,
#                                     n_cores = 0,
#                                     reps_ratio = 0.05, 
#                                     return_plot = F,
#                                     verbose = T) {
#   
#   # require(pbapply)
#   # require(tictoc)
#   # require(pbmcapply)
#   # require(parallelDist)
#   # require(igraph)
#   
#   if (n_cores == 0) {
#     message('Using all cores possible: ', parallel::detectCores() - 1)
#     n_cores = parallel::detectCores() - 1
#   } else {
#     message('Using ', n_cores, ' cores.')
#   }
#   
#   involved_samples = unique(Reduce(lapply(NMFs, function(x) {
#     x$subset_cols
#   }), f = union))
#   involved_genes = unique(Reduce(lapply(NMFs, function(x) {
#     x$subset_rows
#   }), f = union))
#   
#   
#   if (verbose)
#     message('Scaling')
#   
#   normh_lst = pbmcapply::pbmclapply(mc.cores = n_cores, 1:length(NMFs), function(i_rep) {
#     x = NMFs[[i_rep]]
#     w_mt = t(x$nmf_obj@h)
#     norm_term = apply(w_mt, 2, norm, type = '2')
#     norm_w_mt = sweep(w_mt, 2, norm_term, FUN = '/')
#     
#     miss_samples = setdiff(involved_samples, rownames(norm_w_mt))
#     norm_w_mt_mend = rbind(norm_w_mt, matrix(
#       NA,
#       nrow = length(miss_samples),
#       ncol = ncol(norm_w_mt)
#     ))
#     rownames(norm_w_mt_mend) = c(rownames(norm_w_mt), miss_samples)
#     norm_w_mt_mend = norm_w_mt_mend[involved_samples, ]
#     colnames(norm_w_mt_mend) = paste0("Rep_", i_rep, "_", colnames(norm_w_mt))
#     return(t(norm_w_mt_mend))
#   })
#   
#   normh_all_mt = do.call(rbind, normh_lst)
#   
#   # Calculate distance with 10% NA
#   if (verbose)
#     message('Calculating distance, this could take a while...')
#   if (verbose)
#     tictoc::tic()
#   # norm_dist_obj = stats::dist(normh_all_mt)
#   # norm_dist_mt = as.matrix(norm_dist_obj)
#   norm_dist_mt = as.matrix(parallelDist::parDist(normh_all_mt, threads = n_cores))
#   if (verbose)
#     tictoc::toc()
#   
#   # norm_dist_mtmt = as.matrix(norm_dist_mt)
#   
#   if(return_plot) {
#     p1 = norm_dist_mt[, sample(colnames(norm_dist_mt), 1)] %>% enframe() %>%
#       ggplot(aes(x = reorder(name, value), y = value)) +
#       geom_point() +
#       geom_vline(xintercept = n_reps) +
#       theme_classic(base_size = 16) +
#       theme(axis.text.x = element_blank(),
#             axis.ticks.x = element_blank()) +
#       labs(x = 'Ordered index', y = 'Euclidean distance')
#   }
#   
#   diag(norm_dist_mt) = Inf
#   
#   if (verbose) message('Calculating KNN')
#   
#   normh_knn_lst = pbmcapply::pbmclapply(rownames(norm_dist_mt), mc.cores = n_cores, function(i_row) {
#     x = norm_dist_mt[i_row,]
#     return(sort(x, decreasing = F)[1:n_reps])
#   })
#   names(normh_knn_lst) = rownames(norm_dist_mt)
#   
#   median_knn_dists = unlist(lapply(normh_knn_lst, median))
#   
#   
#   
#   neighbor_dist_cutoff = tryCatch({
#     knn_success(median_knn_dists)
#   }, error = function(e){
#     knn_fail(median_knn_dists)
#   })
#   
#   if (return_plot) {
#     p2 = median_knn_dists %>% enframe() %>%
#       mutate(cluster = factor(value > neighbor_dist_cutoff)) %>%
#       ggplot(aes(x = value, fill = cluster)) +
#       geom_histogram(bins = 50) +
#       theme_classic(base_size = 18) +
#       geom_vline(xintercept = neighbor_dist_cutoff, linetype = 2) +
#       labs(
#         title = 'Distribution of the median distances',
#         fill = 'Kmeans cluster',
#         x = 'Median distance of NN',
#         y = 'Number of programs'
#       )
#     #
#   }
#   
#   if (verbose) message('Filtering KNN')
#   # Change KNN to adjacency table
#   normh_flt_knn_dist_lst = pbmcapply::pbmclapply(rownames(norm_dist_mt), mc.cores = n_cores, function(i_row) {
#     x = norm_dist_mt[i_row,]
#     x[x > neighbor_dist_cutoff] = Inf
#     x[!names(x) %in% names(head(sort(x, decreasing = F), 100))] = Inf
#     return(x)
#   })
#   names(normh_flt_knn_dist_lst) = rownames(norm_dist_mt)
#   normh_flt_knn_dist_mt = do.call(rbind, normh_flt_knn_dist_lst)
#   normh_flt_knn_weights_mt = 1 / normh_flt_knn_dist_mt
#   
#   NotZero_flags = `!=`(rowSums(normh_flt_knn_weights_mt), 0)
#   NotZero_connected_nodes =  names(NotZero_flags)[NotZero_flags]
#   normh_flt_knn_weights_mt = normh_flt_knn_weights_mt[NotZero_connected_nodes, NotZero_connected_nodes]
#   
#   
#   if (verbose)
#     message('Identifying clusters')
#   leiden_obj = leiden::leiden(normh_flt_knn_weights_mt, resolution_parameter = 1)
#   
#   clusters_flags = table(leiden_obj) >= round(n_reps * reps_ratio)
#   chosen_clusts = as.numeric(names(clusters_flags)[clusters_flags])
#   
#   flt_leiden_ind = leiden_obj
#   flt_leiden_ind[!flt_leiden_ind %in% chosen_clusts] = 0
#   flt_leiden_ind_nozero = flt_leiden_ind[flt_leiden_ind != 0]
#   
#   if (verbose)
#     message('Merging programs')
#   
#   merged_program_mt =
#     pbapply::pbapply(
#       normh_all_mt[rownames(normh_flt_knn_weights_mt)[flt_leiden_ind != 0], ], 
#       2, 
#       function(x) {
#         tapply(x, INDEX = flt_leiden_ind_nozero, median, na.rm = T)
#       })
#   rownames(merged_program_mt) = sprintf("MergedProgram_%.2d", 1:length(unique(flt_leiden_ind_nozero)))
#   
#   # fill NA with 0 because this gene is completely lost in a program cluster
#   merged_program_mt[is.na(merged_program_mt)] = 0
#   
#   if(return_plot){
#     return(list(merged_programs = merged_program_mt, plots = list(p1, p2)))
#   } else {
#     return(list(merged_programs = merged_program_mt))
#   }
# }

knn_fail = function(median_knn_dists) {
  warning('Cannot find cutoff using K-means.')
  warning(
    'The range of the distribution of distance is:',
    min(median_knn_dists),
    " ",
    max(median_knn_dists)
  )
  stopifnot(min(median_knn_dists) != max(median_knn_dists))
  return(median(median_knn_dists))
}

knn_success = function(median_knn_dists) {
  kmeans_obj = kmeans(median_knn_dists, 2)
  
  Kmeans_higher_group =
    names(kmeans_obj$centers[, 1])[kmeans_obj$centers[, 1] == max(kmeans_obj$centers)]
  
  return(median(c(
    max(median_knn_dists[kmeans_obj$cluster != Kmeans_higher_group]), min(median_knn_dists[kmeans_obj$cluster == Kmeans_higher_group])
  )))
}

# #' Identify consensus programs.ver2
# #' 
# #' A wrapper function to merge similar programs into one
# #' @param NMFs List of Results from repetitive NMF runs.
# #' @param merge_data A named list return by `preprocessing`
# #' @param n_cores Integer. Number of cores to use, recommended to use multicore. Default = 0, all cores possible.
# #' @param n_reps Integer. Number of replication for each K.
# #' @param reps_ratio Double between \(0, 1\). Keeping clusters that have at least `round(n_reps * reps_ratio)` programs.
# #' @param return_plot Logical. If return the plots for checking.
# #' @param verbose Logical. If print the logs.
# #' @return A named list with consensus programs
# #' @export
# #' 
# identifyConsensusProgram2 = function(NMFs,
#                                      merge_data,
#                                      n_reps,
#                                      n_cores = 0,
#                                      reps_ratio = 0.05,
#                                      return_plot = F,
#                                      verbose = T) {
#   
#   # require(pbapply)
#   # require(tictoc)
#   # require(pbmcapply)
#   # require(parallelDist)
#   # require(igraph)
#   
#   if (n_cores == 0) {
#     message('Using all cores possible: ', parallel::detectCores() - 1)
#     n_cores = parallel::detectCores() - 1
#   } else {
#     message('Using ', n_cores, ' cores.')
#   }
#   
#   involved_samples = unique(Reduce(lapply(NMFs, function(x) {
#     x$subset_cols
#   }), f = union))
#   involved_genes = unique(Reduce(lapply(NMFs, function(x) {
#     x$subset_rows
#   }), f = union))
#   
#   
#   if (verbose) message('Merging H matrices')
#   normh_lst = pbmcapply::pbmclapply(mc.cores = n_cores, 1:length(NMFs), function(i_rep) {
#     x = NMFs[[i_rep]]
#     norm_h_mt = t(x$nmf_obj@h)
#     # norm_term = apply(h_mt, 2, norm, type = '2')
#     # norm_h_mt = sweep(h_mt, 2, norm_term, FUN = '/')
#     
#     miss_samples = setdiff(involved_samples, rownames(norm_h_mt))
#     norm_h_mt_mend = rbind(norm_h_mt, matrix(
#       NA,
#       nrow = length(miss_samples),
#       ncol = ncol(norm_h_mt)
#     ))
#     rownames(norm_h_mt_mend) = c(rownames(norm_h_mt), miss_samples)
#     norm_h_mt_mend = norm_h_mt_mend[involved_samples,]
#     colnames(norm_h_mt_mend) = paste0("Rep_", i_rep, "_", colnames(norm_h_mt))
#     return(t(norm_h_mt_mend))
#   })
#   normh_all_mt = do.call(rbind, normh_lst)
#   
#   if (verbose) message('Merging W matrices')
#   normw_lst = pbmcapply::pbmclapply(mc.cores = n_cores, 1:length(NMFs), function(i_rep) {
#     x = NMFs[[i_rep]]
#     norm_w_mt = x$nmf_obj@w
#     # norm_term = apply(w_mt, 2, norm, type = '2')
#     # norm_w_mt = sweep(w_mt, 2, norm_term, FUN = '/')
#     
#     miss_genes = setdiff(involved_genes, rownames(norm_w_mt))
#     norm_w_mt_mend = rbind(norm_w_mt, matrix(
#       NA,
#       nrow = length(miss_genes),
#       ncol = ncol(norm_w_mt)
#     ))
#     rownames(norm_w_mt_mend) = c(rownames(norm_w_mt), miss_genes)
#     norm_w_mt_mend = norm_w_mt_mend[involved_genes,]
#     colnames(norm_w_mt_mend) = paste0("Rep_", i_rep, "_", colnames(norm_w_mt))
#     return(t(norm_w_mt_mend))
#   })
#   normw_all_mt = do.call(rbind, normw_lst)
#   
#   # Calculate distance with 10% NA
#   if (verbose)
#     message('Calculating distance, this could take a while...')
#   if (verbose)
#     tictoc::tic()
#   # norm_dist_obj = stats::dist(normh_all_mt)
#   # norm_dist_mt = as.matrix(norm_dist_obj)
#   norm_dist_mt = as.matrix(parallelDist::parDist(normh_all_mt, threads = n_cores))
#   if (verbose)
#     tictoc::toc()
#   
#   # norm_dist_mtmt = as.matrix(norm_dist_mt)
#   
#   if (return_plot) {
#     p1 = norm_dist_mt[, sample(colnames(norm_dist_mt), 1)] %>% enframe() %>%
#       ggplot(aes(x = reorder(name, value), y = value)) +
#       geom_point() +
#       geom_vline(xintercept = n_reps) +
#       theme_classic(base_size = 16) +
#       theme(axis.text.x = element_blank(),
#             axis.ticks.x = element_blank()) +
#       labs(x = 'Ordered index', y = 'Euclidean distance')
#   }
#   
#   diag(norm_dist_mt) = Inf
#   
#   if (verbose)
#     message('Calculating KNN')
#   
#   normh_knn_lst = pbmcapply::pbmclapply(rownames(norm_dist_mt), mc.cores = n_cores, function(i_row) {
#     x = norm_dist_mt[i_row, ]
#     return(sort(x, decreasing = F)[1:n_reps])
#   })
#   names(normh_knn_lst) = rownames(norm_dist_mt)
#   
#   median_knn_dists = unlist(lapply(normh_knn_lst, median))
#   
#   
#   
#   neighbor_dist_cutoff = tryCatch({
#     knn_success(median_knn_dists)
#   }, error = function(e) {
#     knn_fail(median_knn_dists)
#   })
#   
#   if (return_plot) {
#     p2 = median_knn_dists %>% enframe() %>%
#       mutate(cluster = factor(value > neighbor_dist_cutoff)) %>%
#       ggplot(aes(x = value, fill = cluster)) +
#       geom_histogram(bins = 50) +
#       theme_classic(base_size = 18) +
#       geom_vline(xintercept = neighbor_dist_cutoff, linetype = 2) +
#       labs(
#         title = 'Distribution of the median distances',
#         fill = 'Kmeans cluster',
#         x = 'Median distance of NN',
#         y = 'Number of programs'
#       )
#     #
#   }
#   
#   if (verbose)
#     message('Filtering KNN')
#   # Change KNN to adjacency table
#   normh_flt_knn_dist_lst = pbmcapply::pbmclapply(rownames(norm_dist_mt), mc.cores = n_cores, function(i_row) {
#     x = norm_dist_mt[i_row, ]
#     x[x > neighbor_dist_cutoff] = Inf
#     x[!names(x) %in% names(head(sort(x, decreasing = F), 100))] = Inf
#     return(x)
#   })
#   names(normh_flt_knn_dist_lst) = rownames(norm_dist_mt)
#   normh_flt_knn_dist_mt = do.call(rbind, normh_flt_knn_dist_lst)
#   normh_flt_knn_weights_mt = 1 / normh_flt_knn_dist_mt
#   
#   NotZero_flags = `!=`(rowSums(normh_flt_knn_weights_mt), 0)
#   NotZero_connected_nodes =  names(NotZero_flags)[NotZero_flags]
#   normh_flt_knn_weights_mt = normh_flt_knn_weights_mt[NotZero_connected_nodes, NotZero_connected_nodes]
#   
#   
#   if (verbose) message('Identifying clusters')
#   leiden_obj = leiden::leiden(normh_flt_knn_weights_mt, resolution_parameter = 1)
#   
#   clusters_flags = table(leiden_obj) >= round(n_reps * reps_ratio)
#   chosen_clusts = as.numeric(names(clusters_flags)[clusters_flags])
#   
#   flt_leiden_ind = leiden_obj
#   flt_leiden_ind[!flt_leiden_ind %in% chosen_clusts] = 0
#   flt_leiden_ind_nozero = flt_leiden_ind[flt_leiden_ind != 0]
#   
#   if (verbose) message('Merging programs')
#   
#   consensus_h_mt =
#     pbapply::pbapply(normh_all_mt[rownames(normh_flt_knn_weights_mt)[flt_leiden_ind != 0],],
#                      2,
#                      function(x) {
#                        tapply(x, INDEX = flt_leiden_ind_nozero, median, na.rm = T)
#                      })
#   rownames(consensus_h_mt) = sprintf("MergedProgram_%.2d", 1:length(unique(flt_leiden_ind_nozero)))
#   
#   consensus_w_mt =
#     pbapply::pbapply(normw_all_mt[rownames(normh_flt_knn_weights_mt)[flt_leiden_ind != 0],],
#                      2,
#                      function(x) {
#                        tapply(x, INDEX = flt_leiden_ind_nozero, median, na.rm = T)
#                      })
#   rownames(consensus_h_mt) = sprintf("MergedProgram_%.2d", 1:length(unique(flt_leiden_ind_nozero)))
#   
#   # fill NA with 0 because this gene is completely lost in a program cluster
#   consensus_h_mt[is.na(consensus_h_mt)] = 0
#   consensus_h_mt = apply(consensus_h_mt, 1, function(x){x/sum(x)})
#   
#   consensus_w_mt[is.na(consensus_w_mt)] = 0
#   consensus_w_mt = apply(consensus_w_mt, 1, function(x){x/sum(x)})
#   
#   # consensus_w_mt = RcppML::predict.nmf(w = t(consensus_h_mt), data = Matrix::t(merge_data$merge_MM))
#   # consensus_diag = rowSums(consensus_w_mt)
#   # consensus_w_mt = apply(consensus_w_mt, 1, function(x){x/sum(x)})
#   # consensus_h_mt = t(consensus_h_mt)
#   colnames(consensus_w_mt) = colnames(consensus_h_mt)
#   # names(consensus_diag) = colnames(consensus_w_mt)
#   
#   if(return_plot){
#     return(list(
#       consensus_W = consensus_w_mt, 
#       # consensus_d = consensus_diag, 
#       consensus_H = t(consensus_h_mt), 
#       plots = list(p1, p2)))
#   } else {
#     return(list(
#       consensus_W = consensus_w_mt, 
#       # consensus_d = consensus_diag, 
#       consensus_H = consensus_h_mt
#     ))
#   }
# }


#' Identify consensus programs.ver3
#' 
#' A wrapper function to merge similar programs into one
#' @param NMFs List of Results from repetitive NMF runs.
#' @param merge_data A named list return by `preprocesssing`
#' @param n_cores Integer. Number of cores to use, recommended to use multicore. Default = 0, all cores possible.
#' @param n_reps Integer. Number of replication for each K.
#' @param reps_ratio Double between \(0, 1\). Keeping clusters that have at least `round(n_reps * reps_ratio)` programs.
#' @param return_plot Logical. If return the plots for checking.
#' @param verbose Logical. If print the logs.
#' @return A named list with consensus programs
#' @export
#' 
identifyConsensusProgram3 = function(NMFs,
                                     merge_data,
                                     n_reps,
                                     n_cores = 0,
                                     reps_ratio = 0.05,
                                     return_plot = F,
                                     verbose = T) {
  
  # require(pbapply)
  # require(tictoc)
  # require(pbmcapply)
  # require(parallelDist)
  # require(igraph)
  
  if (n_cores == 0) {
    message('Using all cores possible: ', parallel::detectCores() - 1)
    n_cores = parallel::detectCores() - 1
  } else {
    message('Using ', n_cores, ' cores.')
  }
  
  involved_samples = unique(Reduce(lapply(NMFs, function(x) {
    x$subset_cols
  }), f = union))
  involved_genes = unique(Reduce(lapply(NMFs, function(x) {
    x$subset_rows
  }), f = union))
  
  
  if (verbose) message('Merging H matrices')
  normh_lst = pbmcapply::pbmclapply(mc.cores = n_cores, 1:length(NMFs), function(i_rep) {
    x = NMFs[[i_rep]]
    # norm_h_mt = t(x$nmf_obj@h)
    norm_h_mt = t(as.matrix(Matrix::Diagonal(x = NMFs[[i_rep]]$nmf_obj@d) %*% NMFs[[i_rep]]$nmf_obj@h))
    # norm_term = apply(h_mt, 2, norm, type = '2')
    # norm_h_mt = sweep(h_mt, 2, norm_term, FUN = '/')
    
    miss_samples = setdiff(involved_samples, rownames(norm_h_mt))
    norm_h_mt_mend = rbind(norm_h_mt, matrix(
      NA,
      nrow = length(miss_samples),
      ncol = ncol(norm_h_mt)
    ))
    rownames(norm_h_mt_mend) = c(rownames(norm_h_mt), miss_samples)
    norm_h_mt_mend = norm_h_mt_mend[involved_samples,]
    colnames(norm_h_mt_mend) = paste0("Rep_", i_rep, "_", colnames(norm_h_mt))
    return(t(norm_h_mt_mend))
  })
  normh_all_mt = do.call(rbind, normh_lst)
  
  if (verbose) message('Merging W matrices')
  normw_lst = pbmcapply::pbmclapply(mc.cores = n_cores, 1:length(NMFs), function(i_rep) {
    x = NMFs[[i_rep]]
    norm_w_mt = x$nmf_obj@w
    # norm_term = apply(w_mt, 2, norm, type = '2')
    # norm_w_mt = sweep(w_mt, 2, norm_term, FUN = '/')
    
    miss_genes = setdiff(involved_genes, rownames(norm_w_mt))
    norm_w_mt_mend = rbind(norm_w_mt, matrix(
      NA,
      nrow = length(miss_genes),
      ncol = ncol(norm_w_mt)
    ))
    rownames(norm_w_mt_mend) = c(rownames(norm_w_mt), miss_genes)
    norm_w_mt_mend = norm_w_mt_mend[involved_genes,]
    colnames(norm_w_mt_mend) = paste0("Rep_", i_rep, "_", colnames(norm_w_mt))
    return(t(norm_w_mt_mend))
  })
  normw_all_mt = do.call(rbind, normw_lst)
  
  d_lst = pbapply::pblapply(X = 1:length(NMFs), function(i_rep){
    tmp_vec = setNames(NMFs[[i_rep]]$nmf_obj@d, 
                       nm = paste0('Rep_', i_rep, '_', colnames(NMFs[[i_rep]]$nmf_obj@w)))
    return(tmp_vec)
  })
  d_all_vec = do.call(c, d_lst)
  
  
  # Calculate distance with 10% NA
  if (verbose)
    message('Calculating distance, this could take a while...')
  if (verbose)
    tictoc::tic()
  # norm_dist_obj = stats::dist(normh_all_mt)
  # norm_dist_mt = as.matrix(norm_dist_obj)
  norm_dist_mt = as.matrix(parallelDist::parDist(normh_all_mt, threads = n_cores))
  if (verbose)
    tictoc::toc()
  
  # norm_dist_mtmt = as.matrix(norm_dist_mt)
  
  if (return_plot) {
    p1 = norm_dist_mt[, sample(colnames(norm_dist_mt), 1)] %>% tibble::enframe() %>%
      ggplot(aes(x = reorder(name, value), y = value)) +
      geom_point() +
      geom_vline(xintercept = n_reps) +
      theme_classic(base_size = 16) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      labs(x = 'Ordered index', y = 'Euclidean distance')
  }
  
  diag(norm_dist_mt) = Inf
  
  if (verbose)
    message('Calculating AKNN')
  
  normh_knn_lst = pbmcapply::pbmclapply(rownames(norm_dist_mt), mc.cores = n_cores, function(i_row) {
    x = norm_dist_mt[i_row, ]
    return(sort(x, decreasing = F)[1:n_reps])
  })
  names(normh_knn_lst) = rownames(norm_dist_mt)
  
  median_knn_dists = unlist(lapply(normh_knn_lst, median))
  
  
  
  neighbor_dist_cutoff = tryCatch({
    knn_success(median_knn_dists)
  }, error = function(e) {
    knn_fail(median_knn_dists)
  })
  
  if (return_plot) {
    p2 = median_knn_dists %>% tibble::enframe() %>%
      mutate(cluster = factor(value > neighbor_dist_cutoff)) %>%
      ggplot(aes(x = value, fill = cluster)) +
      geom_histogram(bins = 50) +
      theme_classic(base_size = 18) +
      geom_vline(xintercept = neighbor_dist_cutoff, linetype = 2) +
      labs(
        title = 'Distribution of the median distances',
        fill = 'Kmeans cluster',
        x = 'Median distance of NN',
        y = 'Number of programs'
      )
    #
  }
  
  if (verbose)
    message('Filtering KNN')
  # Change KNN to adjacency table
  normh_flt_knn_dist_lst = pbmcapply::pbmclapply(rownames(norm_dist_mt), mc.cores = n_cores, function(i_row) {
    x = norm_dist_mt[i_row, ]
    x[x > neighbor_dist_cutoff] = Inf
    x[!names(x) %in% names(head(sort(x, decreasing = F), 100))] = Inf
    return(x)
  })
  names(normh_flt_knn_dist_lst) = rownames(norm_dist_mt)
  normh_flt_knn_dist_mt = do.call(rbind, normh_flt_knn_dist_lst)
  normh_flt_knn_weights_mt = 1 / normh_flt_knn_dist_mt
  
  NotZero_flags = `!=`(rowSums(normh_flt_knn_weights_mt), 0)
  NotZero_connected_nodes =  names(NotZero_flags)[NotZero_flags]
  normh_flt_knn_weights_mt = normh_flt_knn_weights_mt[NotZero_connected_nodes, NotZero_connected_nodes]
  
  
  if (verbose) message('Identifying clusters')
  leiden_obj = leiden::leiden(normh_flt_knn_weights_mt, resolution_parameter = 1)
  
  clusters_flags = table(leiden_obj) >= round(n_reps * reps_ratio)
  chosen_clusts = as.numeric(names(clusters_flags)[clusters_flags])
  
  flt_leiden_ind = leiden_obj
  flt_leiden_ind[!flt_leiden_ind %in% chosen_clusts] = 0
  flt_leiden_ind_nozero = flt_leiden_ind[flt_leiden_ind != 0]
  
  if (verbose) message('Merging programs')
  
  consensus_H_mt =
    pbapply::pbapply(normh_all_mt[rownames(normh_flt_knn_weights_mt)[flt_leiden_ind != 0],],
                     2,
                     function(x) {
                       tapply(x, INDEX = flt_leiden_ind_nozero, median, na.rm = T)
                     })
  rownames(consensus_H_mt) = sprintf("cProgram %.2d", 1:length(unique(flt_leiden_ind_nozero)))
  consensus_H_mt[is.na(consensus_H_mt)] = 0
  consensus_d_vec = rowSums(consensus_H_mt)
  # consensus_h_mt2 = apply(consensus_h_mt, 1, function(x){x/sum(x)})
  consensus_h_mt = consensus_H_mt/rowSums(consensus_H_mt)
  
  consensus_w_mt =
    pbapply::pbapply(normw_all_mt[rownames(normh_flt_knn_weights_mt)[flt_leiden_ind != 0],],
                     2,
                     function(x) {
                       tapply(x, INDEX = flt_leiden_ind_nozero, median, na.rm = T)
                     })
  consensus_w_mt[is.na(consensus_w_mt)] = 0
  consensus_w_mt = apply(consensus_w_mt, 1, function(x){x/sum(x)})
  colnames(consensus_w_mt) = rownames(consensus_h_mt)
  
  # consensus_d_vec =
  #   pbapply::pbtapply(
  #     X = d_all_vec[rownames(normh_flt_knn_weights_mt)[flt_leiden_ind != 0]], 
  #     INDEX = flt_leiden_ind_nozero,
  #     FUN = median, na.rm = T)
  # consensus_d_vec = setNames(as.numeric(consensus_d_vec), nm = rownames(consensus_h_mt))
  
  if(return_plot){
    return(list(
      consensus_W = consensus_w_mt, 
      consensus_d = consensus_d_vec, 
      consensus_H = consensus_h_mt, 
      plots = list(p1, p2)))
  } else {
    return(list(
      consensus_W = consensus_w_mt, 
      consensus_d = consensus_d_vec, 
      consensus_H = consensus_h_mt
      ))
  }
}

#' Process the consensus programs
#' 
#' A wrapper function to process the consensus programs.
#' @param merged_program_mt A matrix
#' @param sources A named vector with the sources information
#' @param program_cutoff Double, default = 0.1, between (0-1]
#' @return A processed matrix
#' @export
#' 
processConsensusProgram2 = function(merged_program_mt,
                                   sources,
                                   program_cutoff = 0.1) {
  
  consensus_H_mt = apply(merged_program_mt, 1, function(x) { x / sum(x) }) # normalize by each program
  
  sources = sources[rownames(consensus_H_mt)]
  i_cat = unique(sources)[1]
  cut_mt_normbySample = apply(consensus_H_mt, 2, function(x) {
    y = x
    y[y < (max(y) * program_cutoff)] = 0 # For each program, samples with prog activity less than 10% of the maximum activity are set to 0
    y[is.na(y)] = 0
    y[sources == i_cat] = (y[sources == i_cat]) / (max(y[sources == i_cat]) + .Machine$double.xmin)
    y[sources != i_cat] = (y[sources != i_cat]) / (max(y[sources != i_cat]) + .Machine$double.xmin)
    return(y)
  }) 
  cut_mt_normbySample = t(cut_mt_normbySample)
  
  return(cut_mt_normbySample)
}

#' Process the consensus programs
#' 
#' A wrapper function to process the consensus programs.
#' @param merged_program_mt A matrix
#' @param seurat_obj A seurat object
#' @param program_cutoff Double, default = 0.1, between (0-1]
#' @return A processed matrix
#' @export
#' 
processConsensusProgram = function(merged_program_mt,
                                   seurat_obj,
                                   program_cutoff = 0.1) {
  
  consensus_H_mt = apply(merged_program_mt, 1, function(x) { x / sum(x) }) # normalize by each program
  
  source_vec = setNames(nm = seurat_obj$name, seurat_obj$source)
  source_vec = source_vec[rownames(consensus_H_mt)]
  cut_mt_normbySample = apply(consensus_H_mt, 2, function(x) {
    y = x
    y[y < (max(y) * program_cutoff)] = 0 # For each program, samples with prog activity less than 10% of the maximum activity are set to 0
    y[source_vec == 'bulk'] = (y[source_vec == 'bulk']) / (max(y[source_vec == 'bulk']) + 1e-32)
    y[source_vec != 'bulk'] = (y[source_vec != 'bulk']) / (max(y[source_vec != 'bulk']) + 1e-32)
    return(y)
  }) 
  cut_mt_normbySample = t(cut_mt_normbySample)
  
  return(cut_mt_normbySample)
}


#' Check and plot the consensus programs
#' 
#' @param seurat_obj A seurat object
#' @param program_mt A program matrix
#' @param program_cutoff Double, default = 0.1, between (0-1]
#' @return A named list
#' @export
#' 
checkConsensusProgram = function(seurat_obj, program_mt, program_cutoff = 0.1, raster = F) {
  
  # require(ggplot2)
  # require(dplyr)
  # require(Seurat) # Seurat dependencies should be removed in the future.
  
  plt_tbl = as.data.frame(t(program_mt), stringAsfactor = F) %>% 
    mutate(.before = "cProgram 01", cellID = colnames(program_mt))
  
  seurat_obj@meta.data = left_join(seurat_obj@meta.data, plt_tbl, by = c("name" = "cellID"))
  rownames(seurat_obj@meta.data) = seurat_obj@meta.data$name
  
  
  plt_lst = pbapply::pblapply(1:nrow(program_mt), function(i_program) {
    p1 = Seurat::FeaturePlot(
      seurat_obj,
      reduction = "tsne",
      features = sprintf('cProgram %.2d', i_program),
      order = T,
      raster = T,
      pt.size = 3
    )
    # hist_plt_lst2 =
    tmp_plt = seurat_obj@meta.data[, c('name',
                                       'type',
                                       'celltype',
                                       'source',
                                       sprintf('cProgram %.2d', i_program))]
    colnames(tmp_plt) = c('name', 'type', 'celltype', 'source', 'ProgramValue')
    
    if(!requireNamespace("ggrastr", quietly = TRUE)){
      warning('ggrastr not installed. Setting `raster` to FALSE')
      raster = F
    }
    
    if (raster) {
      p2 =
        tmp_plt %>% ggplot(aes(
          x = reorder(name, ProgramValue),
          y = ProgramValue,
          color = type
        )) +
        ggrastr::rasterise(geom_point(
          alpha = .5,
          shape = 'o',
          size = 6
        )) +
        geom_hline(
          yintercept = max(tmp_plt$ProgramValue, na.rm = T) * program_cutoff,
          linetype = 2
        ) +
        theme_classic(base_size = 16) +
        theme(
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()
        ) +
        labs(title = sprintf('cProgram %.2d', i_program))
    } else {
      p2 =
        tmp_plt %>% ggplot(aes(
          x = reorder(name, ProgramValue),
          y = ProgramValue,
          color = type
        )) +
        geom_point(
          alpha = .5,
          shape = 'o',
          size = 6
        ) +
        geom_hline(
          yintercept = max(tmp_plt$ProgramValue, na.rm = T) * program_cutoff,
          linetype = 2
        ) +
        theme_classic(base_size = 16) +
        theme(
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()
        ) +
        labs(title = sprintf('cProgram %.2d', i_program))
    }
    return(list(p1 = p1, p2 = p2))
  })
  
  return(list(plt_lst = plt_lst, sce_plt = seurat_obj))
}

# Correlative graph ----

# #' Check consensus program values.ver1
# #' 
# #' @param program_mt Matrix, processed
# #' @param sources Named strings, c\('sample_name' = 'bulk|singlecell'\)
# #' @param celltypes Named vector, should be a cell-id-named vector of cell types
# #' @param nonzero_ratio Value between 0-1, c\('sample_name' = 'bulk|singlecell'\)
# #' @return A named list with the correlative information between celltypes.
# #' @export
# #' 
# calcCorrelative = function(program_mt,
#                            sources,
#                            celltypes,
#                            nonzero_ratio = 0.1,
#                            rho_cutoff = 0,
#                            qval_cutoff = 0.1) {
#   
#   # require(dplyr)
#   
#   
#   bulksID = names(sources)[sources == 'bulk']
#   bulk_cut_mt = program_mt[, bulksID]
#   
#   progs2run_idx = rowSums(bulk_cut_mt > 0) > round(nonzero_ratio * length(bulksID))
#   progs2run = names(progs2run_idx)[progs2run_idx]
#   
#   bulk_cut_mt = program_mt[progs2run, bulksID]
#   
#   test_combn_df = data.frame(t(combn(rownames(bulk_cut_mt), 2)),
#                              rho = Inf,
#                              pval = Inf)
#   colnames(test_combn_df)[1:2] = c('Prog1', 'Prog2')
#   
#   message('Identifying correlative programs')
#   test_lst = pbapply::pblapply(
#     1:nrow(test_combn_df),
#     FUN = function(i_row) {
#       x_label = test_combn_df[i_row, 'Prog1', drop = T]
#       y_label = test_combn_df[i_row, 'Prog2', drop = T]
#       x = bulk_cut_mt[x_label, ]
#       y = bulk_cut_mt[y_label, ]
#       # t_idx = x != 0 | y != 0
#       # x = x[t_idx]
#       # y = y[t_idx]
#       tryCatch({
#         suppressWarnings({
#           test_res = cor.test(x, y, method = 'spearman', exact = F)
#         })
#       }, error = function(e) {
#         warning(e)
#         test_res = list()
#       })
#       tmp_df = data.frame(
#         Prog1 = test_combn_df[i_row, 'Prog1', drop = T],
#         Prog2 = test_combn_df[i_row, 'Prog2', drop = T],
#         rho = test_res$estimate,
#         pval = test_res$p.value
#       )
#     }
#   )
#   test_combn_df = do.call(rbind, test_lst)
#   rownames(test_combn_df) = NULL
#   
#   # Mapping correlative programs to cell types
#   scsID = names(sources)[sources != 'bulk']
#   sc_cut_mt = program_mt[, scsID]
#   
#   sctype_prog_mt = sc_cut_mt %>% 
#     apply(1, function(x){tapply(x, celltypes[scsID], mean)}) %>% 
#     apply(1, function(x){x/(sum(x) + 1e-32)}) %>% 
#     apply(1, function(x){x/(max(x) + 1e-32)}) %>% 
#     t()
#   
#   correlative_progs_tbl = test_combn_df %>% dplyr::filter(rho > rho_cutoff, pval < pval_cutoff)
#   
#   shared_progs = with(correlative_progs_tbl, unique(c(Prog1, Prog2)))
#   sctype_prog_mt = sctype_prog_mt[shared_progs, ]
#   
#   message('Mapping correlative programs to cell types...')
#   sharedprogs2celltype_lst = pbapply::pblapply(1:nrow(correlative_progs_tbl), function(i_row) {
#     # correlative_progs_tbl[i_row, ]
#     Prog1_vec = sctype_prog_mt[correlative_progs_tbl[i_row, 'Prog1', drop = T], ]
#     Prog2_vec = sctype_prog_mt[correlative_progs_tbl[i_row, 'Prog2', drop = T], ]
#     tmp_outer_mt = outer(Prog1_vec,
#                          Prog2_vec)
#     tmp_outer_mt = tmp_outer_mt + t(tmp_outer_mt)
#     return(tmp_outer_mt)
#   })
#   sharedprogs2celltype_mt = Reduce('+', sharedprogs2celltype_lst)
# 
#   
#   return(list(test_df = test_combn_df, 
#               all_bulk_mt = program_mt[, bulksID],
#               bulk_mt = bulk_cut_mt, 
#               sctype_prog_mt = sctype_prog_mt, 
#               correlative_mt = sharedprogs2celltype_mt))
# }
# 
# 
# #' Check consensus program values.ver2
# #' 
# #' @param program_mt Matrix, processed
# #' @param sources Named strings, c\('sample_name' = 'bulk|singlecell'\)
# #' @param celltypes Named vector, should be a cell-id-named vector of cell types
# #' @param nonzero_ratio Value between 0-1, c\('sample_name' = 'bulk|singlecell'\)
# #' @param rho_cutoff The remaining program pairs must have spearman correlation rho larger than this. Default = 0.
# #' @param qval_cutoff The remaining program pairs must have qvalue smaller than this. Default = 0.05. BH method.
# #' @return A named list with the correlative information between celltypes.
# #' @export
# #' 
# calcCorrelative2 = function(program_mt,
#                             sources,
#                             celltypes,
#                             nonzero_ratio = 0.1,
#                             rho_cutoff = 0,
#                             qval_cutoff = 0.05) {
#   
#   # require(dplyr)
#   
#   bulksID = names(sources)[sources == 'bulk']
#   bulk_cut_mt = program_mt[, bulksID]
#   
#   progs2run_idx = rowSums(bulk_cut_mt > 0) > round(nonzero_ratio * length(bulksID))
#   progs2run = names(progs2run_idx)[progs2run_idx]
#   
#   bulk_cut_mt = program_mt[progs2run, bulksID]
#   
#   test_combn_df = data.frame(t(combn(rownames(bulk_cut_mt), 2)),
#                              rho = Inf,
#                              pval = Inf)
#   colnames(test_combn_df)[1:2] = c('Prog1', 'Prog2')
#   
#   message('Identifying correlative programs')
#   
#   test_lst = pbapply::pblapply(
#     1:nrow(test_combn_df),
#     FUN = function(i_row) {
#       x_label = test_combn_df[i_row, 'Prog1', drop = T]
#       y_label = test_combn_df[i_row, 'Prog2', drop = T]
#       x = bulk_cut_mt[x_label, ]
#       y = bulk_cut_mt[y_label, ]
#       # t_idx = x != 0 | y != 0
#       # x = x[t_idx]
#       # y = y[t_idx]
#       tryCatch({
#         suppressWarnings({
#           test_res = cor.test(x, y, method = 'spearman')
#         })
#       }, error = function(e) {
#         warning(e)
#         test_res = list()
#       })
#       tmp_df = data.frame(
#         Prog1 = test_combn_df[i_row, 'Prog1', drop = T],
#         Prog2 = test_combn_df[i_row, 'Prog2', drop = T],
#         rho = test_res$estimate,
#         pval = test_res$p.value
#       )
#     }
#   )
#   test_combn_df = do.call(rbind, test_lst)
#   rownames(test_combn_df) = NULL
#   test_combn_df$qval = p.adjust(test_combn_df$pval, method = 'BH')
#   
#   # Mapping correlative programs to cell types
#   scsID = names(sources)[sources != 'bulk']
#   sc_cut_mt = program_mt[, scsID]
#   
#   sctype_prog_mt = sc_cut_mt %>% 
#     apply(1, function(x){tapply(x, celltypes[scsID], mean)}) %>% 
#     apply(1, function(x){x/(sum(x) + 1e-32)}) %>% 
#     apply(1, function(x){x/(max(x) + 1e-32)}) %>% 
#     t()
#   
#   correlative_progs_tbl = test_combn_df %>% dplyr::filter(rho > rho_cutoff, qval < qval_cutoff)
#   
#   shared_progs = with(correlative_progs_tbl, unique(c(Prog1, Prog2)))
#   sctype_prog_mt = sctype_prog_mt[shared_progs, ]
#   
#   message('Mapping correlative programs to cell types...')
#   sharedprogs2celltype_lst = pbapply::pblapply(1:nrow(correlative_progs_tbl), function(i_row) {
#     # correlative_progs_tbl[i_row, ]
#     Prog1_vec = sctype_prog_mt[correlative_progs_tbl[i_row, 'Prog1', drop = T], ]
#     Prog2_vec = sctype_prog_mt[correlative_progs_tbl[i_row, 'Prog2', drop = T], ]
#     tmp_outer_mt = outer(Prog1_vec,
#                          Prog2_vec)
#     tmp_outer_mt = tmp_outer_mt + t(tmp_outer_mt)
#     return(tmp_outer_mt)
#   })
#   sharedprogs2celltype_mt = Reduce('+', sharedprogs2celltype_lst)
#   
#   
#   return(list(test_df = test_combn_df, 
#               all_bulk_mt = program_mt[, bulksID],
#               bulk_mt = bulk_cut_mt, 
#               sctype_prog_mt = sctype_prog_mt, 
#               correlative_mt = sharedprogs2celltype_mt))
# }

#' Check consensus program values.ver3
#' 
#' @param program_mt Matrix, processed
#' @param sources Named strings, c\('sample_name' = 'bulk|singlecell'\)
#' @param celltypes Named vector, should be a cell-id-named vector of cell types
#' @param nonzero_ratio Value between 0-1, c\('sample_name' = 'bulk|singlecell'\)
#' @param rho_cutoff The remaining program pairs must have spearman correlation rho larger than this. Default = 0.
#' @param qval_cutoff The remaining program pairs must have qvalue smaller than this. Default = 0.05. BH method.
#' @param n_permtest Number of permutation rounds to perform, default = 2000
#' @param n_cores Number of cores to use for permutation, default = the maximum available cores - 1
#' @return A named list with the correlative information between celltypes.
#' @export
#' 
calcCorrelative3 = function(program_mt,
                            sources,
                            celltypes,
                            nonzero_ratio = 0.1,
                            rho_cutoff = 0,
                            n_permtest = 2000, 
                            n_cores = parallel::detectCores() - 1,
                            qval_cutoff = 0.05) {
  
  # require(dplyr)
  
  bulksID = names(sources)[sources == 'bulk']
  bulk_cut_mt = program_mt[, bulksID]
  
  progs2run_idx = rowSums(bulk_cut_mt > 0) > round(nonzero_ratio * length(bulksID))
  progs2run = names(progs2run_idx)[progs2run_idx]
  
  bulk_cut_mt = program_mt[progs2run, bulksID]
  
  test_combn_df = data.frame(t(combn(rownames(bulk_cut_mt), 2)), rho = Inf, pval = Inf)
  colnames(test_combn_df)[1:2] = c('Prog1', 'Prog2')
  
  message('Identifying correlative programs')
  
  test_lst = pbapply::pblapply(
    1:nrow(test_combn_df),
    FUN = function(i_row) {
      x_label = test_combn_df[i_row, 'Prog1', drop = T]
      y_label = test_combn_df[i_row, 'Prog2', drop = T]
      x = bulk_cut_mt[x_label, ]
      y = bulk_cut_mt[y_label, ]
      # t_idx = x != 0 | y != 0
      # x = x[t_idx]
      # y = y[t_idx]
      tryCatch({
        suppressWarnings({
          test_res = cor.test(x, y, method = 'spearman')
        })
      }, error = function(e) {
        warning(e)
        test_res = list()
      })
      tmp_df = data.frame(
        Prog1 = test_combn_df[i_row, 'Prog1', drop = T],
        Prog2 = test_combn_df[i_row, 'Prog2', drop = T],
        rho = test_res$estimate,
        pval = test_res$p.value
      )
    }
  )
  test_combn_df = do.call(rbind, test_lst)
  rownames(test_combn_df) = NULL
  test_combn_df$qval = p.adjust(test_combn_df$pval, method = 'BH')
  
  correlative_progs_tbl = test_combn_df %>% dplyr::filter(rho > rho_cutoff, qval < qval_cutoff)
  
  # Mapping correlative programs to cell types
  scsID = names(sources)[sources != 'bulk']
  sc_cut_mt = program_mt[, scsID]
  
  sctype_prog_mt = sc_cut_mt %>% 
    apply(1, function(x){tapply(x, celltypes[scsID], mean)}) %>% 
    apply(1, function(x){x/(sum(x) + 1e-32)}) %>% 
    apply(1, function(x){x/(max(x) + 1e-32)}) %>% 
    t()
  
  
  if(!nrow(correlative_progs_tbl)){
    warning('No programs are significantly correlated!')
    return(
      list(
        test_df = test_combn_df,
        correlative_progs_tbl = correlative_progs_tbl, 
        all_bulk_mt = program_mt[, bulksID],
        bulk_mt = bulk_cut_mt,
        sctype_prog_mt = sctype_prog_mt,
        correlative_mt = NA,
        pval_df = NA, 
        sharedprogs2celltype_lst = NA, 
        perm_lst = NA 
      )
    )
  }
  
  shared_progs = with(correlative_progs_tbl, unique(c(Prog1, Prog2)))
  sctype_prog_mt = sctype_prog_mt[shared_progs, ]
  
  message('Mapping correlative programs to cell types...')
  sharedprogs2celltype_lst = pbapply::pblapply(1:nrow(correlative_progs_tbl), function(i_row) {
    # correlative_progs_tbl[i_row, ]
    Prog1_vec = sctype_prog_mt[correlative_progs_tbl$Prog1[i_row], ]
    Prog2_vec = sctype_prog_mt[correlative_progs_tbl$Prog2[i_row], ]
    
    tmp_pairmin_mt1 = outer(Prog1_vec, Prog1_vec, '-')
    # tmp_pairmin_mt1 = abs(tmp_pairmin_mt1)
    tmp_pairmin_mt1[tmp_pairmin_mt1<0] = 0
    
    tmp_pairmin_mt2 = outer(Prog2_vec, Prog2_vec, '-')
    # tmp_pairmin_mt2 = abs(tmp_pairmin_mt2)
    tmp_pairmin_mt2[tmp_pairmin_mt2<0] = 0
      
    # tmp_outer_mt = tmp_pairmin_mt1 * tmp_pairmin_mt2
    # tmp_outer_mt = outer(Prog1_vec, Prog2_vec, '*')
    tmp_outer_mt = tmp_pairmin_mt1 * t(tmp_pairmin_mt2)
    tmp_outer_mt = (tmp_outer_mt + t(tmp_outer_mt)) * correlative_progs_tbl$rho[i_row]/ sum(correlative_progs_tbl$rho)
    return(tmp_outer_mt)
  })
  sharedprogs2celltype_mt = Reduce('+', sharedprogs2celltype_lst) 
  # rm(sharedprogs2celltype_lst)
  invisible(gc())
  
  if(n_permtest != 0){
    message('Begin permutation')
    perm_lst = pbmcapply::pbmclapply(
      X = 1:n_permtest,
      mc.cores = n_cores,
      FUN = function(i_perm) {
        tmp_sharedprogs2celltype_lst = lapply(1:nrow(correlative_progs_tbl), function(i_row) {
          
          # tmp_celltypes = setNames(
          #   nm = names(celltypes), 
          #   sample(celltypes))
          # 
          # tmp_sctype_prog_mt = sc_cut_mt %>% 
          #   apply(1, function(x){tapply(x, tmp_celltypes[scsID], mean)}) %>% 
          #   apply(1, function(x){x/(sum(x) + 1e-32)}) %>% 
          #   apply(1, function(x){x/(max(x) + 1e-32)}) %>% 
          #   t()
          
          tmp_correlative_progs_tbl = correlative_progs_tbl
          # tmp_correlative_progs_tbl$Prog1 = sample(tmp_correlative_progs_tbl$Prog1)
          # tmp_correlative_progs_tbl$Prog2 = sample(tmp_correlative_progs_tbl$Prog2)
          
          tmp_sctype_prog_mt = sctype_prog_mt
          tmp_sctype_prog_mt = tmp_sctype_prog_mt[, sample(colnames(tmp_sctype_prog_mt))]
          colnames(tmp_sctype_prog_mt) = colnames(sctype_prog_mt)
          # tmp_sctype_prog_mt = tmp_sctype_prog_mt[sample(rownames(tmp_sctype_prog_mt)), ]
          # rownames(tmp_sctype_prog_mt) = rownames(sctype_prog_mt)
          # 
          Prog1_vec = tmp_sctype_prog_mt[tmp_correlative_progs_tbl$Prog1[i_row], ]
          Prog2_vec = tmp_sctype_prog_mt[tmp_correlative_progs_tbl$Prog2[i_row], ]
          # # correlative_progs_tbl[i_row, ]
          tmp_pairmin_mt1 = outer(Prog1_vec, Prog1_vec, '-')
          # tmp_pairmin_mt1 = abs(tmp_pairmin_mt1)
          tmp_pairmin_mt1[tmp_pairmin_mt1<0] = 0
          
          tmp_pairmin_mt2 = outer(Prog2_vec, Prog2_vec, '-')
          # tmp_pairmin_mt2 = abs(tmp_pairmin_mt2)
          tmp_pairmin_mt2[tmp_pairmin_mt2<0] = 0
          
          # tmp_outer_mt = tmp_pairmin_mt1 * tmp_pairmin_mt2
          # tmp_outer_mt = outer(Prog1_vec, Prog2_vec, '*')
          tmp_outer_mt = tmp_pairmin_mt1 * t(tmp_pairmin_mt2)
          tmp_outer_mt = (tmp_outer_mt + t(tmp_outer_mt)) * tmp_correlative_progs_tbl$rho[i_row]/sum(tmp_correlative_progs_tbl$rho)
          return(tmp_outer_mt)
        })
        tmp_sharedprogs2celltype_mt = Reduce('+', tmp_sharedprogs2celltype_lst)
        rm(tmp_sharedprogs2celltype_lst)
        invisible(gc())
        return(tmp_sharedprogs2celltype_mt)
      }
    )
    message('Calculate emperical p values.')
    perm_mt = do.call(pbapply::pblapply(perm_lst, function(x_mt) {
      return(x_mt[1:length(x_mt)])
    }), what = rbind)
    obs_vec = sharedprogs2celltype_mt[1:length(sharedprogs2celltype_mt)]
    
    
    succ_vec = unlist(pbapply::pblapply(1:length(sharedprogs2celltype_mt), function(i_pair) {
      return(sum(perm_mt[, i_pair] >= obs_vec[i_pair]) + 1)
    }))
    succ_mt = matrix(
      succ_vec,
      nrow = nrow(sharedprogs2celltype_mt),
      dimnames = dimnames(sharedprogs2celltype_mt)
    )
    diag(succ_mt) = NA
    # cluster_pairs = paste(rownames(succ_mt)[row(succ_mt)], colnames(succ_mt)[col(succ_mt)], sep = '---')
    succ_df = data.frame(
      row_type = rownames(succ_mt)[row(succ_mt)],
      col_type = colnames(succ_mt)[col(succ_mt)],
      succ_count = succ_mt[1:length(succ_mt)],
      n_trials = n_permtest + 1
    )
    succ_df$EmpiricalP = succ_df$succ_count / succ_df$n_trials
    binomtest_lst = lapply(succ_df$succ_count, function(succ_count){
      if(is.na(succ_count)){
        return(list(conf.int = c(NA, NA)))
      }
      return(binom.test(x = succ_count, n = n_permtest + 1, p = 0.5))
    })
    succ_df$`lb95%CI(p)` = unlist(lapply(binomtest_lst, function(x){x$conf.int[1]}))
    succ_df$`ub95%CI(p)` = unlist(lapply(binomtest_lst, function(x){x$conf.int[2]}))
    
    # pval_vec = unlist(pbapply::pblapply(1:length(sharedprogs2celltype_mt), function(i_pair) {
    #   return((sum(perm_mt[, i_pair] >= obs_vec[i_pair]) + 1) / (nrow(perm_mt) + 1))
    # }))
    # 
    # pval_mt = matrix(
    #   pval_vec,
    #   nrow = nrow(sharedprogs2celltype_mt),
    #   dimnames = dimnames(sharedprogs2celltype_mt)
    # )
    # diag(pval_mt) = NA
    
    return(
      list(
        test_df = test_combn_df,
        correlative_progs_tbl = correlative_progs_tbl, 
        all_bulk_mt = program_mt[, bulksID],
        bulk_mt = bulk_cut_mt,
        sctype_prog_mt = sctype_prog_mt,
        correlative_mt = sharedprogs2celltype_mt,
        pval_df = succ_df, 
        # pval_mt = pval_mt,
        sharedprogs2celltype_lst = sharedprogs2celltype_lst, 
        perm_lst = perm_lst
      )
    )
  }
  
  return(list(test_df = test_combn_df, 
              correlative_progs_tbl = correlative_progs_tbl, 
              all_bulk_mt = program_mt[, bulksID],
              bulk_mt = bulk_cut_mt, 
              sctype_prog_mt = sctype_prog_mt, 
              correlative_mt = sharedprogs2celltype_mt, 
              sharedprogs2celltype_lst = sharedprogs2celltype_lst
              ))
  
}




#' Check and plot the correlation between two programs
#' 
#' @param x_label Name of a program, x
#' @param y_label Name of a program, y
#' @param test_obj A test object
#' @param plt_obj An object for ploting generated by `checkConsensusProgram`
#' @param ... Options pass to Seurat::FeaturePlot
#' @return ggplot objects
#' @export
#' 
checkCorrelative = function(x_label, y_label, test_obj, plt_obj, ...){
  # require(dplyr)
  # require(ggplot2)
  # require(patchwork)
  
  x = test_obj$all_bulk_mt[x_label,]
  y = test_obj$all_bulk_mt[y_label,]
  
  test_res = cor.test(x, y, method = 'spearman')
  rho = test_res$estimate
  p_value = test_res$p.value
  p1 = data.frame(x, y) %>% ggplot(aes(x = x, y = y)) +
    geom_point() +
    geom_smooth(method = 'glm', formula = y ~ x) + 
    theme_classic(base_size = 20) +
    labs(
      x = x_label,
      y = y_label,
      title = sprintf('Spearman cor = %.4f\nP value = %.4f', rho, p_value)
    )
  
  p_x = Seurat::FeaturePlot(
    plt_obj$sce_plt,
    reduction = "tsne",
    features = sprintf('%s', x_label),
    ...
  )
  # p1
  p_y = Seurat::FeaturePlot(
    plt_obj$sce_plt,
    reduction = "tsne",
    features = sprintf('%s', y_label),
    ...
  )
  
  p_tsne = Seurat::DimPlot( # Should be replaced in the future
    plt_obj$sce_plt,
    reduction = "tsne",
    group.by = c('celltype'),
    cols = color_celltype,
    label = T, 
    ...
  )
  
  p_all = p1 + p_tsne + p_x + p_y + patchwork::plot_layout(ncol = 2, nrow = 2) 
  print(p_all)
  invisible(p_all)
}

#' Check and plot the correlation between two programs.ver2
#' 
#' @param x_label Name of a program, x
#' @param y_label Name of a program, y
#' @param test_obj A test object
#' @param plt_obj An object for ploting generated by `checkConsensusProgram`
#' @return ggplot objects
#' @export
#' 
checkCorrelative2 = function(x_label, y_label, test_obj, plt_obj, ...){
  # require(dplyr)
  # require(ggplot2)
  # require(patchwork)
  
  x = test_obj$all_bulk_mt[x_label,]
  y = test_obj$all_bulk_mt[y_label,]
  
  test_res = cor.test(x, y, method = 'spearman')
  rho = test_res$estimate
  p_value = test_res$p.value
  p1 = data.frame(x, y) %>% ggplot(aes(x = x, y = y)) +
    geom_point() +
    geom_smooth(method = 'glm', formula = y ~ x) + 
    theme_classic(base_size = 20) +
    labs(
      x = x_label,
      y = y_label,
      title = sprintf('Spearman cor = %.4f\nP value = %.4f', rho, p_value)
    )
  
  seu_tsne_df = as.data.frame(plt_obj$sce_plt@reductions[['tsne']]@cell.embeddings)
  seu_tsne_df[['celltype']] = plt_obj$sce_plt[['celltype', drop = T]]
  seu_tsne_df[['x_label']] = plt_obj$sce_plt[[x_label, drop = T]]
  seu_tsne_df[['y_label']] = plt_obj$sce_plt[[y_label, drop = T]]
  
  low_bound = min(seu_tsne_df$x_label[seu_tsne_df$x_label > 0])
  p_x = ggplot(seu_tsne_df %>% arrange(x_label), aes(x = tSNE_1, y = tSNE_2, color = x_label)) +
    geom_point() +
    scale_color_gradient2(midpoint = low_bound, high = '#6B6BFF', mid = '#FCFCFF', low = '#B8B8B8')+
    labs(color = x_label) + 
    NULL# ; p_x
  # p1
  low_bound = min(seu_tsne_df$y_label[seu_tsne_df$y_label > 0])
  p_y = ggplot(seu_tsne_df %>% arrange(y_label), aes(x = tSNE_1, y = tSNE_2, color = y_label)) +
    geom_point() +
    scale_color_gradient2(midpoint = low_bound, high = '#6B6BFF', mid = '#FCFCFF', low = '#B8B8B8')+
    labs(color = y_label) + 
    NULL# ; p_x
  
  p_tsne = ggplot(seu_tsne_df, aes(x = tSNE_1, y = tSNE_2, color = celltype)) +
    geom_point() +
    scale_color_manual(values = color_celltype)+
    labs(color = 'Type')
    NULL
  
  p_all = p1 + p_tsne + p_x + p_y + patchwork::plot_layout(ncol = 2, nrow = 2) 
  print(p_all)
  invisible(p_all)
}




#' Plot related genes of a given program
#' 
#' @param seurat_obj A seurat object that holds the expresssion information and tSNE coordinates.
#' @param nmf_obj cNMF object from `identifyConsensusProgram`
#' @param genes The name(s) of genes to plot 
#' @param top_n_genes Top n genes to use, if `genes` is not specified, default = 1000.
#' @param reduction String, the name of the reduction to use, default tsne.
#' @param normalize_by A named group vector. Default NULL. When specified, the values would be normalize by the maximum number according to group.
#' @param programID Name of the program to check, default cProgram 01.
#' @param ... Arguments for geom_point
#' @return A named list of ggplot objects
#' 
#' @importFrom rlang !! :=
#' @export
#' 
plotProgramGenes = function(seurat_obj,
                            nmf_obj,
                            genes = NULL,
                            top_n_genes = 1000,
                            reduction = 'tsne',
                            normalize_by = NULL,
                            programID = 'cProgram 01',
                            raster = T,
                            ...) {
  if (is.null(genes)) {
    genes = names(tail(sort(nmf_obj$consensus_W[, programID]), n = top_n_genes))
  }
  
  if(!requireNamespace("ggrastr", quietly = TRUE)){
    warning('ggrastr not installed. Setting `raster` to FALSE')
    raster = F
  }
  
  chosen_features = intersect(genes, rownames(seurat_obj@assays$RNA@data))
  plt_df = as.data.frame(seurat_obj@reductions[[reduction]]@cell.embeddings) %>%
    mutate(.before = 1,
           cellID = rownames(seurat_obj@reductions[[reduction]]@cell.embeddings))
  
  # Plot correlation with program value
  plt_df = plt_df %>% mutate(!!as.symbol(programID) := nmf_obj$consensus_H[programID,][cellID])
  
  meanScores_vec =
    apply(Matrix::t(seurat_obj@assays$RNA@data[chosen_features,]), 2, function(x) {
      (x - min(x)) / (max(x) - min(x))
    }) %>% rowMeans()
  
  plt_df = plt_df %>% mutate(geneMeanScore = meanScores_vec[cellID])
  
  if (!is.null(normalize_by)) {
    plt_df[['group']] = normalize_by[plt_df[['cellID']]]
    plt_df = plt_df %>%
      group_by(group) %>%
      mutate(norm_score = (geneMeanScore - min(geneMeanScore)) / (max(geneMeanScore) -
                                                                    min(geneMeanScore)))
    
    if(!requireNamespace("ggrastr", quietly = TRUE)){
      warning('ggrastr not installed. Setting `raster` to FALSE')
      raster = F
    }
    if (raster) {
      p_umap = plt_df %>% arrange(norm_score) %>%
        ggplot(aes_string(
          x = colnames(plt_df)[2],
          y = colnames(plt_df)[3],
          color = 'norm_score'
        )) +
        ggrastr::rasterise(geom_point(...), dpi = 300) +
        scale_color_gradient2(
          low = RColorBrewer::brewer.pal(3, 'YlGnBu')[1],
          mid = RColorBrewer::brewer.pal(3, 'YlGnBu')[2],
          midpoint = with(plt_df, (max(norm_score) - min(norm_score)) / 2 + min(norm_score)),
          high = RColorBrewer::brewer.pal(3, 'YlGnBu')[3]
        ) +
        labs(color = 'Mean expression\n(Scaled)',
             title = paste0(
               'cProgram ',
               stringr::str_extract(pattern = '(\\d+)$', programID)
             )) +
        theme(plot.background = element_blank())
      
    } else {
      p_umap = plt_df %>% arrange(norm_score) %>%
        ggplot(aes_string(
          x = colnames(plt_df)[2],
          y = colnames(plt_df)[3],
          color = 'norm_score'
        )) +
        geom_point(...) +
        scale_color_gradient2(
          low = RColorBrewer::brewer.pal(3, 'YlGnBu')[1],
          mid = RColorBrewer::brewer.pal(3, 'YlGnBu')[2],
          midpoint = with(plt_df, (max(norm_score) - min(norm_score)) / 2 + min(norm_score)),
          high = RColorBrewer::brewer.pal(3, 'YlGnBu')[3]
        ) +
        labs(color = 'Mean expression\n(Scaled)',
             title = paste0(
               'cProgram ',
               str_extract(pattern = '(\\d+)$', programID)
             )) +
        theme(plot.background = element_blank())
      
    }
    
    plt_df = plt_df %>% group_by(group) %>% mutate(norm_program = !!as.symbol(programID) / max(!!as.symbol(programID)))
    tmp_vec = plt_df[['norm_program']]
    if (raster) {
      p_umap2 = plt_df %>% arrange(norm_program) %>%
        ggplot(aes_string(
          x = colnames(plt_df)[2],
          y = colnames(plt_df)[3],
          color = 'norm_program'
        )) +
        # geom_point(...) +
        # ggrastr::rasterise(geom_point(...), dpi = 75) +
        ggrastr::rasterise(geom_point(...), dpi = 300) +
        scale_color_gradient2(
          low = RColorBrewer::brewer.pal(3, 'YlGnBu')[1],
          mid = RColorBrewer::brewer.pal(3, 'YlGnBu')[2],
          midpoint = (max(tmp_vec) - min(tmp_vec)) / 2 + min(tmp_vec),
          high = RColorBrewer::brewer.pal(3, 'YlGnBu')[3]
        ) +
        labs(color = 'Program activity\n(Scaled)',
             title = paste0(
               'cProgram ',
               str_extract(pattern = '(\\d+)$', programID)
             )) +
        theme(plot.background = element_blank())
      
    } else {
      p_umap2 = plt_df %>% arrange(norm_program) %>%
        ggplot(aes_string(
          x = colnames(plt_df)[2],
          y = colnames(plt_df)[3],
          color = 'norm_program'
        )) +
        # geom_point(...) +
        # ggrastr::rasterise(geom_point(...), dpi = 75) +
        geom_point(...) +
        scale_color_gradient2(
          low = RColorBrewer::brewer.pal(3, 'YlGnBu')[1],
          mid = RColorBrewer::brewer.pal(3, 'YlGnBu')[2],
          midpoint = (max(tmp_vec) - min(tmp_vec)) / 2 + min(tmp_vec),
          high = RColorBrewer::brewer.pal(3, 'YlGnBu')[3]
        ) +
        labs(color = 'Program activity\n(Scaled)',
             title = paste0(
               'cProgram ',
               str_extract(pattern = '(\\d+)$', programID)
             )) +
        theme(plot.background = element_blank())
      
    }
    
  } else {
    if (raster) {
      p_umap = plt_df %>% arrange(geneMeanScore) %>%
        ggplot(aes_string(
          x = colnames(plt_df)[2],
          y = colnames(plt_df)[3],
          color = 'geneMeanScore'
        )) +
        # geom_point(...) +
        # ggrastr::rasterise(geom_point(...), dpi = 75) +
        ggrastr::rasterise(geom_point(...), dpi = 300) +
        scale_color_gradient2(
          low = RColorBrewer::brewer.pal(3, 'YlGnBu')[1],
          mid = RColorBrewer::brewer.pal(3, 'YlGnBu')[2],
          midpoint = with(plt_df, (
            max(geneMeanScore) - min(geneMeanScore)
          ) / 2 + min(geneMeanScore)),
          high = RColorBrewer::brewer.pal(3, 'YlGnBu')[3]
        ) +
        labs(color = 'Mean expression\n(Scaled)',
             title = paste0(
               'cProgram ',
               str_extract(pattern = '(\\d+)$', programID)
             )) +
        theme(plot.background = element_blank())
      
    } else {
      p_umap = plt_df %>% arrange(geneMeanScore) %>%
        ggplot(aes_string(
          x = colnames(plt_df)[2],
          y = colnames(plt_df)[3],
          color = 'geneMeanScore'
        )) +
        # geom_point(...) +
        # ggrastr::rasterise(geom_point(...), dpi = 75) +
        geom_point(...) +
        scale_color_gradient2(
          low = RColorBrewer::brewer.pal(3, 'YlGnBu')[1],
          mid = RColorBrewer::brewer.pal(3, 'YlGnBu')[2],
          midpoint = with(plt_df, (
            max(geneMeanScore) - min(geneMeanScore)
          ) / 2 + min(geneMeanScore)),
          high = RColorBrewer::brewer.pal(3, 'YlGnBu')[3]
        ) +
        labs(color = 'Mean expression\n(Scaled)',
             title = paste0(
               'cProgram ',
               str_extract(pattern = '(\\d+)$', programID)
             )) +
        theme(plot.background = element_blank())
    }
    
    tmp_vec = plt_df[[programID]]
    if (raster) {
      p_umap2 = plt_df %>% arrange({
        {
          programID
        }
      }) %>%
        ggplot(aes_string(
          x = colnames(plt_df)[2],
          y = colnames(plt_df)[3],
          color = programID
        )) +
        # geom_point(...) +
        # ggrastr::rasterise(geom_point(...), dpi = 75) +
        ggrastr::rasterise(geom_point(...), dpi = 300) +
        scale_color_gradient2(
          low = RColorBrewer::brewer.pal(3, 'YlGnBu')[1],
          mid = RColorBrewer::brewer.pal(3, 'YlGnBu')[2],
          midpoint = (max(tmp_vec) - min(tmp_vec)) / 2 + min(tmp_vec),
          high = RColorBrewer::brewer.pal(3, 'YlGnBu')[3]
        ) +
        labs(color = 'Program activity\n(Scaled)',
             title = paste0(
               'cProgram ',
               str_extract(pattern = '(\\d+)$', programID)
             )) +
        theme(plot.background = element_blank())
      
    } else {
      p_umap2 = plt_df %>% arrange({
        {
          programID
        }
      }) %>%
        ggplot(aes_string(
          x = colnames(plt_df)[2],
          y = colnames(plt_df)[3],
          color = programID
        )) +
        geom_point(...) +
        scale_color_gradient2(
          low = RColorBrewer::brewer.pal(3, 'YlGnBu')[1],
          mid = RColorBrewer::brewer.pal(3, 'YlGnBu')[2],
          midpoint = (max(tmp_vec) - min(tmp_vec)) / 2 + min(tmp_vec),
          high = RColorBrewer::brewer.pal(3, 'YlGnBu')[3]
        ) +
        labs(color = 'Program activity\n(Scaled)',
             title = paste0(
               'cProgram ',
               str_extract(pattern = '(\\d+)$', programID)
             )) +
        theme(plot.background = element_blank())
    }
    
  }
  
  if (!is.null(normalize_by)) {
    cor_res = cor.test(plt_df[[programID]], plt_df$norm_score, method = 'pearson')
    p_cor = plt_df %>%
      ggplot(aes_string(x = sprintf("`%s`", programID), y = "norm_score")) +
      geom_point(aes(color = group)) +
      geom_smooth(aes(color = group),
                  method = 'lm',
                  formula = y ~ x) +
      labs(
        title =
          ifelse(
            cor_res$p.value == 0,
            sprintf('Cor = %s\nP.val < 2.20e-16', round(cor_res$estimate, 4)),
            sprintf(
              'Cor = %s\nP.val = %s',
              round(cor_res$estimate, 4),
              formatC(cor_res$p.value, digits = 2, format = 'e')
            )
          ),
        x = paste0(
          'Program activity (cProgram ',
          str_extract(pattern = '(\\d+)$', programID),
          ')'
        ),
        y = 'Mean expression(Scaled)'
      ) +
      theme(plot.background = element_blank())
  } else {
    cor_res = cor.test(plt_df[[programID]], plt_df$geneMeanScore, method = 'pearson')
    p_cor = plt_df %>%
      ggplot(aes_string(x = sprintf("`%s`", programID), y = "geneMeanScore")) +
      geom_point() +
      geom_smooth(method = 'lm', formula = y ~ x) +
      labs(
        title =
          ifelse(
            cor_res$p.value == 0,
            sprintf('Cor = %s\nP.val < 2.20e-16', round(cor_res$estimate, 4)),
            sprintf(
              'Cor = %s\nP.val = %s',
              round(cor_res$estimate, 4),
              formatC(cor_res$p.value, digits = 2, format = 'e')
            )
          ),
        x = paste0(
          'Program activity (cProgram ',
          str_extract(pattern = '(\\d+)$', programID),
          ')'
        ),
        y = 'Mean expression'
      ) +
      theme(plot.background = element_blank())
    
  }
  
  p_umapall = patchwork::wrap_plots(p_umap, p_umap2, p_cor, design = "AB\nCC") +
    patchwork::plot_annotation(title = paste0(
      'Correlation between gene expression and cProgram ',
      str_extract(pattern = '(\\d+)$', programID)
    ))
  
  # if (plot_each) {
  #   return(list(
  #     genes = genes,
  #     umaps = plt_lst,
  #     all = p_umapall
  #   ))
  # }
  return(
    list(
      genes = genes,
      p_umap = p_umap,
      p_umap2 = p_umap2,
      p_cor = p_cor,
      all = p_umapall,
      plt_df = plt_df
    )
  )
}


#' Extract and scale the data for ploting
#' 
#' @param seurat_obj A seurat object that holds the expresssion information and tSNE coordinates.
#' @param nmf_obj cNMF object from `identifyConsensusProgram`
#' @param genes The name(s) of genes to plot 
#' @param cells The name(s) of cells to plot 
#' @param top_n_genes Top n genes to use, if `genes` is not specified, default = 1000.
#' @param reduction String, the name of the reduction to use, default tsne.
#' @param programID Name of the program to check, default cProgram 01.
#' @param ... Arguments for geom_point
#' @return A tbl/df
#' 
#' @importFrom rlang !! :=
#' @export
#' 
extractProgramData = function(seurat_obj,
                            nmf_obj,
                            genes = NULL,
                            cells = NULL,
                            top_n_genes = 1000,
                            reduction = 'tsne',
                            programID = 'cProgram 01',
                            ...) {
  if (is.null(genes)) {
    genes = names(tail(sort(nmf_obj$consensus_W[, programID]), n = top_n_genes))
  }
  
  chosen_features = intersect(genes, rownames(seurat_obj@assays$RNA@data))
  plt_df = as.data.frame(seurat_obj@reductions[[reduction]]@cell.embeddings) %>%
    mutate(.before = 1,
           cellID = rownames(seurat_obj@reductions[[reduction]]@cell.embeddings))
  
  if(!is.null(cells)){
    plt_df = plt_df %>% dplyr::filter(cellID %in% cells)
  }
  
  # Plot correlation with program value
  origin_score = nmf_obj$consensus_H[programID,][plt_df$cellID]
  scaled_program_score = (origin_score - min(origin_score))/(max(origin_score) - min(origin_score))
  names(scaled_program_score) = names(origin_score)
  plt_df = plt_df %>% 
    mutate(target_program = programID, 
           program_score = scaled_program_score)
  
  meanScores_vec =
    apply(Matrix::t(seurat_obj@assays$RNA@data[chosen_features, plt_df$cellID]), 2, function(x) {
      (x - min(x)) / (max(x) - min(x))
    }) %>% rowMeans()
  
  plt_df = plt_df %>% mutate(topGeneMeanScore = meanScores_vec[cellID]) %>% arrange(topGeneMeanScore)
  
  
  return(plt_df)
}


#' Plot the enrichment for celltypes DEGs using the related genes of a given program
#' 
#' @param nmf_obj cNMF object from `identifyConsensusProgram`
#' @param programID Name of the program to check, default cProgram 01.
#' @param TERM2GENE A data.frame pass to `clusterProfiler::GSEA`
#' @param all_types A string of all the cell types.
#' @param genes A ranked and named genes vector. Default NULL.
#' @param top_n_genes Top n genes for the enricher, if `genes` is not specified. Default = 1000.
#' @param maxGSSize Default Inf
#' @param verbose Logical, whether to print log.
#' @param GSEA_args Other options that passed to `clusterProfiler::GSEA`, default empty.
#' @param enricher_args Other options that passed to `clusterProfiler::enricher`, default empty.
#' @param gseaplot_args Other options that passed to `clusterProfiler::gseaplot`, default empty.
#' @return A named list of ggplot objects
#' 
#' @importFrom clusterProfiler GSEA enricher
#' @export
#' 
plotEnrichCelltypes = function(
  nmf_obj,
  programID = 'cProgram 01',
  TERM2GENE, 
  all_types,
  genes = NULL,
  top_n_genes = 1000, 
  max_frac = NULL,
  maxGSSize = Inf, 
  verbose = TRUE, 
  GSEA_args = list(), 
  enricher_args = list(), 
  gseaplot_args = list()
  ) {
  
  all_genes = rownames(nmf_obj$consensus_W)
  if (is.null(genes)) {
    genes = sort(nmf_obj$consensus_W[, programID], decreasing = T)
  }
  
  
  gsea_res = do.call(what = clusterProfiler::GSEA, args = c(list(
    genes,
    TERM2GENE = TERM2GENE,
    scoreType = 'pos',
    verbose = verbose,
    maxGSSize = Inf
  ), GSEA_args))
  gsea_df = as.data.frame(gsea_res)
  
  gsea_plt_lst = list()
  for(i_type in gsea_df$ID) {
    gsea_plt_lst[[i_type]] = do.call(what = clusterProfiler::gseaplot,
                                     args = c(list(
                                       gsea_res,
                                       geneSetID = i_type,
                                       title = sprintf('%s: %s\nQval = %s',
                                                       programID,
                                                       i_type,
                                                       formatC(gsea_df[which(gsea_df$ID == i_type), 'qvalues'], digits = 2)) 
                                     ),
                                     gseaplot_args))
  }
  
  if(!is.null(top_n_genes)) {
    enr_res = do.call(what = clusterProfiler::enricher,
                      args = c(
                        list(
                          names(head(genes, top_n_genes)),
                          TERM2GENE = TERM2GENE,
                          universe = all_genes,
                          maxGSSize = Inf
                        ),
                        enricher_args
                      ))
  } else if (!is.null(max_frac)){
    frac_cutoff = max(genes) * max_frac
    enr_res = do.call(what = clusterProfiler::enricher,
                      args = c(
                        list(
                          names(genes)[genes > frac_cutoff],
                          TERM2GENE = TERM2GENE,
                          universe = all_genes,
                          maxGSSize = Inf
                        ),
                        enricher_args
                      ))
  } else {
    stop('Neither max_frac or top_n_genes is specified.')
  }
  
  enr_df = as.data.frame(enr_res)
  # gsea_res %>% as.data.frame() %>% View()
  
  return(list(
    gsea_obj = gsea_res, 
    enrich_obj = enr_res,
    gsea_plots = gsea_plt_lst
  ))
}

#' Get the lower triangle of a symmetric matrix
#'
#' @param x_mt A symmetric matrix
#' @param keep_diag Logical. Keep diagonal values or not.
#' @return A data.frame
#' @export 
getLowerTri2df = function(x_mt, keep_diag = T){
  stopifnot(all(rownames(x_mt) == colnames(x_mt)))
  
  clusters = rownames(x_mt)
  row_tags = clusters[row(x_mt)]
  col_tags = clusters[col(x_mt)]
  low_idx = lower.tri(x_mt, diag = keep_diag)
  
  out_df = data.frame(row_tag = row_tags[low_idx], col_tag = col_tags[low_idx], 
                      value = x_mt[low_idx])
  return(out_df)
}

