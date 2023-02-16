#' @import utils
#' @import stats
#' @import ggplot2
#' @import patchwork
#' @import dplyr 
#' @import tictoc
NULL

# Context-dependent ====

#' Encode context
#'
#' This function remove NA, and assign numbers to each 
#' context, in order to perform multi-variable linear regression
#'
#' @param annotation named vector, indicate the context of each bulk sample
#' @param context_num_vec named numeric vector, indicate the number assign to each context, to perform multi-variable linear regression
#' @export
#' @return a tibble
#' 
encodeContext <- function(annotation, 
                          context_num = NULL) {
  
  if(is.null(context_num)){
    categories = unique(annotation)[!is.na(unique(annotation))]
    context_num = setNames(nm = categories,
                           as.numeric(factor(unique(MSI_annotation)[!is.na(unique(MSI_annotation))])))
    
  }
  # remove NA and useless context annotation
  annotation <- annotation[!is.na(annotation)]
  annotation <-
    annotation[which(annotation %in% names(context_num))]
  
  annotation_tbl <-
    cbind(names(annotation), annotation) %>% as_tibble()
  colnames(annotation_tbl) <- c("label", "context")
  
  
  # assign values to context
  context_num_vec <- c()
  
  for (i in c(1:dim(annotation_tbl)[1])) {
    context_num_vec <- c(context_num_vec,
                         context_num[annotation_tbl[[i, 2]]])
  }
  
  assertthat::assert_that(sum(names(context_num_vec) == annotation_tbl[, 2]) == length(context_num_vec))
  
  annotation_tbl <- cbind(annotation_tbl,
                          "num" = context_num_vec)
  
  return(annotation_tbl)
  
}


#' Detect different co-varying pattern in contexts
#' This function performs multi-variable linear regression on context and programs, then select
#' program pairs that are significantly influenced by context
#'
#' @param bulk_pg_dH_mt matrix, with columns as programs, rows as bulk samples, calculated as `d` %*% `H`, where `d` and `H` are from `cNMF`
#' @param context_num_tbl tibble, with colnames `label`, `context`, and `num`, indicating the sample name, context, and context number respectively
#' @param signif_cutoff Cutoff to determine cProgram pairs with significant context-factor
#'
#' @return A list contains the result of multi-linear regreassion
#' @export
detectDiffContext <- function(bulk_pg_dH_mt,
                     context_num_tbl, 
                     signif_cutoff = 0.05) {
  
  # perform regression ----
  num_programs <- dim(bulk_pg_dH_mt)[2]
  num_comparison <- choose(num_programs, 2)
  
  # combine with bulk_program_matrix
  bulk_pg_dH_mt <- as.data.frame(bulk_pg_dH_mt) %>% rownames_to_column(var = "label")
  
  combined_mt <- dplyr::left_join(context_num_tbl,
                                  bulk_pg_dH_mt, by = 'label')
  context_num <- context_num_tbl[,3]
  
  glm_res <- list()
  # glm_res_summary <- list()
  prog_names = grep(pattern = '\\d+', colnames(combined_mt), value = T)
  combn_mt1 = t(combn(prog_names, 2))
  combn_mt2 = combn_mt1[, rev(1:ncol(combn_mt1))]
  combn_mt = rbind(combn_mt1, combn_mt2)
  test_df = as.data.frame(combn_mt, stringsAsFactors = F)
  colnames(test_df) = c('Prog1', 'Prog2')
  test_df[['b1_estimate']] = NA
  test_df[['b2_estimate']] = NA
  test_df[['b1_pval']] = NA
  test_df[['b2_pval']] = NA
  for(i_row in 1:nrow(test_df)){
    pg_1 = test_df[i_row, 1]
    pg_2 = test_df[i_row, 2]
      
    program_1 = combined_mt[, pg_1, drop = T]
    program_2 = combined_mt[, pg_2, drop = T]
    
    pg_1_str = pg_1
    pg_2_str = pg_2
    
    pg_pair = sprintf('%s:%s', pg_1_str, pg_2_str)
    
    glm_res[[pg_pair]] <- glm(program_1 ~ program_2 + context_num)
    glm_summary <- summary(glm_res[[pg_pair]])
    
    test_df[i_row, 'b1_estimate'] = glm_summary[['coefficients']][2,1]
    test_df[i_row, 'b2_estimate'] = glm_summary[['coefficients']][3,1]
    test_df[i_row, 'b1_pval'] = glm_summary[['coefficients']][2,4]
    test_df[i_row, 'b2_pval'] = glm_summary[['coefficients']][3,4]
  }
  
  flt_test_df = test_df %>% dplyr::filter(b1_estimate > 0, b1_pval < signif_cutoff, b2_pval < signif_cutoff)
  
  return(list(
    all_test_res = test_df, 
    flt_test_res = flt_test_df, 
    model = glm_res
    ))
  
}

#' Visualize different co-varying pattern in contexts
#' This function visualize the results of program pairs that are influenced significantly by the context, for each program pair, it includes the model and the summary
#'
#' @param bulk_pg_dH_mt a matrix, with columns as programs, rows as bulk samples, calculated as `d` %*% `H`, where `d` and `H` are from `cNMF`
#' @param test_df a data.frame, with columns `Prog1` and `Prog2`
#' @param context_num_tbl a tibble, result of `preprocess`
#' @param cor_method A correlation method passed to `cor.test`
#' @param both Default = TRUE, return both signifcant.
#' @param max_row Maximum number of rows in test_df to plot.
#' @param plt_dir a string specify the ploting directory, default NULL
#'
#' @export
#' @return A list of plots
#' 
#' 
visualizeDiffContext <- function(bulk_pg_dH_mt,
                        test_df = NULL,
                        context_num_tbl,
                        cor_method = 'spearman',
                        both = T,
                        max_row = 30, 
                        plt_dir = NULL) {
  
  # test_df = detection_obj$flt_test_res
  context_vec = with(context_num_tbl, {setNames(nm = label, context)})
  context_vec = context_vec[rownames(bulk_pg_dH_mt)]
  if(is.null(test_df)){
    message('Test table not provided. Calculating de novo.')
    
    mt2test = t(combn(colnames(bulk_pg_dH_mt), 2))
    tmp_test_lst = list()
    k = 1
    for(i_row in 1:nrow(mt2test)){
      prog1 = mt2test[i_row, 1]
      prog2 = mt2test[i_row, 2]
      
      for(i_context in unique(context_vec)){
        x_vec = bulk_pg_dH_mt[context_vec == i_context, prog1] 
        y_vec = bulk_pg_dH_mt[context_vec == i_context, prog2]
        tmp_test_res = cor.test(x_vec, y_vec, method = 'spearman')
        tmp_test_lst[[k]] = tibble(
          prog1 = prog1, prog2 = prog2, 
          context = i_context, rho = tmp_test_res$estimate, pval = tmp_test_res$p.value)
        k = k+1
      }
    }
    ttl_test_tbl = bind_rows(tmp_test_lst)
    if(both){
      test_df = ttl_test_tbl %>% dplyr::group_by(prog1, prog2) %>% 
        dplyr::summarise(anyPos = sum(rho > 0), anySignif = sum(pval < 0.05), cor_diff = max(rho) - min(rho)) %>% 
        dplyr::filter(anySignif > 1, anyPos > 0) %>% arrange(desc(cor_diff))
    } else {
      test_df = ttl_test_tbl %>% dplyr::group_by(prog1, prog2) %>% 
        dplyr::summarise(anyPos = sum(rho > 0), anySignif = sum(pval < 0.05), cor_diff = max(rho) - min(rho)) %>% 
        dplyr::filter(anySignif > 0, anyPos > 0) %>% arrange(desc(cor_diff))
    }
    test_df = test_df %>% dplyr::arrange()
  } else {
    stopifnot(all(test_df$Prog1 %in% colnames(bulk_pg_dH_mt)))
    stopifnot(all(test_df$Prog2 %in% colnames(bulk_pg_dH_mt)))
  }
  
  test_df = test_df %>% head(max_row)
  
  
  bulk_pg_dH_mt <- as.data.frame(bulk_pg_dH_mt) %>% rownames_to_column(var = "label")
  
  combined_mt <- dplyr::left_join(context_num_tbl,
                                  bulk_pg_dH_mt, by = 'label')
  
  # pg_pairs <- names(detection_obj)
  
  wrap <- list()
  
  for (i_row in 1:nrow(test_df)) {
    
    pg_1 = test_df[i_row, 1, drop = T]
    pg_2 = test_df[i_row, 2, drop = T]
    pg_pair = sprintf('%s:%s', pg_1, pg_2)
    
    program_1 <- combined_mt[,pg_1] %>% unlist()
    program_2 <- combined_mt[,pg_2] %>% unlist()
    
    res <- data.frame(program_1, program_2, context_vec)
    
    # Calculate cor & p val
    title_vec = c()
    for(cat2test in unique(context_vec)){
      tmp_cor = cor.test(res$program_1[res$context_vec == cat2test], res$program_2[res$context_vec == cat2test], method = cor_method)
      tmp_title_str = sprintf('%s: Cor = %.2f, P.val = %.4f', cat2test, tmp_cor$estimate, tmp_cor$p.value)
      title_vec = c(title_vec, tmp_title_str)
    }
    
    prog <- res %>%
      ggplot2::ggplot(aes(x = program_1, y = program_2, color = context_vec)) +
      geom_point(alpha = 0.7) +
      geom_smooth(aes(color = context_vec), method = "lm", formula = y ~ x) +
      labs(
        x = pg_1,
        y = pg_2, 
        title = paste0(title_vec, collapse = '\n'), 
        color = 'Context'
      )
    
    wrap[[pg_pair]] <- prog
    
    if(!is.null(plt_dir)){
      dir.create(plt_dir, showWarnings = F, mode = '0755')
      ggsave(paste0(plt_dir, "/",pg_1, "_", pg_2, ".pdf"),
             width = 8,
             height = 6,
             device = "pdf",
             limitsize = FALSE)
      
    }
    
  }
  
  return(wrap)
  
}

# #' Calculate correlative pattern consensus program values in a specific context
# #' 
# #' @param program_mt Matrix, processed
# #' @param sources Named strings, c\('sample_name' = 'bulk|singlecell'\)
# #' @param celltypes Named vector, should be a cell-id-named vector of cell types
# #' @param nonzero_ratio Value between 0-1, c\('sample_name' = 'bulk|singlecell'\)
# #' @param rho_cutoff The remaining program pairs must have spearman correlation rho larger than this. Default = 0.
# #' @param qval_cutoff The remaining program pairs must have qvalue smaller than this. Default = 0.05. BH method.
# #' @param n_permtest Number of permutation rounds to perform, default = 2000
# #' @param n_cores Number of cores to use for permutation, default = the maximum available cores - 1
# #' 
# #' @return A named list with the correlative information between celltypes.
# #' @export
# calcCorrelative4 <- function(program_mt,
#                              sources,
#                              celltypes,
#                              contexts = NULL,
#                              progs,
#                              .context = NULL,
#                              nonzero_ratio = 0.1,
#                              rho_cutoff = 0,
#                              n_permtest = 2000,
#                              n_cores = parallel::detectCores() - 1,
#                              qval_cutoff = 0.05) {
#   
#   
#   if (.context != "all") {
#     if (is.null(context_mt)) {
#       stop("The matrix to describe the context of bulk should be supplied")
#     }
#     
#     # drop other context bulks ----
#     all_bulks <- annotation_vec[annotation_vec == "bulk"]
#     
#     ## find index of other context in the merged matrix
#     all_colnames <- program_mt %>% colnames()
#     context_bulk_names <- context_mt %>% dplyr::filter(!is.na(context), context == .context) %>% .[, 1] %>% unlist()
#     
#     other_bulk_names <- all_bulks[which(!names(all_bulks) %in% context_bulk_names)] %>% names()
#     
#     context_bulk_names_idx <- which(all_colnames %in% context_bulk_names)
#     other_bulk_names_idx <- which(all_colnames %in% other_bulk_names)
#     
#     sig_pg <- c()
#     
#     for (pg_pair in progs) {
#       pgs <- pg_pair %>% strsplit(split = "_") %>% .[[1]]
#       sig_pg <- c(sig_pg, pgs)
#     }
#     
#     sig_pg <- sig_pg %>% unique() %>% sort()
#     
#     sig_pg_str <- paste0("MergedProgram_", sig_pg)
#     
#     ## drop other context bulks
#     sig_program_mt <- program_mt[sig_pg_str, -other_bulk_names_idx]
#     
#     # select source ----
#     sources <- sources[(sig_program_mt %>% colnames())]
#     # checkmate::assert(sum(colnames(program_mt) == names(sources)) == length(sources))
#   }
#   
#   return(calcCorrelative3(sig_program_mt,
#                           sources,
#                           celltypes,
#                           nonzero_ratio,
#                           rho_cutoff,
#                           n_permtest,
#                           n_cores,
#                           qval_cutoff))
#   
# }

#' Identify correlative pattern.ver4
#' 
#' @param program_mt Matrix, processed
#' @param sources Named strings, c\('sample_name' = 'bulk|singlecell'\)
#' @param celltypes Named vector, should be a cell-id-named vector of cell types
#' @param nonzero_ratio Value between 0-1, c\('sample_name' = 'bulk|singlecell'\)
#' @param rho_cutoff The remaining program pairs must have spearman correlation rho larger than this. Default = 0.
#' @param qval_cutoff The remaining program pairs must have qvalue smaller than this. Default = 0.05. BH method.
#' @param n_permtest Number of permutation rounds to perform, default = 2000
#' @param n_cores Number of cores to use for permutation, default = the maximum available cores - 1
#' @param contexts Named vector, specify the context information of each bulk sample. If not specified, all bulk samples are considered to be from the same conditions.
#' @param progs Character vector, names of the list that derived from `detectDiffContext`, indicating program pairs that are significantly influenced by the context. Must be provided when running in context mode.
#' @param .context Character, default NULL. Specify one of the contexts to run.
#' 
#' @return A named list with the correlative information between celltypes.
#' @export
#' 
calcCorrelative4 = function(program_mt,
                            sources,
                            celltypes,
                            nonzero_ratio = 0.1,
                            rho_cutoff = 0,
                            n_permtest = 2000, 
                            n_cores = parallel::detectCores() - 1,
                            qval_cutoff = 0.05, 
                            contexts = NULL,
                            progs = NULL,
                            # prog_pairs = NULL,
                            .context = NULL, 
                            verbose = T) {
  
  if(!is.null(contexts)){
    contexts_cat = unique(contexts)
    if(verbose) message(sprintf('%d conditions found.', length(contexts_cat)))
    
    if(is.null(progs)){
      warning('Progs not specified. Running with all programs.')
      progs = rownames(program_mt)
    }
    
    if(is.null(.context)){
      return_lst = list()
      for(.context in contexts_cat){
        if(verbose) message(sprintf('Now running %s', .context))
        
        # Drop bulks from other contexts
        bulksID = names(sources)[sources == 'bulk']
        bulksID = bulksID[bulksID %in% names(contexts)]
        other_bulks_idx = which(colnames(program_mt) %in% bulksID[contexts != .context])
        i_program_mt = program_mt[progs,-other_bulks_idx]
        other_bulks_idx = which(names(sources) %in% bulksID[contexts != .context])
        i_sources = sources[-other_bulks_idx]
        
        res_obj = calcCorrelative4(
          program_mt = i_program_mt,
          sources = i_sources,
          celltypes = celltypes,
          nonzero_ratio = nonzero_ratio,
          rho_cutoff = rho_cutoff,
          n_permtest = n_permtest,
          n_cores = n_cores,
          qval_cutoff = qval_cutoff, 
          verbose = verbose
        )
        return_lst[[.context]] = res_obj
      }
      return(return_lst)
    } else {
      
      # Drop bulks from other contexts
      bulksID = names(sources)[sources == 'bulk']
      bulksID = bulksID[bulksID %in% names(contexts)]
      other_bulks_idx = which(colnames(program_mt) %in% bulksID[contexts != .context])
      i_program_mt = program_mt[progs, -other_bulks_idx]
      other_bulks_idx = which(names(sources) %in% bulksID[contexts != .context])
      i_sources = sources[-other_bulks_idx]
      
      res_obj = calcCorrelative4(
        program_mt = i_program_mt, 
        sources = i_sources, 
        celltypes = celltypes,
        nonzero_ratio = nonzero_ratio, 
        rho_cutoff = rho_cutoff, 
        n_permtest = n_permtest, 
        n_cores = n_cores, 
        qval_cutoff = qval_cutoff, 
        verbose = verbose
      )
      return_lst = list()
      return_lst[[.context]] = res_obj
      return(return_lst)
    }
    
  } 
  # recurse to normal mode
  
  bulksID = names(sources)[sources == 'bulk']
  bulk_cut_mt = program_mt[, bulksID]
  
  progs2run_idx = rowSums(bulk_cut_mt > 0) > round(nonzero_ratio * length(bulksID))
  progs2run = names(progs2run_idx)[progs2run_idx]
  
  bulk_cut_mt = program_mt[progs2run, bulksID]
  
  scsID = names(sources)[sources != 'bulk']
  sc_cut_mt = program_mt[, scsID]
  
  # 
  test_combn_df = data.frame(t(combn(rownames(bulk_cut_mt), 2)), rho = Inf, pval = Inf)
  colnames(test_combn_df)[1:2] = c('Prog1', 'Prog2')
  
  if(verbose) message('Identifying correlative programs')
  
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
  
  # Mapping correlative programs to cell types
  
  sctype_prog_mt = sc_cut_mt %>% 
    apply(1, function(x){tapply(x, celltypes[scsID], mean)}) %>% 
    apply(1, function(x){x/(sum(x) + 1e-32)}) %>% 
    apply(1, function(x){x/(max(x) + 1e-32)}) %>% 
    t()
  
  correlative_progs_tbl = test_combn_df %>% dplyr::filter(rho > rho_cutoff, qval < qval_cutoff)
  
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
  
  if(verbose) message('Mapping correlative programs to cell types...')
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
    if(verbose) message('Begin permutation')
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
    if(verbose) message('Calculate emperical p values.')
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


#' Check and plot the correlation between two programs within specific context
#' 
#' @param x_label Name of a program, x
#' @param y_label Name of a program, y
#' @param test_obj_context A test object with multiple contexts
#' @param plt_obj An object for ploting generated by `checkConsensusProgram`
#' @param ... Options pass to Seurat::FeaturePlot
#' @return ggplot objects
#' @export
#' 
checkCorrelative_context = function(x_label, y_label, test_obj_context, plt_obj, ...){
  # require(dplyr)
  # require(ggplot2)
  # require(patchwork)
  
  reduction_df = plt_obj$sce_plt@reductions$tsne@cell.embeddings %>% as.data.frame()
  all_bulks = plt_obj$sce_plt$name[plt_obj$sce_plt$source == 'bulk']
  reduction_df = reduction_df[all_bulks, ]
  
  args = list(...)
  args$size = ifelse(is.null(args$pt.size), 1, args$pt.size)
  
  plt_lst = list()
  for(i_context in names(test_obj_context)){
    x = test_obj_context[[i_context]]$all_bulk_mt[x_label, ]
    y = test_obj_context[[i_context]]$all_bulk_mt[y_label, ]
    
    test_res = cor.test(x, y, method = 'spearman')
    rho = test_res$estimate
    p_value = test_res$p.value
    p1 = data.frame(x, y) %>% ggplot(aes(x = x, y = y)) +
      geom_point(size= args$size) +
      geom_smooth(method = 'glm', formula = y ~ x) +
      theme_classic(base_size = 20) +
      labs(
        x = x_label,
        y = y_label,
        title = sprintf('Cor = %.4f\nP value = %.4f', rho, p_value)
      )
    
    # bulks2exclude = 
    
    reduction_df$label = rownames(reduction_df) %in% colnames(test_obj_context[[i_context]]$all_bulk_mt)
    reduction_df$value = NA
    reduction_df$value[reduction_df$label] = plt_obj$sce_plt[[x_label, drop = T]][rownames(reduction_df)[reduction_df$label]]
    p_x = ggplot(reduction_df[reduction_df$label, ]) + 
      geom_point(aes(x = tSNE_1, y = tSNE_2, color = rank(value)), size = args$size) +
      labs(color = sprintf('Rank(%s)', x_label)) + 
      scale_color_gradient(low = "#1946FA", high = "#FF2929", na.value = "#CCCCCC") + 
      theme(axis.text = element_blank(), 
            axis.ticks = element_blank(), 
            axis.line = element_blank(), 
            axis.title = element_blank())
    
    reduction_df$value = NA
    reduction_df$value[reduction_df$label] = plt_obj$sce_plt[[y_label, drop = T]][rownames(reduction_df)[reduction_df$label]]
    p_y = ggplot(reduction_df[reduction_df$label, ]) + 
      geom_point(aes(x = tSNE_1, y = tSNE_2, color = rank(value)), size = args$size) +
      labs(color = sprintf('Rank(%s)', y_label)) + 
      scale_color_gradient(low = "#1946FA", high = "#FF2929", na.value = "#CCCCCC") + 
      theme(axis.text = element_blank(), 
            axis.ticks = element_blank(), 
            axis.line = element_blank(), 
            axis.title = element_blank())
    
    plt_lst[[i_context]] = p1+p_x+p_y + patchwork::plot_layout(ncol = 3) +
      patchwork::plot_annotation(title = i_context, theme = theme_classic(base_size = 20))
  }
  
  
  p_x = Seurat::FeaturePlot(plt_obj$sce_plt,
                            reduction = "tsne",
                            features = sprintf('%s', x_label),
                            ...)
  # p1
  p_y = Seurat::FeaturePlot(plt_obj$sce_plt,
                            reduction = "tsne",
                            features = sprintf('%s', y_label),
                            ...)
  
  p_tsne = Seurat::DimPlot(
    # Should be replaced in the future
    plt_obj$sce_plt,
    reduction = "tsne",
    group.by = c('celltype'),
    cols = color_celltype,
    label = T,
    ...
  )
  
  p_all = p_tsne + p_x + p_y + patchwork::plot_layout(ncol = 3, nrow = 1)
  plt_lst[['all']] = p_all
  
  # print(p_all)
  # invisible(p_all)
  return(plt_lst)
}

