#### Prune the candidate network by elastic-net regression ####
#
# Step 5, the final step, of the GRN reconstruction pipeline. For every
# candidate edge target in the base network (3_preparing_base_network.R),
# fits an elastic-net regression of that target's expression or accessibility
# against all of its candidate regulators (predictor weights set from
# TOBIAS footprinting score for TF-to-peak edges, or from a power-law
# distance decay for enhancer-to-promoter edges), separately for the target's
# aCRE and tCRE measurement where both exist. Regulators with a non-zero
# fitted coefficient become edges in the final network; every retained edge
# additionally carries:
#   - a permutation-test p-value (`randomize = TRUE`): the target is
#     reshuffled across metacells `num_perms` times, the same model refit
#     each time, and the edge's p-value is the fraction of permutations whose
#     coefficient magnitude matches or exceeds the real one (this is what
#     downstream scripts filter on as `adjpvalue`, after BH correction);
#   - `dev_ratio`, the model's deviance ratio (fraction of null deviance
#     explained), used downstream as a per-edge model-fit-quality filter;
#   - `deltaR2` and an OLS p-value, from a post-hoc ordinary-least-squares
#     refit on only the lasso-selected predictors — this gives each retained
#     predictor's unique contribution to R^2, which the elastic-net fit alone
#     does not provide.
#
# The `use_hic` / `use_footprinting` control-panel flags must match the
# base-network and dataset variant read; `use_enhancer_tCRE = FALSE`
# additionally excludes tCRE as a predictor/target for enhancers, producing
# the "no_enhancer_tCRE" variant used to ask what modelling enhancer
# transcriptional activity contributes beyond accessibility alone (see
# Data_and_code/tCRE_topology/). This script is by far the most
# computationally expensive step of the pipeline — with `randomize = TRUE`
# and the full metacell set, expect the two model-fitting loops below (one
# per edge class) together to need on the order of a day of wall time even
# with substantial parallelism.
#
# Reads (see ../config.R to set raw_data_folder and primary_data_folder):
#   [primary_data_folder]/3_preparing_base_network/base_network(_no_hic)
#     (_no_footprinting).rds        (from 3_preparing_base_network.R)
#   [primary_data_folder]/4_preparing_dataset/dataset(_no_hic)
#     (_no_footprinting).rds        (from 4_preparing_dataset.R)
#   [raw_data_folder]/pseudotime.scRNA.tsv
#
# Writes into output/5_lasso_models/:
#   lasso_network[_no_hic][_no_footprinting][_with_randomization]
#     [_no_enhancer_tCRE].rds
#     one row per retained edge; columns include source, target, source_type,
#     target_type, source_gene, target_gene, link_type, link_score, estimate
#     (the fitted coefficient), pvalue, adjpvalue, dev_ratio, deltaR2,
#     ols_pvalue, ols_adjpvalue, source_measurement_type,
#     target_measurement_type
#   ..._sample.csv   (1000-row random sample, for quick review)

#### set up ####

rm(list = ls())
library(tidyverse)
library(foreach)
library(doParallel)
library(glmnet)

# ── CONTROL PANEL ─────────────────────────────────────────────────────────────
source('../config.R')  # defines raw_data_folder, primary_data_folder

use_hic <- TRUE
use_footprinting <- TRUE
randomize <- TRUE
use_enhancer_tCRE <- TRUE  # if FALSE: enhancer tCRE excluded as predictor in enhancer-->promoter models and as target in all models
num_perms <- ifelse(randomize, 100, 1)
n.cores <- 90 # parallel::detectCores() - 2
alpha <- 0.95
nfolds <- 10
CRE_types <- c('aCRE', 'tCRE')
power_law <- 0.75 # as in https://www.cell.com/molecular-cell/pdfExtended/S1097-2765(18)30547-1

# set random seed for reproducibility
set.seed(42)

# load data
base_network_file <- file.path(primary_data_folder, '3_preparing_base_network',
                               paste0('base_network',
                                     ifelse(use_hic, '', '_no_hic'),
                                     ifelse(use_footprinting, '', '_no_footprinting'),
                                     '.rds'))
dataset_file <- file.path(primary_data_folder, '4_preparing_dataset',
                          paste0('dataset',
                                ifelse(use_hic, '', '_no_hic'),
                                ifelse(use_footprinting, '', '_no_footprinting'),
                                '.rds'))
pseudo_time_file <- file.path(raw_data_folder, 'pseudotime.scRNA.tsv')
results_folder <- file.path(primary_data_folder, '5_lasso_models')
dir.create(results_folder, showWarnings = FALSE, recursive = TRUE)

# parallel computation environment
my.cluster <- parallel::makeCluster(n.cores)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()

# lading data
base_network <- readRDS(base_network_file)
dataset <- readRDS(file = dataset_file)

# loading pseudotime
pseudotime <- read.table(pseudo_time_file, row.names = 1, header = TRUE)

# ordering dataset according to pseudotime
pseudotime <- pseudotime[colnames(dataset), ]
pseudotime <- pseudotime[order(pseudotime$pseudotime, decreasing = FALSE), ]
dataset <- dataset[, rownames(pseudotime)]

# creating the subdivision in folds for the CV
folds_ids <- rep(1:nfolds, ceiling(ncol(dataset)/nfolds))
folds_ids <- folds_ids[1:ncol(dataset)]

#### Helper function for permutation testing ####

perform_permutation_test <- function(tmp_res, current_predictors, target_values,
                                    alpha, folds_ids, current_weights,
                                    randomize, num_perms, i, j, n_targets) {

  # if randomize, we compute actual p-values and permuted dev ratio
  if(randomize){

    num_success <- 0
    tmp_res$pvalue <- 0
    tmp_res$dev_ratio_mean_permuted <- NA
    tmp_res$dev_ratio_std_permuted <- NA
    dev_ratios <- rep(NA, num_perms)

    k <- 1
    for(k in 1:num_perms){

      # different seed for each target, yet replicable
      set.seed(k * (i + ((j-1) * n_targets)))

      target_values_k <- sample(target_values,
                                size = length(target_values), replace = FALSE)

      cv_res_k <- tryCatch({glmnet::cv.glmnet(x = current_predictors, y = target_values_k,
                                              alpha = alpha, foldid = folds_ids,
                                              penalty.factor = current_weights)},
                           error = function(e){NULL})

      if(!is.null(cv_res_k)){

        # a single failed permutation should not abort the whole loop
        m_k <- tryCatch({glmnet::glmnet(x = current_predictors, y = target_values_k,
                                        alpha = alpha, lambda = cv_res_k$lambda.min,
                                        penalty.factor = current_weights)},
                        error = function(e){NULL})

        if(!is.null(m_k)){

          num_success <- num_success + 1

          tmp_res_k <- data.frame(estimate = as.matrix(coef(m_k))[, 1])
          tmp_res_k <- tmp_res_k[-which(rownames(tmp_res_k) == '(Intercept)'), , drop = FALSE]
          tmp_res$pvalue <- tmp_res$pvalue + ifelse(abs(tmp_res_k$estimate) > abs(tmp_res$estimate), 1, 0)
          dev_ratios[k] <- m_k$dev.ratio

        }

      }

    }

    tmp_res$pvalue <- (tmp_res$pvalue + 1) / (num_success + 1) # avoiding 0 p-values
    sel <- which(tmp_res$estimate != 0)
    tmp_res$adjpvalue <- NA
    if(length(sel) > 0){
      tmp_res$adjpvalue[sel] <- p.adjust(tmp_res$pvalue[sel], method = 'fdr') # adjusting p-values
    }else{
      tmp_res$adjpvalue <- tmp_res$pvalue
    }
    tmp_res$dev_ratio_mean_permuted <- mean(dev_ratios, na.rm = TRUE)
    tmp_res$dev_ratio_std_permuted <- sd(dev_ratios, na.rm = TRUE)

  }

  return(tmp_res)
}

#### Helper function for fitting elastic net model with OLS post-hoc analysis ####

fit_elastic_net_model <- function(current_predictors, target_values,
                                  alpha, folds_ids, current_weights,
                                  randomize, num_perms, i, j, n_targets) {

  tmp_res <- tryCatch({

    cv_res <- glmnet::cv.glmnet(x = current_predictors, y = target_values,
                                alpha = alpha, foldid = folds_ids,
                                penalty.factor = current_weights)
    m <- glmnet::glmnet(x = current_predictors, y = target_values,
                        alpha = alpha, lambda = cv_res$lambda.min,
                        penalty.factor = current_weights)

    tmp_res <- data.frame(estimate = as.matrix(coef(m))[, 1])
    tmp_res <- tmp_res[-which(rownames(tmp_res) == '(Intercept)'), , drop = FALSE]
    tmp_res$pvalue <- ifelse(tmp_res$estimate == 0, 1, 0)
    tmp_res$dev_ratio <- m$dev.ratio
    tmp_res$deltaR2 <- 0            # unique contribution to deltaR2
    tmp_res$ols_pvalue <- NA_real_  # two-sided OLS p-values
    tmp_res$ols_adjpvalue <- NA_real_

    # indices of selected predictors (non-zero lasso/EN coefficients)
    sel <- which(tmp_res$estimate != 0)

    # if any predictor selected, try OLS
    if(length(sel) > 0) {

      # refit OLS on selected predictors only (fast base R)
      Xsel <- current_predictors[, sel, drop = FALSE]
      y    <- target_values

      # guards
      n <- nrow(Xsel)
      s <- ncol(Xsel)
      ysd <- stats::sd(y)

      if(is.finite(ysd) && ysd > 0 && s < (n - 1)) {

        fit_ols <- tryCatch(
          stats::lm.fit(x = cbind(1, Xsel), y = y),
          error = function(e) NULL
        )

        if(!is.null(fit_ols)) {

          # compute R2
          rss <- sum(fit_ols$residuals^2)
          tss <- sum((y - mean(y))^2)
          R2  <- if(tss > 0) 1 - (rss / tss) else 0

          # residual df and MSE
          df_res <- n - (s + 1)
          if(df_res > 0 && is.finite(R2)) {

            mse <- rss / df_res

            # compute t-stats for selected predictors
            Rmat <- fit_ols$qr$qr
            r    <- fit_ols$qr$rank

            if(is.finite(r) && r == (s + 1)) { # full rank (incl intercept)

              XtX_inv <- tryCatch(
                chol2inv(Rmat[1:r, 1:r, drop = FALSE]),
                error = function(e) NULL
              )

              if(!is.null(XtX_inv)) {

                se   <- sqrt(diag(XtX_inv) * mse)
                beta <- fit_ols$coefficients

                # exclude intercept
                beta_j <- beta[-1]
                se_j   <- se[-1]
                tvals  <- beta_j / se_j
                t2     <- tvals^2

                # deltaR2 for selected predictors (fast formula)
                deltaR2_sel <- as.numeric(t2) * (1 - R2) / df_res
                deltaR2_sel[!is.finite(deltaR2_sel)] <- NA_real_
                tmp_res$deltaR2[sel] <- deltaR2_sel

                # two-sided OLS p-values
                ols_p_sel <- 2 * stats::pt(-abs(tvals), df = df_res)
                ols_p_sel[!is.finite(ols_p_sel)] <- NA_real_
                tmp_res$ols_pvalue[sel] <- as.numeric(ols_p_sel)
                tmp_res$ols_adjpvalue[sel] <- p.adjust(as.numeric(ols_p_sel), method = 'fdr')
              }

            } else {
              # rank-deficient: avoid misleading inference
              tmp_res$deltaR2[sel] <- NA_real_
              tmp_res$ols_pvalue[sel] <- NA_real_
              tmp_res$ols_adjpvalue[sel] <- NA_real_
            }
          }
        }
      } else {
        # too many selected vars or degenerate target
        tmp_res$deltaR2[sel] <- NA_real_
        tmp_res$ols_pvalue[sel] <- NA_real_
        tmp_res$ols_adjpvalue[sel] <- NA_real_
      }
    }

    # if randomize, we compute actual p-values and permuted dev ratio
    tmp_res <- perform_permutation_test(tmp_res, current_predictors, target_values,
                                       alpha, folds_ids, current_weights,
                                       randomize, num_perms, i, j, n_targets)

    # return
    tmp_res

  }, error = function(e){

    # if the model fails, we select no predictor
    NULL # we check for this immediately after

  })

  return(tmp_res)
}

#### TF to peaks models ####

# considering only connections from footprinting
idx <- base_network$link_type == 'TF_peak'
tf_peak_network <- base_network[idx, ]

# set parallel-safe RNG
RNGkind("L'Ecuyer-CMRG")
set.seed(42)

# looping over targets
targets <- unique(tf_peak_network$target)
n_targets <- length(targets)

tmp <- foreach(i = 1:n_targets, .packages = c('glmnet', 'tidyverse')) %dopar% {

  # current target
  target <- targets[i]
  if(i %% 1000 == 0){
    print(target)
    sink(paste0('tf_peak_', i, '_', n_targets, '.csv'))
    print('done!')
    sink()
  }

  # general tryCatch
  to_return <- tryCatch({

    # selecting base network
    current_network <- tf_peak_network[tf_peak_network$target == target, ]
    current_network <- unique(current_network)

    # if no elements in the current network, then return
    if(nrow(current_network) == 0){
      return(paste0('No predictors at iteration ', i))
    }

    # summarizing over link score, if needed
    current_network <- current_network %>% group_by(across(c(-link_score))) %>%
      summarise(link_score = max(link_score)) %>% as.data.frame()

    # check if multiple target gene in case of promoters
    if(all(current_network$target_type == 'Promoter') && length(unique(current_network$target_gene)) > 1){
      return(paste0('The same promoter serving multiple gene at iteration ', i))
    }

    # selecting the current predictors
    current_network$source_measurement <- paste0(current_network$source, '_Expression')
    current_predictors <- t(dataset[current_network$source_measurement, , drop = FALSE])

    # scaling, see https://glmnet.stanford.edu/articles/glmnet.html
    current_predictors <- scale(current_predictors, center = TRUE, scale = TRUE)

    # predictors weights... Inf would led a variable to be excluded
    current_weights <- 1 / current_network$link_score

    # looping between the aCRE and tCRE of the target
    current_CRE_types <- if(!use_enhancer_tCRE && all(current_network$target_type == 'Enhancer')) 'aCRE' else CRE_types
    model_res <- vector('list', length(current_CRE_types))
    names(model_res) <- current_CRE_types
    j <- 1
    for(j in 1:length(current_CRE_types)){

      # selecting the values of the target... if not present, skip
      target_measurement <- paste0(target, '_', current_CRE_types[j])
      if(!(target_measurement %in% rownames(dataset))){
        next() # model_res contents remain NULL
      }
      target_values <- as.numeric(dataset[target_measurement, ])

      # if anything to predict, otherwise we forgo the model
      if(sd(target_values) > 0){

        # scaling, see https://glmnet.stanford.edu/articles/glmnet.html
        target_values <- as.numeric(scale(target_values, center = TRUE, scale = TRUE))

        tmp_res <- fit_elastic_net_model(current_predictors, target_values,
                                        alpha, folds_ids, current_weights,
                                        randomize, num_perms, i, j, n_targets)

        # adding to the network... if successful
        if(!is.null(tmp_res) && all(current_network$source_measurement == rownames(tmp_res))){

          # aligning the results
          rownames(tmp_res) <- rownames(current_network)
          tmp <- cbind(current_network, tmp_res)
          tmp$source_measurement <- NULL

          # important! type of measurements
          tmp$source_measurement_type <- 'Expression'
          tmp$target_measurement_type <- names(model_res)[j]

          # important! removing zero coeff predictors
          tmp <- tmp[tmp$estimate != 0, ]

          # adding to model_res
          if(dim(tmp)[1] > 0){
            model_res[[j]] <- tmp
          } # otherwise, model_res contents remain NULL

        } # otherwise, model_res contents remain NULL

      } # if no variation in the target, model_res contents remain NULL

    }

    # return object
    network_to_return <- do.call(rbind, model_res)
    rownames(network_to_return) <- NULL

    # returning
    if(!is.null(network_to_return)){
      return(network_to_return)
    }else{
      return(paste0('No model fitted at iteration ', i))
    }

  }, error = function(e){
    return(paste0('Error at iteration ', i, '\n', e))
  })

}

# removing errors
errors <- which(sapply(tmp, class) == 'character')
if(length(errors) > 0){
  tmp <- tmp[-errors]
}

# binding
tf_peak_network <- do.call('rbind', tmp)

#### enhancer to promoter ####

# re-starting the cluster
stopCluster(my.cluster)
my.cluster <- parallel::makeCluster(n.cores)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()

# considering only connections from enhancer to promoter
idx <- base_network$link_type == 'Enhancer_promoter'
enhancer_promoter_network <- base_network[idx, ]

# set parallel-safe RNG
RNGkind("L'Ecuyer-CMRG")
set.seed(42)

# looping over targets
targets <- unique(enhancer_promoter_network$target)
n_targets <- length(targets)
tmp <- foreach(i = 1:n_targets, .packages = c('glmnet', 'tidyverse')) %dopar% {

  # current target
  target <- targets[i]
  if(i %% 1000 == 0){
    print(target)
    sink(paste0('enhancer_promoter_', i, '_', n_targets, '.csv'))
    print('done!')
    sink()
  }

  # general tryCatch
  to_return <- tryCatch({

    # selecting base network
    current_network <- enhancer_promoter_network[enhancer_promoter_network$target == target, ]
    current_network <- unique(current_network)

    # avoiding loops: peaks flagged as both source and target of the same edge
    current_network <- current_network[current_network$source != current_network$target, ]

    # if no elements in the current network, then return
    if(nrow(current_network) == 0){
      return(paste0('No predictors at iteration ', i))
    }

    # ensuring uniquenss
    current_network <- current_network %>% group_by(across(c(-link_score))) %>%
      summarise(link_score = max(link_score)) %>% as.data.frame()

    # check if multiple target gene in case of promoters
    if(all(current_network$target_type == 'Promoter') && length(unique(current_network$target_gene)) > 1){
      return(paste0('The same promoter serving multiple gene at iteration ', i))
    }

    # selecting the current predictors... we use a trick to cover both aCRE and tCRE
    current_network$source_measurement <- paste0(current_network$source, '_aCRE')
    if(use_enhancer_tCRE){
      current_network_2 <- current_network
      current_network_2$source_measurement <- paste0(current_network_2$source, '_tCRE')
      current_network <- rbind(current_network, current_network_2)
      rm(current_network_2)
    }
    current_network <- current_network[current_network$source_measurement %in% rownames(dataset), ]
    current_predictors <- t(dataset[current_network$source_measurement, , drop = FALSE])

    # scaling, see https://glmnet.stanford.edu/articles/glmnet.html
    current_predictors <- scale(current_predictors, center = TRUE, scale = TRUE)

    # predictors weights... note that it is in Kb
    current_weights <- (1 - (current_network$link_score / 1000) ^ (-power_law))
    current_weights <- pmax(current_weights, 0) # to ensure positive weights

    # looping between the aCRE and tCRE of the target
    model_res <- vector('list', length(CRE_types))
    names(model_res) <- CRE_types
    for(j in 1:length(CRE_types)){

      # selecting the values of the target... if not present, skip
      target_measurement <- paste0(target, '_', CRE_types[j])
      if(!(target_measurement %in% rownames(dataset))){
        next() # model_res contents remain NULL
      }
      target_values <- as.numeric(dataset[target_measurement, ])

      # if anything to predict, otherwise we forgo the model
      if(sd(target_values) > 0){

        # scaling, see https://glmnet.stanford.edu/articles/glmnet.html
        target_values <- as.numeric(scale(target_values, center = TRUE, scale = TRUE))

        tmp_res <- fit_elastic_net_model(current_predictors, target_values,
                                        alpha, folds_ids, current_weights,
                                        randomize, num_perms, i, j, n_targets)

        # adding to the network... if successful
        if(!is.null(tmp_res) && all(current_network$source_measurement == rownames(tmp_res))){

          # aligning the results
          rownames(tmp_res) <- rownames(current_network)
          tmp <- cbind(current_network, tmp_res)

          # important! type of measurements
          tmp$source_measurement_type <- sub(".*_(aCRE|tCRE)$", "\\1",
                                             tmp$source_measurement)
          tmp$target_measurement_type <- names(model_res)[j]
          tmp$source_measurement <- NULL

          # important! removing zero coeff predictors
          tmp <- tmp[tmp$estimate != 0, ]

          # adding to model_res
          if(dim(tmp)[1] > 0){
            model_res[[j]] <- tmp
          } # otherwise, model_res contents remain NULL

        } # otherwise, model_res contents remain NULL

      } # if no variation in the target, model_res contents remain NULL

    }

    # return object
    network_to_return <- do.call(rbind, model_res)
    rownames(network_to_return) <- NULL

    # returning
    if(!is.null(network_to_return)){
      return(network_to_return)
    }else{
      return(paste0('No model fitted at iteration ', i))
    }

  }, error = function(e){
    return(paste0('Error at iteration ', i, '\n', e))
  })

}

# removing errors
errors <- which(sapply(tmp, class) == 'character')
if(length(errors) > 0){
  tmp <- tmp[-errors]
}

# binding
enhancer_promoter_network <- do.call('rbind', tmp)

#### saving ####

# binding the networks together
lasso_network <- rbind(tf_peak_network, enhancer_promoter_network)

# writing
idx <- sample(1:dim(lasso_network)[1], 1000)
write.csv(lasso_network[idx, ], row.names = FALSE,
          file.path(results_folder, paste0('lasso_network',
                                           ifelse(use_hic, '', '_no_hic'),
                                           ifelse(use_footprinting, '', '_no_footprinting'),
                                           ifelse(randomize, '_with_randomization', ''),
                                           ifelse(use_enhancer_tCRE, '', '_no_enhancer_tCRE'),
                                           '_sample.csv')))
saveRDS(lasso_network, file.path(results_folder,
                                 paste0('lasso_network',
                                        ifelse(use_hic, '', '_no_hic'),
                                        ifelse(use_footprinting, '', '_no_footprinting'),
                                        ifelse(randomize, '_with_randomization', ''),
                                        ifelse(use_enhancer_tCRE, '', '_no_enhancer_tCRE'),
                                        '.rds')))

message('Lasso network edges retained: ', nrow(lasso_network))

# stopping the cluster
stopCluster(my.cluster)
