## Analysis functions

#' creates job list for all models
#'
#' @export
#' @param model_type 'w' or 'ne'
#' @param data data
#' @param cv_cluster if cross-validation to be used, what to group by
#' @param random_effects logical
#' @param fixed_effects list of trait combinations
#' @param trait_param null or both
#' @importFrom gtools combinations
#' @return jobs list
#'
create_jobs <- function (model_type,
                         data,
                         cv_cluster = NULL,
                         random_effects = TRUE,
                         fixed_effects = NULL,
                         trait_param = NULL) {

  if (is.null(fixed_effects)) {
    param_formula <- list('~ 1')
  }

  if (!is.null(fixed_effects)) {
    param_formula <- lapply(fixed_effects, mypaste)
  }

  # for rows that are FE and traits on one only
  if (!is.null(cv_cluster) & is.null(trait_param)) {

    clusters <- unique(data[, cv_cluster])
    tmp <- lapply(param_formula, expand_models, clusters)
    tmp_expanded <- unlist(tmp, recursive = FALSE)

    jobs <- lapply(tmp_expanded, model_details,
                   model_type = model_type,
                   random_effects = random_effects,
                   cv = TRUE)
  }

  # not FE traits on one only
  if (is.null(cv_cluster) & is.null(trait_param)) {

    jobs <- lapply(param_formula, model_details,
                   model_type = model_type,
                   random_effects = random_effects,
                   cv = FALSE)
  }

  # for rows that are FE traits on both
  if (!is.null(fixed_effects) & !is.null(trait_param)) {
    unlist_formulas <- unlist(param_formula)
    x <- gtools::combinations(n = length(unlist_formulas), r = 2, v = unlist_formulas, repeats.allowed = TRUE)
    ext_param_formulas <- tapply(x, rep(1:nrow(x), ncol(x)), function(i) i)
  }

  if (!is.null(cv_cluster) & !is.null(trait_param)) {

    clusters <- unique(data[, cv_cluster])
    tmp <- lapply(ext_param_formulas, expand_models, clusters, both = TRUE)
    tmp_expanded <- unlist(tmp, recursive = FALSE)

    jobs <- lapply(tmp_expanded, model_details,
                   model_type = model_type,
                   random_effects = random_effects,
                   cv = TRUE)
  }

  if (is.null(cv_cluster) & !is.null(trait_param)) {

    jobs <- lapply(ext_param_formulas, model_details,
                   model_type = model_type,
                   random_effects = random_effects,
                   cv = FALSE)
  }


  return(jobs)
}


#' puts the parameter formulas in correct format
#'
#' @param vector trait list to be modified into formula form
#' @return trait formula in proper form
#'
mypaste <- function (vector) {

  vectornew <- c('1', vector)
  vectornew <- paste(vectornew, sep = '', collapse = ' + ')
  vectornew <- paste('~ ', vectornew, sep = ' ')
  vectornew

}


#' put all model details in correct format
#'
#' @param param_formula formula of fixed effect to be estimated
#' @param model_type 'w' or 'ne'
#' @param random_effects logical
#' @param cv logical
#' @return model dataframe in proper form
#'
model_details <- function (param_formula, model_type, random_effects, cv) {

  my_df <- data.frame(param_formula,
                      model_type = model_type,
                      random_effects = random_effects,
                      cross_validation = cv)

}


#' replicate model rows for each cross-validation cluster
#'
#' @param input fixed effect model
#' @param clusters list of clusters for cross validation
#' @param both whether traits on one or both traits
#' @return model dataframe expanded for cross val levels
#'
expand_models <- function (input, clusters, both = FALSE) {

  # repeat each model n = number of species times

  if (!isTRUE(both)) {
    expand <- expand.grid(cluster = clusters,
                          param_formula = input)

    expand$param_formula <- as.character(expand$param_formula)
  }


  if (isTRUE(both)) {
    expand <- expand.grid(cluster = clusters,
                          alpha_formula = input[1],
                          beta_formula = input[2])

    expand$alpha_formula <- as.character(expand$alpha_formula)
    expand$beta_formula <- as.character(expand$beta_formula)
  }


  list <- split(expand, seq(nrow(expand)))

  return(list)
}


#' runs all models from jobs list
#'
#' @export
#' @param jobs_list list of jobs to be run
#' @param data data
#' @param initial_mass initial mass variable
#' @param removal_mass removal mass variable
#' @param time time variable
#' @param group group for random effects
#' @param n_cores number of cores to run model on
#' @param save_fit logical
#' @param trait_param alpha or beta - traits on which one?
#' @return all models outputs
#'
run_models <- function (jobs_list, data, initial_mass,
                        removal_mass, time, group,
                        n_cores, save_fit = FALSE,
                        trait_param) {

  doMC::registerDoMC(n_cores)

  `%dopar%` <- foreach::`%dopar%`

  # so my decaymod outputs the fit. i want to append it to the list.
  model_output <- foreach::foreach(i = 1:length(jobs_list)) %dopar% {
    output_i <- evaluate_decaymod(jobs_list[[i]], data, initial_mass,
                                  removal_mass, time, group, save_fit,
                                  trait_param)
    return(output_i)
  }

  # negative log likelihood results
  mod_specs <- dplyr::bind_rows(lapply(1:length(model_output), function(x) {
    return(model_output[[x]]$mod_specs)
  }))

  diag <- dplyr::bind_rows(lapply(1:length(model_output), function(x) {
    return(model_output[[x]]$diag)
  }))

  if (isTRUE(save_fit)) {
    model_results <- list(mod_specs = mod_specs,
                          diagnostics = diag,
                          fit = model_output[[1]]$fit)
  }

  if (!isTRUE(save_fit)) {
    model_results <- list(mod_specs = mod_specs,
                          diagnostics = diag)
  }

  return(model_results)
}


#' call the function to fit the model, extracts the fit, diagnostics, and negative loglikelihood
#'
#' @param df job list item as dataframe row
#' @param data data
#' @param initial_mass initial mass variable
#' @param removal_mass removal mass variable
#' @param time time variable
#' @param group group for random effects
#' @param save_fit logical
#' @return evaluate model for individual list item
#'
evaluate_decaymod <- function (df, data, initial_mass, removal_mass,
                               time, group, save_fit, trait_param) {

  if ('cluster' %in% colnames(df)) {
    cross_validation <- TRUE
    group_id <- df$cluster
  } else {
    cross_validation <- FALSE
    group_id <- NA
  }

  fit <- decaymod(data = data,
                  initial_mass = initial_mass,
                  removal_mass = removal_mass,
                  time = time,
                  group = group,
                  model_type = df$model_type,
                  random_effects = df$random_effects,
                  cross_validation = cross_validation,
                  group_id = group_id,
                  param_formula = df$param_formula,
                  trait_param = trait_param)

  # diagnostics
  fit_summary <- rstan::summary(fit)$summary
  abs_rhat <- max(abs(fit_summary[,'Rhat'] - 1))
  neff_min <- min(fit_summary[,'n_eff'])
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  sum_div <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
  max_treedepth <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
  diagnostics <- data.frame(chain = 1:4,
                            abs_rhat = abs_rhat,
                            neff_min = neff_min,
                            sum_div = sum_div,
                            max_treedepth = max_treedepth,
                            stringsAsFactors = FALSE)

  if (isTRUE(cross_validation)) {
    neg_loglik_df <- as.data.frame(fit, 'neg_loglik')
    df$neg_loglik <- mean(neg_loglik_df$neg_loglik, na.rm = TRUE)
  }

  # merge df and diagnostics.
  df_expanded <- df[rep(row.names(df), 4),]
  df_diag <- cbind(df_expanded, diagnostics)

  if (isTRUE(save_fit)) {
    fit_list <- list(mod_specs = df,
                     diag = df_diag,
                     fit = fit)
  }

  if (!isTRUE(save_fit)) {
    fit_list <- list(mod_specs = df,
                     diag = df_diag)
  }

  return(fit_list)
}


#' function to fit the model
#'
#' @param data data
#' @param initial_mass initial mass variable
#' @param removal_mass removal mass variable
#' @param time time variable
#' @param group group for random effects
#' @param model_type 'w' or 'ne'
#' @param random_effects logical
#' @param cross_validation logical
#' @param group_id random effects grouping
#' @param param_formula formula for parameters i model
#' @return model output
#'
decaymod <- function (data, initial_mass, removal_mass, time, group,
                      model_type, random_effects,
                      cross_validation, group_id = NA,
                      param_formula, trait_param) {

  if (isTRUE(cross_validation)) {
    train <- data[data[, group] != group_id, ]
    test <- data[data[, group] == group_id, ]

    mT_test <- test[, removal_mass]
    m0_test <- test[, initial_mass]
    time_test <- test[, time]
    X_test <- unique(model.matrix(as.formula(param_formula), test))
    X_null_test <- unique(model.matrix(as.formula('~ 1'), test))
  }

  if (!isTRUE(cross_validation)) {
    train <- data
  }

  # number of levels
  # assign level to each group in the training dataset (as now has one fewer item)
  train$group_level <- as.factor(as.numeric(as.factor(as.character(train[, group]))))

  J <- length(unique(train$group_level))

  # assign level to each group in the training dataset (as now has one fewer item)
  sp <- as.numeric(as.factor(train$group_level))

  mT <- train[, removal_mass]
  m0 <- train[, initial_mass]
  time <- train[, time]

  X_null <- as.matrix(model.matrix(as.formula('~ 1'), train)[1:J,])

  if (param_formula == '~1' | param_formula == '~ 1') {
    X <- X_null
  } else {
    X <- unique(model.matrix(as.formula(as.character(param_formula)), train))
  }

  sim_df <- data.frame(m0_sim = rep(log(4100), 5800),
                       time_sim = rep(seq(0, 0.7, length.out = 200), 29),
                       sp_sim = rep(1:29, each = 200))

  if (model_type == 'ne' & isTRUE(cross_validation) & !isTRUE(random_effects)) {
    fit <- ne_CV_noRE_stan(mT, mT_test, m0, m0_test, time,
                           time_test, sp, J, X, X_test)
  }

  if (model_type == 'ne' & isTRUE(cross_validation) & isTRUE(random_effects)) {
    fit <- ne_CV_RE_stan(mT, mT_test, m0, m0_test, time,
                         time_test, sp, J, X, X_test)
  }

  if (model_type == 'ne' & !isTRUE(cross_validation) & !isTRUE(random_effects)) {
    fit <- ne_noCV_noRE_stan(mT, m0, time, sp, J, X)
  }

  if (model_type == 'ne' & !isTRUE(cross_validation) & isTRUE(random_effects)) {
    fit <- ne_noCV_RE_stan(mT, m0, time, sp, J, X,
                           m0_sim = sim_df$m0_sim,
                           time_sim = sim_df$time_sim,
                           sp_sim = sim_df$sp_sim)
  }

  # w, cv, no re
  if (model_type == 'w' & isTRUE(cross_validation) & !isTRUE(random_effects)
      & trait_param == 'alpha') {
    fit <- w_CV_noRE_stan(mT, mT_test, m0, m0_test, time,
                          time_test, sp, J,
                          X_alpha = X, X_alpha_test = X_test,
                          X_beta = X_null, X_beta_test = X_null_test)
  }

  if (model_type == 'w' & isTRUE(cross_validation) & !isTRUE(random_effects)
      & trait_param == 'beta') {
    fit <- w_CV_noRE_stan(mT, mT_test, m0, m0_test, time,
                          time_test, sp, J,
                          X_alpha = X_null, X_alpha_test = X_null_test,
                          X_beta = X, X_beta_test = X_test)
  }

  if (model_type == 'w' & isTRUE(cross_validation) & !isTRUE(random_effects)
      & trait_param == 'both') {
    fit <- w_CV_noRE_stan(mT, mT_test, m0, m0_test, time,
                          time_test, sp, J,
                          X_alpha = X, X_alpha_test = X_test,
                          X_beta = X, X_beta_test = X_test)
  }

  # w, cv, re
  if (model_type == 'w' & isTRUE(cross_validation) & isTRUE(random_effects)
      & trait_param == 'alpha') {
    fit <- w_CV_RE_stan(mT, mT_test, m0, m0_test, time,
                        time_test, sp, J,
                        X_alpha = X, X_alpha_test = X_test,
                        X_beta = X_null, X_beta_test = X_null_test)
  }

  if (model_type == 'w' & isTRUE(cross_validation) & isTRUE(random_effects)
      & trait_param == 'beta') {
    fit <- w_CV_RE_stan(mT, mT_test, m0, m0_test, time,
                        time_test, sp, J,
                        X_alpha = X_null, X_alpha_test = X_null_test,
                        X_beta = X, X_beta_test = X_test)
  }

  if (model_type == 'w' & isTRUE(cross_validation) & isTRUE(random_effects)
      & trait_param == 'both') {
    fit <- w_CV_RE_stan(mT, mT_test, m0, m0_test, time,
                        time_test, sp, J,
                        X_alpha = X, X_alpha_test = X_test,
                        X_beta = X, X_beta_test = X_test)
  }

  # w, no cv, no re
  if (model_type == 'w' & !isTRUE(cross_validation) & !isTRUE(random_effects)
      & trait_param == 'alpha') {
    fit <- w_noCV_noRE_stan(mT, m0, time, sp, J, X_alpha = X, X_beta = X_null)
  }

  if (model_type == 'w' & !isTRUE(cross_validation) & !isTRUE(random_effects)
      & trait_param == 'beta') {
    fit <- w_noCV_noRE_stan(mT, m0, time, sp, J, X_alpha = X_null, X_beta = X)
  }

  if (model_type == 'w' & !isTRUE(cross_validation) & !isTRUE(random_effects)
      & trait_param == 'both') {
    fit <- w_noCV_noRE_stan(mT, m0, time, sp, J, X_alpha = X, X_beta = X)
  }

  # w, no cv, re
  if (model_type == 'w' & !isTRUE(cross_validation) & isTRUE(random_effects)
      & trait_param == 'alpha') {
    fit <- w_noCV_RE_stan(mT, m0, time, sp, J, X_alpha = X,
                          X_beta = X_null,
                          m0_sim = sim_df$m0_sim,
                          time_sim = sim_df$time_sim,
                          sp_sim = sim_df$sp_sim)
  }

  if (model_type == 'w' & !isTRUE(cross_validation) & isTRUE(random_effects)
      & trait_param == 'beta') {
    fit <- w_noCV_RE_stan(mT, m0, time, sp, J, X_alpha = X_null,
                          X_beta = X,
                          m0_sim = sim_df$m0_sim,
                          time_sim = sim_df$time_sim,
                          sp_sim = sim_df$sp_sim)
  }

  if (model_type == 'w' & !isTRUE(cross_validation) & isTRUE(random_effects)
      & trait_param == 'both') {
    fit <- w_noCV_RE_stan(mT, m0, time, sp, J, X_alpha = X,
                          X_beta = X,
                          m0_sim = sim_df$m0_sim,
                          time_sim = sim_df$time_sim,
                          sp_sim = sim_df$sp_sim)
  }

  fit
}
