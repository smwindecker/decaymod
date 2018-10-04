# try and find old versions of the create_jobs

# and model calls.
# Want them to save .rds for each model - putting each cv iteration into a list

# Create job list of all model iterations to be run
#
# @param models dataframe listing model types and specifications
# @param data dataframe containing data for models
# @param cv_cluster optional columns name to specify cv clusters
# @param fixed_effects_list optional list of fixed effects - corresponding to column names in data
# @return list of all model iterations
# traits <- c('N', 'C', 'SLA', 'DMC', 'HC', 'CL', 'LG')
# model_df <- data.frame(model_type = c('w', 'w', 'ne'),
#                        fixed_effects = c(TRUE, FALSE, FALSE),
#                        random_effects = c(FALSE, TRUE, FALSE))
# jobs_list <- create_jobs(models = model_df,
#                          data = df,
#                          cv_cluster = 'species_code',
#                          fixed_effects_list = traits)

# want this to create all my jobs.therefore doesn't need a cv cluster or fixed effects

create_jobs <- function (models, data, cv_cluster = NULL, fixed_effects_list = NULL) {

  if (is.null(fixed_effects_list) & nrow(models[isTRUE(models$fixed_effects), ]) > 0) {
    stop('fixed effects reported but no list provided')
  }

  if (!is.null(fixed_effects_list)) {

    # create list of traits as expressions with tilda and intercept terms
    formulas <- c('~ 1', sprintf("%s + %s", '~ 1', fixed_effects_list))

    # prep for negative exponential model
    n <- length(formulas)
    ne_formulas <- data.frame(model_type = rep('ne', n),
                              fixed_effects = rep('FE', n),
                              formula_k = formulas,
                              formula_alpha = NA,
                              formula_beta = NA,
                              stringsAsFactors = FALSE)
    models_ne <- merge(models, ne_formulas, by = c('model_type', 'fixed_effects'), all.x = TRUE)

    # prep two parameter merge
    formula_grid <- expand.grid(formulas, formulas, stringsAsFactors = FALSE)
    g <- nrow(formula_grid)

    # prep for weibull model
    w_formulas <- data.frame(model_type = rep('w', g),
                             fixed_effects = rep('FE', g),
                             formula_k = NA,
                             formula_alpha = formula_grid[,1],
                             formula_beta = formula_grid[,2],
                             stringsAsFactors = FALSE)

    models <- merge(models_ne, w_formulas, by = c('model_type', 'fixed_effects',
                                                  'formula_k', 'formula_alpha', 'formula_beta'),
                    all.y = TRUE)

  }

  if (is.null(fixed_effects_list)) {
    models$formula_k <- NA
    models$formula_alpha <- NA
    models$formula_beta <- NA
  }

  # give formulas to non fixed effects models
  models$formula_k[models$model_type == 'ne' & is.na(models$formula_k)] <- '~ 1'
  models$formula_alpha[models$model_type == 'w' & is.na(models$formula_alpha)] <- '~ 1'
  models$formula_beta[models$model_type == 'w' & is.na(models$formula_beta)] <- '~ 1'

  # number the models
  models$model <- seq(1:nrow(models))

  models <- models[, - which(names(models) %in% 'fixed_effects')]

  # change dataframe to list
  jobs <- split(models, seq(nrow(models)))

  # add filename attribute -- what's the job ID in this case?
  # jobs$filename <- NULL
  # jobs$filename <- sprintf('R/stanOUTPUT/%s.rds', jobs$job_id)

  # if cv cluster is provided, these jobs are cv. and therefore need to be expanded
  if (!is.null(cv_cluster)) {
    clusters <- unique(data[, cv_cluster])
    tmp <- lapply(jobs, expand_models, clusters)
    jobs <- unlist(tmp, recursive = FALSE)
  }

  jobs

}

expand_models <- function (input, clusters) {

  # repeat each model n = number of species times
  expand <- expand.grid(cv_cluster = clusters,
                        model = input$model)
  expand$model <- as.character(expand$model)

  # merge with formula dataframe so also have the parameter functions listed
  job <- merge(expand, input, by = 'model')

  list <- split(job, seq(nrow(job)))

  return(list)

}

## Model analysis

# runs all models from jobs list
run_models <- function (jobs_list, data, initial_mass, removal_mass, time, group, n_cores) {

  registerDoMC(n_cores)

  # so my decaymod outputs the fit. i want to append it to the list.
  model_output <- foreach::foreach(i = 1:length(jobs_list)) %dopar% {
    output_i <- evaluate_decaymod(jobs_list[[i]], data, initial_mass, removal_mass, time, group)
    return(output_i)
  }

  return(model_output)
}

# call the function to fit the model, extracts the fit, diagnostics, and negative loglikelihood
evaluate_decaymod <- function (df, data, initial_mass, removal_mass, time, group) {

  if ('cv_cluster' %in% colnames(df)) {
    cross_validation <- TRUE
    group_id <- df$cv_cluster
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
                  formula_k = df$formula_k,
                  formula_alpha = df$formula_alpha,
                  formula_beta = df$formula_beta)

  # diagnostics
  fit_summary <- rstan::summary(fit)$summary
  abs_rhat <- max(abs(fit_summary[,'Rhat'] - 1))
  neff_min <- min(fit_summary[,'n_eff'])
  sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
  sum_div <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
  max_treedepth <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
  diagnostics <- data.frame(abs_rhat = abs_rhat,
                            neff_min = neff_min,
                            sum_div = sum_div,
                            max_treedepth = max_treedepth,
                            stringsAsFactors = FALSE)

  if (isTRUE(cross_validation)) {
    neg_loglik_df <- as.data.frame(fit, 'neg_loglik')
    df$neg_loglik <- mean(neg_loglik_df$neg_loglik, na.rm = TRUE)
  }

  fit_list <- list(mod_specs = df,
                   fit = fit,
                   diagnostics = diagnostics)

  fit_list
}

# conducts the model fit itself
decaymod <- function (data, initial_mass, removal_mass, time, group,
                      model_type, random_effects, cross_validation = FALSE, group_id = NA,
                      formula_k = NA, formula_alpha = NA, formula_beta = NA) {

  if (isTRUE(cross_validation)) {
    train <- data[data[, group] != group_id, ]
    test <- data[data[, group] == group_id, ]

    # assign level to each group in the training dataset (as now has one fewer item)
    train$group_level <- as.factor(as.numeric(as.factor(as.character(train[, group]))))

    # number of levels
    n_train_levels <- length(unique(train$group_level))

    if (model_type == 'ne') {
      X <- make_matrix(model_type = 'ne', X_type = 'train', n_train_levels. = n_train_levels,
                       train. = train, formula_k. = formula_k)
      X_test <- make_matrix(model_type = 'ne', X_type = 'test', n_train_levels. = n_train_levels,
                            test. = test, formula_k. = formula_k)

      # for model without random effects
      if (!isTRUE(random_effects)) {
        fit <- ne_CV_noRE_stan(mT = train[, removal_mass],
                               mT_test = test[, removal_mass],
                               m0 = train[, initial_mass],
                               m0_test = test[, initial_mass],
                               time = train[, time],
                               time_test = test[, time],
                               sp = as.numeric(as.factor(train$group_level)),
                               J = n_train_levels,
                               X = X,
                               X_test = X_test)
      }

      # for model with random effects
      if (isTRUE(random_effects)) {
        fit <- ne_CV_RE_stan(mT = train[, removal_mass],
                             mT_test = test[, removal_mass],
                             m0 = train[, initial_mass],
                             m0_test = test[, initial_mass],
                             time = train[, time],
                             time_test = test[, time],
                             sp = as.numeric(as.factor(train$group_level)),
                             J = n_train_levels,
                             X = X,
                             X_test = X_test)
      }
    }

    if (model_type == 'w') {
      X_alpha <- make_matrix(model_type = 'w', X_type = 'train', n_train_levels. = n_train_levels,
                             train. = train, formula_alpha. = formula_alpha)
      X_beta <- make_matrix(model_type = 'w', X_type = 'train', n_train_levels. = n_train_levels,
                            train. = train, formula_beta. = formula_beta)
      X_alpha_test <- make_matrix(model_type = 'w', X_type = 'test', n_train_levels. = n_train_levels,
                                  test. = test, formula_alpha. = formula_alpha)
      X_beta_test <- make_matrix(model_type = 'w', X_type = 'test', n_train_levels. = n_train_levels,
                                 test. = test, formula_beta. = formula_beta)

      # for weibull model without random effects
      if (!isTRUE(random_effects)) {
        fit <- w_CV_noRE_stan(mT = train[, removal_mass],
                              mT_test = test[, removal_mass],
                              m0 = train[, initial_mass],
                              m0_test = test[, initial_mass],
                              time = train[, time],
                              time_test = test[, time],
                              sp = as.numeric(as.factor(train$group_level)),
                              J = n_train_levels,
                              X_alpha = X_alpha,
                              X_alpha_test = X_alpha_test,
                              X_beta = X_beta,
                              X_beta_test = X_beta_test)
      }

      # for weibull model with random effects
      if (isTRUE(random_effects)) {
        fit <- w_CV_RE_stan(mT = train[, removal_mass],
                            mT_test = test[, removal_mass],
                            m0 = train[, initial_mass],
                            m0_test = test[, initial_mass],
                            time = train[, time],
                            time_test = test[, time],
                            sp = as.numeric(as.factor(train$group_level)),
                            J = n_train_levels,
                            X_alpha = X_alpha,
                            X_alpha_test = X_alpha_test,
                            X_beta = X_beta,
                            X_beta_test = X_beta_test)
      }
    }
  }

  if(!isTRUE(cross_validation)) {

    train <- data

    # assign level to each group in the training dataset (as now has one fewer item)
    train$group_level <- as.factor(as.numeric(as.factor(as.character(train[, group]))))

    # number of levels
    n_train_levels <- length(unique(train$group_level))

    if (model_type == 'ne') {
      X <- make_matrix(model_type = 'ne', X_type = 'train', n_train_levels. = n_train_levels,
                       train. = train, formula_k. = formula_k)

      if (isTRUE(random_effects)) {
        sim_df <- data.frame(m0_sim = rep(log(4100), 5800),
                             time_sim = rep(seq(0, 0.7, length.out = 200), 29),
                             sp_sim = rep(1:29, each = 200))

        fit <- ne_noCV_RE_stan(mT = train[, removal_mass],
                               m0 = train[, initial_mass],
                               time = train[, time],
                               sp = as.numeric(as.factor(train$group_level)),
                               J = n_train_levels,
                               X = X,
                               m0_sim = sim_df$m0_sim,
                               time_sim = sim_df$time_sim,
                               sp_sim = sim_df$sp_sim)
      }

      if (!isTRUE(random_effects)) {
        fit <- ne_noCV_noRE_stan(mT = train[, removal_mass],
                                 m0 = train[, initial_mass],
                                 time = train[, time],
                                 sp = as.numeric(as.factor(train$group_level)),
                                 J = n_train_levels,
                                 X = X)
      }
    }

    if (model_type == 'w') {
      X_alpha <- make_matrix(model_type = 'w', X_type = 'train', n_train_levels. = n_train_levels,
                             train. = train, formula_alpha. = formula_alpha)
      X_beta <- make_matrix(model_type = 'w', X_type = 'train', n_train_levels. = n_train_levels,
                            train. = train, formula_beta. = formula_beta)

      if (isTRUE(random_effects)) {
        sim_df <- data.frame(m0_sim = rep(log(4100), 5800),
                             time_sim = rep(seq(0, 0.7, length.out = 200), 29),
                             sp_sim = rep(1:29, each = 200))

        fit <- w_noCV_RE_stan(mT = train[, removal_mass],
                              m0 = train[, initial_mass],
                              time = train[, time],
                              sp = as.numeric(as.factor(train$group_level)),
                              J = n_train_levels,
                              X_alpha = X_alpha,
                              X_beta = X_beta,
                              m0_sim = sim_df$m0_sim,
                              time_sim = sim_df$time_sim,
                              sp_sim = sim_df$sp_sim)
      }

      if (!isTRUE(random_effects)) {
        fit <- w_noCV_noRE_stan(mT = train[, removal_mass],
                                m0 = train[, initial_mass],
                                time = train[, time],
                                sp = as.numeric(as.factor(train$group_level)),
                                J = n_train_levels,
                                X_alpha = X_alpha,
                                X_beta = X_beta)
      }
    }
  }

  fit
}

# creates model matrix
make_matrix <- function (model_type, X_type, n_train_levels. = n_train_levels,
                         train. = NULL, test. = NULL, formula_k. = NULL,
                         formula_alpha. = NULL, formula_beta. = NULL) {

  # negative exponential model matrices
  if (X_type == 'train') {
    if (model_type == 'ne') {
      if (formula_k. == '~1' | formula_k. == '~ 1') {
        mat <- as.matrix(model.matrix(as.formula(formula_k.), train.)[1:n_train_levels.,])
      }
      else {
        mat <- unique(model.matrix(as.formula(formula_k.), train.))
      }
    }

    if (model_type == 'w' & !is.null(formula_alpha.)) {

      # create alpha model matrix for training data
      if (formula_alpha. == '~1' | formula_alpha. == '~ 1') {
        mat <- as.matrix(model.matrix(as.formula(formula_alpha.), train.)[1:n_train_levels., ])
      }
      else {
        mat <- unique(model.matrix(as.formula(formula_alpha.), train.))
      }
    }

    if (model_type == 'w' & !is.null(formula_beta.)) {

      # create beta model matrix for training data
      if (formula_beta. == '~1' | formula_beta. == '~ 1') {
        mat <- as.matrix(model.matrix(as.formula(formula_beta.), train.)[1:n_train_levels., ])
      }
      else {
        mat <- unique(model.matrix(as.formula(formula_beta.), train.))
      }
    }
  }

  if (X_type == 'test') {
    if (model_type == 'ne') mat <- unique(model.matrix(as.formula(formula_k.), test.))
    if (model_type == 'w') {
      if (!is.null(formula_alpha.)) mat <- unique(model.matrix(as.formula(formula_alpha.), test.))
      if (!is.null(formula_beta.)) mat <- unique(model.matrix(as.formula(formula_beta.), test.))
    }
  }

  mat
}

