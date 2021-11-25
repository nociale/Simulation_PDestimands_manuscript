create_data_ice <- function(data, imp_drug, imp_placebo) {

    levels_visit <- levels(data$visit)
    J <- length(levels_visit)

    data_hasICE <- data %>%
        group_by(patnum) %>%
        mutate(hasICE = any(!trt_stop_visit %in% c(levels_visit[J], "Inf"))) %>%
        filter(hasICE) %>%
        dplyr::select(-hasICE)

    data_ice <- data_hasICE %>%
        group_by(patnum) %>%
        summarise("strategy" = ifelse(group[1] == "Intervention", imp_drug, imp_placebo),
                  "visit" = ifelse(trt_stop_visit[1] == "0", "1", levels_visit[which(trt_stop_visit[1] == levels_visit) + 1]))

    return(data_ice)
}

impute_analyze_pool <- function(draws_obj, data_ice, references, vars, visit_analysis) {

    res_imp <- impute(
        draws = draws_obj,
        update_strategy = data_ice,
        references = references
    )

    vars$covariates <- "x_bl"
    analyze_res <- analyse(
        imputations = res_imp,
        fun = ancova,
        vars = vars,
        visits = visit_analysis
    )

    if(class(draws_obj)[2] == "condmean") {
        res <- pool(
            analyze_res,
            alternative = "two.sided",
            type = "normal"
        )
    } else {
        res <- pool(
            analyze_res,
            alternative = "two.sided"
        )
    }

    return(res)
}

run_simul <- function(data) {

    # set "outcome" later
    vars <- set_vars(
        subjid = "patnum",
        visit = "visit",
        group = "group",
        covariates = c("x_bl*visit","group*visit"),
        strategy = "strategy"
    )

    n_imputations_bayes <- 100

    references <- c("Control" = "Control", "Intervention" = "Control")

    # create data_ice
    data_ice_JR <- create_data_ice(data, "JR", "JR")
    data_ice_CIR <- create_data_ice(data, "CIR", "CIR")

    levels_visit <- levels(data[[vars$visit]])
    J <- length(levels_visit)

    method <- method_bayes(
        n_samples = n_imputations_bayes,
        burn_between = 50,
        verbose = FALSE
    )

    ########################### HYBRID ESTIMAND ##############################

    print("Running Bayesian multiple imputation (mixed estimand)")

    vars$outcome <- "y"

    ########################### CIR-BASED ESTIMATOR

    draws_obj <- tryCatch_WEM(
        expr = draws(
            data,
            data_ice = data_ice_CIR,
            vars,
            method = method
        ),
        NULL
    )

    if(!is.null(draws_obj$value)) {
        res_mixed_CIR <- impute_analyze_pool(
            draws_obj$value,
            data_ice_CIR,
            references,
            vars,
            visit_analysis = levels_visit[J]
        )
        res_mixed_CIR$failed <- FALSE
        res_mixed_CIR$error <- NULL
    } else {
        res_mixed_CIR <- list(
            failed = TRUE,
            error = draws_obj$error
        )
    }

    ########################### JR-BASED ESTIMATOR
    if(!is.null(draws_obj$value)) {
        res_mixed_JR <- impute_analyze_pool(
            draws_obj$value,
            data_ice_JR,
            references,
            vars,
            visit_analysis = levels_visit[J]
        )
        res_mixed_JR$failed <- FALSE
        res_mixed_JR$error <- NULL
    } else {
        res_mixed_JR <- list(
            failed = TRUE,
            error = draws_obj$error
        )
    }


    ########################### time_from_ice-BASED ESTIMATOR
    vars_time_var <- vars
    vars_time_var$covariates <- c("x_bl*visit", "group*visit", "time_from_trt_stop*group")

    draws_obj <- tryCatch_WEM(
        expr = draws(
            data,
            data_ice = NULL,
            vars_time_var,
            method = method
        ),
        NULL
    )

    if(!is.null(draws_obj$value)) {
        res_mixed_time_from_trt_stop <- impute_analyze_pool(
            draws_obj$value,
            NULL,
            references,
            vars_time_var,
            visit_analysis = levels_visit[J]
        )
        res_mixed_time_from_trt_stop$failed <- FALSE
        res_mixed_time_from_trt_stop$error <- NULL
    } else {
        res_mixed_time_from_trt_stop <- list(
            failed = TRUE,
            error = draws_obj$error
        )
    }

    ########################### is_post_ice-BASED ESTIMATOR
    vars_time_var$covariates <- c("x_bl*visit", "group*visit", "is_post_trt_stop*group")

    draws_obj <- tryCatch_WEM(
        expr = draws(
            data,
            data_ice = NULL,
            vars_time_var,
            method = method
        ),
        NULL
    )

    if(!is.null(draws_obj$value)) {
        res_mixed_is_post_trt_stop <- impute_analyze_pool(
            draws_obj$value,
            NULL,
            references,
            vars_time_var,
            visit_analysis = levels_visit[J]
        )
        res_mixed_is_post_trt_stop$failed <- FALSE
        res_mixed_is_post_trt_stop$error <- NULL
    } else {
        res_mixed_is_post_trt_stop <- list(
            failed = TRUE,
            error = draws_obj$error
        )
    }

    ########################### MMRM ESTIMATOR

    # repeat draws to include post treatment disc data in the base imputation model fit
    draws_obj <- tryCatch_WEM(
        expr = draws(
            data,
            data_ice = NULL,
            vars,
            method = method
        ),
        NULL
    )

    if(!is.null(draws_obj$value)) {
        res_mixed_MAR <- impute_analyze_pool(
            draws_obj$value,
            NULL,
            references,
            vars,
            visit_analysis = levels_visit[J]
        )
        res_mixed_MAR$failed <- FALSE
        res_mixed_MAR$error <- NULL
    } else {
        res_mixed_MAR <- list(
            failed = TRUE,
            error = draws_obj$error
        )
    }

    ########################### TREATMENT POLICY ESTIMAND ##############################

    print("Running Bayesian multiple imputation (treatment policy estimand)")

    vars$outcome <- "y_dropout"
    ########################### MMRM ESTIMATOR

    draws_obj <- tryCatch_WEM(
        expr = draws(
            data,
            data_ice = NULL,
            vars,
            method = method
        ),
        NULL
    )

    if(!is.null(draws_obj$value)) {
        res_treatment_policy_MAR <- impute_analyze_pool(
            draws_obj$value,
            NULL,
            references,
            vars,
            visit_analysis = levels_visit[J]
        )
        res_treatment_policy_MAR$failed <- FALSE
        res_treatment_policy_MAR$error <- NULL
    } else {
        res_treatment_policy_MAR <- list(
            failed = TRUE,
            error = draws_obj$error
        )
    }


    ################# is_post_start_dopa, is_post_trt_stop
    vars_time_var <- vars
    vars_time_var$covariates <- c("x_bl*visit", "group*visit", "is_post_start_dopa*group", "is_post_trt_stop*group")

    draws_obj <- tryCatch_WEM(
        expr = draws(
            data,
            data_ice = NULL,
            vars_time_var,
            method = method
        ),
        NULL
    )

    if(!is.null(draws_obj$value)) {
        res_treatment_policy_is_post_trt_stop <- impute_analyze_pool(
            draws_obj$value,
            NULL,
            references,
            vars_time_var,
            visit_analysis = levels_visit[J]
        )
        res_treatment_policy_is_post_trt_stop$failed <- FALSE
        res_treatment_policy_is_post_trt_stop$error <- NULL
    } else {
        res_treatment_policy_is_post_trt_stop <- list(
            failed = TRUE,
            error = draws_obj$error
        )
    }

    ################# is_post_start_dopa, time_from_trt_stop
    vars_time_var <- vars
    vars_time_var$covariates <- c("x_bl*visit", "group*visit", "is_post_start_dopa*group", "time_from_trt_stop*group")

    draws_obj <- tryCatch_WEM(
        expr = draws(
            data,
            data_ice = NULL,
            vars_time_var,
            method = method
        ),
        NULL
    )

    if(!is.null(draws_obj$value)) {
        res_treatment_policy_time_from_trt_stop <- impute_analyze_pool(
            draws_obj$value,
            NULL,
            references,
            vars_time_var,
            visit_analysis = levels_visit[J]
        )
        res_treatment_policy_time_from_trt_stop$failed <- FALSE
        res_treatment_policy_time_from_trt_stop$error <- NULL
    } else {
        res_treatment_policy_time_from_trt_stop <- list(
            failed = TRUE,
            error = draws_obj$error
        )
    }



    ########################### HYPOTHETICAL ESTIMAND ##############################

    print("Running Bayesian multiple imputation (hypothetical estimand)")
    vars$outcome <- "y_MAR"
    ########################### MMRM ESTIMATOR

    draws_obj <- tryCatch_WEM(
        expr = draws(
            data,
            data_ice = NULL,
            vars,
            method = method
        ),
        NULL
    )

    if(!is.null(draws_obj$value)) {
        res_hypothetical_MAR <- impute_analyze_pool(
            draws_obj$value,
            NULL,
            references,
            vars,
            visit_analysis = levels_visit[J]
        )
        res_hypothetical_MAR$failed <- FALSE
        res_hypothetical_MAR$error <- NULL
    } else {
        res_hypothetical_MAR <- list(
            failed = TRUE,
            error = draws_obj$error
        )
    }

    res_obj <- list(
        "estimators_mixed" = list(
            res_CIR = res_mixed_CIR,
            res_JTR = res_mixed_JR,
            res_MAR = res_mixed_MAR,
            res_is_post_trt_stop = res_mixed_is_post_trt_stop,
            res_time_from_trt_stop = res_mixed_time_from_trt_stop
        ),
        "estimators_treatment_policy" = list(
            res_MAR = res_treatment_policy_MAR,
            res_is_post_trt_stop = res_treatment_policy_is_post_trt_stop,
            res_time_from_trt_stop = res_treatment_policy_time_from_trt_stop
        ),
        "estimators_hypothetical" = list(
            res_MAR = res_hypothetical_MAR
        )
    )

    return(res_obj)
}
