dir.create("Simulation/Results", showWarnings = FALSE)

compute_power_type1error <- function(pvalues, alpha = 0.05) {
    pvalues <- pvalues[!is.na(pvalues)]
    return((sum(pvalues <= 0.05))/length(pvalues))
}

compute_rmse <- function(a, b) {
    return(sqrt(mean((a - b)^2, na.rm = TRUE)))
}

extract_statistics <- function(res_simul, estimand, methods, statistics, truth = NULL) {

    stats_data <- list()
    for(i in 1:length(methods)) {
        stats_data[[i]] <- lapply(
            res_simul,
            function(x) x[[estimand]][[methods[i]]]$pars$trt_6[statistics]
        )
        stats_data[[i]] <- as.data.frame(t(sapply(
            stats_data[[i]],
            function(x) {
                if(!is.null(x)) {
                    return(unlist(x))
                } else{
                    ret <- rep(NA, length(statistics))
                    names(ret) <- names(statistics)
                    return(ret)
                }
            }
        )))
    }
    names(stats_data) <- methods
    stats_data$truth <- truth

    return(stats_data)

}

summary_stats <- function(res_simul, n, H0, truth = NULL, remove_wrong_est = TRUE) {

    methods <- list("estimators_mixed" =
                        c("res_CIR", "res_JTR", "res_MAR", "res_is_post_trt_stop", "res_time_from_trt_stop"),
                    "estimators_treatment_policy" =
                        c("res_MAR", "res_is_post_trt_stop", "res_time_from_trt_stop"),
                    "estimators_hypothetical" =
                        c("res_MAR")
    )

    statistics <- c("est", "se", "pvalue")

    estimands <- names(methods)
    stats_data <- sapply(
        estimands,
        function(x) extract_statistics(res_simul, x, methods[[x]], statistics, truth[[x]]),
        simplify = FALSE,
        USE.NAMES = TRUE
    )

    summary_stats <- sapply(
        estimands,
        function(x) {

            sapply(
                methods[[x]],
                function(method) {

                    y <- stats_data[[x]][[method]]

                    # some estimates are around +-Inf. Remove them from the analysis
                    if(remove_wrong_est) {
                        wrong_est <- y$est < -10e5 | y$est > 10e5
                    } else {
                        wrong_est <- rep(FALSE, length(y$est))
                    }

                    ret_obj <- list("mean_est" = mean(y$est[!wrong_est], na.rm = TRUE),
                                    "sd_est" = sd(y$est[!wrong_est], na.rm = TRUE),
                                    "mean_se" = mean(y$se[!wrong_est], na.rm = TRUE),
                                    "type1error" = compute_power_type1error(y$pvalue[!wrong_est])
                    )

                    if(remove_wrong_est) {
                        ret_obj$wrong_est <- y$est[which(wrong_est)]
                        names(ret_obj$wrong_est) <- which(wrong_est)
                    }

                    # adjust name of type1error/power according to H0
                    if(!H0) names(ret_obj)[which(names(ret_obj) == "type1error")] = "power"

                    if(!is.null(truth)) ret_obj$rmse <- compute_rmse(y$est[!wrong_est], stats_data[[x]]$truth)

                    return(ret_obj)
                },
                simplify = FALSE,
                USE.NAMES = TRUE)
        },
        simplify = FALSE,
        USE.NAMES = TRUE)

    return(summary_stats)

}

# working directory: "Parkinson_disease_simul/Analysis_for_manuscript/Code"
res_simul <- readRDS("Simulation/Results/results.rds")
res_truth <- readRDS("True_values/Results/res_table.rds")
res_failures <- readRDS("Simulation/Results/res_failures.rds")

n_visits <- nrow(res_truth)
truth <- list(
    estimators_mixed = -res_truth$diff_mixed[n_visits],
    estimators_treatment_policy = -res_truth$diff_tp[n_visits],
    estimators_hypothetical = -res_truth$diff_hypoth[n_visits]
)

prob_drop_out <- as.numeric(names(res_simul))
N <- as.numeric(names(res_simul[[1]]))

res <- lapply(1:length(prob_drop_out), function(p) lapply(1:length(N), function(i) summary_stats(res_simul[[p]][[i]], N[i], H0 = FALSE, truth = truth, remove_wrong_est = TRUE)))
names(res) <- prob_drop_out
res <- lapply(res, function(x) {names(x) <- N; return(x)})

# add to res_failures very large estimates
res_failures <- lapply(1:length(prob_drop_out), function(p) lapply(1:length(N), function(n) lapply(names(truth), function(e) {
    curr_fail <- res_failures[[p]][[n]][[e]]
    curr_wrong <- res[[p]][[n]][[e]]
    for(method in names(curr_fail)) {
        curr_fail[[method]] <- curr_fail[[method]] + length(curr_wrong[[method]]$wrong_est)
    }
    return(curr_fail)
}
)))
names(res_failures) <- prob_drop_out
res_failures <- lapply(res_failures, function(x) {names(x) <- N; return(x)})
res_failures <- lapply(res_failures, function(p) lapply(p, function(n) { names(n) <- names(truth); return(n)}))

saveRDS(res, file = "Simulation/Results/table_res.rds")
saveRDS(res_failures, file = "Simulation/Results/res_failures_final.rds")
