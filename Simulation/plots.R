library(ggplot2)

dir.create("Simulation/Plots", showWarnings = FALSE)

data <- readRDS("Simulation/Results/table_res.rds")
res_truth <- readRDS("True_values/Results/res_table.rds")
res_failures <- readRDS("Simulation/Results/res_failures_final.rds")

n_visits <- nrow(res_truth)
truth <- c(
    mixed = res_truth$diff_mixed[n_visits],
    treatment_policy = res_truth$diff_tp[n_visits],
    hypothetical = res_truth$diff_hypoth[n_visits]
)

prob_drop_out <- as.numeric(names(data))
N <- as.numeric(names(data[[1]]))
estimators <- sapply(data[[1]][[1]], names)
estimators$estimators_treatment_policy <- c("res_MAR", "res_is_post_trt_stop_dopa", "res_time_from_trt_stop_dopa")
estimands <- gsub("estimators_", "", names(estimators))
n_estimands <- length(estimands)
n_estimators_per_estimand <- sapply(data[[1]][[1]], length)
n_estimators <- sum(n_estimators_per_estimand)

data_power <- unlist(lapply(data, function(p) lapply(p, function(x) lapply(x, function(y) lapply(y, function(z) z$power)))))
data_failures <- unlist(res_failures)
data_est <- abs(unlist(lapply(data, function(p) lapply(p, function(x) lapply(x, function(y) lapply(y, function(z) z$mean_est))))))
data_truth <- rep(rep(unlist(sapply(1:n_estimands, function(i) rep(truth[i], n_estimators_per_estimand[i]))), length(N)), length(prob_drop_out))
data_plot <- data.frame(
    "prob_drop_out" = rep(prob_drop_out, each = n_estimators*length(N)),
    "sample_size" = rep(rep(N, each = n_estimators),length(prob_drop_out)),
    "estimand" = rep(rep(unlist(sapply(1:n_estimands, function(i) rep(estimands[i], n_estimators_per_estimand[i]))), length(N)),length(prob_drop_out)),
    "estimator" = rep(rep(gsub("res_", "", unlist(estimators)), length(N)),length(prob_drop_out)),
    "power" = data_power,
    "failures" = data_failures,
    "mean_est" = data_est,
    "truth" = data_truth,
    "bias" = data_est - data_truth,
    "mean_se" = unlist(lapply(data, function(p) lapply(p, function(x) lapply(x, function(y) lapply(y, function(z) z$mean_se))))),
    "rmse" = unlist(lapply(data, function(p) lapply(p, function(x) lapply(x, function(y) lapply(y, function(z) z$rmse)))))
)

for(dropout in prob_drop_out) {

    curr_data <- subset(data_plot, data_plot$prob_drop_out == dropout)

    png(filename = paste0("Simulation/Plots/power_", dropout, ".png"),
        width = 1100,
        height = 500)
    p <- ggplot(
        data = curr_data,
        mapping = aes(x = sample_size, y = power, group = estimator, color = estimator)
    ) +
        geom_line() +
        geom_point() +
        facet_wrap(facets = ~ estimand)
    print(p)
    dev.off()

    ############ plot number of failures for each scenario varying the sample size

    png(filename = paste0("Simulation/Plots/failures_", dropout, ".png"),
        width = 1100,
        height = 500)
    p <- ggplot(
        data = curr_data,
        mapping = aes(x = sample_size, y = failures, group = estimator, color = estimator)
    ) +
        geom_line() +
        geom_point() +
        facet_wrap(facets = ~ estimand)
    print(p)
    dev.off()

    for(n in N) {

        curr_data <- subset(data_plot, data_plot$sample_size == n & data_plot$prob_drop_out == dropout)

        curr_data$sample_size <- paste0("sample size: ", curr_data$sample_size)
        png(filename = paste0("Simulation/Plots/bias_",n, "_", dropout, ".png"),
            width = 1100,
            height = 500)
        p <- ggplot(
            data = curr_data,
            mapping = aes(x = estimator, y = bias, group = estimator, fill = estimator)
        ) +
            geom_col() +
            theme(axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) +
            facet_wrap(facets = ~ estimand*sample_size)
        print(p)
        dev.off()

        png(filename = paste0("Simulation/Plots/se_",n, "_", dropout, ".png"),
            width = 1100,
            height = 500)
        p <- ggplot(
            data = curr_data,
            mapping = aes(x = estimator, y = mean_se, group = estimator, fill = estimator)
        ) +
            geom_col() +
            theme(axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) +
            facet_wrap(facets = ~ estimand*sample_size)
        print(p)
        dev.off()

        png(filename = paste0("Simulation/Plots/rmse_",n, "_", dropout, ".png"),
            width = 1100,
            height = 500)
        p <- ggplot(
            data = curr_data,
            mapping = aes(x = estimator, y = rmse, group = estimator, fill = estimator)
        ) +
            geom_col() +
            theme(axis.text.x=element_blank(),
                  axis.ticks.x=element_blank()) +
            facet_wrap(facets = ~ estimand*sample_size)
        print(p)
        dev.off()
    }

}
