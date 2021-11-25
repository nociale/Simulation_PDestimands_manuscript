# working directory: "Parkinson_disease_simul/Analysis_for_manuscript/Code"

dir.create("True_values", showWarnings = FALSE)
dir.create("True_values/Plots", showWarnings = FALSE)
dir.create("True_values/Plots/mean_trajectories", showWarnings = FALSE)
dir.create("True_values/Results", showWarnings = FALSE)

source("simul_functions/simul_functions.R")
# "TRUE" SIMULATION RESULTS
# simulate data with large number of subjects per arm and recover the "true" values
# under treatment policy and "CIR" based assumptions

library(dplyr)
library(MASS)
library(ggplot2)

#----- Parameter which are identical across scenarios

set.seed(1000)
n1 <- n2 <- 3E6 # sample size per group

time <- seq(0,12,by=2)
J <- length(time)

# create Sigma
sd_intercept <- 10
sd_slope <- 5
cor_slope_inter <- 0.5

sd_error <- 6

covRE <- matrix(c(sd_intercept^2,cor_slope_inter*sd_intercept*sd_slope,
                  cor_slope_inter*sd_intercept*sd_slope,sd_slope^2),ncol=2)
Sigma <- cbind(1,time/12)%*%covRE%*%rbind(1,time/12)+diag(sd_error^2,nrow=length(time))

# mean trajectory control
mu2 <- 30+10/12*time

# random drop if subjects initiates dopaminergic treatment
dopa_init_change <- rbeta(n1+n2,shape1=2,shape2=1.5)*25-25
dopa_slope <- 0

# probability of drop-out
prob_drop_out <- 0.5

# mean trajectory intervention: 40% relative reduction
mu1 <- 30+6/12*time

# change in slope after treatment stop
stop_slope2 <- 10/12
stop_slope1 <- 10/12

res <- list()
scenarios <- list()
methods <- list()

for (scenario in 1:2){

    if (scenario==1){
        scenarioName <- "High_unequal"
        method = "CIR"

        # parameters for stopping treatment
        model_stop_trt_intervention <-        ~1 + I(visit >= 0)
        model_coef_stop_trt_intervention <-  c(0, log(0.03/(1-0.03)))
        model_stop_trt_control <-        ~1 + I(visit >= 0)
        model_coef_stop_trt_control <-  c(0,  log(0.02/(1-0.02)))

        # parameters for start of dopaminergic treatment
        model_start_dopa <-      ~1
        model_coef_start_dopa <- -1E6
    }

    if (scenario==2){
        scenarioName <- "High_unequal"
        method = "T.P"

        # parameters for stopping treatment
        model_stop_trt_intervention <-        ~1 + I(visit >= 0)
        model_coef_stop_trt_intervention <-  c(0, log(0.03/(1-0.03)))
        model_stop_trt_control <-        ~1 + I(visit >= 0)
        model_coef_stop_trt_control <-  c(0,  log(0.02/(1-0.02)))

        # parameters for start of dopaminergic treatment
        model_start_dopa <-        ~1+I(visit==0)+I(visit%in%c(1,2))+  I(visit>2)+      I((x-30)/10)
        model_coef_start_dopa <-  c(0,      -1E6,   log(0.025/0.975),log(0.075/(1-0.075)), log(1.5))
    }

    # simulate data

    data <- simple_trajectories(mu2, Sigma, mu1, Sigma, n2, n1, time)

    data <- simulate_ices_and_adjust_trajectories(data,
                                                  mu2,
                                                  mu1,
                                                  model_stop_trt_intervention,
                                                  model_coef_stop_trt_intervention,
                                                  model_stop_trt_control,
                                                  model_coef_stop_trt_control,
                                                  model_start_dopa,
                                                  model_coef_start_dopa,
                                                  prob_drop_out,
                                                  stop_slope2,
                                                  stop_slope1,
                                                  dopa_init_change,
                                                  dopa_slope)


    if(scenario == 2) {
        prob_obs_post_trt_stop <- data %>%
            group_by(patnum, group) %>%
            summarise(has_obs_post_trt_stop = any(visit > trt_stop_visit & visit <= drop_out_visit & visit <= dopa_start_visit)) %>%
            group_by(group) %>%
            summarize(prob_obs_post_trt_stop = sum(has_obs_post_trt_stop)/length(has_obs_post_trt_stop))
    }

    d_tab <- data %>%
        mutate(y_change=y-x_bl, y_change_mar = x - x_bl) %>%
        group_by(group, time) %>%
        summarise(mean_trajectory = mean(y_change),
                  mean_trajectory_mar = mean(y_change_mar),
                  .groups = "drop")

    res <- rbind(res, as.data.frame(d_tab))
    scenarios <- c(scenarios, scenarioName)
    methods <- c(methods, method)
}

scenarios <- unlist(scenarios)
res$scenario <- rep(scenarios, each = 2*J)
methods <- unlist(methods)
res$arm <- rep(methods, each = 2*J)

res$is_dopa <- as.factor(c(rep("No", nrow(res)/2), rep("Yes", nrow(res)/2)))
#res$mean_trajectory[res$group == "Control" & res$is_dopa == "No"] <- mu2 - mu2[1]

res_intervention_mar <- res[res$group == "Control" & res$is_dopa == "No",]
res_intervention_mar$group <- "Intervention"
res_intervention_mar$mean_trajectory <- res$mean_trajectory_mar[res$group == "Intervention" & res$is_dopa == "No"]
res_intervention_mar$arm <- "Interv. (Hypoth)"
res <- rbind(res, res_intervention_mar)

res$arm[res$group == "Control" & res$is_dopa == "No"] <- "Control"
res$arm[res$group == "Control" & res$is_dopa == "Yes"] <- "Control (T.P)"
res$arm[res$arm == "CIR"] <- "Interv. (Mixed est)"
res$arm[res$arm == "T.P"] <- "Interv. (T.P)"

res$arm <- as.factor(res$arm)

cols <- c(RColorBrewer::brewer.pal(n = 4, name = "Set1"),
          RColorBrewer::brewer.pal(n = 4, name = "Blues"))

res$group <- NULL
colnames(res)[colnames(res) == "arm"] <- "group"
labels <- rep("", nrow(res))
labels[seq(J, nrow(res), by = J)] <- round(res$mean_trajectory[seq(J, nrow(res), by = J)],2)

for(scenario in scenarios) {

    data_scenario <- res[res$scenario == scenario,]
    labels_scenario <- labels[res$scenario == scenario]
    data_no_dopa <- data_scenario[data_scenario$is_dopa == "No",]
    labels_no_dopa <- labels_scenario[data_scenario$is_dopa == "No"]

    png(filename = paste0("True_values/Plots/mean_trajectories/",scenario,".png"),
        width = 500*(1+1/3),
        height = 400*(1+1/3),
        bg = "transparent")
    p <- ggplot(data_no_dopa, aes(x = time, y = mean_trajectory, group = group, linetype = group)) +
        geom_line(aes(col = group), size = 1.5) + geom_point(size = 3) +
        geom_text(label = labels_no_dopa, size=4.5, x = rep(time, 3), y = data_no_dopa$mean_trajectory, vjust = 2) +
        scale_color_manual(values=cols[c(1,4,8)]) +
        scale_linetype_manual(values = c("solid","solid","solid")) +
        scale_x_continuous(breaks=seq(0,12,by=2)) +
        theme(legend.justification = "bottom") +
        labs(title="Simulation setting", subtitle="Mean UPDRS trajectories", y="Mean change from baseline",
             x="Time (months)", caption = paste("Obtained with simulation with large sample size (1'000'000 patients per arm)", sep = ""))

    print(p)
    dev.off()

    png(filename = paste0("True_values/Plots/mean_trajectories/",scenario,"_dopa.png"),
        width = 500*(1+1/3),
        height = 400*(1+1/3),
        bg = "transparent")
    p <- ggplot(data_scenario, aes(x = time, y = mean_trajectory, group = group, linetype = group)) +
        geom_line(aes(col = group), size = 1.5) + geom_point(size = 3) +
        geom_text(label = labels_scenario, size=4.5, x = rep(time, 5), y = data_scenario$mean_trajectory, vjust = 2) +
        scale_color_manual(values=cols[c(1,1,4,8,3)]) +
        scale_linetype_manual(values = c("solid","dashed","solid","solid", "dashed")) +
        scale_x_continuous(breaks=seq(0,12,by=2)) +
        labs(title="Simulation setting", subtitle="Mean UPDRS trajectories", y="Mean change from baseline",
             x="Time (months)", caption = paste("Obtained with simulation with large sample size (1'000'000 patients per arm)", sep = "")) +
        theme(legend.justification = "bottom") +
        guides(fill = guide_legend(keywidth = 1, keyheight = 1),
               linetype=guide_legend(keywidth = 1, keyheight = 1),
               colour=guide_legend(keywidth = 3, keyheight = 1))
    print(p)
    dev.off()

}

data <- data.frame("group" = rep(c("Control", "Interv. (Hypoth)"), each = J),
                   "mean_trajectory" = c(mu2-mu2[1], mu1-mu1[1]),
                   "visit" = time)

labels <- rep("", nrow(data))
labels[seq(J, nrow(data), by = J)] <- round(data$mean_trajectory[seq(J, nrow(data), by = J)],2)

png(filename = paste0("True_values/Plots/mean_trajectories/initial_plot.png"),
    width = 500*(1+1/3),
    height = 400*(1+1/3),
    bg = "transparent")
p <- ggplot(data, aes(x = visit, y = mean_trajectory, group = group, linetype = group)) +
    geom_line(aes(col = group), size = 1.5) + geom_point(size = 3) +
    geom_text(label = labels, size=4.5, x = rep(time,2), y = data$mean_trajectory, vjust = 2) +
    scale_color_manual(values=cols[c(1,4)]) +
    scale_linetype_manual(values = c("solid","solid")) +
    scale_x_continuous(breaks=seq(0,12,by=2)) +
    theme(legend.justification = "bottom") +
    labs(title="Simulation setting", subtitle="Mean UPDRS trajectories", y="Mean change from baseline",
         x="Time (months)", caption = paste("Obtained with simulation with large sample size (1'000'000 patients per arm)", sep = ""))
print(p)
dev.off()

res_table <- data.frame(
    "time" = time,
    "control" = res$mean_trajectory[res$group == "Control"],
    "interv_hypoth" = res$mean_trajectory[res$group == "Interv. (Hypoth)"],
    "diff_hypoth" = res$mean_trajectory[res$group == "Control"] - res$mean_trajectory[res$group == "Interv. (Hypoth)"],
    "interv_mixed" = res$mean_trajectory[res$group == "Interv. (Mixed est)"],
    "diff_mixed" = res$mean_trajectory[res$group == "Control"] - res$mean_trajectory[res$group == "Interv. (Mixed est)"],
    "control_tp" = res$mean_trajectory[res$group == "Control (T.P)"],
    "interv_tp" = res$mean_trajectory[res$group == "Interv. (T.P)"],
    "diff_tp" = res$mean_trajectory[res$group == "Control (T.P)"] - res$mean_trajectory[res$group == "Interv. (T.P)"]
)

N <- seq(from = 75, to = 300, by = 25)
prob_failure_mmrm <- c()

for(n in N) {
    probs <- c(
        pbinom(q = 0, p = prob_obs_post_trt_stop$prob_obs_post_trt_stop, size = n),
        prod(pbinom(q = 0, p = prob_obs_post_trt_stop$prob_obs_post_trt_stop, size = n))
    )
    prob_failure_mmrm <- c(prob_failure_mmrm, as.numeric(c(1,1,-2)%*%probs))
}

names(prob_failure_mmrm) <- N
prob_failure <- data.frame(
    sample_size = N,
    prob_nopostdiscdata_in_at_least_one_arm = prob_failure_mmrm
)

png(filename = paste0("True_values/Plots/prob_nopostdiscdata_in_at_least_one_arm.png"),
    width = 500*(1+1/3),
    height = 400*(1+1/3)
)
p <- ggplot(data = prob_failure, aes(x = sample_size, y = prob_nopostdiscdata_in_at_least_one_arm*100)) +
    geom_line() +
    geom_point()
print(p)
dev.off()

saveRDS(prob_obs_post_trt_stop, "True_values/Results/prob_obs_post_trt_stop.rds")
saveRDS(prob_failure_mmrm, "True_values/Results/prob_nopostdiscdata_in_at_least_one_arm.rds")
saveRDS(res_table, "True_values/Results/res_table.rds")
