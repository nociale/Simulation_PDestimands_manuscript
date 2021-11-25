include_compliance_indicator <- function(data) {

    data <- data %>%
        group_by(patnum) %>%
        mutate(time_from_trt_stop = cumsum(visit > trt_stop_visit),
               time_from_dopa = cumsum(visit > dopa_start_visit),
               is_post_trt_stop = ifelse(time_from_trt_stop > 0, 1, 0),
               is_post_start_dopa = ifelse(time_from_dopa > 0, 1, 0)
        )

    return(data)
}

convert_vars <- function(data, numeric_vars, factor_vars) {

    data[,factor_vars] <- as.data.frame(lapply(data[,factor_vars], as.factor))
    data[,numeric_vars] <- as.data.frame(lapply(data[,numeric_vars], as.numeric))

    return(data)
}

compute_change_bl <- function(data) {

    # remove baseline
    data <- data[data$visit != "0",]
    data$visit <- factor(data$visit, levels = levels(data$visit)[-1])

    # compute change from baseline
    data[,c("y", "y_dropout", "y_MAR")] <- data[,c("y", "y_dropout", "y_MAR")] - data$x_bl

    return(data)
}

generate_data <- function(i, N, prob_drop_out) {

    # set-up simulation
    n1 <- n2 <- N

    time <- seq(0,12,by=2)

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

    # mean trajectory intervention
    mu1 <- 30+6/12*time

    # slope after treatment stop
    stop_slope2 <- 10/12
    stop_slope1 <- 10/12

    # parameters for stopping treatment
    model_stop_trt_intervention <-        ~1 + I(visit >= 0)
    model_coef_stop_trt_intervention <-  c(0, log(0.03/(1-0.03)))
    model_stop_trt_control <-        ~1 + I(visit >= 0)
    model_coef_stop_trt_control <-  c(0,  log(0.02/(1-0.02)))

    # parameters for start of dopaminergic treatment
    model_start_dopa <-        ~1+I(visit==0)+I(visit%in%c(1,2))+  I(visit>2)+      I((x-30)/10)
    model_coef_start_dopa <-  c(0,      -1E6,   log(0.025/0.975),log(0.075/(1-0.075)), log(1.5))

    # Simulate trajectories
    data <- simulate_ices_and_adjust_trajectories(data = simple_trajectories(mu2 = mu2, Sigma1 = Sigma,
                                                                             mu1 = mu1, Sigma2 = Sigma,
                                                                             n1 = n1, n2 = n2, time = time),
                                                  mu2 = mu2,
                                                  mu1 = mu1,
                                                  model_stop_trt_intervention = model_stop_trt_intervention,
                                                  model_coef_stop_trt_intervention = model_coef_stop_trt_intervention,
                                                  model_stop_trt_control = model_stop_trt_control,
                                                  model_coef_stop_trt_control = model_coef_stop_trt_control,
                                                  model_start_dopa = model_start_dopa,
                                                  model_coef_start_dopa = model_coef_start_dopa,
                                                  prob_drop_out = prob_drop_out,
                                                  stop_slope2 = stop_slope2,
                                                  stop_slope1 = stop_slope1,
                                                  dopa_init_change = dopa_init_change,
                                                  dopa_slope = dopa_slope)


    # data for pure treatment policy estimand
    data$y_dropout <- do.call(c, lapply(split(data[,c("visit","drop_out_visit","y")], data$patnum),
                                        function(x)  {
                                            x$y[x$visit > x$drop_out_visit] = NA; return(x$y)
                                        }))

    # data for hybrid estimand
    data$y <- do.call(c, lapply(split(data[,c("visit","dopa_start_visit","drop_out_visit","y")], data$patnum),
                                function(x)  {
                                    x$y[x$visit > x$dopa_start_visit | x$visit > x$drop_out_visit] = NA; return(x$y)
                                }))


    # under MAR and JTR we should consider as missing those values in the intervention arm between treatment stop and start of dopa
    # data for pure hypothetical estimand
    data$y_MAR <- data$y
    data$y_MAR <- do.call(c, lapply(split(data[,c("visit","trt_stop_visit","y")], data$patnum),
                                    function(x)  {
                                        x$y[x$visit > x$trt_stop_visit] = NA; return(x$y)
                                    }))


    data <- include_compliance_indicator(data)
    factor_vars <- c("patnum", "group", "visit")
    numeric_vars <- c("x_bl", "y", "y_dropout", "y_MAR", "is_post_trt_stop", "time_from_trt_stop", "is_post_start_dopa")
    data <- convert_vars(
        data = data,
        numeric_vars = numeric_vars,
        factor_vars = factor_vars
    )

    data <- compute_change_bl(data)

    selected_vars <- c("patnum", "group", "visit", "x_bl", "is_post_trt_stop", "time_from_trt_stop", "is_post_start_dopa", "trt_stop_visit", "y", "y_MAR", "y_dropout")
    data <- as.data.frame(data[,selected_vars])

    saveRDS(data, paste0("Simulation/in/data_", i,"_", N, "_", prob_drop_out*100,".rds"))

}
