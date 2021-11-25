# ------ Simulate endpoint trajectories without intercurrent events ----

# mu2, Sigma2: mean and covariance matrix of control group [including baseline]
# mu1, Sigma1: mean and covariance matrix of intervention group 
# n1, n2     : Sample size of control (n2) and intervention (n1) group 
# time       : Vector of time points [must be of same length mu2]

simple_trajectories <- function(mu2, Sigma2, mu1 = mu2, Sigma1 = Sigma2, n2, n1, time = NULL){

  m <- length(mu2) # number of visits per subject (including baseline)
  if (is.null(time)) time <- 0:(m-1)
  # Control group
  r2 <- data.frame("patnum" = rep(1:n2,each=m),
                   "group"    = "Control",
                   "visit"  = rep(0:(m-1),n2),
                   "time"   = rep(time,n2),
                   "x_bl"   =  NA,
                   "x"      = c(t(mvrnorm(n=n2,mu=mu2,Sigma=Sigma2))))
  # Intervention group
  r1 <- data.frame("patnum" = n2+rep(1:n1,each=m),
                   "group"    = "Intervention",
                   "visit"  = rep(0:(m-1),n1),
                   "time"   = rep(time,n1),
                   "x_bl"   = NA,
                   "x"      = c(t(mvrnorm(n=n1,mu=mu1,Sigma=Sigma1))))
  # Pool both groups and add baseline
  r <- rbind(r2,r1)
  r$x_bl <- rep(r$x[time==0],each=m)
  r$group <- factor(r$group,levels=c("Control","Intervention"))
  return(r)
}

# ------ Simulate time of treatment stop and initiation of dopaminergic treatment ------

# - probability of treatment stop at each visit is simulated according to a logistic model 
#   (arguments model and model_coef)
# - number of visits from treatment stop to initation of dopaminergic treatment is simulated 
#   according to a simple geometric distribution with probability prob_dopa 
#   (default for prob_dopa is 1 indicating initiation of dopaminergic treatment at the same visit as 
#    treatment stop)

simulate_ices_and_adjust_trajectories <- function(data,
                                                  mu2,
                                                  mu1,
                                                  model_stop_trt_intervention,
                                                  model_coef_stop_trt_intervention,
                                                  model_stop_trt_control,
                                                  model_coef_stop_trt_control,
                                                  model_start_dopa,
                                                  model_coef_start_dopa,
                                                  prob_drop_out = 0.5,
                                                  stop_slope2 = 0,
                                                  stop_slope1 = 0,
                                                  dopa_init_change = 0,
                                                  dopa_slope = 0){
  
  # extend numeric slope changes to vectors
  visits <- levels(as.factor(data$visit))
  J = length(visits)
  if (length(stop_slope2)==1) { stop_slope2 <- rep(stop_slope2,J-1) }
  if (length(stop_slope1)==1) { stop_slope1 <- rep(stop_slope1,J-1) }
  if (length(dopa_slope)==1)  { dopa_slope <- rep(dopa_slope,J-1) }
  if ((length(stop_slope2)!=J-1)|(length(stop_slope1)!=J-1)|
      (length(dopa_slope)!=J-1)) { error("Argument stop_slope1 or 2 or dopa_slope_change has wrong length")}
  
  
  data$stop_trt <- 0
  data$drop_out <- 0
  data$start_dopa <- 0
  data$trt_stop_visit <- Inf
  data$dopa_start_visit <- Inf
  data$drop_out_visit <- Inf
  
  data$y <- data$x
  
  visits = as.numeric(visits)
  
  for(j in visits) {
    data_j <- data[data$visit == j,] # data related to current visit
    
    if(j > 0) {
      data_previous <- data[data$visit == j-1,] # data related to previous visit
    } else {
      data_previous <- data_j
    }
    
    time_increment <- data_j$time - data_previous$time
    
    is_intervention <- data_j$group == "Intervention"
    is_control <- data_j$group == "Control"
    
    
    #### adjust trajectory due to discontinuation from randomized treatment ####
    
    curr_index_change <- which(visits == j-1) # j-1 because if j = 1 we want this index to be 1 etc..
    
    intervention_stop <- which(is_intervention & data_previous$trt_stop_visit <= j-1)
    control_stop <- which(is_control & data_previous$trt_stop_visit <= j-1)
    
    data_j$y[intervention_stop] <-
      data_previous$y[intervention_stop] + (data_j$x[intervention_stop] - mu1[visits == j]) - (data_previous$x[intervention_stop] - mu1[visits == j-1]) + time_increment[intervention_stop]*stop_slope1[curr_index_change]
    
    data_j$y[control_stop] <- 
      data_previous$y[control_stop] + (data_j$x[control_stop] - mu2[visits == j]) - (data_previous$x[control_stop] - mu2[visits == j-1]) + time_increment[control_stop]*stop_slope2[curr_index_change]
    
    
    #### adjust trajectories due to dopaminergic treatment ####
    
    start_dopa_previous_intervention <- which(data_previous$dopa_start_visit == j-1 & is_intervention)
    data_j$y[start_dopa_previous_intervention] <-
      data_previous$y[start_dopa_previous_intervention] + (data_j$x[start_dopa_previous_intervention] - mu1[visits == j]) - (data_previous$x[start_dopa_previous_intervention] - mu1[visits == j-1]) + dopa_init_change[start_dopa_previous_intervention]
    
    start_dopa_previous_control <- which(data_previous$dopa_start_visit == j-1 & is_control)
    data_j$y[start_dopa_previous_control] <-
      data_previous$y[start_dopa_previous_control] + (data_j$x[start_dopa_previous_control] - mu2[visits == j]) - (data_previous$x[start_dopa_previous_control] - mu2[visits == j-1]) + dopa_init_change[start_dopa_previous_control]
    
    if(j > 1) {
      
      start_dopa_previous2_intervention <- which(data_previous$dopa_start_visit <= j-2 & is_intervention)
      data_j$y[start_dopa_previous2_intervention] <-
        data_previous$y[start_dopa_previous2_intervention] + (data_j$x[start_dopa_previous2_intervention] - mu1[visits == j]) - (data_previous$x[start_dopa_previous2_intervention] - mu1[visits == j-1]) + time_increment[start_dopa_previous2_intervention]*dopa_slope[curr_index_change]
      
      start_dopa_previous2_control <- which(data_previous$dopa_start_visit <= j-2 & is_control)
      data_j$y[start_dopa_previous2_control] <-
        data_previous$y[start_dopa_previous2_control] + (data_j$x[start_dopa_previous2_control] - mu2[visits == j]) - (data_previous$x[start_dopa_previous2_control] - mu2[visits == j-1]) + time_increment[start_dopa_previous2_control]*dopa_slope[curr_index_change]
      
    }
    
    
    #### simulate discontinuation from randomized treatment ####
    
    no_trt_stop_previous <- data_previous$stop_trt == 0
    data_j$stop_trt <- data_previous$stop_trt # those who stopped the treatment previously will have stopped the treatment currently
    
    # intervention arm
    if(sum(is_intervention & no_trt_stop_previous) > 0) {
      lp <- c(model.matrix(model_stop_trt_intervention,data_j[is_intervention & no_trt_stop_previous, ])%*%model_coef_stop_trt_intervention) # linear predictor
      prob_stop_trt <- binomial(link = "logit")$linkinv(lp)
      data_j$stop_trt[is_intervention & no_trt_stop_previous] <- rbinom(n=length(lp),size=1,prob=prob_stop_trt)
    }
    
    # control arm
    if(sum(is_control & no_trt_stop_previous) > 0) {
      lp <- c(model.matrix(model_stop_trt_control,data_j[is_control & no_trt_stop_previous, ])%*%model_coef_stop_trt_control) # linear predictor
      prob_stop_trt <- binomial(link = "logit")$linkinv(lp)
      data_j$stop_trt[is_control & no_trt_stop_previous] <- rbinom(n=length(lp),size=1,prob=prob_stop_trt)
    }
    
    # update trt_stop_visit
    data_j$trt_stop_visit <- data_previous$trt_stop_visit
    data_j$trt_stop_visit[no_trt_stop_previous & data_j$stop_trt == 1] <- j # set to current visit the visit of trt stop for those patients who stop the treatment at this visit
    
    
    #### simulate drop-out from the study ####
    
    data_j$drop_out <- data_previous$drop_out
    data_j$drop_out[data_j$trt_stop_visit == j] <- rbinom(n=sum(data_j$trt_stop_visit == j),size=1,prob=prob_drop_out)
    
    data_j$drop_out_visit <- data_previous$drop_out_visit
    data_j$drop_out_visit[data_j$trt_stop_visit == j & data_j$drop_out == 1] <- j
    
    
    #### simulate start of rescue medication ####
    no_start_dopa_previous <- data_previous$start_dopa == 0
    data_j$start_dopa <- data_previous$start_dopa # those who started dopa previously will be on dopa currently
    
    if (sum(no_start_dopa_previous) > 0) {
      lp <- c(model.matrix(model_start_dopa, data_j[no_start_dopa_previous, ])%*%model_coef_start_dopa) # linear predictor
      prob_start_dopa <- binomial(link = "logit")$linkinv(lp)
      data_j$start_dopa[no_start_dopa_previous] <- rbinom(n=length(lp),size=1,prob=prob_start_dopa)
    }
    
    data_j$dopa_start_visit <- data_previous$dopa_start_visit
    data_j$dopa_start_visit[no_start_dopa_previous & data_j$start_dopa==1] <- j # set to current visit the visit of trt stop for those patients who start dopa at this visit
    
    
    #### merge into data ####
    data[data$visit == j,] <- data_j
  }
  
  data$stop_trt <- NULL # drop variable again ("relevant" variable is trt_stop_visit)
  data$start_dopa <- NULL # drop variable again ("relevant" variable is dopa_start_visit)
  data$drop_out <- NULL # drop variable again ("relevant" variable is drop_out_visit)
  
  data <- data %>%
    group_by(patnum) %>%
    mutate(trt_stop_visit = trt_stop_visit[J],
           dopa_start_visit = dopa_start_visit[J],
           drop_out_visit = drop_out_visit[J]
           #,y = ifelse(visit > drop_out_visit, NA, y)
    )
  
  # return dataset
  return(data)
  
}
