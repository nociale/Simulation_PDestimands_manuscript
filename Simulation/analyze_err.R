dir.create("Simulation/Results", showWarnings = FALSE)

# N <- seq(from = 75, to = 300, by = 25)
# dir.create("Simulation/err_files", showWarnings = FALSE)
# err_files <- list.files(path = "Simulation/out", pattern = ".err")
# out_files <- list.files(path = "Simulation/out", pattern = ".out")
# rds_files <- list.files(path = "Simulation/out", pattern = ".rds")
# N_sim <- length(err_files)
#
# rds_files_mod <- rds_files
# for(n in N) {
#     rds_files_mod <- gsub(paste0(n, "_.rds"), "", rds_files_mod)
# }
#
# failed_jobs <- c()
# for(i in seq.int(N_sim)) {
#
#     if(length(rds_files_mod[grep(paste0("_",i,"_"), rds_files_mod)]) < length(N)) {
#         failed_jobs <- c(failed_jobs, i)
#         file.copy(from = paste0("Simulation/out/", err_files[grep(paste0("_",i,".err"), err_files)]), to = "Simulation/err_files")
#         file.copy(from = paste0("Simulation/out/", out_files[grep(paste0("_",i,".out"), out_files)]), to = "Simulation/err_files")
#     }
# }
#
# # reproduce an example
# input_files <- list.files(path = "Simulation/in")
# data_failed <- input_files[grep(paste0("_",failed_jobs[2],"_"), input_files)][1]
# data_failed <- readRDS(paste0("Simulation/in/",data_failed))
#
# res <- run_simul(data_failed)


results <- readRDS("Simulation/Results/results.rds")
failed <- lapply(results, function(p) lapply(p, function(x) sapply(x, function(y) sapply(y, function(z)  sapply(z, function(w) w$failed)))))
failed <- lapply(failed, function(p) lapply(p, function(x) apply(x, 1, function(x) rowSums(as.data.frame(x)))))

error_messages <- lapply(results,
                         function(p)
                             lapply(p,
                                    function(x)
                                            sapply(x,
                                                   function(y)
                                                       as.data.frame(lapply(y,
                                                              function(z)
                                                                  lapply(z, function(w)
                                                                  {
                                                                  if(w$failed) {
                                                                      return(w$error)
                                                                  } else{
                                                                      return(NA)
                                                                  }
                                                              }
                                                       )
                                            ))
                             )
)
)

error_messages <- lapply(error_messages, function(p) lapply(p, function(y) apply(y, 1, function(z) {z <- z[!is.na(z)]; return(unlist(z))})))
saveRDS(failed, file = "Simulation/Results/res_failures.rds")
saveRDS(error_messages, file = "Simulation/Results/error_messages.rds")
