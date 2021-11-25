dir.create("Simulation/Results", showWarnings = FALSE)
N <- seq(from = 75, to = 300, by = 25)
prob_drop_out <- c(25, 50)
dir.create("Simulation/Results", showWarnings = FALSE)
all <- list.files(path='Simulation/out', pattern = '*.rds$', full.names=TRUE)
results <- sapply(prob_drop_out, function(p) sapply(N, function(n) lapply(all[grep(paste0('_', n, "_", p, '_.rds'), all)], readRDS), simplify = FALSE), simplify = FALSE)
names(results) <- prob_drop_out
results <- lapply(results, function(x) {names(x) <- N; return(x)})
saveRDS(results, file = "Simulation/Results/results.rds")
