library(MASS)
library(rbmi)
library(dplyr)
source("Simulation/run_simul.R")
source("Simulation/generate_data.R")
source("simul_functions/simul_functions.R")
source("Simulation/tryCatch-WEM.R")

# handle input/output
dir.create("Simulation/in", showWarnings = FALSE)
dir.create('Simulation/out', showWarnings = FALSE)

i <- Sys.getenv('LSB_JOBINDEX')
j <- Sys.getenv('LSB_JOBID')

N <- seq(from = 75, to = 300, by = 25)
prob_drop_out <- c(0.25)

for(n in N) {
    generate_data(i, n, prob_drop_out)
    result <- run_simul(readRDS(paste0("Simulation/in/data_", i,"_", n, "_", prob_drop_out*100,".rds")))
    saveRDS(result, paste('Simulation/out/out', j, i, n, prob_drop_out*100, '.rds', sep = '_'))
}
