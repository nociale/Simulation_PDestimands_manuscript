#BSUB -J pd_simul[1-10000]
#BSUB -q preempt
#BSUB -n 1
#BSUB -M 1000
#BSUB -W 03:00
singularity exec https://ross.science.roche.com/singularity-cbs/githubroche/rinstallation/rockerplus_4.1.1.sif Rscript --vanilla Simulation/run_dropout50.R
