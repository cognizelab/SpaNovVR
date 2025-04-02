#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=64G
#SBATCH --partition=blade
#SBATCH --job-name=calc_ts
#SBATCH --error=log/job.%J.err
#SBATCH --output=log/job.%J.out
#SBATCH --mail-user=jaquent@fudan.edu.cn
#SBATCH --mail-type=ALL #Send all job status to the email address abov

# Load modules
module load R/4.2.0 

# Run the script
Rscript --vanilla Get_HC_2_cortex_conn_matrices.R
