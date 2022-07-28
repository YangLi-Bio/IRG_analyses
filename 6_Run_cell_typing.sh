#!/bin/bash
#SBATCH --time=7:50:59
#SBATCH --output=6_Cell_typing.out
#SBATCH --account=PCON0022
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=100GB
#SBATCH --gpus-per-node=1

set -e

module load R/4.1.0-gnu9.1

cd /fs/ess/PCON0022/liyang/STREAM/case_1_IRG/Codes
Rscript 6_Cell_typing.R
