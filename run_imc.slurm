#!/bin/bash
#SBATCH --job-name=imc_scportrait
#SBATCH --output=logs/imc_scportrait_%j.out
#SBATCH --error=logs/imc_scportrait_%j.err
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=24:00:00

set -euo pipefail

module load anaconda/3 || true
source activate scportrait || conda activate scportrait

python C:/GitHub/scPortrait_to_IMC/imc_to_single_cells.py \
  --channels-dir /scratch/$USER/imc_run/images \
  --mask-dir /scratch/$USER/imc_run/masks \
  --projects-root /scratch/$USER/imc_projects \
  --config C:/GitHub/scPortrait_to_IMC/config_imc.yml \
  --mask-expand-px 2 \
  --overwrite
