#!/bin/bash
#SBATCH --job-name=rfd4_s5_if
#SBATCH --partition=gpuhe.4h
#SBATCH --gpus=nvidia_geforce_rtx_4090:1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=03:00:00
#SBATCH --output=/cluster/scratch/gomest/hackaton/enzyme-design-hackathon-2026/rfd4/enzdes_tutorial/slurm-%j.out

eval "$(conda shell.bash hook)"
conda activate pd_hackaton

export HF_HOME="/cluster/scratch/gomest/.cache/huggingface"
export HF_HUB_OFFLINE=1

BASE_DIR="/cluster/scratch/gomest/hackaton/enzyme-design-hackathon-2026/rfd4/enzdes_tutorial"
BEST_CIF=$(ls "${BASE_DIR}/runs/scaffold_n25/best/rank1_"*.cif* | head -n 1)
OUT_DIR="${BASE_DIR}/runs/inverse_folding_n10_ll2.0"

rfd4-generate \
  inference_experiment=inference_on_user_inputs \
  generation=inverse_folding \
  inputs="${BEST_CIF}" \
  num_replicates=10 \
  batch_size=1 \
  generation.model.local_latents.simulation_step_params.tsr_k=2.0 \
  metrics_tags="[design,3d_coordinate_conditioning]" \
  out_dir="${OUT_DIR}" \
  run_name_=s3invfold
