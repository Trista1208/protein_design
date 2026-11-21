#!/bin/bash
#SBATCH --job-name=rfd4_s2
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

THEOZYME="/cluster/scratch/gomest/hackaton/enzyme-design-hackathon-2026/rfd4/enzdes_tutorial/M0732_1xs1_theozyme.cif"
OUT_DIR="/cluster/scratch/gomest/hackaton/enzyme-design-hackathon-2026/rfd4/enzdes_tutorial/runs/scaffold_n25"

rfd4-generate \
  inference_experiment=inference_on_user_inputs \
  inputs="${THEOZYME}" \
  num_replicates=25 \
  batch_size=1 \
  generation.n_recycle=2 \
  metrics_tags="[design,unindexed,3d_coordinate_conditioning]" \
  out_dir="${OUT_DIR}" \
  run_name_=s2design
