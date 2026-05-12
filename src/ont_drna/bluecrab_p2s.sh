#!/bin/bash
#SBATCH -p gpu
#SBATCH -c 4
#SBATCH -t 16:00:00
#SBATCH --gpus-per-node=2g.20gb:1
#SBATCH --mem=64G
#SBATCH --gres=tmpspace:64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<email>
#SBATCH --error=<data_filepath>blow5/bluecrab_p2s_%j.err
#SBATCH --output=<data_filepath>blow5/bluecrab_p2s_%j.out
set -euo pipefail

#pod5_dir=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <scripts_path>pod5_dirs_input.txt)
#pod5_dirname=$(basename "${pod5_dir}")
#output_dir="<data_filepath>blow5/${pod5_dirname}/"

pod5_dir=`realpath $1`
output_dir=`realpath $2`

#pod5_dir="<scratch_data_path>pod5_pilot2_fus_stressed/"
pod5_dirname=$(basename "${pod5_dir}")
#output_dir="<data_filepath>pilot2_fus_stressed_complete/blow5/"

mkdir -p ${output_dir}

echo "Start bluecrab p2s ${pod5_dirname}"


set +u; env - APPTAINER_TMPDIR="$TMPDIR" APPTAINER_CACHEDIR=<singularity_cache_path> \
singularity exec  -B $TMPDIR:$TMPDIR -B /hpc:/hpc  docker://ghcr.io/kennalab/blue-crab:0.5.0 \
blue-crab p2s \
    -c zlib \
    ${pod5_dir} -d ${output_dir} \
    --threads 4

set +u; env - APPTAINER_TMPDIR="$TMPDIR" APPTAINER_CACHEDIR=<singularity_cache_path> \
singularity exec  -B $TMPDIR:$TMPDIR -B /hpc:/hpc  docker://ghcr.io/kennalab/blue-crab:0.5.0 \
blue-crab --version > ${output_dir}/versions.yml

echo "Completed bluecrab ${pod5_dirname}"
