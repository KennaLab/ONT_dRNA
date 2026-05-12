#!/bin/bash
#SBATCH --job-name=modkit_pileup
#SBATCH --output=<data_filepath>modkit_pileup_%j.out
#SBATCH --error=<data_filepath>modkit_pileup_%j.err
#SBATCH -c 8
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --gres=tmpspace:32G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<email>
set -euo pipefail

input_bam=`realpath $1`
output_dir=`realpath $2`
prefix=$3
ref=`realpath $4`

echo "${input_bam}"
echo "${output_dir}"
echo "${prefix}"
echo "${ref}"

mkdir -p $output_dir && cd $output_dir

set +u; env - APPTAINER_TMPDIR="$TMPDIR" APPTAINER_CACHEDIR=<singularity_cache_path> \
singularity exec  -B $TMPDIR:$TMPDIR -B /hpc:/hpc  docker://quay.io/biocontainers/ont-modkit:0.6.1--hcdda2d0_0 \
modkit pileup \
  --modified-bases m6A \
  --prefix "${prefix}" \
  --threads 8 \
  --log-filepath "${output_dir}/${prefix}.modkit_pileup.log" \
  --suppress-progress \
  --reference "${ref}" \
  "${input_bam}" \
  "${output_dir}/${prefix}.pileup.bed"

set +u; env - APPTAINER_TMPDIR="$TMPDIR" APPTAINER_CACHEDIR=<singularity_cache_path> \
singularity exec  -B $TMPDIR:$TMPDIR -B /hpc:/hpc  docker://quay.io/biocontainers/ont-modkit:0.6.1--hcdda2d0_0 \
modkit --version > "${output_dir}/versions.txt"
