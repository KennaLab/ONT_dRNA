#!/bin/bash
#SBATCH -c 12
#SBATCH -t 08:00:00
#SBATCH --mem=64G
#SBATCH --gres=tmpspace:64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<email>
#SBATCH --error=<data_filepath>slow5tools_merge_%j.err
#SBATCH --output=<data_filepath>slow5tools_merge_%j.out
set -euo pipefail

blow5_dir=`realpath $1`
output_dir=`realpath $2`
prefix=$3
#blow5_dir=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <scripts_path>blow5_dirs_input.txt)
#blow5_dirname=$(basename "${blow5_dir}")
#output_dir="<data_filepath>blow5/merged/${blow5_dirname}/"

echo "Input blow5 in ${blow5_dir}:"
echo "$(ls -1 ${blow5_dir}/*.blow5)"
echo "Output directory ${output_dir}"

echo "Input blow5 in ${blow5_dir}:" 1>&2
echo "$(ls -1 ${blow5_dir}/*.blow5)" 1>&2
echo "Output directory ${output_dir}" 1>&2

mkdir -p ${output_dir}

echo "Start slow5tools merge ${prefix}"
echo "Start slow5tools merge ${prefix}" 1>&2

set +u; env - APPTAINER_TMPDIR="$TMPDIR" APPTAINER_CACHEDIR=<singularity_cache_path> \
singularity exec  -B $TMPDIR:$TMPDIR -B /hpc:/hpc  docker://quay.io/biocontainers/slow5tools:1.3.0--hee927d3_2 \
    slow5tools merge \
    ${blow5_dir} -o ${output_dir}/${prefix}.blow5 \
    --threads 12

echo "Completed slow5tools merge ${prefix}"
echo "Completed slow5tools merge ${prefix}" 1>&2

set +u; env - APPTAINER_TMPDIR="$TMPDIR" APPTAINER_CACHEDIR=<singularity_cache_path> \
singularity exec  -B $TMPDIR:$TMPDIR -B /hpc:/hpc  docker://quay.io/biocontainers/slow5tools:1.3.0--hee927d3_2 \
slow5tools --version > ${output_dir}/versions.yml
