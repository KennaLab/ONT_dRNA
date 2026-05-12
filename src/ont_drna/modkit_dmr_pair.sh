#!/bin/bash
#SBATCH --job-name=modkit_dmr_pair
#SBATCH --output=<data_filepath>modkit_dmr_pair_%j.out
#SBATCH --error=<data_filepath>modkit_dmr_pair_%j.err
#SBATCH -c 8
#SBATCH --time=04:00:00
#SBATCH --mem=24G
#SBATCH --gres=tmpspace:24G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<email>
set -euo pipefail

custom_params=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <scripts_path>modkit_dmr_input.txt)
IFS=',' read -r -a array_params <<< "$custom_params"
type="${array_params[0]}"
prefix="${array_params[1]}"
modkit_params="${array_params[@]:2}"

#type="pair"
#prefix="stressed_vs_control"
#input_dir="<processed_output_path>"
output_dir="<processed_output_path>modkit_dmr_${type}.${prefix}/"
mkdir -p $output_dir && cd $output_dir

echo "Start modkit dmr ${type} comparing the condition ${prefix}."
echo "Following parameters are added: ${modkit_params}"
set +u; env - APPTAINER_TMPDIR="$TMPDIR" APPTAINER_CACHEDIR=<singularity_cache_path> \
singularity exec  -B $TMPDIR:$TMPDIR -B /hpc:/hpc  docker://quay.io/biocontainers/ont-modkit:0.6.1--hcdda2d0_0 \
modkit dmr pair \
  ${modkit_params} \
  -o ${output_dir}/${prefix}.dmr_${type}.bed \
  --header \
  --ref <references_path>epi2me_alignment_ref_genome_gencode_v49/GRCh38.p14.genome.fa \
  --base A \
  -t 8 \
  --force \
  --log-filepath "${output_dir}/${prefix}.${type}.log"

#  --segment "${output_dir}/${prefix}.single_base_segments.${type}.bed" \
#  --fine-grained \
