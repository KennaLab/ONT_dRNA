#!/bin/bash
#SBATCH --job-name=methylartist_segmeth
#SBATCH --output=<data_filepath>methylartist_segmeth_%j.out
#SBATCH --error=<data_filepath>methylartist_segmeth_%j.err
#SBATCH -c 8
#SBATCH --time=04:00:00
#SBATCH --mem=24G
#SBATCH --gres=tmpspace:24G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<email>
set -euo pipefail

output_dir="<processed_output_path>methylartist_segment_genes/"
mkdir -p $output_dir && cd $output_dir

set +u; env - APPTAINER_TMPDIR="$TMPDIR" APPTAINER_CACHEDIR=<singularity_cache_path> \
singularity  exec -B $TMPDIR:$TMPDIR -B /hpc:/hpc  docker://quay.io/biocontainers/methylartist:1.5.3--pyhdfd78af_0 \
methylartist segmeth \
  -q 9 \
  -i <references_path>gencode_human_v49/gencode.v49_genes.bed \
  -b <processed_output_path>modkit_pileup_pilot1_ico_control/pilot1_ico_control.pileup.bed.gz,<processed_output_path>modkit_pileup_pilot1_ico_stressed/pilot1_ico_stressed.pileup.bed.gz,<processed_output_path>modkit_pileup_pilot2_co_control/pilot2_co_control.pileup.bed.gz,<processed_output_path>modkit_pileup_pilot2_co_stressed/pilot2_co_stressed.pileup.bed.gz,<processed_output_path>modkit_pileup_pilot2_fus_control/pilot2_fus_control.pileup.bed.gz,<processed_output_path>modkit_pileup_pilot2_fus_stressed/pilot2_fus_stressed.pileup.bed.gz \
  --bedmethyl \
  -o pilot1_and_pilot2.genetypes.segmeth.tsv \
  --predict_dmr \
  -p 8
