#!/bin/bash
#SBATCH -c 16
#SBATCH -t 24:00:00
#SBATCH --mem=64G
#SBATCH --gres=tmpspace:64G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<email>
#SBATCH --error=<data_filepath>slurm-m6anet_%j.err
#SBATCH --output=<data_filepath>slurm-m6anet_%j.out
set -euo pipefail

prefix=$1
blow5=`realpath $2`
reads_fastq=`realpath $3`
bam=`realpath $4`
fasta=`realpath $5`
event_align=`realpath $6`
output_dir=`realpath $7`

#prefix="pilot2_fus_control_subset"
#output_dir="<data_filepath>m6anet/pilot2_fus_control_subset_v2/"
#blow5="<data_filepath>input_nfcore_nanoseq/pilot2_fus_control_subset/blow5/PBG84424_f1d488ee_43206fa6_0.blow5"
#reads_fastq="<data_filepath>input_nfcore_nanoseq/pilot2_fus_control_subset/fastq/PBG84424_pass_f1d488ee_43206fa6_0.fastq.gz"
#bam="<data_filepath>input_nfcore_nanoseq/pilot2_fus_control_subset/bam/fus_control_subset_R1.sorted.bam"
#fasta="<ref_genomes_path>GRCh38_gencode_v22_CTAT_lib_Mar012021/GRCh38_gencode_v22_CTAT_lib_Mar012021.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa"
#event_align="<data_filepath>f5c/pilot2_fus_control_subset/pilot2_fus_control_subset.events.tsv"
blow5_name=$(basename "${blow5}" .blow5)

mkdir -p "${output_dir}/dataprep_${blow5_name}/" "${output_dir}/inference_${blow5_name}/"

echo "Start m6anet dataprep ${prefix} ${blow5_name}"
echo "Start m6anet dataprep ${prefix} ${blow5_name}" 1>&2
set +u; env - APPTAINER_TMPDIR="$TMPDIR" APPTAINER_CACHEDIR=<singularity_cache_path> \
singularity exec  -B $TMPDIR:$TMPDIR -B /hpc:/hpc  docker://quay.io/biocontainers/m6anet:2.1.0--pyhdfd78af_0 \
m6anet dataprep --eventalign ${event_align} \
--out_dir "${output_dir}/dataprep_${blow5_name}/" \
--n_processes 16
echo "Completed m6anet dataprep ${prefix} ${blow5_name}"
echo "Completed m6anet dataprep ${prefix} ${blow5_name}" 1>&2

echo "Start m6anet inference ${prefix}  ${blow5_name}"
echo "Start m6anet inference ${prefix}  ${blow5_name}" 1>&2
set +u; env - APPTAINER_TMPDIR="$TMPDIR" APPTAINER_CACHEDIR=<singularity_cache_path> \
singularity exec  -B $TMPDIR:$TMPDIR -B /hpc:/hpc  docker://quay.io/biocontainers/m6anet:2.1.0--pyhdfd78af_0 \
m6anet inference \
--input_dir "${output_dir}/dataprep_${blow5_name}/" \
--out_dir "${output_dir}/inference_${blow5_name}/" \
--n_processes 16 \
--pretrained_model HEK293T_RNA004 \
--num_iterations 1000
echo "Completed m6anet inference ${prefix} ${blow5_name}"
echo "Completed m6anet inference ${prefix} ${blow5_name}" 1>&2

set +u; env - APPTAINER_TMPDIR="$TMPDIR" APPTAINER_CACHEDIR=<singularity_cache_path> \
singularity exec  -B $TMPDIR:$TMPDIR -B /hpc:/hpc  docker://quay.io/biocontainers/m6anet:2.1.0--pyhdfd78af_0 \
m6anet --version > ${output_dir}/versions.yml

echo "Completed m6anet ${prefix} ${blow5_name}"
echo "Completed m6anet ${prefix} ${blow5_name}" 1>&2
