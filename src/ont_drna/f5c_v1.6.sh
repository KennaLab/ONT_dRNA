#!/bin/bash
#SBATCH -p gpu
#SBATCH -c 16
#SBATCH -t 05:00:00
#SBATCH --gpus-per-node=7g.79gb:1
#SBATCH --mem=36G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<email>
#SBATCH --error=<data_filepath>f5c_v1.6_%j.err
#SBATCH --output=<data_filepath>f5c_v1.6_%j.out
set -euo pipefail

blow5_file=`realpath $1`
reads_fastq=`realpath $2`
bam=`realpath $3`
fasta=`realpath $4`
output_dir=`realpath $5`

#blow5_file=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <data_filepath>blow5/input_pod5_pilot2_fus_control.txt)
blow5_name=$(basename "${blow5_file}" .blow5)
#output_dir="<data_filepath>f5c/pilot2_fus_control/"
#reads_fastq="<data_filepath>input_nfcore_nanoseq/pilot2_fus_control/fastq/PBG84424_pass_f1d488ee_43206fa6_0.fastq.gz"
#bam="<data_filepath>input_nfcore_nanoseq/pilot2_fus_control/bam/fus_control_R1.sorted.bam"
mkdir -p ${output_dir}

cd ${output_dir}

if [ ! -f "${blow5_file}.idx" ]; then
    echo "Start f5c index ${blow5_name}"
    echo "Start f5c index ${blow5_name}" 1>&2
    set +u; env - APPTAINER_TMPDIR="$TMPDIR" APPTAINER_CACHEDIR=<singularity_cache_path> \
    singularity exec  -B $TMPDIR:$TMPDIR -B /hpc:/hpc  docker://quay.io/biocontainers/f5c:1.6--hee927d3_0 \
        f5c index \
        --slow5 ${blow5_file} \
        ${reads_fastq} \
        --threads 16
    echo "Completed f5c index ${blow5_name}"
    echo "Completed f5c index ${blow5_name}" 1>&2
fi

echo "Start f5c eventalign ${blow5_name}"
echo "Start f5c eventalign ${blow5_name}" 1>&2

set +u; env - APPTAINER_TMPDIR="$TMPDIR" APPTAINER_CACHEDIR=<singularity_cache_path> \
singularity exec  -B $TMPDIR:$TMPDIR -B /hpc:/hpc  docker://quay.io/biocontainers/f5c:1.6--hee927d3_0 \
    f5c eventalign \
        --signal-index \
        --scale-events \
        --min-mapq 0 \
        --rna \
        --pore rna004 \
        --bam ${bam} \
        --genome ${fasta} \
        --reads ${reads_fastq} \
        --slow ${blow5_file} \
        --m6anet \
        --summary ${output_dir}/${blow5_name}.summary.txt \
        --threads 16 \
        > ${output_dir}/${blow5_name}.events.tsv

echo "Completed f5c eventalign ${blow5_name}"
echo "Completed f5c eventalign ${blow5_name}" 1>&2

set +u; env - APPTAINER_TMPDIR="$TMPDIR" APPTAINER_CACHEDIR=<singularity_cache_path> \
singularity exec  -B $TMPDIR:$TMPDIR -B /hpc:/hpc  docker://quay.io/biocontainers/f5c:1.6--hee927d3_0 \
f5c --version > ${output_dir}/${blow5_name}.versions.yml
