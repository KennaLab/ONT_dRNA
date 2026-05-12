#!/bin/bash
#SBATCH -c 4
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<email>
#SBATCH --error=<data_filepath>fast5/slurm-pod5_to_fast5_%j.err
#SBATCH --output=<data_filepath>fast5/slurm-pod5_to_fast5_%j.out

set -euo pipefail

# Set input and output dirs
line_input_output=$(sed "${SLURM_ARRAY_TASK_ID}q;d" <data_filepath>convert_pod5_input_output.csv)

input_pod5="$(cut -d';' -f1 <<< ${line_input_output})"
output_dir="$(cut -d';' -f2 <<< ${line_input_output})"

echo "Input ${input_pod5}"
echo "Output ${output_dir}"

# Create environment in $TMPDIR
mkdir -p $TMPDIR/pod5_file_format_${SLURM_ARRAY_TASK_ID}/
cd $TMPDIR/pod5_file_format_${SLURM_ARRAY_TASK_ID}/
python3.12 -m venv venv
. venv/bin/activate
pip install --upgrade pip
pip install pod5==0.3.35


# Convert pod5 to fast5
pod5 convert to_fast5 "${input_pod5}" --output "${output_dir}" --threads 4 --force-overwrite

# Exit environment
deactivate
rm -rf $TMPDIR/pod5_file_format_${SLURM_ARRAY_TASK_ID}/

echo "Completed pod5 convert to fast5 ${input_pod5}."
