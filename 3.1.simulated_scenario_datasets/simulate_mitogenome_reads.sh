#!/bin/bash -le

#SBATCH -A naiss2025-22-1155
#SBATCH --job-name=gargammel
#SBATCH --partition=shared
#SBATCH --cpus-per-task=8
#SBATCH --time=1:00:00
#SBATCH --output=slurm_output/slurm_%A_%a.out
#SBATCH --error=slurm_output/slurm_%A_%a.out


################################################################################
# Generate simulated aDNA test data using gargammel
# Creates 20 samples with 100K reads each (runs as parallel array jobs)
# Ancient samples include DNA damage simulation
################################################################################

### CONFIGURATION
num_reads=100000
damage="-damage 0.03,0.4,0.01,0.3"
length_command="--loc 3.5 --scale 0.35 --minsize 30 --maxsize 100 -rl 150"
gargammel_input="/cfs/klemming/projects/snic/sllstore2017093/dnaharvester/2_test_run_sim/1_data/gargammel_input_3/"
output_dir="/cfs/klemming/projects/snic/sllstore2017093/dnaharvester/2_test_run_sim/1_data/simulated_mitogenome_reads"
output_filename="simulated_data_samplesheet.csv"



source ~/.bashrc
load_conda
conda activate gargammel

mkdir -p "${output_dir}"
cd "${output_dir}"
### Use SLURM array task ID for unique sample names
output_prefix="mammoth_mitogenome"
mkdir -p ${output_prefix}
gargammel -n ${num_reads} --comp 0,0,1 ${length_command} ${damage} \
    -o ${output_dir}/${output_prefix} ${gargammel_input}
rm ${output_dir}/*fa.gz
