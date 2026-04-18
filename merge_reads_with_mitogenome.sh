#!/bin/bash -le

#SBATCH -A naiss2025-22-1155
#SBATCH --job-name=merge_mito
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --output=slurm_output/slurm_merge_mito_%j.out
#SBATCH --error=slurm_output/slurm_merge_mito_%j.out

################################################################################
# Merge simulated reads (with bacteria) with mitogenome reads
#
# Mitogenome reads per library for each sample:
# sim-1: 100 reads per library (300 total)
# sim-2: 1000 reads per library (3000 total)
# sim-3: 2000 reads per library (6000 total)
# sim-4: 5000 reads per library (15000 total)
# sim-5: 10000 reads per library (30000 total)
################################################################################

### CONFIGURATION
original_dir="/cfs/klemming/projects/snic/sllstore2017093/dnaharvester/2_test_run_sim/1_data/simulated_dataset_dnaharvester_with_bacteria"
mito_dir="/cfs/klemming/projects/snic/sllstore2017093/dnaharvester/2_test_run_sim/1_data/simulated_mitogenome_reads"
output_dir="/cfs/klemming/projects/snic/sllstore2017093/dnaharvester/2_test_run_sim/1_data/simulated_dataset_dnaharvester_with_bacteria_and_mitogenome"
output_filename="simulated_data_with_bacteria_and_mitogenome_samplesheet.csv"
temp_dir="${output_dir}/temp"

ml seqtk

mkdir -p "${output_dir}"
mkdir -p "${temp_dir}"

# Mitogenome read files
mito_r1="${mito_dir}/mammoth_mitogenome_s1.fq.gz"
mito_r2="${mito_dir}/mammoth_mitogenome_s2.fq.gz"

# Define mitogenome reads per library for each sample
declare -A mito_per_lib
mito_per_lib[1]=100     # sim-1: 100 reads per library
mito_per_lib[2]=1000    # sim-2: 1000 reads per library
mito_per_lib[3]=2000    # sim-3: 2000 reads per library
mito_per_lib[4]=5000    # sim-4: 5000 reads per library
mito_per_lib[5]=10000   # sim-5: 10000 reads per library

# Function to randomly sample N paired reads from mitogenome files
sample_mito_reads() {
    local n_reads=$1
    local out_r1=$2
    local out_r2=$3
    local seed=$4

    # Use seqtk for efficient random sampling
    seqtk sample -s${seed} "${mito_r1}" ${n_reads} | gzip > "${out_r1}"
    seqtk sample -s${seed} "${mito_r2}" ${n_reads} | gzip > "${out_r2}"
}

# Process each sample (sim-1 to sim-5)
for sample_num in 1 2 3 4 5; do
    sample="sim-${sample_num}"
    mkdir -p "${output_dir}/${sample}"

    n_mito=${mito_per_lib[${sample_num}]}
    total_mito=$((n_mito * 3))

    echo "Processing ${sample}: ${n_mito} mitogenome reads per library (${total_mito} total)"

    for lib_num in 1 2 3; do
        echo "  Processing ${sample} library ${lib_num}..."

        # Original files (with bacteria)
        orig_r1="${original_dir}/${sample}/${sample}_${lib_num}_s1.fq.gz"
        orig_r2="${original_dir}/${sample}/${sample}_${lib_num}_s2.fq.gz"

        # Temporary mitogenome sample files
        temp_mito_r1="${temp_dir}/${sample}_${lib_num}_mito_r1.fq.gz"
        temp_mito_r2="${temp_dir}/${sample}_${lib_num}_mito_r2.fq.gz"

        # Output files
        out_r1="${output_dir}/${sample}/${sample}_${lib_num}_s1.fq.gz"
        out_r2="${output_dir}/${sample}/${sample}_${lib_num}_s2.fq.gz"

        # Sample mitogenome reads with unique seed per sample/library
        seed=$((sample_num * 100 + lib_num))
        sample_mito_reads ${n_mito} "${temp_mito_r1}" "${temp_mito_r2}" ${seed}

        # Merge original + mitogenome reads
        zcat "${orig_r1}" "${temp_mito_r1}" | gzip > "${out_r1}"
        zcat "${orig_r2}" "${temp_mito_r2}" | gzip > "${out_r2}"

        # Clean up temp mitogenome files
        rm -f "${temp_mito_r1}" "${temp_mito_r2}"
    done
done

# Clean up temp directory
rm -rf "${temp_dir}"

echo "Merging complete!"

### GENERATE SAMPLESHEET
echo "sample_id,library_id,lane,sample_type,library_type,flowcell_id,seq_platform,fastq_1,fastq_2" > \
    "${output_dir}/${output_filename}"

# sim-1: Modern (3 libraries)
echo "sim-1,lib-1,L001,modern,double-stranded,FC12345,Illumina,${output_dir}/sim-1/sim-1_1_s1.fq.gz,${output_dir}/sim-1/sim-1_1_s2.fq.gz" >> "${output_dir}/${output_filename}"
echo "sim-1,lib-1,L002,modern,double-stranded,FC12345,Illumina,${output_dir}/sim-1/sim-1_2_s1.fq.gz,${output_dir}/sim-1/sim-1_2_s2.fq.gz" >> "${output_dir}/${output_filename}"
echo "sim-1,lib-2,L001,modern,double-stranded,FC12345,Illumina,${output_dir}/sim-1/sim-1_3_s1.fq.gz,${output_dir}/sim-1/sim-1_3_s2.fq.gz" >> "${output_dir}/${output_filename}"

# sim-2: Ancient (3 libraries)
echo "sim-2,lib-1,L001,ancient,double-stranded,FC12345,Illumina,${output_dir}/sim-2/sim-2_1_s1.fq.gz,${output_dir}/sim-2/sim-2_1_s2.fq.gz" >> "${output_dir}/${output_filename}"
echo "sim-2,lib-1,L002,ancient,double-stranded,FC12345,Illumina,${output_dir}/sim-2/sim-2_2_s1.fq.gz,${output_dir}/sim-2/sim-2_2_s2.fq.gz" >> "${output_dir}/${output_filename}"
echo "sim-2,lib-2,L001,ancient,double-stranded,FC12345,Illumina,${output_dir}/sim-2/sim-2_3_s1.fq.gz,${output_dir}/sim-2/sim-2_3_s2.fq.gz" >> "${output_dir}/${output_filename}"

# sim-3: Ancient (3 libraries)
echo "sim-3,lib-1,L001,ancient,double-stranded,FC12345,Illumina,${output_dir}/sim-3/sim-3_1_s1.fq.gz,${output_dir}/sim-3/sim-3_1_s2.fq.gz" >> "${output_dir}/${output_filename}"
echo "sim-3,lib-1,L002,ancient,double-stranded,FC12345,Illumina,${output_dir}/sim-3/sim-3_2_s1.fq.gz,${output_dir}/sim-3/sim-3_2_s2.fq.gz" >> "${output_dir}/${output_filename}"
echo "sim-3,lib-2,L001,ancient,double-stranded,FC12345,Illumina,${output_dir}/sim-3/sim-3_3_s1.fq.gz,${output_dir}/sim-3/sim-3_3_s2.fq.gz" >> "${output_dir}/${output_filename}"

# sim-4: Ancient (3 libraries)
echo "sim-4,lib-1,L001,ancient,double-stranded,FC12345,Illumina,${output_dir}/sim-4/sim-4_1_s1.fq.gz,${output_dir}/sim-4/sim-4_1_s2.fq.gz" >> "${output_dir}/${output_filename}"
echo "sim-4,lib-1,L002,ancient,double-stranded,FC12345,Illumina,${output_dir}/sim-4/sim-4_2_s1.fq.gz,${output_dir}/sim-4/sim-4_2_s2.fq.gz" >> "${output_dir}/${output_filename}"
echo "sim-4,lib-2,L001,ancient,double-stranded,FC12345,Illumina,${output_dir}/sim-4/sim-4_3_s1.fq.gz,${output_dir}/sim-4/sim-4_3_s2.fq.gz" >> "${output_dir}/${output_filename}"

# sim-5: Ancient (3 libraries)
echo "sim-5,lib-1,L001,ancient,double-stranded,FC12345,Illumina,${output_dir}/sim-5/sim-5_1_s1.fq.gz,${output_dir}/sim-5/sim-5_1_s2.fq.gz" >> "${output_dir}/${output_filename}"
echo "sim-5,lib-1,L002,ancient,double-stranded,FC12345,Illumina,${output_dir}/sim-5/sim-5_2_s1.fq.gz,${output_dir}/sim-5/sim-5_2_s2.fq.gz" >> "${output_dir}/${output_filename}"
echo "sim-5,lib-2,L001,ancient,double-stranded,FC12345,Illumina,${output_dir}/sim-5/sim-5_3_s1.fq.gz,${output_dir}/sim-5/sim-5_3_s2.fq.gz" >> "${output_dir}/${output_filename}"

echo "Samplesheet created: ${output_dir}/${output_filename}"
