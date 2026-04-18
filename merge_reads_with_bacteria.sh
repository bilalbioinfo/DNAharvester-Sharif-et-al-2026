#!/bin/bash -le

#SBATCH -A naiss2025-22-1155
#SBATCH --job-name=merge_reads
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --output=slurm_output/slurm_merge_%j.out
#SBATCH --error=slurm_output/slurm_merge_%j.out

################################################################################
# Merge simulated reads with bacterial reads (randomly sampled from pool)
#
# Target bacterial reads per sample:
# Sim-1: 45,000 (15,000 per library)
# Sim-2: 50,000 (16,667 per library)
# Sim-3: 895,000 (298,334 per library)
# Sim-4: 499,000 (166,334 per library)
# Sim-5: 500,000 (166,667 per library)
################################################################################

### CONFIGURATION
original_dir="/cfs/klemming/projects/snic/sllstore2017093/dnaharvester/2_test_run_sim/1_data/simulated_dataset_dnaharvester"
bacterial_dir="/cfs/klemming/projects/snic/sllstore2017093/dnaharvester/2_test_run_sim/1_data/simulated_bacterial_reads"
output_dir="/cfs/klemming/projects/snic/sllstore2017093/dnaharvester/2_test_run_sim/1_data/simulated_dataset_dnaharvester_with_bacteria"
output_filename="simulated_data_with_bacteria_samplesheet.csv"
temp_dir="${output_dir}/temp"

ml seqtk

mkdir -p "${output_dir}"
mkdir -p "${temp_dir}"
cd "${output_dir}"

pooled_r1="/cfs/klemming/projects/snic/sllstore2017093/dnaharvester/2_test_run_sim/1_data/simulated_bacterial_reads/gtdb_GCF_reads_s1.fq.gz"
pooled_r2="/cfs/klemming/projects/snic/sllstore2017093/dnaharvester/2_test_run_sim/1_data/simulated_bacterial_reads/gtdb_GCF_reads_s2.fq.gz"

# Function to randomly sample N paired reads from pooled files
sample_bacterial_reads() {
    local n_reads=$1
    local out_r1=$2
    local out_r2=$3
    local seed=$4

    # Use seqtk for efficient random sampling
    seqtk sample -s${seed} "${pooled_r1}" ${n_reads} | gzip > "${out_r1}"
    seqtk sample -s${seed} "${pooled_r2}" ${n_reads} | gzip > "${out_r2}"
}

# Define bacterial reads per sample (total, will be divided by 3 for libraries)
declare -A bact_reads
bact_reads[1]=45000    # sim-1: 45k total
bact_reads[2]=50000    # sim-2: 50k total
bact_reads[3]=895000   # sim-3: 895k total
bact_reads[4]=499000   # sim-4: 499k total
bact_reads[5]=500000   # sim-5: 500k total

# Process each sample
for sample_num in 1 2 3 4 5; do
    sample="sim-${sample_num}"
    mkdir -p "${output_dir}/${sample}"

    total_bact=${bact_reads[${sample_num}]}
    # Divide into 3 libraries (last one gets remainder)
    bact_per_lib=$((total_bact / 3))
    bact_lib3=$((total_bact - 2 * bact_per_lib))

    echo "Processing ${sample}: ${total_bact} bacterial reads (${bact_per_lib}, ${bact_per_lib}, ${bact_lib3} per library)"

    for lib_num in 1 2 3; do
        echo "  Processing ${sample} library ${lib_num}..."

        # Determine bacterial reads for this library
        if [ ${lib_num} -eq 3 ]; then
            n_bact=${bact_lib3}
        else
            n_bact=${bact_per_lib}
        fi

        # Original files
        orig_r1="${original_dir}/${sample}/${sample}_${lib_num}_s1.fq.gz"
        orig_r2="${original_dir}/${sample}/${sample}_${lib_num}_s2.fq.gz"

        # Temporary bacterial sample files
        temp_bact_r1="${temp_dir}/${sample}_${lib_num}_bact_r1.fq.gz"
        temp_bact_r2="${temp_dir}/${sample}_${lib_num}_bact_r2.fq.gz"

        # Output files
        out_r1="${output_dir}/${sample}/${sample}_${lib_num}_s1.fq.gz"
        out_r2="${output_dir}/${sample}/${sample}_${lib_num}_s2.fq.gz"

        # Sample bacterial reads with unique seed per sample/library
        seed=$((sample_num * 100 + lib_num))
        sample_bacterial_reads ${n_bact} "${temp_bact_r1}" "${temp_bact_r2}" ${seed}

        # Merge original + bacterial reads
        zcat "${orig_r1}" "${temp_bact_r1}" | gzip > "${out_r1}"
        zcat "${orig_r2}" "${temp_bact_r2}" | gzip > "${out_r2}"

        # Clean up temp bacterial files
        rm -f "${temp_bact_r1}" "${temp_bact_r2}"
    done
done

# Clean up pooled files
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
