#!/bin/bash

################################################################################
# Generate simulated aDNA test data using gargammel
# Creates 5 samples with varying Mammoth/Human/Yersinia compositions
# Ancient samples include DNA damage simulation
#
# Target reads per sample (before adding bacteria):
# Sim-1: Mammoth 900k, Yersinia 10k, Human 45k = 955k (+ 45k bacteria = 1M)
# Sim-2: Mammoth 50k, Yersinia 50k, Human 850k = 950k (+ 50k bacteria = 1M)
# Sim-3: Mammoth 10k, Yersinia 5k, Human 90k = 105k (+ 895k bacteria = 1M)
# Sim-4: Mammoth 1k, Yersinia 1k, Human 499k = 501k (+ 499k bacteria = 1M)
# Sim-5: Mammoth 0, Yersinia 0, Human 500k = 500k (+ 500k bacteria = 1M)
################################################################################

### CONFIGURATION
damage="-damage 0.03,0.4,0.01,0.3"
length_command_ancient="--loc 3.5 --scale 0.35 --minsize 25 --maxsize 100 -rl 150"
length_command_historical="--loc 4.6 --scale 0.5 --minsize 25 --maxsize 300 -rl 150"
length_command_modern="--loc 5.0 --scale 0.5 --minsize 25 --maxsize 300 -rl 150"
gargammel_input="/cfs/klemming/projects/snic/sllstore2017093/dnaharvester/2_test_run_sim/1_data/gargammel_input_1/"
output_dir="/cfs/klemming/projects/snic/sllstore2017093/dnaharvester/2_test_run_sim/1_data/simulated_dataset_dnaharvester/"
output_filename="simulated_data_samplesheet.csv"

mkdir -p "${output_dir}"
cd "${output_dir}"

## sim-1: Modern, 955k reads (Mammoth 900k, Yersinia 10k, Human 45k)
## comp order: Yersinia, Human, Mammoth
## Proportions: 10k/955k=0.01047, 45k/955k=0.04712, 900k/955k=0.94241
(
    output_prefix="sim-1"
    num_reads=318334  # 955k / 3 libraries
    mkdir -p ${output_prefix}
    for i in {1..3}; do
        gargammel -n ${num_reads} --comp 0.01047,0.04712,0.94241 ${length_command_modern} \
            -o ${output_prefix}/${output_prefix}_${i} ${gargammel_input}
    done
    rm -f ${output_prefix}/*fa.gz
) &

## sim-2: Historical, 950k reads (Mammoth 50k, Yersinia 50k, Human 850k)
## Proportions: 50k/950k=0.0526, 850k/950k=0.8948, 50k/950k=0.0526 (sum=1.0)
(
    output_prefix="sim-2"
    num_reads=316667  # 950k / 3 libraries
    mkdir -p ${output_prefix}
    for i in {1..3}; do
        gargammel -n ${num_reads} --comp 0.0526,0.8948,0.0526 ${damage} ${length_command_historical} \
            -o ${output_prefix}/${output_prefix}_${i} ${gargammel_input}
    done
    rm -f ${output_prefix}/*fa.gz
) &

### sim-3: Ancient, 105k reads (Mammoth 10k, Yersinia 5k, Human 90k)
### Proportions: 5k/105k=0.04762, 90k/105k=0.85714, 10k/105k=0.09524
(
    output_prefix="sim-3"
    num_reads=35000  # 105k / 3 libraries
    mkdir -p ${output_prefix}
    for i in {1..3}; do
        gargammel -n ${num_reads} --comp 0.04762,0.85714,0.09524 ${damage} ${length_command_ancient} \
            -o ${output_prefix}/${output_prefix}_${i} ${gargammel_input}
    done
    rm -f ${output_prefix}/*fa.gz
) &

### sim-4: Ancient, 501k reads (Mammoth 1k, Yersinia 1k, Human 499k)
### Proportions: 1k/501k=0.001996, 499k/501k=0.996008, 1k/501k=0.001996
(
    output_prefix="sim-4"
    num_reads=167000  # 501k / 3 libraries
    mkdir -p ${output_prefix}
    for i in {1..3}; do
        gargammel -n ${num_reads} --comp 0.001996,0.996008,0.001996 ${damage} ${length_command_ancient} \
            -o ${output_prefix}/${output_prefix}_${i} ${gargammel_input}
    done
    rm -f ${output_prefix}/*fa.gz
) &

### sim-5: Ancient, 500k reads (Mammoth 0, Yersinia 0, Human 500k)
### Proportions: 0, 1.0, 0
(
    output_prefix="sim-5"
    num_reads=166667  # 500k / 3 libraries
    mkdir -p ${output_prefix}
    for i in {1..3}; do
        gargammel -n ${num_reads} --comp 0,1.0,0 ${damage} ${length_command_ancient} \
            -o ${output_prefix}/${output_prefix}_${i} ${gargammel_input}
    done
    rm -f ${output_prefix}/*fa.gz
) &

wait
echo "Simulated data generation complete."

### GENERATE SAMPLESHEET
echo "sample_id,library_id,lane,sample_type,library_type,flowcell_id,seq_platform,fastq_1,fastq_2" > \
    "${output_dir}${output_filename}"

# sim-1: Modern (3 libraries)
echo "sim-1,lib-1,L001,modern,double-stranded,FC12345,Illumina,${PWD}/sim-1/sim-1_1_s1.fq.gz,${PWD}/sim-1/sim-1_1_s2.fq.gz" >> "${output_dir}${output_filename}"
echo "sim-1,lib-1,L002,modern,double-stranded,FC12345,Illumina,${PWD}/sim-1/sim-1_2_s1.fq.gz,${PWD}/sim-1/sim-1_2_s2.fq.gz" >> "${output_dir}${output_filename}"
echo "sim-1,lib-2,L001,modern,double-stranded,FC12345,Illumina,${PWD}/sim-1/sim-1_3_s1.fq.gz,${PWD}/sim-1/sim-1_3_s2.fq.gz" >> "${output_dir}${output_filename}"

# sim-2: Ancient (3 libraries)
echo "sim-2,lib-1,L001,ancient,double-stranded,FC12345,Illumina,${PWD}/sim-2/sim-2_1_s1.fq.gz,${PWD}/sim-2/sim-2_1_s2.fq.gz" >> "${output_dir}${output_filename}"
echo "sim-2,lib-1,L002,ancient,double-stranded,FC12345,Illumina,${PWD}/sim-2/sim-2_2_s1.fq.gz,${PWD}/sim-2/sim-2_2_s2.fq.gz" >> "${output_dir}${output_filename}"
echo "sim-2,lib-2,L001,ancient,double-stranded,FC12345,Illumina,${PWD}/sim-2/sim-2_3_s1.fq.gz,${PWD}/sim-2/sim-2_3_s2.fq.gz" >> "${output_dir}${output_filename}"

# sim-3: Ancient (3 libraries)
echo "sim-3,lib-1,L001,ancient,double-stranded,FC12345,Illumina,${PWD}/sim-3/sim-3_1_s1.fq.gz,${PWD}/sim-3/sim-3_1_s2.fq.gz" >> "${output_dir}${output_filename}"
echo "sim-3,lib-1,L002,ancient,double-stranded,FC12345,Illumina,${PWD}/sim-3/sim-3_2_s1.fq.gz,${PWD}/sim-3/sim-3_2_s2.fq.gz" >> "${output_dir}${output_filename}"
echo "sim-3,lib-2,L001,ancient,double-stranded,FC12345,Illumina,${PWD}/sim-3/sim-3_3_s1.fq.gz,${PWD}/sim-3/sim-3_3_s2.fq.gz" >> "${output_dir}${output_filename}"

# sim-4: Ancient (3 libraries)
echo "sim-4,lib-1,L001,ancient,double-stranded,FC12345,Illumina,${PWD}/sim-4/sim-4_1_s1.fq.gz,${PWD}/sim-4/sim-4_1_s2.fq.gz" >> "${output_dir}${output_filename}"
echo "sim-4,lib-1,L002,ancient,double-stranded,FC12345,Illumina,${PWD}/sim-4/sim-4_2_s1.fq.gz,${PWD}/sim-4/sim-4_2_s2.fq.gz" >> "${output_dir}${output_filename}"
echo "sim-4,lib-2,L001,ancient,double-stranded,FC12345,Illumina,${PWD}/sim-4/sim-4_3_s1.fq.gz,${PWD}/sim-4/sim-4_3_s2.fq.gz" >> "${output_dir}${output_filename}"

# sim-5: Ancient (3 libraries)
echo "sim-5,lib-1,L001,ancient,double-stranded,FC12345,Illumina,${PWD}/sim-5/sim-5_1_s1.fq.gz,${PWD}/sim-5/sim-5_1_s2.fq.gz" >> "${output_dir}${output_filename}"
echo "sim-5,lib-1,L002,ancient,double-stranded,FC12345,Illumina,${PWD}/sim-5/sim-5_2_s1.fq.gz,${PWD}/sim-5/sim-5_2_s2.fq.gz" >> "${output_dir}${output_filename}"
echo "sim-5,lib-2,L001,ancient,double-stranded,FC12345,Illumina,${PWD}/sim-5/sim-5_3_s1.fq.gz,${PWD}/sim-5/sim-5_3_s2.fq.gz" >> "${output_dir}${output_filename}"

echo "Test data generation completed!"
echo "Samplesheet: ${output_dir}${output_filename}"














