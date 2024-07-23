#!/bin/bash

parent_dir="/home/projects2/metagnm_asm/paper/03simulations/ancientGut/simul/samples"
#"/home/projects/metagnm_asm/paper/data/ancientHorse/simul/samples"
output_base="/home/projects2/metagnm_asm/paper/03simulations/ancientGut/simul/reads"
sizeDist="/home/projects2/metagnm_asm/paper/genreads/meta200/Vi33.19.insize"
dmid="/home/projects2/metagnm_asm/10_gnm/sim/damage/dmid"
dhigh="/home/projects2/metagnm_asm/10_gnm/sim/damage/dhigh"

# Iterate over each subdirectory in the parent directory
for subdir in ${parent_dir}/*; do
    if [ -d "${subdir}" ]; then
        base=$(basename ${subdir})

        # Update the paths based on the current subdirectory
        input="${subdir}/input"
        cat "${subdir}/input/bact/sum.txt"
        size=$(cat "${subdir}/input/bact/sum.txt") #| tr -cd '0-9')
        # cov3=$(( (3 * $size) / 50 ))
        # cov5=$(( (5 * $size) / 50 ))
        # cov10=$(( (10 * $size) / 50 ))

        # Iterate over both dmid and dhigh profiles
        for dprof in ${dmid} ${dhigh}; do
            profile_name=$(basename ${dprof})

            # Iterate over all coverage values
            cov_values=(3 5 10)
            for cov_multiplier in "${cov_values[@]}"; do
                current_cov=$(( ($cov_multiplier * $size) / 50 ))
                outfile="${output_base}/${base}_${profile_name}_cov${cov_multiplier}/out/${base}_${profile_name}_cov${cov_multiplier}"
                log="${output_base}/${base}_${profile_name}_cov${cov_multiplier}/out/${base}_${profile_name}_cov${cov_multiplier}.log"

                # Run gargammel.pl for the current subdirectory with the provided parameters
                mkdir -p "${output_base}/${base}_${profile_name}_cov${cov_multiplier}/out"
                #nice -n19 gargammel.pl --comp 1,0,0 -n ${current_cov} -s ${sizeDist} -matfile ${dprof} -o ${outfile} ${input} &> ${log}
                nice -n19 /home/ctools/leehom-1.2.17 --ancientdna -t 42 -fq1 "${outfile}_s1.fq.gz" -fq2 "${outfile}_s2.fq.gz" -fqo "${outfile}_fin"
            done
        done
    fi
done








# input="/home/projects/metagnm_asm/paper/genreads/meta100/input"
# outfile="/home/projects/metagnm_asm/paper/genreads/meta100/out/meta100"
# log="/home/projects/metagnm_asm/paper/genreads/meta100/out/meta100.log"

# # Parameters
# dmid="/home/projects/metagnm_asm/paper/damage/dmid/dmid"
# dhigh="path/to/dhigh"

# sizeDist="/home/projects/metagnm_asm/paper/genreads/meta200/Vi33.19.insize"

# cov="path/to/count.txt"

# gargammel.pl --comp 1,0,0 -n ${cov} -s ${sizeDist} -matfile ${dprof} -o ${outfile} ${input} &> ${log}
