rule assembly_mapping_evaluation_carpedeam:
    input:
        contigs=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.carpedeam.{{config}}.fasta',
        reads=f'{config["rdir"]}/assembly-input-2/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.fastq.gz',
    output:
        bam_sorted=f'{config["rdir"]}/assembly-mapping/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mapping.carpedeam.{{config}}.sorted.bam',
        flagstat=f'{config["rdir"]}/assembly-mapping/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mapping.carpedeam.{{config}}.stats.txt',
    wildcard_constraints:
        assm_input_1="\w+",
        assm_input_2="\w+",
        config="\w+",
    params:
        sam=f'{config["rdir"]}/assembly-mapping/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mapping.carpedeam.{{config}}.sam',
        bam=f'{config["rdir"]}/assembly-mapping/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mapping.carpedeam.{{config}}.bam',
        bowtie2_bin=config["bowtie2_bin"],
        bowtie2_build_bin=config["bowtie2_build_bin"],
        samtools_bin=config["samtools_bin"],
        wdir=config["wdir"],
        rdir=config["rdir"] + "/assembly-mapping",
        flags=config["bowtie2_map_mode"],
        bt_index=f'{config["rdir"]}/assembly-mapping/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mapping.carpedeam.{{config}}',
    threads: config["bowtie2_threads"]
    log:
        mapping_log=f'{config["rdir"]}/assembly-mapping/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mapping.carpedeam.{{config}}.log',
    benchmark:
        f'{config["rdir"]}/benchmarks/assembly-mapping/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mapping.carpedeam.{{config}}.bmk'
    conda:
        "/miniconda3/envs/bowtie2"
    message:
        """--- bowtie2 mappping of carpedeam assemblies. """
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        N=$(grep -c '>' {input.contigs} || [[ $? == 1 ]])

        if [[ ${{N}} -eq 0 ]]; then
            touch {output.bam_sorted}
            touch {output.flagstat}
            exit 0
        fi

        {params.bowtie2_build_bin} \
            --large-index \
            --threads 16 \
            {input.contigs} \
            {params.bt_index}

        {params.bowtie2_bin} \
            {params.flags} \
            -x {params.bt_index} \
            -U {input.reads} \
            -S {params.sam} >> {log.mapping_log} 2>&1

        {params.samtools_bin} view -bS {params.sam} > {params.bam}
        {params.samtools_bin} sort {params.bam} -o {output.bam_sorted}
        {params.samtools_bin} index {output.bam_sorted}
        {params.samtools_bin} flagstat {output.bam_sorted} > {output.flagstat}

        rm -rf {params.sam}
        rm -rf {params.bam}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """

