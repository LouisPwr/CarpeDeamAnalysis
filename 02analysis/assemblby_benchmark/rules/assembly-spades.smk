def get_sp_conf(wc):
    return config["spades_config"][wc.config]

rule assembly_spades:
    input:
        reads=f'{config["rdir"]}/assembly-input-2/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.fastq.gz',
        reads1=(
        lambda wildcards: f'{config["sdir"]}/{sample_table_read.fw_reads[wildcards.smp]}'
        ),
        reads2=(
        lambda wildcards: f'{config["sdir"]}/{sample_table_read.rev_reads[wildcards.smp]}'
        )
    output:
        contigs=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.spades.{{config}}.fasta',
        config=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.spades.{{config}}.conf',
    threads: config["spades_threads"]
    params:
        assm_input_1="{assm_input_1}",
        assm_input_2="{assm_input_2}",
        spades_bin=config["spades_bin"],
        seqkit_bin=config["seqkit_bin"],
        tmp_dir=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.spades.{{config}}_tmp',
        min_contig_length=config["spades_min_contig_length"],
        spades_parms=lambda wc: get_sp_conf(wc),
        config="{config}",
        wdir=config["wdir"],
        rdir=config["rdir"] + "/assembly",
        pattern=".+",
        subst="{smp}_mh_\{nr\}",
        name="{smp}_spades--{assm_input_1}-{assm_input_2}.{config}",
    log:
        config["rdir"]
        + "/logs/assembly/{smp}.{assm_input_1}-{assm_input_2}.assm.spades.{config}.log",    
    benchmark:
        (
            config["rdir"]
            + "/benchmarks/assembly/{smp}.{assm_input_1}-{assm_input_2}.assm.spades.{config}.bmk"
        ) 
    resources:
        time="1000:00:00",
        nodes=1,
        ntasks_per_node=1,
        cpus_per_task=32,
        mem_gb=400,
        partition="debug"
    message:
        """--- Assembly with metaSPADES"""
    shell:
        """
        set -x
        if [ -d {params.tmp_dir} ]; then 
            rm -Rf {params.tmp_dir}
        fi
        CONFIG="{params.spades_parms}"
        echo "${{CONFIG}}" > {output.config}
        {params.spades_bin} \
            -1 {input.reads1} \
            -2 {input.reads2} \
            --merged {input.reads} \
            {params.spades_parms} \
            -o {params.tmp_dir} \
            -t {threads} \
             >> {log} 2>&1

        if [ -s {params.tmp_dir}/contigs.fasta ]; then
            {params.seqkit_bin} seq \
                -j {threads} \
                -m {params.min_contig_length} \
                {params.tmp_dir}/contigs.fasta \
                > {params.tmp_dir}/contigs.{params.min_contig_length}.fasta

            {params.seqkit_bin} replace \
                -j {threads} \
                -p {params.pattern} \
                -r {params.subst} \
                -o {output.contigs} \
                --nr-width 12 \
                {params.tmp_dir}/contigs.{params.min_contig_length}.fasta \
                >> {log} 2>&1
        else
            touch {output.contigs}
        fi

        rm -rf {params.tmp_dir}
        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
