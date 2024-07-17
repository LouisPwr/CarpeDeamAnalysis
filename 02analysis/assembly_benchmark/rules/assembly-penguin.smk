def get_pg_conf(wc):
    return config["penguin_config"][wc.config]


rule assembly_penguin:
    input:
        reads=f'{config["rdir"]}/assembly-input-2/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.fastq.gz',
    output:
        contigs=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.penguin.{{config}}.fasta',
        config=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.penguin.{{config}}.conf',
    threads: config["penguin_threads"]
    params:
        assm_input_1="{assm_input_1}",
        assm_input_2="{assm_input_2}",
        penguin_bin=config["penguin_bin"],
        seqkit_bin=config["seqkit_bin"],
        tmp_dir=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.penguin.{{config}}_tmp',
        contigs=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.penguin.{{config}}.tmp.fasta',
        penguin_parms=lambda wc: get_pg_conf(wc),
        rdir=config["rdir"] + "/assembly",
        wdir=config["wdir"],
        pattern=".+",
        subst="{smp}_pg_\{nr\}",
        name="{smp}_penguin--{assm_input_1}-{assm_input_2}.{config}",
    log:
        config["rdir"]
        + "/logs/assembly/{smp}.{assm_input_1}-{assm_input_2}.assm.penguin.{config}.log",
    benchmark:
        (
            config["rdir"]
            + "/benchmarks/assembly/{smp}.{assm_input_1}-{assm_input_2}.assm.penguin.{config}.bmk"
        )
    resources:
        time="1000:00:00",
        nodes=1,
        ntasks_per_node=1,
        cpus_per_task=32,
        mem_gb=400,
        partition="debug"
    message:
        """--- Assembling reads with penguin."""
    shell:
        """
        set -x

        mkdir -p {params.tmp_dir}
        rm -rf {params.tmp_dir}/*
        rm -rf {params.contigs}

        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        CONFIG="{params.penguin_parms}"
        echo "${{CONFIG}}" > {output.config}

        {params.penguin_bin} guided_nuclassemble \
            {input.reads} {params.contigs} {params.tmp_dir} \
            {params.penguin_parms} \
            --threads {threads} \
            >> {log} 2>&1

        {params.seqkit_bin} replace \
            -j {threads} \
            -p {params.pattern}  \
            -r {params.subst} \
            -o {output.contigs} \
            --nr-width 12 \
            {params.contigs} \
            >> {log} 2>&1


        rm -rf {params.tmp_dir}/ {params.contigs}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
