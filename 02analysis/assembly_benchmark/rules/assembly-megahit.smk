def get_mh_conf(wc):
    return config["megahit_config"][wc.config]


rule assembly_megahit:
    input:
        reads=f'{config["rdir"]}/assembly-input-2/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.fastq.gz',
    output:
        contigs=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.megahit.{{config}}.fasta',
        config=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.megahit.{{config}}.conf',
    wildcard_constraints:
        assm_input_1="\w+",
        assm_input_2="\w+",
        config="\w+",
    threads: config["megahit_threads"]
    params:
        assm_input_1="{assm_input_1}",
        assm_input_2="{assm_input_2}",
        megahit_bin=config["megahit_bin"],
        megahit_tk_bin=config["megahit_tk_bin"],
        seqkit_bin=config["seqkit_bin"],
        tmp_dir=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.megahit.{{config}}_tmp',
        min_contig_length=config["megahit_min_contig_length"],
        megahit_parms=lambda wc: get_mh_conf(wc),
        config="{config}",
        wdir=config["wdir"],
        rdir=config["rdir"] + "/assembly",
        pattern=".+",
        subst="{smp}_mh_\{nr\}",
        name="{smp}_megahit--{assm_input_1}-{assm_input_2}.{config}",
    log:
        config["rdir"]
        + "/logs/assembly/{smp}.{assm_input_1}-{assm_input_2}.assm.megahit.{config}.log",
    benchmark:
        (
            config["rdir"]
            + "/benchmarks/assembly/{smp}.{assm_input_1}-{assm_input_2}.assm.megahit.{config}.bmk"
        )
    resources:
        time="1000:00:00",
        nodes=1,
        ntasks_per_node=1,
        cpus_per_task=32,
        mem_gb=400,
        partition="debug"
    message:
        """--- Assembling reads with Megahit."""
    shell:
        """
        set -x
        if [ -d {params.tmp_dir} ]; then 
            rm -Rf {params.tmp_dir}
        fi
        CONFIG="{params.megahit_parms}"
        echo "${{CONFIG}}" > {output.config}
        {params.megahit_bin} \
            -r {input} \
            {params.megahit_parms} \
            --min-contig-len {params.min_contig_length} \
            -o {params.tmp_dir} \
            --num-cpu-threads {threads} \
             >> {log} 2>&1

        if [ -s {params.tmp_dir}/final.contigs.fa ]; then
            {params.seqkit_bin} replace \
                -j {threads} \
                -p {params.pattern} \
                -r {params.subst} \
                -o {output.contigs} \
                --nr-width 12 \
                {params.tmp_dir}/final.contigs.fa \
                >> {log} 2>&1
        else
            touch {output.contigs}
        fi

        rm -rf {params.tmp_dir}
        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
