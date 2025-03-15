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
    threads: config["threads_32"]
    params:
        assm_input_1="{assm_input_1}",
        assm_input_2="{assm_input_2}",
        mmseqs_bin=config["mmseqs_bin"],
        megahit_bin=config["megahit_bin"],
        megahit_tk_bin=config["megahit_tk_bin"],
        seqkit_bin=config["seqkit_bin"],
        stmp_dir=config["rdir"] + "/assembly",
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
        
        nice -n 19  {params.megahit_bin} \
            -r {input} \
            {params.megahit_parms} \
            --min-contig-len {params.min_contig_length} \
            -o {params.tmp_dir} \
            --num-cpu-threads {threads} \
             >> {log} 2>&1

        #N=$(grep -c '>' {params.tmp_dir}/final.contigs.fa)
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

        # cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}
        # {params.mmseqs_bin} easy-cluster \
        #     {output.contigs} {params.rdir}/{params.name}_clu \
        #     {params.tmp_dir} --min-seq-id 0.9 --cov-mode 1 -c 0.9 >> {log} 2>&1
        rm -rf {params.tmp_dir}
        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
