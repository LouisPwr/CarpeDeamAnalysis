rule assm_input_2:
    input:
        reads=config["rdir"] + "/assembly-input-1/{smp}.{assm_input_1}.fastq.gz",
        stats=config["rdir"] + "/stats/{smp}.{assm_input_1}.stats.txt",
    output:
        reads=f'{config["rdir"]}/assembly-input-2/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.fastq.gz',
        stats=f'{config["rdir"]}/stats/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.stats.txt',
    wildcard_constraints:
        assm_input_1="\w+",
        assm_input_2="\w+",
    threads: config["threads_32"]
    params:
        assm_input="{assm_input_2}",
        seqkit_bin=config["seqkit_bin"],
        vsearch_bin=config["vsearch_bin"],
        derep_minlen=config["read_minlen"],
        derep_parms=config["derep_parms"],
        bbnorm_bin=config["bbnorm_bin"],
        bbnorm_parms=config["bbnorm_parms"],
        memory=config["bbnorm_memory"],
    log:
        config["rdir"]
        + "/logs/assembly-input-2/{smp}.{assm_input_1}-{assm_input_2}.assembly-input-2.log",
    benchmark:
        (
            config["rdir"]
            + "/benchmarks/assembly-input-2/{smp}.{assm_input_1}-{assm_input_2}.assembly-input-2.bmk"
        )
    message:
        """--- Dereplicate reads"""
    shell:
        """
        if [[ {params.assm_input} == "derep" ]]; then
            # Derep and rename reads
            {params.vsearch_bin} \
                --fastx_uniques \
                {input.reads} \
                --minseqlength {params.derep_minlen} \
                {params.derep_parms} \
            | gzip > {output.reads} 

            {params.seqkit_bin} stats -j {threads} -T -a --quiet \
                {output.reads} -o {output.stats}
        elif [[ {params.assm_input} == "norm" ]]; then
            {params.bbnorm_bin} \
                -Xmx{params.memory} \
                in={input.reads} \
                out={output.reads} \
                {params.bbnorm_parms}

        {params.seqkit_bin} stats -j {threads} -T -a --quiet \
            {output.reads} > {output.stats}
        else
            ln -s {input.reads} {output.reads}
            ln -s {input.stats} {output.stats}
        fi
        """
