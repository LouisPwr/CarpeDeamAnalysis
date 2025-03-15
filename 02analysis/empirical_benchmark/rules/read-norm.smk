rule norm_reads:
    input:
        reads=config["rdir"] + "/read-extension/{smp}.extended.fastq.gz",
    output:
        norm_reads=config["rdir"] + "/read-norm/{smp}.fq.gz",
        norm_stats=config["rdir"] + "/stats/{smp}.stats-norm.txt",
    threads: config["norm_threads"]
    params:
        seqkit_bin=config["seqkit_bin"],
        bbnorm_bin=config["bbnorm_bin"],
        bbnorm_parms=config["bbnorm_parms"],
        memory=config["bbnorm_memory"],
        clumpify_bin=config["clumpify_bin"],
        clumpify_parms=config["clumpify_parms"],
        norm_reads=config["rdir"] + "/read-norm/{smp}.tmp.fq.gz",
    conda:
        "../envs/qc.yaml"
    log:
        config["rdir"] + "/logs/read-norm/{smp}.read-norm.log",
    benchmark:
        config["rdir"] + "/benchmarks/read-norm/{smp}.read-norm.bmk"
    message:
        """--- Normalize reads"""
    shell:
        """
        # Derep and rename reads
        {params.bbnorm_bin} \
            -Xmx{params.memory} \
            in={input.reads} \
            out={output.norm_reads} \
            {params.bbnorm_parms}

        # {params.clumpify_bin} \
        #     -Xmx{params.memory} \
        #     in={params.norm_reads} \
        #     out={output.norm_reads} \
        #     {params.clumpify_parms}

        {params.seqkit_bin} stats -j {threads} -T -a --quiet \
            {output.norm_reads} > {output.norm_stats}

        rm -rf {params.norm_reads}
        """
