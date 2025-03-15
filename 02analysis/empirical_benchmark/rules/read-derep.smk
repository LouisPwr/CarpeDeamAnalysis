rule derep_reads:
    input:
        reads=config["rdir"] + "/read-extension/{smp}.extended.fastq.gz",
    output:
        derep_reads=config["rdir"] + "/read-derep/{smp}.fa.gz",
        derep_stats=config["rdir"] + "/stats/{smp}.stats-derep.txt",
    threads: config["read_derep_threads"]
    params:
        seqkit_bin=config["seqkit_bin"],
        vsearch_bin=config["vsearch_bin"],
        derep_minlen=config["read_minlen"],
        derep_parms=config["derep_parms"],
    conda:
        "../envs/qc.yaml"
    log:
        config["rdir"] + "/logs/read-derep/{smp}.read-derep.log",
    benchmark:
        config["rdir"] + "/benchmarks/read-derep/{smp}.read-derep.bmk"
    message:
        """--- Dereplicate reads"""
    shell:
        """
        # Derep and rename reads
        {params.vsearch_bin} \
            --fastx_uniques \
            {input.reads} \
            --minseqlength {params.derep_minlen} \
            {params.derep_parms} \
        | sed -e 's|;size=|__|' \
        | gzip > {output.derep_reads} 

        {params.seqkit_bin} stats -j {threads} -T -a --quiet \
            {output.derep_reads} -o {output.derep_stats}
        """
