rule sanitize_and_pair_reads:
    input:
        fw_reads=lambda wildcards: config["sdir"] + "/" + sample_table_read.fw_reads[wildcards.smp],
        rev_reads=lambda wildcards: config["sdir"] + "/" + sample_table_read.rev_reads[wildcards.smp]
    output:
        fw_paired=f'{config["rdir"]}/reads_paired/{{smp}}_R1.sana.paired.fq.gz',
        rev_paired=f'{config["rdir"]}/reads_paired/{{smp}}_R2.sana.paired.fq.gz',
    threads: config["threads_32"]
    params:
        sana_r1=f'{config["rdir"]}/reads_paired/{{smp}}_R1.sana.fq.gz',
        sana_r2=f'{config["rdir"]}/reads_paired/{{smp}}_R2.sana.fq.gz',
        wdir=config["wdir"],
        outdir=f'{config["rdir"]}/reads_paired'
    log:
        config["rdir"] + "/logs/reads_paired/{smp}.pair.log"
    benchmark:
        config["rdir"] + "/benchmarks/reads_paired/{smp}.pair.bmk"
    message:
        """--- Sanitize and Pair Reads. ---"""
    shell:
        """
        set -x

        # Sanitize the forward reads
        seqkit sana {input.fw_reads} -j {threads} -o {params.sana_r1} >> {log} 2>&1
        # Sanitize the reverse reads
        seqkit sana {input.rev_reads} -j {threads} -o {params.sana_r2} >> {log} 2>&1

        # Pair the sanitized reads
        seqkit pair -1 {params.sana_r1} -2 {params.sana_r2} \
                    -o {params.outdir} -u -j {threads} >> {log} 2>&1

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
