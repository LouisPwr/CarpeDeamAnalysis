rule merge_reads:
    input:
        fw_reads=lambda wildcards: f'{config["rdir"]}/reads_paired/{wildcards.smp}_R1.sana.paired.fq.gz',
        rev_reads=lambda wildcards: f'{config["rdir"]}/reads_paired/{wildcards.smp}_R2.sana.paired.fq.gz'
    output:
        reads_out=f'{config["rdir"]}/reads_merged/{{smp}}.merge.fq.gz',
    threads: config["threads_32"]
    params:
        leehom_bin=config["leehom_bin"],
        reads_prefix=f'{config["rdir"]}/reads_merged/{{smp}}.merge',
        reads_prefix_tmp=f'{config["rdir"]}/reads_merged/{{smp}}.merge.tmp',
        rdir=config["rdir"] + "/reads_merged",
        wdir=config["wdir"],
        reads_inter=f'{config["rdir"]}/reads_merged/{{smp}}.merge.fq.gz',
    log:
        config["rdir"] + "/logs/reads_merged/{smp}.merge.log",
    benchmark:
        config["rdir"] + "/benchmarks/assembly/{smp}.merge.bmk"
    message:
        """--- Merge Reads. ---"""
    shell:
        """
        set -x

        nice -n 19 {params.leehom_bin} --ancientdna -t {threads} \
            -fq1 {input.fw_reads} -fq2 {input.rev_reads} \
            -fqo {params.reads_prefix_tmp} >> {log} 2>&1

        seqkit seq -m 35 -Q 20 {params.reads_prefix_tmp}.fq.gz | gzip > {params.reads_prefix}.fq.gz

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
