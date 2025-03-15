rule assembly_stats:
    input:
        contigs_refined=config["rdir"] + "/assembly-refined/{smp}.assm.refined.fasta",
        contigs_combined=(
            config["rdir"] + "/assembly-combined/{smp}.assm.combined.fasta"
        ),
        contigs_extended=(
            config["rdir"] + "/assembly-extended/{smp}.assm.megahit-extended.fasta"
        ),
    output:
        stats_refined=config["rdir"] + "/stats/{smp}.assembly-refined.txt",
        stats_combined=config["rdir"] + "/stats/{smp}.assembly-combined.txt",
        stats_extended=config["rdir"] + "/stats/{smp}.assembly-extended.txt",
    threads: config["threads_8"]
    params:
        seqkit_bin=config["seqkit_bin"],
        label="{smp}",
    log:
        config["rdir"] + "/logs/stats/{smp}.assembly-stats.log",
    benchmark:
        config["rdir"] + "/benchmarks/stats/{smp}.assembly-stats.bmk"
    message:
        """--- Get assembly refined stats"""
    shell:
        """
        {params.seqkit_bin} fx2tab -j {threads} -l -n --quiet {input.contigs_refined} \
            | awk -vS={params.label} '{{print S"\t"$1"\t"$2}}'  > {output.stats_refined}
        {params.seqkit_bin} fx2tab -j {threads} -l -n --quiet {input.contigs_combined} \
            | awk -vS={params.label} '{{print S"\t"$1"\t"$2}}'  > {output.stats_combined}
        {params.seqkit_bin} fx2tab -j {threads} -l -n --quiet {input.contigs_extended} \
            | awk -vS={params.label} '{{print S"\t"$1"\t"$2}}'  > {output.stats_extended}
        """
