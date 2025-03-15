rule map_reads_assembly:
    input:
        reads=config["rdir"] + "/read-renamed/{smp}.fq.gz",
        contigs=config["rdir"] + "/assembly-refined/{smp}.assm.refined.fasta",
        stats=config["rdir"] + "/stats/{smp}.stats-initial.txt",
    output:
        bam=config["rdir"] + "/map-assembly/{smp}.assm.refined.sorted.bam",
        bam_dedup=(
            config["rdir"] + "/map-assembly/{smp}.assm.refined.markdup.sorted.bam"
        ),
        bam_dedup_metrics=(
            config["rdir"] + "/map-assembly/{smp}.assm.refined.markdup.metrics"
        ),
    threads: config["threads_32"]
    params:
        bowtie2_bin=config["bowtie2_bin"],
        bowtie2_build_bin=config["bowtie2_build_bin"],
        bowtie2_db=config["rdir"] + "/map-assembly/{smp}-db",
        bowtie2_parms=config["bowtie2_parms"],
        samtools_bin=config["samtools_bin"],
        samtools_view_parms=config["samtools_view_parms"],
        samtools_sort_parms=config["samtools_sort_parms"],
        picard_bin=config["picard_bin"],
        java_opts=config["java_opts"],
        bam_dedup_tmp=(
            config["rdir"] + "/map-assembly/{smp}.assm.refined.markdup.tmp.bam"
        ),
        bam_dedup_sorted=(
            config["rdir"] + "/map-assembly/{smp}.assm.refined.markdup.sorted.bam"
        ),
        rdir=config["rdir"] + "/map-assembly",
        wdir=config["wdir"],
        name="{smp}_map-assembly",
    log:
        config["rdir"] + "/logs/map-assembly/{smp}.assm.map-assembly.log",
    benchmark:
        config["rdir"] + "/benchmarks/map-assembly/{smp}.map-assembly.bmk"
    conda:
        "../envs/mapping.yaml"
    resources:
        time="1000:00:00",
        nodes=1,
        ntasks_per_node=1,
        cpus_per_task=16,
        mem_gb=100,
        partition="debug"
    message:
        """--- Mapping reads to assembly."""
    shell:
        """
        set -x

        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        IS_FASTQ=$(awk 'NR==2{{print $2}}' {input.stats})

        {params.bowtie2_build_bin} {input.contigs} {params.bowtie2_db}  \
             >> {log} 2>&1
        if [ ${{IS_FASTQ}} = "FASTQ" ]; then
            {params.bowtie2_bin} -p {threads} -x {params.bowtie2_db} -q -U {input.reads} {params.bowtie2_parms} \
                | samtools view {params.samtools_view_parms} \
                | samtools sort -@ {threads} {params.samtools_sort_parms} > {output.bam}
        else
            {params.bowtie2_bin} -p {threads} -x {params.bowtie2_db} -f -U {input.reads} {params.bowtie2_parms} \
                | samtools view {params.samtools_view_parms} \
                | samtools sort -@ {threads} {params.samtools_sort_parms} > {output.bam}
        fi

        {params.picard_bin} MarkDuplicates \
            -INPUT {output.bam} \
            -OUTPUT {params.bam_dedup_tmp} \
            -METRICS_FILE {output.bam_dedup_metrics} \
            -AS TRUE \
            -VALIDATION_STRINGENCY LENIENT \
            -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 \
            -REMOVE_DUPLICATES TRUE  \
             >> {log} 2>&1

        {params.samtools_bin} sort -@ {threads} {params.samtools_sort_parms} {params.bam_dedup_tmp} > {output.bam_dedup}
        {params.samtools_bin} index {output.bam_dedup}

        rm {params.bam_dedup_tmp}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
