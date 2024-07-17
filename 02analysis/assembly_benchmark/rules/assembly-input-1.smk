rule assm_input_1:
    input:
        reads=config["rdir"] + "/read-renamed/{smp}.fq.gz",
        stats=config["rdir"] + "/stats/{smp}.stats-initial.txt",
    wildcard_constraints:
        assm_input_1="\w+",
    output:
        assm_input_1=config["rdir"] + "/assembly-input-1/{smp}.{assm_input_1}.fastq.gz",
        stats=config["rdir"] + "/stats/{smp}.{assm_input_1}.stats.txt",
    threads: config["assm_input_threads"]
    params:
        assm_input="{assm_input_1}",
    resources:
        time="1000:00:00",
        nodes=1,
        ntasks_per_node=1,
        cpus_per_task=16,
        mem_gb=100,
        partition="debug"
    log:
        config["rdir"]
        + "/logs/assembly-input-1/{smp}.{assm_input_1}.assembly-input-1.log",
    benchmark:
        (
            config["rdir"]
            + "/benchmarks/assembly-input-1/{smp}.{assm_input_1}.assembly-input-1.bmk"
        )
    message:
        """--- ASSEMBLY INPUT 1"""
    shell:
        """
        ln -s {input.reads} {output.assm_input_1}
        ln -s {input.stats} {output.stats}
        """
