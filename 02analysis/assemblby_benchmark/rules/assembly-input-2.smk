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
    threads: config["assm_input_threads"]
    params:
        assm_input="{assm_input_2}",
    resources:
        time="1000:00:00",
        nodes=1,
        ntasks_per_node=1,
        cpus_per_task=16,
        mem_gb=100,
        partition="debug"
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
        ln -s {input.reads} {output.reads}
        ln -s {input.stats} {output.stats}
        """
