rule assm_input_1:
    input:
        reads=config["rdir"] + "/read-renamed/{smp}.fq.gz",
        stats=config["rdir"] + "/stats/{smp}.stats-initial.txt",
    wildcard_constraints:
        #assm_input_1="\w+",
        assm_input_1="\w+",
    output:
        assm_input_1=config["rdir"] + "/assembly-input-1/{smp}.{assm_input_1}.fastq.gz",
        stats=config["rdir"] + "/stats/{smp}.{assm_input_1}.stats.txt",
    threads: config["threads_32"]
    params:
        assm_input="{assm_input_1}",
        seqkit_bin=config["seqkit_bin"],
        tadpole_bin=config["tadpole_bin"],
        loglog_bin=config["loglog_bin"],
        mode=config["tadpole_mode"],
        length=config["tadpole_length"],
        memory=config["tadpole_memory"],
        k=config["tadpole_k"],
        ibb=config["ibb"],
        ecc=config["ecc"],
        ecco=config["ecco"],
        conservative=config["conservative"],
        seed=config["loglog_seed"],
        bits=config["loglog_bits"],
        maxmem=config["maxmem"],
        hashes=config["loglog_hashes"],
        table_cap=config["loglog_table_cap"],
        trimends=config["trimends"],
        prefilter=config["prefilter"],
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
        if [[ {params.assm_input} == "extended" ]]; then
            set -x
            MEM=$({params.loglog_bin} seed={params.seed} k={params.k} in={input.reads} ignorebadquality 2> >(grep Cardinality) \
                | awk -vP={params.table_cap} -vB={params.bits} -vH={params.hashes} '{{print int( (((B*H)/8)*$2)/P )}}')
            echo ${{MEM}}
            if [[ ${{MEM}} -gt {params.maxmem} ]]; then
                MEM={params.maxmem}
            fi

            {params.tadpole_bin} \
                -Xmx{params.memory} \
                k={params.k} \
                in={input.reads} out={output.assm_input_1} mode={params.mode} ibb={params.ibb} prefilter={params.prefilter} \
                el={params.length} er={params.length} threads={threads} overwrite=true trimends={params.trimends} \
                ecc={params.ecc} ecco={params.ecco} filtermem=${{MEM}} conservative={params.conservative} ignorebadquality \
                >> {log} 2>&1

            {params.seqkit_bin} stats -j {threads} -T -a --quiet \
                {output.assm_input_1} -o {output.stats}
        else
            ln -s {input.reads} {output.assm_input_1}
            ln -s {input.stats} {output.stats}
        fi
        """
