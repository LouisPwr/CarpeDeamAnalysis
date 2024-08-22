rule assembly_linclust_evaluation_penguin:
    input:
        #contigs=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.penguin.{{config}}.fasta',
        contigs=f'{config["rdir"]}/assembly-annotation-eval/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.proteins.penguin.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.proteins.penguin.{{config}}.faa',
    output:
        cluster=f'{config["rdir"]}/assembly-linclust/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.penguin.{{config}}_rep_seq.fasta',
    wildcard_constraints:
        assm_input_1="\w+",
        assm_input_2="\w+",
        config="\w+",
    params:
        tmp=f'{config["rdir"]}/assembly-linclust/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.penguin.{{config}}_tmp',
        mmseqs_bin=config["mmseqs_bin"],
        wdir=config["wdir"],
        rdir=config["rdir"] + "/assembly-linclust",
        ref_folder=config["ref_folder"],
        # reference=(
        #     lambda wildcards: f'{config["ref_folder"]}/{sample_table_read.reference[wildcards.smp]}'
        # )
        clusterID=f'{config["rdir"]}/assembly-linclust/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.penguin.{{config}}',
    threads: config["mmseqs_threads"]
    log:
        mmseqs_log=f'{config["rdir"]}/assembly-linclust/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.penguin.{{config}}.log',
    benchmark:
        f'{config["rdir"]}/benchmarks/assembly-linclust/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.penguin.{{config}}.bmk'
    message:
        """--- mmseqs assemblies vs. reference penguin """
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        N=$(grep -c '>' {input.contigs} || [[ $? == 1 ]])

        if [[ ${{N}} -eq 0 ]]; then
            touch {output.cluster}
            exit 0
        fi

        {params.mmseqs_bin} easy-linclust \
            {input.contigs} \
            {params.clusterID} \
            {params.tmp} \
            --min-seq-id 1 --kmer-per-seq 80 --cov-mode 1 -c 0.8 \
            --threads 32 --split-memory-limit 300G >> {log.mmseqs_log} 2>&1


        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
