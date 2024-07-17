rule assembly_mmseqs_evaluation_carpedeam:
    input:
        contigs=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.carpedeam.{{config}}.fasta',
        #contigs=f'{config["rdir"]}/assembly-annotation-eval/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.proteins.carpedeam.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.proteins.carpedeam.{{config}}.faa',
    output:
        tsv=f'{config["rdir"]}/assembly-mmseqs/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mmseqs.carpedeam.{{config}}.tsv',
    wildcard_constraints:
        assm_input_1="\w+",
        assm_input_2="\w+",
        config="\w+",
    params:
        tmp=f'{config["rdir"]}/assembly-mmseqs/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mmseqs.carpedeam.{{config}}_tmp',
        mmseqs_bin=config["mmseqs_bin"],
        wdir=config["wdir"],
        rdir=config["rdir"] + "/assembly-mmseqs",
        ref_folder=config["ref_folder"],
        # reference=(
        #     lambda wildcards: f'{config["ref_folder"]}/{sample_table_read.reference[wildcards.smp]}'
        # )
        reference=config["adir"] + "/all.faa",
    threads: config["mmseqs_threads"]
    log:
        mmseqs_log=f'{config["rdir"]}/assembly-mmseqs/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mmseqs.carpedeam.{{config}}.log',
    benchmark:
        f'{config["rdir"]}/benchmarks/assembly-mmseqs/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mmseqs.carpedeam.{{config}}.bmk'
    message:
        """--- mmseqs assemblies vs. reference carpedeam """
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        N=$(grep -c '>' {input.contigs} || [[ $? == 1 ]])

        if [[ ${{N}} -eq 0 ]]; then
            touch {output.tsv}
            exit 0
        fi

        {params.mmseqs_bin} easy-search \
            {input.contigs} \
            {params.reference} \
            {output.tsv} \
            {params.tmp} \
            --max-seq-len 1000000 --search-type 4 --cov-mode 1 -c 0.5 \
            --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov \
            --threads 32 --split-memory-limit 300G >> {log.mmseqs_log} 2>&1

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
