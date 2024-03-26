rule assembly_evaluation_quast_spades:
    input:
        contigs=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.spades.{{config}}.fasta',
    output:
        report=f'{config["rdir"]}/assembly-evaluation-quast/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.spades.{{config}}.report.tsv',
    wildcard_constraints:
        assm_input_1="\w+",
        assm_input_2="\w+",
        config="\w+",
    threads: config["quast_threads"]
    params:
        quast_bin=config["quast_bin"],
        ref_folder=config["ref_folder"],
        reference=(
            lambda wildcards: f'{config["ref_folder"]}/{sample_table_read.reference[wildcards.smp]}'
        ),
        odir=f'{config["rdir"]}/assembly-evaluation-quast/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.spades.{{config}}',
        wdir=config["wdir"],
        rdir=config["rdir"] + "/assembly-evaluation-quast",
    log:
        f'{config["rdir"]}/logs/assembly-evaluation-quast/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.spades.{{config}}.log',
    benchmark:
        f'{config["rdir"]}/benchmarks/assembly-evaluation-quast/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.spades.{{config}}.bmk'
    conda:
        "/vol/cloud/louis/miniconda3/envs/metaquast"
    message:
        """--- Evaluating spades assemblies. """
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        N=$(grep -c '>' {input.contigs} || [[ $? == 1 ]])

        if [[ ${{N}} -eq 0 ]]; then
            touch {output.report}
            exit 0
        fi

        {params.quast_bin} \
            --threads {threads} \
            -r {reference} \
            --fragmented \
            --min-identity 89.99 \
            --unique-mapping \
            --no-icarus \
            -o {params.odir} \
            {input.contigs}

        cp {params.odir}/combined_reference/transposed_report.tsv {output.report}
        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
