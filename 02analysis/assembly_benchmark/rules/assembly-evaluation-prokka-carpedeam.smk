rule assembly_annotation_carpedeam:
    input:
        contigs=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.carpedeam.{{config}}.fasta',
    output:
        proteinsAA=f'{config["rdir"]}/assembly-annotation-eval/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.proteins.carpedeam.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.proteins.carpedeam.{{config}}.faa',
    wildcard_constraints:
        assm_input_1="\w+",
        assm_input_2="\w+",
        config="\w+",
    threads: config["prokka_threads"]
    params:
        prokka_bin=config["prokka_bin"],
        wdir=config["wdir"],
        rdir=config["rdir"] + "/assembly-annotation-eval",
        prefix=f'{{smp}}.{{assm_input_1}}-{{assm_input_2}}.proteins.carpedeam.{{config}}',
        outDir=f'{config["rdir"]}/assembly-annotation-eval/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.proteins.carpedeam.{{config}}',
    log:
        prokka_log=f'{config["rdir"]}/logs/assembly-annotation-eval/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.proteins.carpedeam.{{config}}.log',
    benchmark:
        f'{config["rdir"]}/benchmarks/assembly-annotation-eval/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.proteins.carpedeam.{{config}}.bmk'
    conda:
        "miniconda3/envs/prokka"
    message:
        """--- prokka protein annotation for carpedeam assemblies. """
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        N=$(grep -c '>' {input.contigs} || [[ $? == 1 ]])

        if [[ ${{N}} -eq 0 ]]; then
            touch {output.proteinsAA}
            exit 0
        fi

        {params.prokka_bin} \
            --fast \
            --metagenome \
            --outdir {params.outDir} \
            --prefix {params.prefix} \
            --cpus {threads} \
            --force \
            {input.contigs}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
