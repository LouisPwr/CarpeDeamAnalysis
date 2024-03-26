rule assembly_genes_evaluation_megahit:
    input:
        contigs=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.megahit.{{config}}.fasta',
    output:
        genesNUC=f'{config["rdir"]}/assembly-genes-eval/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.genes.megahit.{{config}}.fna',
    wildcard_constraints:
        assm_input_1="\w+",
        assm_input_2="\w+",
        config="\w+",
    params:
        genesAA=f'{config["rdir"]}/assembly-genes-eval/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.genes.megahit.{{config}}.faa',
        genesNUC_full=f'{config["rdir"]}/assembly-genes-eval/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.genes.megahit.fullorfs.{{config}}.fna',
        genesGBK=f'{config["rdir"]}/assembly-genes-eval/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.genes.megahit.{{config}}.gbk',
        prodigal_bin=config["prodigal_bin"],
        wdir=config["wdir"],
        rdir=config["rdir"] + "/assembly-genes-eval",
    log:
        prodigal_log=f'{config["rdir"]}/logs/assembly-genes-eval/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.genes.megahit.{{config}}.log',
    benchmark:
        f'{config["rdir"]}/benchmarks/assembly-genes-eval/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.genes.megahit.{{config}}.bmk'
    conda:
        "/miniconda3/envs/taxonomy"
    message:
        """--- prodigal gene prediction for megahit assemblies. """
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        N=$(grep -c '>' {input.contigs} || [[ $? == 1 ]])

        if [[ ${{N}} -eq 0 ]]; then
            touch {params.genesAA}
            touch {output.genesNUC}
            touch {params.genesGBK}
            touch {params.genesNUC_full}

            exit 0
        fi

        {params.prodigal_bin} \
            -i {input.contigs} \
            -o {params.genesGBK} \
            -c \
            -a {params.genesAA} \
            -d {output.genesNUC} \
            -p meta >> {log.prodigal_log} 2>&1

        seqtk seq {output.genesNUC} | grep --no-group-separator -A 1 "partial=00" > {params.genesNUC_full}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
