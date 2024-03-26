rule assembly_easytaxonomy_evaluation_penguin:
    input:
        contigs=f'{config["rdir"]}/assembly-genes-eval/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.genes.penguin.{{config}}.faa',
    output:
        mmseqs_tophit_aln=f'{config["rdir"]}/assembly-easytaxonomy-eval/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.easytaxonomy.penguin.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.easytaxonomy.penguin.{{config}}_tophit_aln'
    wildcard_constraints:
        assm_input_1="\w+",
        assm_input_2="\w+",
        config="\w+",
    threads: config["mmseqs_tax_threads"]
    params:
        mmseqs_bin=config["mmseqs_bin"],
        wdir=config["wdir"],
        rdir=config["rdir"] + "/assembly-easytaxonomy-eval",
        taxonomy_db=config["tax_db"],
        mmseqs_results_easy=f'{config["rdir"]}/assembly-easytaxonomy-eval/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.easytaxonomy.penguin.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.easytaxonomy.penguin.{{config}}'
    log:
        mmseqs_log=f'{config["rdir"]}/logs/assembly-easytaxonomy-eval/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.easytaxonomy.penguin.{{config}}.log'
    benchmark:
        f'{config["rdir"]}/benchmarks/assembly-easytaxonomy-eval/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.easytaxonomy.penguin.{{config}}.bmk'
    conda:
        "miniconda3/envs/taxonomy"
    message:
        """--- penguin assembly taxonomic profiling. """
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        NN=$(grep -c '>' {input.contigs} || [[ $? == 1 ]])

        if [[ ${{NN}} -eq 0 ]]; then
            touch {output.mmseqs_tophit_aln}
            exit 0
        fi
        
        {params.mmseqs_bin} easy-taxonomy {input.contigs} {params.taxonomy_db} {params.mmseqs_results_easy} {params.rdir}/tmp --split-memory-limit 300G >> {log.mmseqs_log} 2>&1

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
