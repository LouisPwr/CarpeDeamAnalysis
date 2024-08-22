rule assembly_skani_evaluation_penguin:
    input:
        contigs=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.penguin.{{config}}.fasta',
    output:
        skani=f'{config["rdir"]}/assembly-skani/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.genes.penguin.{{config}}.ani.tsv',
        searchskani=f'{config["rdir"]}/assembly-skani/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.genes.penguin.{{config}}.search-skani.tsv',
    wildcard_constraints:
        assm_input_1="\w+",
        assm_input_2="\w+",
        config="\w+",
    params:
        skani_bin="/vol/cloud/louis/apps/skani/skani",
        reference=(
            lambda wildcards: f'{config["ref_folder"]}/{sample_table_read.reference[wildcards.smp]}'
        ),
        wdir=config["wdir"],
        rdir=config["rdir"] + "/assembly-skani",
        skaniDB=config["refdir"] + "/all_skani"
    log:
        skani_log=f'{config["rdir"]}/logs/assembly-skani/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.skani.penguin.{{config}}.log',
    benchmark:
        f'{config["rdir"]}/benchmarks/assembly-skani/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.skani.penguin.{{config}}.bmk'
    message:
        """--- skani """
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        N=$(grep -c '>' {input.contigs} || [[ $? == 1 ]])

        if [[ ${{N}} -eq 0 ]]; then
            touch {output.skani}
            touch {output.searchskani}
            exit 0
        fi

        {params.skani_bin} triangle -E -t 16 -i {input.contigs} -o {output.skani} --small-genomes

        {params.skani_bin} search --qi {input.contigs} -t 16 -d {params.skaniDB} -o {output.searchskani}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
