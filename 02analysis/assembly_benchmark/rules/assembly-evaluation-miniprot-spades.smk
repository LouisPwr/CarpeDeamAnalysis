rule assembly_miniprot_evaluation_spades:
    input:
        contigs=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.spades.{{config}}.fasta',
    output:
        paf=f'{config["rdir"]}/assembly-miniprot/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}.paf',
    wildcard_constraints:
        assm_input_1="\w+",
        assm_input_2="\w+",
        config="\w+",
    params:
        miniprot_bin=config["miniprot_bin"],
        samtools_bin=config["samtools_bin"],
        wdir=config["wdir"],
        rdir=config["rdir"] + "/assembly-miniprot",
        reference=config["adir"] + "/all.faa",
        paf_filter=f'{config["rdir"]}/assembly-miniprot/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}.filter.paf',
    threads: config["minimap_threads"]
    log:
        minimap_log=f'{config["rdir"]}/assembly-miniprot/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}.log',
    benchmark:
        f'{config["rdir"]}/benchmarks/assembly-miniprot/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}.bmk'
    message:
        """--- miniprot assemblies vs. reference spades """
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        N=$(grep -c '>' {input.contigs} || [[ $? == 1 ]])

        if [[ ${{N}} -eq 0 ]]; then
            touch {output.paf}
            exit 0
        fi

        {params.miniprot_bin} \
            -S \
            -t 16 \
            {input.contigs} \
            {params.reference} > {output.paf}

        awk '$10 == $11 && ($11 / ($2 * 3)) >= 0.8' {output.paf} > {params.paf_filter}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """

