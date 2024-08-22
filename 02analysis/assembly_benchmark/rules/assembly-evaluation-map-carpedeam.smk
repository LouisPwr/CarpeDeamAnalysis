rule assembly_map_evaluation_carpedeam:
    input:
        #contigs=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.carpedeam.{{config}}.fasta',
        contigs=f'{config["rdir"]}/assembly-linclust/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.carpedeam.{{config}}_rep_seq.fasta',
    output:
        mappingtsv=f'{config["rdir"]}/assembly-map/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.carpedeam.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.carpedeam.{{config}}_all.tsv',
    wildcard_constraints:
        assm_input_1="\w+",
        assm_input_2="\w+",
        config="\w+",
    params:
        tmp=f'{config["rdir"]}/assembly-map/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.carpedeam.{{config}}_tmp',
        mmseqs_bin=config["mmseqs_bin"],
        wdir=config["wdir"],
        rdir=config["rdir"] + "/assembly-map",
        ref_folder=config["ref_folder"],
        reference=config["adir"] + "/all.faa",
        contigsDB=f'{config["rdir"]}/assembly-map/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.carpedeam.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.carpedeam.{{config}}_contigs_db',
        targetDB=f'{config["rdir"]}/assembly-map/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.carpedeam.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.carpedeam.{{config}}_target_db',
        alignment=f'{config["rdir"]}/assembly-map/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.carpedeam.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.carpedeam.{{config}}_aln',
        besthit=f'{config["rdir"]}/assembly-map/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.carpedeam.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.carpedeam.{{config}}_besthit',
        besthittsv=f'{config["rdir"]}/assembly-map/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.carpedeam.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.carpedeam.{{config}}_besthit.tsv',
    threads: config["mmseqs_threads"]
    log:
        mmseqs_log=f'{config["rdir"]}/assembly-map/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.carpedeam.{{config}}.log',
    benchmark:
        f'{config["rdir"]}/benchmarks/assembly-map/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.carpedeam.{{config}}.bmk'
    message:
        """--- mmseqs assemblies vs. reference carpedeam """
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        N=$(grep -c '>' {input.contigs} || [[ $? == 1 ]])

        if [[ ${{N}} -eq 0 ]]; then
            touch {output.mappingtsv}
            exit 0
        fi

        {params.mmseqs_bin} createdb \
            {input.contigs} \
            {params.contigsDB}

        {params.mmseqs_bin} createdb \
            {params.reference} \
            {params.targetDB}

        {params.mmseqs_bin} map \
            {params.contigsDB} \
            {params.targetDB} \
            {params.alignment} \
            {params.tmp} \
            -a

        {params.mmseqs_bin} filterdb \
            {params.alignment} \
            {params.besthit} \
            --extract-lines 1

        {params.mmseqs_bin} convertalis \
            {params.contigsDB} \
            {params.targetDB} \
            {params.alignment} \
            {output.mappingtsv}

        {params.mmseqs_bin} convertalis \
            {params.contigsDB} \
            {params.targetDB} \
            {params.besthit} \
            {params.besthittsv}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
