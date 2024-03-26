rule assembly_minimap_evaluation_penguin:
    input:
        contigs=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.penguin.{{config}}.fasta',
    output:
        sorted_sam=f'{config["rdir"]}/assembly-minimap/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.minimap.penguin.{{config}}.sorted.sam',
        sorted_md_bam=f'{config["rdir"]}/assembly-minimap/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.minimap.penguin.{{config}}.sorted_md.bam',
    wildcard_constraints:
        assm_input_1="\w+",
        assm_input_2="\w+",
        config="\w+",
    params:
        sam=f'{config["rdir"]}/assembly-minimap/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.minimap.penguin.{{config}}.sam',
        minimap_bin=config["minimap2_bin"],
        samtools_bin=config["samtools_bin"],
        wdir=config["wdir"],
        rdir=config["rdir"] + "/assembly-minimap",

        bt_index=f'{config["rdir"]}/assembly-minimap/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.minimap.penguin.{{config}}',

        ref_folder=config["ref_folder"],
        reference=(
            lambda wildcards: f'{config["ref_folder"]}/{sample_table_read.reference[wildcards.smp]}'
        )
    threads: config["minimap_threads"]
    log:
        minimap_log=f'{config["rdir"]}/assembly-minimap/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.minimap.penguin.{{config}}.log',
    benchmark:
        f'{config["rdir"]}/benchmarks/assembly-minimap/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.minimap.penguin.{{config}}.bmk'
    message:
        """--- minimap assemblies vs. reference penguin """
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        N=$(grep -c '>' {input.contigs} || [[ $? == 1 ]])

        if [[ ${{N}} -eq 0 ]]; then
            touch {output.sorted_sam}
            touch {output.sorted_md_bam}
            exit 0
        fi

        {params.minimap_bin} \
            -ax asm20 \
            -t 16 \
            -MD \
            -o {params.sam} \
            {params.reference} \
            {input.contigs}

        {params.samtools_bin} view -bS {params.sam} | samtools sort -o {output.sorted_sam} -
        {params.samtools_bin} calmd -Q -b {output.sorted_sam} {params.reference} > {output.sorted_md_bam}
        {params.samtools_bin} index {output.sorted_md_bam}

        rm -rf {params.sam}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """

