rule assembly_metadmg_evaluation_megahit:
    input:
        bam_sorted=f'{config["rdir"]}/assembly-mapping/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mapping.megahit.{{config}}.sorted.bam',
    output:
        dfit=f'{config["rdir"]}/assembly-metadmg/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.metadmg.megahit.{{config}}.dfit.gz',
    wildcard_constraints:
        assm_input_1="\w+",
        assm_input_2="\w+",
        config="\w+",
    params:
        name=f'{config["rdir"]}/assembly-metadmg/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.metadmg.megahit.{{config}}',
        bam_name_sorted=f'{config["rdir"]}/assembly-metadmg/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.metadmg.megahit.{{config}}.name.sorted.bam',  
        metaDMG_bin=config["mdmg_cpp"],
        samtools_bin=config["samtools_bin"],
        wdir=config["wdir"],
        rdir=config["rdir"] + "/assembly-metadmg",
    threads: config["mdmg_threads"]
    log:
        metadmg_log=f'{config["rdir"]}/assembly-metadmg/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.metadmg.megahit.{{config}}.log',
    benchmark:
        f'{config["rdir"]}/benchmarks/assembly-metadmg/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.metadmg.megahit.{{config}}.bmk'
    message:
        """--- metaDMG local damage estimation of megahit assemblies. """
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        {params.samtools_bin} sort \
            -n \
            -@ 16 \
            {input.bam_sorted} > {params.bam_name_sorted}

        {params.metaDMG_bin} getdamage {params.bam_name_sorted} -p 25 --threads 16 -r 1 -o {params.name}
        {params.metaDMG_bin} dfit {params.name}.bdamage.gz --bam {params.bam_name_sorted} --showfits 2 --nopt 5 --seed 1234 --lib ds --nthreads 16 --out {params.name} >> {log.metadmg_log} 2>&1

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """

