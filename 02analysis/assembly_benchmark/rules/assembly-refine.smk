rule refine_assemblies:
    input:
        contigs=config["rdir"] + "/assembly-combined/{smp}.assm.combined.fasta",
    output:
        contigs=config["rdir"] + "/assembly-refined/{smp}.assm.refined.fasta",
    threads: config["refineC_threads"]
    params:
        seqkit_bin=config["seqkit_bin"],
        refineC_split_parms=config["refineC_split_parms"],
        refineC_merge_parms=config["refineC_merge_parms"],
        refineC_bin=config["refineC_bin"],
        rdir=config["rdir"] + "/assembly-refined",
        wdir=config["wdir"],
        name="{smp}",
    log:
        config["rdir"] + "/logs/assembly-refined/{smp}.assm.refined.log",
    benchmark:
        config["rdir"] + "/benchmarks/assembly-refined/{smp}.assm.refined.bmk"
    conda:
        "../envs/refineC.yaml"
    message:
        """--- Merge megahit and penguin assemblies."""
    shell:
        """
        set -x

        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        rm -rf "{params.name}*"

        if [ ! -s {input.contigs} ]; then
            echo "No assemblies found"
            touch {output.contigs}
            exit 0
        fi

        N=$(grep -c '>' {input.contigs} || [[ $? == 1 ]])

        if [[ ${{N}} -lt 2 ]]; then
            cp {input.contigs} {output.contigs}
            exit 0
        fi

        {params.refineC_bin} split \
            --contigs {input.contigs} \
            --threads {threads} \
            --prefix {params.name} \
            --output {params.name} \
            {params.refineC_split_parms}

        if [[ ! -f {params.name}.split.fasta.gz ]]; then
            cp {input.contigs} {output.contigs}
            exit 0
        fi

        {params.refineC_bin} merge \
            --contigs "{params.name}.split.fasta.gz" \
            --threads {threads} \
            --prefix {params.name} \
            --output {params.name} \
            {params.refineC_merge_parms}

        if [[ -f {params.name}.merged.fasta.gz ]]; then
            zcat {params.name}.merged.fasta.gz > {output.contigs}
        else
            zcat {params.name}.split.fasta.gz > {output.contigs}
        fi

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
