rule assembly_rna_evaluation_carpedeam:
    input:
        proteinsAA=f'{config["rdir"]}/assembly-annotation-eval/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.proteins.carpedeam.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.proteins.carpedeam.{{config}}.faa',
    output:
        tsv=f'{config["rdir"]}/assembly-rna/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.rna.carpedeam.{{config}}.tsv',
    wildcard_constraints:
        assm_input_1="\w+",
        assm_input_2="\w+",
        config="\w+",
    params:
        tmp=f'{config["rdir"]}/assembly-rna/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.rna.carpedeam.{{config}}_tmp',
        mmseqs_bin=config["mmseqs_bin"],
        features=f'{config["rdir"]}/assembly-annotation-eval/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.proteins.carpedeam.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.proteins.carpedeam.{{config}}.ffn',
        rna=f'{config["rdir"]}/assembly-rna/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.rna.carpedeam.{{config}}.fa',
        rna_clust=f'{config["rdir"]}/assembly-rna/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.rna.carpedeam.{{config}}.clust',
        wdir=config["wdir"],
        rdir=config["rdir"] + "/assembly-rna",
        #ref_folder=config["ref_folder"],
        reference=config["rna"],
    threads: config["threads_14"]
    log:
        mmseqs_log=f'{config["rdir"]}/assembly-rna/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.rna.carpedeam.{{config}}.log',
    benchmark:
        f'{config["rdir"]}/benchmarks/assembly-rna/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.rna.carpedeam.{{config}}.bmk'
    message:
        """--- mmseqs rna vs. database """
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        seqtk seq {params.features} | grep -A1 RNA > {params.rna}

        N=$(grep -c '>' {params.rna} || [[ $? == 1 ]])

        if [[ ${{N}} -eq 0 ]]; then
            touch {output.tsv}
            exit 0
        fi

        # {params.mmseqs_bin} easy-search \
        #     {params.rna} \
        #     {params.reference} \
        #     {output.tsv} \
        #     {params.tmp} \
        #     --max-seq-len 1000000 --search-type 3 \
        #     --threads 32 --split-memory-limit 300G >> {log.mmseqs_log} 2>&1

        {params.mmseqs_bin} easy-linclust \
            {params.rna} \
            {params.rna_clust} \
            {params.tmp} \
            -c 1 \
            --min-seq-id 0.999 \
            --threads {threads} --split-memory-limit 300G >> {log.mmseqs_log} 2>&1  

        {params.mmseqs_bin} easy-search \
            {params.rna_clust}_rep_seq.fasta \
            {params.reference} \
            {output.tsv} \
            {params.tmp} \
            --max-seq-len 1000000 --search-type 3 \
            --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov \
            --threads {threads} --split-memory-limit 300G >> {log.mmseqs_log} 2>&1

        # {params.mmseqs_bin} easy-search \
        #     {params.rna} \
        #     {params.reference} \
        #     {output.tsv} \
        #     {params.tmp} \
        #     --max-seq-len 1000000 --search-type 1 \
        #     --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,qcov,tcov \
        #     --threads 32 --split-memory-limit 300G >> {log.mmseqs_log} 2>&1

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
