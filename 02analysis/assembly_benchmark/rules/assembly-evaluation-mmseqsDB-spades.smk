rule assembly_mmseqsDB_evaluation_spades:
    input:
        #contigs=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.spades.{{config}}.fasta',
        #contigs=f'{config["rdir"]}/assembly-linclust/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}_rep_seq.fasta',
        reference=config["adir"] + "/all.faa",
    output:
        mappingtsv=f'{config["rdir"]}/assembly-mmseqsDB/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}_all.tsv',
    wildcard_constraints:
        assm_input_1="\w+",
        assm_input_2="\w+",
        config="\w+",
    params:
        tmp=f'{config["rdir"]}/assembly-mmseqsDB/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}_tmp',
        mmseqs_bin=config["mmseqs_bin"],
        wdir=config["wdir"],
        rdir=config["rdir"] + "/assembly-mmseqsDB",
        ref_folder=config["ref_folder"],
        contigs=f'{config["rdir"]}/assembly-annotation-eval/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.proteins.spades.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.proteins.spades.{{config}}.faa',
        contigsDB=f'{config["rdir"]}/assembly-mmseqsDB/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}_contigs_db',
        referenceID=f'{config["rdir"]}/assembly-mmseqsDB/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}_reference',
        reference_clust=f'{config["rdir"]}/assembly-mmseqsDB/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}_reference_rep_seq.fasta',
        reference_clust_DB=f'{config["rdir"]}/assembly-mmseqsDB/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}_reference_rep_seq',
        alignment=f'{config["rdir"]}/assembly-mmseqsDB/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}_aln',
        besthit=f'{config["rdir"]}/assembly-mmseqsDB/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}_besthit',
        besthittsv=f'{config["rdir"]}/assembly-mmseqsDB/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}_besthit.tsv',
    threads: config["mmseqs_threads"]
    log:
        mmseqs_log=f'{config["rdir"]}/assembly-mmseqsDB/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}.log',
    benchmark:
        f'{config["rdir"]}/benchmarks/assembly-mmseqsDB/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.spades.{{config}}.bmk'
    message:
        """--- mmseqs assemblies vs. reference spades """
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        N=$(grep -c '>' {params.contigs} || [[ $? == 1 ]])

        if [[ ${{N}} -eq 0 ]]; then
            touch {output.mappingtsv}
            exit 0
        fi

        {params.mmseqs_bin} easy-linclust \
            {input.reference} \
            {params.referenceID} \
            {params.tmp} \
            --min-seq-id 1 --kmer-per-seq 80 --cov-mode 1 -c 0.8 \
            --threads 32 --split-memory-limit 300G

        {params.mmseqs_bin} createdb \
            {params.reference_clust} \
            {params.reference_clust_DB}

        {params.mmseqs_bin} createdb \
            {params.contigs} \
            {params.contigsDB}

        {params.mmseqs_bin} map \
            {params.reference_clust_DB} \
            {params.contigsDB} \
            {params.alignment} \
            {params.tmp} \
            -a

        {params.mmseqs_bin} filterdb \
            {params.alignment} \
            {params.besthit} \
            --extract-lines 1

        {params.mmseqs_bin} convertalis \
            {params.reference_clust_DB} \
            {params.contigsDB} \
            {params.alignment} \
            {output.mappingtsv}

        {params.mmseqs_bin} convertalis \
            {params.reference_clust_DB} \
            {params.contigsDB} \
            {params.besthit} \
            {params.besthittsv}

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
