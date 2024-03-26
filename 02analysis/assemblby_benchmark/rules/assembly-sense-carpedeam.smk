rule assembly_precision_sense_carpedeam:
    input:
        contigs=f'{config["rdir"]}/assembly/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.assm.carpedeam.{{config}}.fasta',
    output:
        tsv=f'{config["rdir"]}/assembly-sense/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mmseqs.carpedeam.{{config}}/dummy.txt',
    wildcard_constraints:
        assm_input_1="\w+",
        assm_input_2="\w+",
        config="\w+",
    params:
        tmp=f'{config["rdir"]}/assembly-sense/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mmseqs.carpedeam.{{config}}_tmp',
        result=f'{config["rdir"]}/assembly-sense/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mmseqs.carpedeam.{{config}}',
        mmseqs_bin=config["mmseqs_bin"],
        wdir=config["wdir"],
        rdir=config["rdir"] + "/assembly-sense",
        ref_folder=config["ref_folder"],
        reference=(
            lambda wildcards: f'{config["ref_folder"]}/{sample_table_read.reference[wildcards.smp]}'
        )
    threads: config["mmseqs_threads"]
    log:
        mmseqs_log=f'{config["rdir"]}/assembly-sense/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mmseqs.carpedeam.{{config}}.log',
    benchmark:
        f'{config["rdir"]}/benchmarks/assembly-sense/{{smp}}.{{assm_input_1}}-{{assm_input_2}}.mmseqs.carpedeam.{{config}}.bmk'
    message:
        """--- mmseqs assemblies vs. reference carpedeam """
    shell:
        """
        cd {params.rdir} || {{ echo "Cannot change dir"; exit 1; }}

        N=$(grep -c '>' {input.contigs} || [[ $? == 1 ]])

        if [[ ${{N}} -eq 0 ]]; then
            touch {output.tsv}
            exit 0
        fi

        EVALUATION="eval_script.sh"
        len=500

        tmp=${params.result}/$RANDOM
        mkdir -p $tmp

        mmseqs easy-linclust {params.reference} $tmp/"$(basename "{params.reference}")".CLUST $tmp --min-seq-id 0.99 -c 0.99 --remove-tmp-files 1 --threads 32

        mmseqs createdb {input.contigs} $tmp/"$(basename "{input.contigs}")".db
        mmseqs createdb {params.reference} $tmp/"$(basename "{params.reference}")".db 
        mmseqs createdb $tmp/"$(basename "{params.reference}")".CLUST_rep_seq.fasta $tmp/"$(basename "{params.reference}")".CLUST_rep_seq.db

        sh $EVALUATION $tmp/"$(basename "{input.contigs}")".db $tmp/"$(basename "{params.reference}")".db $tmp/"$(basename "{params.reference}")".CLUST_rep_seq.db {params.result} $len

        touch {output.tsv}

        rm -rf $tmp

        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
