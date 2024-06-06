rule get_cds:
    input:
         "results/megahit_coassembly/final.contigs.fa"
    output:
        protein = "results/pyrodigal/cds_proteins.faa",
        nucl = "results/pyrodigal/cds_nucl.fa",
        gff = "results/pyrodigal/cds.gff"
    threads: 6
    conda: "../envs/pyrodigal.yaml"
    resources:
        mem="30G",
        time="05:00:00"
    shell:
        """
        pyrodigal -i {input} -a {output.protein} -d {output.nucl} -f gff -o {output.gff} -p meta --jobs {threads}
        """

rule diamond_uniref90:
    input:
         "results/pyrodigal/cds_proteins.faa"
    output:
        "results/diamond_contigs/uniref90.tsv"
    params:
        db_loc = "/projects/b1052/databases/uniref90/uniref90_240317"
    threads: 10
    conda: "../envs/diamond.yaml"
    resources:
        mem="50G",
        time="05:00:00"
    shell:
        """
        diamond blastp -d {params.db_loc} -q {input} -o {output} \
        --very-sensitive --iterate --verbose --threads {threads} 
        """