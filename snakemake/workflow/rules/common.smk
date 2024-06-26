### packages ###
from snakemake.utils import validate
import pandas as pd

### read in sample sheet ###
sample_sheet = pd.read_csv(config["sample_sheet"]).set_index("sample_name", drop=False)
sample_sheet.index.names = ["sample_name"]

### get outputs ###
def get_rules(wildcards):
    all_rules = []
    if config["readqc"]:
        all_rules = all_rules + expand(
            "results/fastp_out/{sample}/{sample}.fastp.r1.fastq.gz", sample=sample_sheet["sample_name"])
        all_rules = all_rules + expand(
            "results/fastp_out/{sample}/{sample}.fastp.r2.fastq.gz", sample=sample_sheet["sample_name"])

    if config["readanalysis"]:
        all_rules = all_rules + expand(
            "results/kraken/{sample}/{sample}_kraken.txt", sample=sample_sheet["sample_name"])
        
    if config["coassembly"]:
        all_rules = all_rules + expand(
            "results/megahit_coassembly/final.contigs.fa", sample=sample_sheet["sample_name"])

        all_rules = all_rules + expand(
            "results/quast_coassembly/report.html", sample=sample_sheet["sample_name"])

        all_rules = all_rules + expand(
            "results/coverage_coassembly/{sample}/{sample}.sorted.bam", sample=sample_sheet["sample_name"])

    if config["binning"]:
        all_rules = all_rules + directory(expand(
            "results/metabat_coassembly/{sample}/", sample=sample_sheet["sample_name"]))

        all_rules = all_rules + expand(
            "results/metabat_checkm_coassembly/{sample}/{sample}_checkm_output.txt", sample=sample_sheet["sample_name"])
        
        all_rules = all_rules + directory(expand(
            "results/metabat_filt_coassembly/{sample}/", sample=sample_sheet["sample_name"]))
            
        all_rules = all_rules + directory(expand(
            "results/gtdbtk_coassembly/{sample}/", sample=sample_sheet["sample_name"]))

        all_rules = all_rules + expand(
            "results/prokka_coassembly/{sample}/aggregate.txt", sample=sample_sheet["sample_name"])
    
    if config["contiganalysis"]: 
        all_rules.append("results/pyrodigal/cds_proteins.faa")

        all_rules.append("results/pyrodigal/cds.gff")

        all_rules.append("results/diamond_contigs/uniref90_denit.tsv")

        all_rules = all_rules + expand(
            "results/pileup/{sample}/{sample}.pileup.txt", sample=sample_sheet["sample_name"])

        all_rules = all_rules + expand(
            "results/pileup/{sample}/{sample}.depth.txt", sample=sample_sheet["sample_name"])

    return all_rules


### helper functions ###
# get names of reads, files, etc

def get_read_path(wildcards):
    return sample_sheet.loc[wildcards.sample, ["sample_name", "forward", "reverse"]].dropna()

# raw reads
def get_r1(wildcards):
    tmp = get_read_path(wildcards)
    return tmp["forward"]

def get_r2(wildcards):
    tmp = get_read_path(wildcards)
    return tmp["reverse"]

# trimmed reads
def get_trimmed_r1(wildcards):
    return "results/fastp_out/{sample}/{sample}.fastp.r1.fastq.gz"

def get_trimmed_r2(wildcards):
    return "results/fastp_out/{sample}/{sample}.fastp.r2.fastq.gz"

# all trimmed reads
def get_all_r1(wildcards):
    return expand("results/fastp_out/{sample}/{sample}.fastp.r1.fastq.gz", sample=sample_sheet["sample_name"])

def get_all_r2(wildcards):
    return expand("results/fastp_out/{sample}/{sample}.fastp.r2.fastq.gz", sample=sample_sheet["sample_name"])

def get_all_bams(wildcards):
    return expand("results/megahit_coverage/{sample}/{sample}.sorted.bam", sample=sample_sheet["sample_name"])