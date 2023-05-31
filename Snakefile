import os
import subprocess
import pandas as pd

from glob import glob
from pathlib import Path
from snakemake.logging import logger

#base = "/vol/projects/CIIM/SIcohort/RNAseq"
base = "/vol/projects/nvanunen/code/rnaseq-snakemake"
#indir = base+"/raw"

indir = "/vol/projects/CIIM/SIcohort/RNAseq/raw"
#outdir = base+"/processed"
outdir = "/vol/projects/CIIM/SIcohort/RNAseq/processed"
envs = base+"/envs"
scripts = base+"/scripts"

logdir = outdir+"/logs"
bmdir = outdir+"/benchmarks"

# prepare input
ID,STIM,READ = glob_wildcards(indir+"/{id}-{stim}_{read}.fastq.gz")
logger.info("ID: " + ", ".join(set(ID)))
logger.info("STIM: " + ", ".join(set(STIM)))
logger.info("READ: " + ", ".join(set(READ)))

# get all stimulations
#STIM = [i.split("-")[-1] for i in set(SAMPLE)]
#id = [i.split("-")[0] for i in set(SAMPLE)]
#logger.info("STIMS: " + ", ".join(STIM))

# R1 = forward in 5' end 
# R2 = reverse in 3' end

localrules: all

rule all:
    input:  #expand(outdir+"/fastqc/{id}_{dir}_fastqc.{ext}", id=ID, dir=["R1", "R2"], ext=["html", "zip"]),
            #expand(outdir+"/trimmed/{id}_{dir}.fastq.gz", id=ID, dir=["R1", "R2"])
            expand(outdir+"/qc/{id}-{stim}_{read}_fastqc.{ext}", zip, id=ID, stim=STIM, read=READ, ext=["html", "zip"]),
            #expand(outdir+"/trimmed/{sample}_{sid}_{read}.fastq.gz", zip, sample=SAMPLE, sid=SID, read=READ),
            expand(outdir+"/trimmed-qc/{id}-{stim}_{read}_fastqc.{ext}", zip, id=ID, stim=STIM, read=READ, ext=["html", "zip"]),
            expand(outdir+"/count/{id}-{stim}.count", zip, id=ID, stim=STIM),
            #expand(outdir+"/align/{sample}_{sid}_Aligned.sortedByCoord.out.bam", zip, sample=SAMPLE, sid=SID),
            #expand(outdir+"/count_merge/{stim}.count", zip, stim=STIM),            expand(outdir+"/count_norm/{sample}_{sid}.count", zip, sample=SAMPLE, sid=SID),
            expand(outdir+"/count_norm/{stim}.count", stim=set(STIM))

rule qc:
    input:      indir+"/{id}-{stim}_{read}.fastq.gz"
    output:     expand(outdir+"/qc/{{id}}-{{stim}}_{{read}}_fastqc.{ext}", ext=["html", "zip"])
    params:     out = "qc", threads = 8
    log:        logdir+"/qc/{id}-{stim}_{read}.log"
    benchmark:  bmdir+"/qc/{id}-{stim}_{read}.txt"
    shell: "fastqc {input} -o {outdir}/{params.out} -t {params.threads}"

rule trim:
    input:      R1 = indir+"/{id}-{stim}_R1.fastq.gz",
                R2 = indir+"/{id}-{stim}_R2.fastq.gz"
    output:     R1 = outdir+"/trimmed/{id}-{stim}_R1.fastq.gz",
                R2 = outdir+"/trimmed/{id}-{stim}_R2.fastq.gz"
    log:        logdir+"/trim/{id}-{stim}.log"
    benchmark:  bmdir+"/trim/{id}-{stim}.txt"
    resources:  mem_mb = 10000, runtime = 120
    shell: "fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} --detect_adapter_for_pe"

use rule qc as qc_trim with:
    input:      rules.trim.output
    output:     expand(outdir+"/trimmed-qc/{{id}}-{{stim}}_{{read}}_fastqc.{ext}", ext=["html", "zip"])
    params:     out = "trimmed-qc", threads = 8
	log:    	logdir+"/qc_trim/{id}-{stim}_{read}.log" 
	benchmark: 	bmdir+"/qc_trim/{id}-{stim}_{read}.txt"

#rule multiqc:
#    input:  
#    output: 
#    shell: 
#        "multiqc {input}"

rule index:
    input:      fa = "/vol/projects/CIIM/refs/refdata-gex-GRCh38-2020-A/fasta/genome.fa",
                gtf = "/vol/projects/CIIM/refs/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    output:     outdir+"/ref_index/genomeParameters.txt"
    params:     out = outdir+"/ref_index", threads = 20, sjdbOverhang = 50
    log:    	logdir+"/index/index.log"
    benchmark: 	bmdir+"/index/index.txt"
    resources:  mem_mb = 200000, runtime = 180
    shell: "cd {params.out} && STAR --runMode genomeGenerate --genomeDir {params.out} --genomeFastaFiles {input.fa} \
        --sjdbGTFfile {input.gtf} --sjdbOverhang {params.sjdbOverhang} --runThreadN {params.threads}"

rule align:
    input:      R1 = rules.trim.output.R1,
                R2 = rules.trim.output.R2,
                index = rules.index.output
    output:     outdir+"/align/{id}-{stim}_Aligned.sortedByCoord.out.bam"
    params:     out = outdir+"/align", out_prefix = "{id}-{stim}_",
                refdir = outdir+"/ref_index", threads = 10,
                gtf = "/vol/projects/CIIM/refs/refdata-gex-GRCh38-2020-A/genes/genes.gtf",
                sam_type = "BAM SortedByCoordinate"
    log:    	logdir+"/align/{id}-{stim}.log" 
    benchmark: 	bmdir+"/align/{id}-{stim}.txt"
    resources:  mem_mb = 100000, runtime = 120
    shell: "cd {params.out} && STAR --genomeDir {params.refdir} --readFilesCommand zcat --readFilesIn {input.R1} {input.R2} \
        --outSAMtype {params.sam_type} --outFileNamePrefix {params.out_prefix} --sjdbGTFfile {params.gtf} --runThreadN {params.threads}"
 
rule count:
    input:      rules.align.output
    output:     outdir+"/count/{id}-{stim}.count"
    params:     gtf = "/vol/projects/CIIM/refs/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
    log:   	    logdir+"/count/{id}-{stim}.log"
    benchmark:  bmdir+"/count/{id}-{stim}.txt"
    resources:  mem_mb = 10000, runtime = 120
    conda:      envs+"/htseq.yaml"
    shell: "htseq-count -f bam -r pos -s no -t exon -i gene_id {input} {params.gtf} > {output}" # --additional-attr=gene_name

rule checkpoint:
    """This rule makes sure merge_counts is not run until every sample 
    and stim that exists has gone through the count rule first."""
    input:      expand(rules.count.output, zip, id=ID, stim=STIM)
    output:     outdir+"/count/checkpoint.txt"
    log:        logdir+"/checkpoint/checkpoint.log"
    benchmark:  bmdir+"/checkpoint/checkpoint.txt"
    shell: "touch {output}"

# merge all samples for each stimulation
rule merge_counts:
    input:      #outdir+"/count",
                rules.checkpoint.output #rules.intermediate.output #expand(rules.count.output, zip, ) #expand(outdir+"/count/{id}-{{stim}}.count", id=ID)
    output:     outdir+"/count_merge/{stim}.count"
    params:     indir = outdir+"/count"
    log:        logdir+"/merge_counts/{stim}.log"
    benchmark:  bmdir+"/merge_counts/{stim}.txt"
    run:
        merged_df = pd.DataFrame()
        for file in glob(os.path.join(f"{params.indir}/*-{wildcards.stim}.count")):
            print(f"Merging {file}...")
            sample = Path(file).stem.split("-")[0]
            df = pd.read_csv(file, index_col=0, names=["gene_id", sample], sep="\t")
            merged_df = pd.concat([merged_df, df], axis=1)
        merged_df.to_csv(output[0], sep="\t")

# Define rule for normalization
rule normalise_counts:
    input:      rules.merge_counts.output
    output:     outdir+"/count_norm/{stim}.count"
    log:        logdir+"/normalise_counts/{stim}.log"
    benchmark:  bmdir+"/normalise_counts/{stim}.txt"
    resources:  mem_mb = 10000, runtime = 120
    conda:      envs+"/edgeR.yaml"
    shell: "Rscript {scripts}/normalise_counts.R {input} {output}"

# Define rule for differential expression analysis
rule run_deseq:
    input:      rules.normalise_counts.output
    output:     outdir+"/deseq/{stim}.deseq"
    log:        logdir+"/run_deseq/{stim}.log"
    benchmark:  bmdir+"/run_deseq/{stim}.txt"
    resources:  mem_mb = 10000, runtime = 120
    #conda:      envs+"/.yaml"
    shell: "Rscript {scripts}/run_deseq.R {input} {output}"

# Define rule for generating summary report
#rule generate_report: