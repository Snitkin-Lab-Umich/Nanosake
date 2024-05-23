# Author: Ali Pirani
configfile: "config/config.yaml"

import pandas as pd
import os
import numpy as np

samples_df = pd.read_csv(config["samples"])
BARCODE = list(samples_df['barcode_id'])
SAMPLE = list(samples_df['sample_id'])
PREFIX = config["prefix"]
SHORTREADS = list(samples_df['sample_id'])

samples_df['combination'] = samples_df[['barcode_id', 'sample_id', 'sample_id']].agg('/'.join, axis=1)
COMBINATION = list(samples_df['combination'])

if not os.path.exists("results/"):
    os.system("mkdir %s" % "results/")

if not os.path.exists("results/" + PREFIX):
    os.system("mkdir %s" % "results/" + PREFIX)

rule all:
    input:
        #pycoqc_report = expand("results/{prefix}/pycoqc/{prefix}.html", prefix=PREFIX),
        trimmed = expand("results/{prefix}/filtlong/{combination}.trimmed.fastq.gz", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
	    nanoplot = expand("results/{prefix}/nanoplot/{combination}_preqcNanoPlot-report.html", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        flye_assembly = expand("results/{prefix}/flye/{combination}_flye.fasta", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        flye_circ_assembly = expand("results/{prefix}/flye/{combination}_flye_circ.fasta", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        medaka_out = expand("results/{prefix}/medaka/{combination}_medaka.fasta", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        polypolish = expand("results/{prefix}/polypolish/{combination}_flye_medaka_polypolish.fasta", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        polypolish_unicycler = expand("results/{prefix}/polypolish_unicycler/{combination}_unicycler_polypolish.fasta", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        prokka = expand("results/{prefix}/prokka/{combination}_unicycler.gff", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX,combination=COMBINATION),
        quast_out = expand("results/{prefix}/quast/{combination}_unicycler/report.txt", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        busco_out = expand("results/{prefix}/busco/{combination}.unicycler/busco_unicycler.txt", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),

rule pycoqc:
    input:
        seq_summary = config["long_reads"] + "/sequencing_summary.txt",
    output:
       "results/{prefix}/pycoqc/{prefix}.html",
    log:
        "logs/{prefix}/pycoqc/{prefix}.log"        
    conda:
        "envs/pycoqc.yaml"
    log:
        "logs/{prefix}/pycoqc/{barcode}/{sample}/{sample}.log"  
    shell:
       'pycoQC -f {input.seq_summary} -o {output} &>{log}'

# Deprecated: Porechop is no longer being supported. 
# rule trim_nano_adaptors:
#     input:
#         longreads = config["long_reads"] + "/{barcode}/",
#         #samplename = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
#     output:
#         trimmed = "results/{prefix}/porechop/{barcode}/{sample}/{sample}.trimmed.fastq"
#     log:
#         "logs/{prefix}/porechop/{barcode}/{sample}/{sample}.log"        
#     conda:
#         "envs/porechop.yaml"
#     log:
#         "logs/{prefix}/porechop/{barcode}/{sample}/{sample}.log"
#     shell:
#         'porechop -i {input.longreads} -o {output.trimmed} -t 8 --discard_middle &>{log}'

rule nanopore_filtlong:
    input:
        longreads = config["long_reads"] + "/{barcode}/",
        #samplename = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
    output:
        trimmed = "results/{prefix}/filtlong/{barcode}/{sample}/{sample}.trimmed.fastq.gz"
    log:
        "logs/{prefix}/porechop/{barcode}/{sample}/{sample}.log"        
    conda:
        "envs/filtlong.yaml"
    #singularity:
        #"docker://staphb/filtlong:0.2.1"
    params:
        prefix="{sample}",
    log:
        "logs/{prefix}/porechop/{barcode}/{sample}/{sample}.log"
    shell:
        #'filtlong --min_length 1000 --keep_percent 95 {input.longreads}/*.fastq.gz | gzip > {output.trimmed}'
        'cat {input.longreads}/*.fastq.gz > /tmp/{params.prefix}.gz | filtlong --min_length 1000 --keep_percent 95 /tmp/{params.prefix}.gz | gzip > {output.trimmed} && rm /tmp/{params.prefix}.gz'

rule nanoplot:
    input:
        longreads = config["long_reads"] + "/{barcode}/",
        trimmed = lambda wildcards: expand(str("results/" + f"{wildcards.prefix}" + "/filtlong/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}" + ".trimmed.fastq.gz")),
    output:
        nanoplot_preqc = "results/{prefix}/nanoplot/{barcode}/{sample}/{sample}_preqcNanoPlot-report.html"
    log:
        "logs/{prefix}/nanoplot/{barcode}/{sample}/{sample}.log"
    params:
        outdir="results/{prefix}/nanoplot/{barcode}/{sample}",
        prefix="{sample}",
    # conda:
    #     "envs/nanoplot.yaml"
    singularity:
        "docker://staphb/nanoplot:1.42.0"
    shell:
        "cat {input.longreads}/*.fastq.gz > /tmp/{params.prefix}.gz && NanoPlot -o {params.outdir} -p {params.prefix}_preqc --tsv_stats --info_in_report --N50 --title {params.prefix}_preqc --fastq /tmp/{params.prefix}.gz && NanoPlot -o {params.outdir} -p {params.prefix}_postqc --tsv_stats --info_in_report --N50 --title {params.prefix}_postqc --fastq {input.trimmed} && rm /tmp/{params.prefix}.gz"

rule trimmomatic_pe:
    input:
        r1 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R1.fastq.gz")),
        r2 = lambda wildcards: expand(str(config["short_reads"] + "/" + f"{wildcards.sample}_R2.fastq.gz")),
        
    output:
        r1 = f"results/{{prefix}}/trimmomatic/{{barcode}}/{{sample}}/{{sample}}_R1_paired.fastq.gz",
        r2 = f"results/{{prefix}}/trimmomatic/{{barcode}}/{{sample}}/{{sample}}_R2_paired.fastq.gz", 
        # reads where trimming entirely removed the mate
        r1_unpaired = f"results/{{prefix}}/trimmomatic/{{barcode}}/{{sample}}/{{sample}}_R1_unpaired.fastq.gz",
        r2_unpaired = f"results/{{prefix}}/trimmomatic/{{barcode}}/{{sample}}/{{sample}}_R2_unpaired.fastq.gz",
    params:
        adapter_filepath=config["adapter_file"],
        seed=config["seed_mismatches"],
        palindrome_clip=config["palindrome_clipthreshold"],
        simple_clip=config["simple_clipthreshold"],
        minadapterlength=config["minadapterlength"],
        keep_both_reads=config["keep_both_reads"],
        window_size=config["window_size"],
        window_size_quality=config["window_size_quality"],
        minlength=config["minlength"],
        headcrop_length=config["headcrop_length"],
        threads = config["ncores"],
    log:
        "logs/{prefix}/trimmomatic/{barcode}/{sample}/{sample}.log"
    # conda:
    #     "envs/trimmomatic.yaml"
    singularity:
        "docker://staphb/trimmomatic:0.39"
    shell:
        "trimmomatic PE {input.r1} {input.r2} {output.r1} {output.r1_unpaired} {output.r2} {output.r2_unpaired} -threads {params.threads} ILLUMINACLIP:{params.adapter_filepath}:{params.seed}:{params.palindrome_clip}:{params.simple_clip}:{params.minadapterlength}:{params.keep_both_reads} SLIDINGWINDOW:{params.window_size}:{params.window_size_quality} MINLEN:{params.minlength} HEADCROP:{params.headcrop_length} &>{log}"

rule flye:
    input:
        trimmed = lambda wildcards: expand(str("results/" + f"{wildcards.prefix}" + "/filtlong/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}" + ".trimmed.fastq.gz")),
    output:
        assembly = f"results/{{prefix}}/flye/{{barcode}}/{{sample}}/{{sample}}_flye.fasta",
    params:
        assembly_dir = "results/{prefix}/flye/{barcode}/{sample}/",
        size = config["genome_size"],
        threads = config["threads"],
        flye_options = config["flye_options"],
        prefix = "{sample}",
    log:
        "logs/{prefix}/flye/{barcode}/{sample}/{sample}_flye.log"  
    # conda:
    #     "envs/flye.yaml"
    singularity:
        "docker://staphb/flye:2.9.3"
    shell:
        "flye --nano-hq {input.trimmed} -g {params.size} -o {params.assembly_dir} -t {params.threads} {params.flye_options} && cp {params.assembly_dir}/assembly.fasta {params.assembly_dir}/{params.prefix}_flye.fasta &>{log}"


rule flye_add_circ:
    input:
        flye_assembly = lambda wildcards: expand(str("results/" + f"{wildcards.prefix}" + "/flye/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}" + "_flye.fasta")),
    output:
        assembly = f"results/{{prefix}}/flye/{{barcode}}/{{sample}}/{{sample}}_flye_circ.fasta",
    params:
        assembly_dir = "results/{prefix}/flye/{barcode}/{sample}/",
        size = config["genome_size"],
        threads = config["threads"],
        flye_options = config["flye_options"],
        prefix = "{sample}",
    log:
        "logs/{prefix}/flye/{barcode}/{sample}/{sample}_flye.log"  
    run:
        shell("cp {params.assembly_dir}/{params.prefix}_flye.fasta {params.assembly_dir}/{params.prefix}_flye_circ.fasta")
        assembly_info = pd.read_csv("%s/assembly_info.txt" % params.assembly_dir, sep='\t', header=0)
        assembly_info["circular"] = np.where(assembly_info["circ."] == "Y", "true", "false")
        flye_assembly = "%s" % input.flye_assembly
        flye_assembly_circ = "%s" % output.assembly
        for index, row in assembly_info.iterrows():
            circular = "%s;circular=%s" % (row['#seq_name'], row['circular'])
            shell("sed -i 's/\<%s\>/%s/g' %s" % (row['#seq_name'], circular, flye_assembly_circ))

rule medaka:
    input:
        trimmed = lambda wildcards: expand(str("results/" + f"{wildcards.prefix}" + "/filtlong/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}" + ".trimmed.fastq.gz")),
        flye_assembly = lambda wildcards: expand(str("results/" + f"{wildcards.prefix}" + "/flye/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}" + "_flye_circ.fasta")),
    output:
        medaka_out = f"results/{{prefix}}/medaka/{{barcode}}/{{sample}}/{{sample}}_medaka.fasta"
    params:
        medaka_out_dir = "results/{prefix}/medaka/{barcode}/{sample}",
        threads = config["threads"],
        prefix = f"{{sample}}",
        model = config["medaka_model"],
    #conda:
    #    "envs/medaka.yaml"
    log:
        "logs/{prefix}/medaka/{barcode}/{sample}/{sample}.log"
    shell:
        "module load Bioinformatics medaka bcftools/1.12-g4b275e && medaka_consensus -i {input.trimmed} -d {input.flye_assembly} -o {params.medaka_out_dir} -t {params.threads} -m {params.model} && cp {params.medaka_out_dir}/consensus.fasta {params.medaka_out_dir}/{params.prefix}_medaka.fasta &>{log}"

rule bwaalign:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_paired.fastq.gz"),
        r1_unpaired = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_unpaired.fastq.gz"),
        r2_unpaired = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_unpaired.fastq.gz"),
        medaka_assembly = f"results/{{prefix}}/medaka/{{barcode}}/{{sample}}/{{sample}}_medaka.fasta"
    output:
        samout_1 = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}/{{sample}}_1.sam",
        samout_2 = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}/{{sample}}_2.sam",
        samout_3 = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}/{{sample}}_3.sam",
        samout_4 = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}/{{sample}}_4.sam",
    params:
        threads = config["ncores"],
    # conda:
    #     "envs/bwa.yaml"
    singularity:
        "docker://staphb/bwa:0.7.17"
    shell:
        "bwa index {input.medaka_assembly} && bwa mem -t12 -a {input.medaka_assembly} {input.r1} > {output.samout_1} && bwa mem -t12 -a {input.medaka_assembly} {input.r2} > {output.samout_2} && bwa mem -t12 -a {input.medaka_assembly} {input.r1_unpaired} > {output.samout_3} && bwa mem -t12 -a {input.medaka_assembly} {input.r2_unpaired} > {output.samout_4}"

# Deprecated: This rule was used for pilon 
#rule bwaalign_polypolish:
#    input:
#        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_paired.fastq.gz"),
#        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_paired.fastq.gz"),
        #r1_unpaired = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_unpaired.fastq.gz"),
        #r2_unpaired = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_unpaired.fastq.gz"),
#        polypolish_assembly = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}/{{sample}}_flye_medaka_polypolish.fasta"
#    output:
#       samout = f"results/{{prefix}}/pilon/{{barcode}}/{{sample}}/{{sample}}.sam",
#        bamout = f"results/{{prefix}}/pilon/{{barcode}}/{{sample}}/{{sample}}.bam",
#        bamout_sorted = f"results/{{prefix}}/pilon/{{barcode}}/{{sample}}/{{sample}}_sorted.bam",
#    params:
#        threads = config["ncores"],
    # conda:
    #     "envs/bwa.yaml"
#    singularity:
#        "docker://staphb/bwa:0.7.17"
#    shell:
#        "bwa index {input.polypolish_assembly} && bwa mem -t12 {input.polypolish_assembly} {input.r1} {input.r2} > {output.samout} && samtools view -Sb {output.samout} > {output.bamout} && samtools sort -o {output.bamout_sorted} {output.bamout} && samtools index {output.bamout_sorted}"

rule bwaalign_unicycler:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_paired.fastq.gz"),
        r1_unpaired = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_unpaired.fastq.gz"),
        r2_unpaired = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_unpaired.fastq.gz"),
        unicycler_assembly = f"results/{{prefix}}/unicycler/{{barcode}}/{{sample}}/{{sample}}_unicycler.fasta"
    output:
        samout_1 = f"results/{{prefix}}/polypolish_unicycler/{{barcode}}/{{sample}}/{{sample}}_1.sam",
        samout_2 = f"results/{{prefix}}/polypolish_unicycler/{{barcode}}/{{sample}}/{{sample}}_2.sam",
        samout_3 = f"results/{{prefix}}/polypolish_unicycler/{{barcode}}/{{sample}}/{{sample}}_3.sam",
        samout_4 = f"results/{{prefix}}/polypolish_unicycler/{{barcode}}/{{sample}}/{{sample}}_4.sam",
    params:
        threads = config["ncores"],
    # conda:
    #     "envs/bwa.yaml"
    singularity:
        "docker://staphb/bwa:0.7.17"
    shell:
        "bwa index {input.unicycler_assembly} && bwa mem -t12 -a {input.unicycler_assembly} {input.r1} > {output.samout_1} && bwa mem -t12 -a {input.unicycler_assembly} {input.r2} > {output.samout_2} && bwa mem -t12 -a {input.unicycler_assembly} {input.r1_unpaired} > {output.samout_3} && bwa mem -t12 -a {input.unicycler_assembly} {input.r2_unpaired} > {output.samout_4}"

rule polypolish:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_paired.fastq.gz"),
        medaka_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/medaka/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}" + "_medaka.fasta"),
        samout_1 = lambda wildcards: expand(f"results/{wildcards.prefix}/polypolish/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_1.sam"),
        samout_2 = lambda wildcards: expand(f"results/{wildcards.prefix}/polypolish/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_2.sam"),
    output:
        filtersam1 = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}_filtered_1.sam",
        filtersam2 = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}_filtered_2.sam",
        flye_medaka_polypolish = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}/{{sample}}_flye_medaka_polypolish.fasta",
    params:
        threads = config["ncores"],
    # conda:
    #     "envs/polypolish.yaml"
    singularity:
        "docker://staphb/polypolish:0.6.0"
    shell:
        "polypolish filter --in1 {input.samout_1} --in2 {input.samout_2} --out1 {output.filtersam1} --out2 {output.filtersam2} && polypolish polish {input.medaka_assembly} {output.filtersam1} {output.filtersam2} > {output.flye_medaka_polypolish}"

rule unicycler:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_paired.fastq.gz"),
        trimmed_long = "results/{prefix}/filtlong/{barcode}/{sample}/{sample}.trimmed.fastq.gz",
        r1_unpaired = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_unpaired.fastq.gz"),
        r2_unpaired = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_unpaired.fastq.gz"),
    output:
        unicycler_assembly = f"results/{{prefix}}/unicycler/{{barcode}}/{{sample}}/{{sample}}_unicycler.fasta",
    params:
        unicycler_assembly_out = "results/{prefix}/unicycler/{barcode}/{sample}",
        threads = config["ncores"],
        prefix = f"{{sample}}",
    # conda:
    #     "envs/unicycler.yaml"
    singularity:
        "docker://staphb/unicycler:0.5.0"
    shell:
        "unicycler -1 {input.r1} -2 {input.r2} -s {input.r1_unpaired} -s {input.r2_unpaired} -l {input.trimmed_long} -o {params.unicycler_assembly_out} && cp {params.unicycler_assembly_out}/assembly.fasta {params.unicycler_assembly_out}/{params.prefix}_unicycler.fasta && sed -i 's/length=.*depth=.*circular/circular/g' {params.unicycler_assembly_out}/{params.prefix}_unicycler.fasta && sed -i 's/ circular/;circular/g' {params.unicycler_assembly_out}/{params.prefix}_unicycler.fasta"

rule polypolish_unicycler:
    input:
        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_paired.fastq.gz"),
        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_paired.fastq.gz"),
        unicycler_assembly = f"results/{{prefix}}/unicycler/{{barcode}}/{{sample}}/{{sample}}_unicycler.fasta",
        samout_1 = lambda wildcards: expand(f"results/{wildcards.prefix}/polypolish_unicycler/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_1.sam"),
        samout_2 = lambda wildcards: expand(f"results/{wildcards.prefix}/polypolish_unicycler/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_2.sam"),
    output:
        filtersam1 = f"results/{{prefix}}/polypolish_unicycler/{{barcode}}/{{sample}}/{{sample}}_filtered_1.sam",
        filtersam2 = f"results/{{prefix}}/polypolish_unicycler/{{barcode}}/{{sample}}/{{sample}}_filtered_2.sam",
        unicycler_polypolish = f"results/{{prefix}}/polypolish_unicycler/{{barcode}}/{{sample}}/{{sample}}_unicycler_polypolish.fasta",
    params:
        threads = config["ncores"],
    # conda:
    #     "envs/polypolish.yaml"
    singularity:
        "docker://staphb/polypolish:0.6.0"
    shell:
        "polypolish filter --in1 {input.samout_1} --in2 {input.samout_2} --out1 {output.filtersam1} --out2 {output.filtersam2} && polypolish polish {input.unicycler_assembly} {output.filtersam1} {output.filtersam2} > {output.unicycler_polypolish}"

# Deprecated: This rule was used for pilon
#rule bwaalign_polypolish_unicycler:
#    input:
#        r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_paired.fastq.gz"),
#        r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_paired.fastq.gz"),
#        polypolish_unicycler_assembly = f"results/{{prefix}}/polypolish_unicycler/{{barcode}}/{{sample}}/{{sample}}_unicycler_polypolish.fasta"
#    output:
#        samout = f"results/{{prefix}}/pilon_unicycler/{{barcode}}/{{sample}}/{{sample}}.sam",
#        bamout = f"results/{{prefix}}/pilon_unicycler/{{barcode}}/{{sample}}/{{sample}}.bam",
#        bamout_sorted = f"results/{{prefix}}/pilon_unicycler/{{barcode}}/{{sample}}/{{sample}}_sorted.bam",
#    params:
#        threads = config["ncores"],
    # conda:
    #     "envs/bwa.yaml"
#    singularity:
#        "docker://staphb/bwa:0.7.17"
#    shell:
#        "bwa index {input.polypolish_unicycler_assembly} && bwa mem -t12 {input.polypolish_unicycler_assembly} {input.r1} {input.r2} > {output.samout} && samtools view -Sb {output.samout} > {output.bamout} && samtools sort -o {output.bamout_sorted} {output.bamout} && samtools index {output.bamout_sorted}"

rule prokka:
    input:
        unicycler_polypolish = f"results/{{prefix}}/polypolish_unicycler/{{barcode}}/{{sample}}/{{sample}}_unicycler_polypolish.fasta",
        unicycler_assembly = f"results/{{prefix}}/unicycler/{{barcode}}/{{sample}}/{{sample}}_unicycler.fasta",
        flye_medaka_polypolish = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}/{{sample}}_flye_medaka_polypolish.fasta",
        flye_assembly = f"results/{{prefix}}/flye/{{barcode}}/{{sample}}/{{sample}}_flye_circ.fasta",
    output:
        unicycler_annotation = f"results/{{prefix}}/prokka/{{barcode}}/{{sample}}/{{sample}}_unicycler.gff",
    params:
        threads = config["ncores"],
        prefix = "{sample}",
        options = config["prokka_options"],
        prokka_dir = directory("results/{prefix}/prokka/{barcode}/{sample}/"),
    # conda:
    #     "envs/prokka.yaml"
    singularity:
        "docker://staphb/prokka:1.14.6"
    shell:
        "prokka {params.options} --strain {params.prefix} -outdir {params.prokka_dir} -prefix {params.prefix}_unicycler {input.unicycler_assembly} && prokka {params.options} --strain {params.prefix} -outdir {params.prokka_dir} -prefix {params.prefix}_flye {input.flye_assembly} && prokka {params.options} --strain {params.prefix} -outdir {params.prokka_dir} -prefix {params.prefix}_flye_medaka_polypolish {input.flye_medaka_polypolish} && prokka {params.options} --strain {params.prefix} -outdir {params.prokka_dir} -prefix {params.prefix}_unicycler_polypolish {input.unicycler_polypolish}"

rule quast:
    input:
        unicycler_polypolish = f"results/{{prefix}}/polypolish_unicycler/{{barcode}}/{{sample}}/{{sample}}_unicycler_polypolish.fasta",
        unicycler_assembly = f"results/{{prefix}}/unicycler/{{barcode}}/{{sample}}/{{sample}}_unicycler.fasta",
        flye_medaka_polypolish = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}/{{sample}}_flye_medaka_polypolish.fasta",
        flye_assembly = f"results/{{prefix}}/flye/{{barcode}}/{{sample}}/{{sample}}_flye_circ.fasta",
    output:
        quast_out = f"results/{{prefix}}/quast/{{barcode}}/{{sample}}/{{sample}}_unicycler/report.txt",
    params:
        threads = config["ncores"],
        quast_dir = directory("results/{prefix}/quast/{barcode}/{sample}/{sample}"),
    # conda:
    #     "envs/quast.yaml"
    singularity:
        "docker://staphb/quast:5.0.2"
    shell:
        "quast.py {input.unicycler_assembly} -o {params.quast_dir}_unicycler && quast.py {input.unicycler_polypolish} -o {params.quast_dir}_unicycler_polypolish && quast.py {input.flye_medaka_polypolish} -o {params.quast_dir}_flye_medaka_polypolish && quast.py {input.flye_assembly} -o {params.quast_dir}_flye"

rule busco:
    input:
        unicycler_polypolish = f"results/{{prefix}}/polypolish_unicycler/{{barcode}}/{{sample}}/{{sample}}_unicycler_polypolish.fasta",
        unicycler = f"results/{{prefix}}/unicycler/{{barcode}}/{{sample}}/{{sample}}_unicycler.fasta",
        flye_medaka_polypolish = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}/{{sample}}_flye_medaka_polypolish.fasta",
        flye_assembly = f"results/{{prefix}}/flye/{{barcode}}/{{sample}}/{{sample}}_flye.fasta",
    output:
        busco_out = f"results/{{prefix}}/busco/{{barcode}}/{{sample}}/{{sample}}.unicycler/busco_unicycler.txt",
    params:
        busco_outpath = f"results/{{prefix}}/busco/{{barcode}}/{{sample}}/{{sample}}",
        unicycler_busco_out = f"short_summary.specific.bacteria_odb10.{{sample}}.unicycler.txt",
        unicycler_polypolish_busco_out = f"short_summary.specific.bacteria_odb10.{{sample}}.unicycler_polypolish.txt",
        flye_medaka_polypolish_busco_out = f"short_summary.specific.bacteria_odb10.{{sample}}.flye_medaka_polypolish.txt",
        medaka_assembly_busco_out = f"short_summary.specific.bacteria_odb10.{{sample}}.medaka_assembly.txt",
        flye_assembly_busco_out = f"short_summary.specific.bacteria_odb10.{{sample}}.flye_assembly.txt",
        flye_medaka_polypolish_pilon_busco_out = f"short_summary.specific.bacteria_odb10.{{sample}}.flye_medaka_polypolish_pilon.txt",
        unicycler_polypolish_pilon_busco_out = f"short_summary.specific.bacteria_odb10.{{sample}}.unicycler_polypolish_pilon.txt",
        threads = config["ncores"],
    # conda:
    #     "envs/busco.yaml"
    singularity:
        "docker://staphb/busco:latest"
    shell:
        "busco -f -i {input.unicycler} -m genome -l bacteria_odb10 -o {params.busco_outpath}.unicycler && cp {params.busco_outpath}.unicycler/{params.unicycler_busco_out} {params.busco_outpath}.unicycler/busco_unicycler.txt && busco -f -i {input.unicycler_polypolish} -m genome -l bacteria_odb10 -o {params.busco_outpath}.unicycler_polypolish && cp {params.busco_outpath}.unicycler_polypolish/{params.unicycler_polypolish_busco_out} {params.busco_outpath}.unicycler_polypolish/busco_unicycler_polypolish.txt && busco -f -i {input.flye_medaka_polypolish} -m genome -l bacteria_odb10 -o {params.busco_outpath}.flye_medaka_polypolish && cp {params.busco_outpath}.flye_medaka_polypolish/{params.flye_medaka_polypolish_busco_out} {params.busco_outpath}.flye_medaka_polypolish/busco_flye_medaka_polypolish.txt && busco -f -i {input.flye_assembly} -m genome -l bacteria_odb10 -o {params.busco_outpath}.flye_assembly && cp {params.busco_outpath}.flye_assembly/{params.flye_assembly_busco_out} {params.busco_outpath}.flye_assembly/busco_flye_assembly.txt"

rule multiqc:
    input:
        quast_out = f"results/{{prefix}}/quast/",
        busco_dir_out = f"results/{{prefix}}/busco/",
        unicycler_annotation = f"results/{{prefix}}/prokka/",
    output:
        multiqc_out = f"results/{{prefix}}/multiqc/multiqc_report.html",
    params:
        multiqc_out_dir = directory("results/{prefix}/multiqc/"),
        prokka_dir_out = directory("results/{prefix}/prokka"),
        quast_dir_out = directory("results/{prefix}/quast"),
        #pycoqc_dir_out = directory("results/{prefix}/pycoqc"),
        busco_dir_out = directory("results/{prefix}/busco"),
        threads = config["ncores"],
    priority: 100000    
    # conda:
    #     "envs/multiqc.yaml"
    singularity:
        "docker://staphb/multiqc:1.8"
    shell:
        #{params.pycoqc_dir_out}
        "multiqc -o {params.multiqc_out_dir} -f {params.prokka_dir_out} {params.busco_dir_out} {params.quast_dir_out}"

