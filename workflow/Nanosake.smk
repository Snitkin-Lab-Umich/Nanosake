# Author: Ali Pirani and Dhatri Badri
configfile: "config/config.yaml"

import pandas as pd
import os
import numpy as np

samples_df = pd.read_csv(config["samples"])
BARCODE = list(samples_df['barcode_id'])
SAMPLE = list(samples_df['sample_id'])
PREFIX = config["prefix"]

samples_df['combination'] = samples_df[['barcode_id', 'sample_id', 'sample_id']].agg('/'.join, axis=1)
COMBINATION = list(samples_df['combination'])

if not os.path.exists("results/" + PREFIX):
    os.system("mkdir %s" % "results/" + PREFIX)

rule all:
    input:
        flye_assembly = expand("results/{prefix}/flye/{combination}_flye.fasta", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        flye_circ_assembly = expand("results/{prefix}/flye/{combination}_flye_circ.fasta", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        medaka_out = expand("results/{prefix}/medaka/{combination}_medaka.fasta", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        polypolish = expand("results/{prefix}/polypolish/{combination}_flye_medaka_polypolish.fasta", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        flye_medaka_polypolish_wo_circ = expand("results/{prefix}/polypolish/{combination}_flye_medaka_polypolish_wo_circ.fasta", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        prokka = expand("results/{prefix}/prokka/{combination}_flye_medaka_polypolish.gff", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX,combination=COMBINATION),
        quast_out = expand("results/{prefix}/quast/{combination}_flye_medaka_polypolish/report.txt", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        busco_out = expand("results/{prefix}/busco/{combination}.flye_medaka_polypolish/busco_flye_medaka_polypolish.txt", barcode=BARCODE, sample=SAMPLE, prefix=PREFIX, combination=COMBINATION),
        
rule flye:
    input:
        trimmed = lambda wildcards: expand(str(config["long_reads"] + "/"  + f"{wildcards.barcode}"  + ".trimmed.fastq.gz")), 
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
    #conda:
    #    "envs/flye.yaml"
    singularity:
        "docker://staphb/flye:2.9.5"
    #envmodules:
    #    "Bioinformatics",
    #    "flye"
    shell:
        "module load Bioinformatics flye && flye --nano-hq {input.trimmed} -g {params.size} -o {params.assembly_dir} -t {params.threads} {params.flye_options} && cp {params.assembly_dir}/assembly.fasta {params.assembly_dir}/{params.prefix}_flye.fasta &>{log}"

rule flye_add_circ:
    input:
        flye_assembly = "results/{prefix}/flye/{barcode}/{sample}/{sample}_flye.fasta",
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
        # Load and process the assembly info
        assembly_info = pd.read_csv("%s/assembly_info.txt" % params.assembly_dir, sep='\t', header=0)
        assembly_info["circular"] = np.where(assembly_info["circ."] == "Y", "true", "false")
        # Build a replacement dictionary: {original_seq_name: new_seq_name}
        replacements = {
            row['#seq_name']: f"hybrid_{params.prefix}_{row['#seq_name']};circular={row['circular']}"
            for _, row in assembly_info.iterrows()
        }
        with open(f"{params.assembly_dir}/{params.prefix}_flye.fasta", "r") as infile, open(f"{params.assembly_dir}/{params.prefix}_flye_circ.fasta", "w") as outfile:
            for line in infile:
                if line.startswith(">"):
                    seq_name = line[1:].strip()
                    new_name = replacements.get(seq_name, seq_name)
                    outfile.write(f">{new_name}\n")
                else:
                    outfile.write(line)
        
rule medaka:
    input:
        trimmed = lambda wildcards: expand(str(config["long_reads"] + "/"  + f"{wildcards.barcode}" + ".trimmed.fastq.gz")), #+ "/" + f"{wildcards.barcode}" 
        #trimmed = lambda wildcards: expand(str("results/" + f"{wildcards.prefix}" + "/filtlong/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}" + ".trimmed.fastq.gz")),
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
    singularity:
        "docker://staphb/medaka:1.2.0"
    #envmodules:
    #    "Bioinformatics",
    #    "medaka",
    #    "bcftools"
    log:
        "logs/{prefix}/medaka/{barcode}/{sample}/{sample}.log"
    shell:
        """
        medaka_consensus -i {input.trimmed} -d {input.flye_assembly} -o {params.medaka_out_dir} -t {params.threads} -m {params.model} && 
        cp {params.medaka_out_dir}/consensus.fasta {params.medaka_out_dir}/{params.prefix}_medaka.fasta &>{log}
        """

rule bwaalign:
    input:
        #r1 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_paired.fastq.gz"),
        #r2 = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_paired.fastq.gz"),
        r1 = config["trimmed_reads"] + "/{sample}"  + "_R1_trim_paired.fastq.gz",
        r2 = config["trimmed_reads"] + "/{sample}"  + "_R2_trim_paired.fastq.gz",
        #r1_unpaired = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R1_unpaired.fastq.gz"),
        #r2_unpaired = lambda wildcards: expand(f"results/{wildcards.prefix}/trimmomatic/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_R2_unpaired.fastq.gz"),
        medaka_assembly = f"results/{{prefix}}/medaka/{{barcode}}/{{sample}}/{{sample}}_medaka.fasta"
    output:
        samout_1 = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}/{{sample}}_1.sam",
        samout_2 = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}/{{sample}}_2.sam",
        #samout_3 = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}/{{sample}}_3.sam",
        #samout_4 = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}/{{sample}}_4.sam",
    params:
        threads = config["ncores"],
    #conda:
    #    "envs/bwa.yaml"
    singularity:
        "docker://staphb/bwa:0.7.19"
    #envmodules:
    #    "Bioinformatics",
    #    "bwa"
    shell:
        "bwa index {input.medaka_assembly} && bwa mem -t12 -a {input.medaka_assembly} {input.r1} > {output.samout_1} && bwa mem -t12 -a {input.medaka_assembly} {input.r2} > {output.samout_2}"

rule polypolish:
    input:
        r1 = config["trimmed_reads"] + "/{sample}"  + "_R1_trim_paired.fastq.gz",
        r2 = config["trimmed_reads"] + "/{sample}"  + "_R2_trim_paired.fastq.gz",
        medaka_assembly = lambda wildcards: expand(f"results/{wildcards.prefix}/medaka/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}" + "_medaka.fasta"),
        samout_1 = lambda wildcards: expand(f"results/{wildcards.prefix}/polypolish/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_1.sam"),
        samout_2 = lambda wildcards: expand(f"results/{wildcards.prefix}/polypolish/" + f"{wildcards.barcode}/{wildcards.sample}/{wildcards.sample}_2.sam"),
    output:
        filtersam1 = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}_filtered_1.sam",
        filtersam2 = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}_filtered_2.sam",
        flye_medaka_polypolish = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}/{{sample}}_flye_medaka_polypolish.fasta",
        flye_medaka_polypolish_wo_circ = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}/{{sample}}_flye_medaka_polypolish_wo_circ.fasta",
    params:
        threads = config["ncores"],
    #conda:
    #    "envs/polypolish.yaml"
    singularity:
        "docker://staphb/polypolish:0.6.0"
    #envmodules:
    #    "Bioinformatics",
    #    "polypolish"
    shell:
        """
        polypolish filter --in1 {input.samout_1} --in2 {input.samout_2} --out1 {output.filtersam1} --out2 {output.filtersam2} && \
        polypolish polish {input.medaka_assembly} {output.filtersam1} {output.filtersam2} > {output.flye_medaka_polypolish} && \
        cp {output.flye_medaka_polypolish} {output.flye_medaka_polypolish_wo_circ} && \
        sed -i 's/;.*//' {output.flye_medaka_polypolish_wo_circ}
        """
     
rule prokka:
    input:
        flye_medaka_polypolish_wo_circ = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}/{{sample}}_flye_medaka_polypolish_wo_circ.fasta",
    output:
        flye_medaka_polypolish_annotation = f"results/{{prefix}}/prokka/{{barcode}}/{{sample}}/{{sample}}_flye_medaka_polypolish.gff",
    params:
        threads = config["ncores"],
        prefix = "{sample}",
        options = config["prokka_options"],
        prokka_dir = directory("results/{prefix}/prokka/{barcode}/{sample}/"),
    #conda:
    #    "envs/prokka.yaml"
    singularity:
        "docker://staphb/prokka:1.14.6"
    #envmodules:
    #    "Bioinformatics",
    #    "prokka"
    shell:
        """
        prokka {params.options} --strain {params.prefix} -outdir {params.prokka_dir} -prefix {params.prefix}_flye_medaka_polypolish {input.flye_medaka_polypolish_wo_circ} 
        """

rule quast:
    input:
        flye_assembly = f"results/{{prefix}}/flye/{{barcode}}/{{sample}}/{{sample}}_flye_circ.fasta",
        medaka_out = f"results/{{prefix}}/medaka/{{barcode}}/{{sample}}/{{sample}}_medaka.fasta",
        flye_medaka_polypolish = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}/{{sample}}_flye_medaka_polypolish.fasta"
    output:
        quast_out_flye = f"results/{{prefix}}/quast/{{barcode}}/{{sample}}/{{sample}}_flye/report.txt",
        quast_out_medaka = f"results/{{prefix}}/quast/{{barcode}}/{{sample}}/{{sample}}_medaka/report.txt",
        quast_out_medaka_polypolish = f"results/{{prefix}}/quast/{{barcode}}/{{sample}}/{{sample}}_flye_medaka_polypolish/report.txt",
    params:
        threads = config["ncores"],
        quast_dir = directory("results/{prefix}/quast/{barcode}/{sample}/{sample}"),
    #conda:
    #    "envs/quast.yaml"
    singularity:
        "docker://staphb/quast:5.3.0"
    #envmodules:
    #    "Bioinformatics",
    #    "quast"
    shell:
       """
       quast.py {input.flye_assembly} -o {params.quast_dir}_flye &&
       quast.py {input.medaka_out} -o {params.quast_dir}_medaka &&
       quast.py {input.flye_medaka_polypolish} -o {params.quast_dir}_flye_medaka_polypolish
       """

rule busco:
    input:
        flye_medaka_polypolish = f"results/{{prefix}}/polypolish/{{barcode}}/{{sample}}/{{sample}}_flye_medaka_polypolish.fasta",
    output:
        busco_out = f"results/{{prefix}}/busco/{{barcode}}/{{sample}}/{{sample}}.flye_medaka_polypolish/busco_flye_medaka_polypolish.txt", 
    params:
        busco_outpath = f"results/{{prefix}}/busco/{{barcode}}/{{sample}}/{{sample}}",
        flye_medaka_polypolish_busco_out = f"short_summary.specific.bacteria_odb12.{{sample}}.flye_medaka_polypolish.txt",
        flye_assembly_busco_out = f"short_summary.specific.bacteria_odb12.{{sample}}.flye_assembly.txt",
        threads = config["ncores"],
    #conda:
    #    "envs/busco.yaml"
    singularity:
        "docker://staphb/busco:5.8.2-prok-bacteria_odb12_2024-11-14"
    #envmodules:
    #    "Bioinformatics",
    #    "busco"
    shell:
        """
        set -e

        # Run BUSCO on polypolish assembly with retry
        for i in {{1..2}}; do
            echo "Attempt $i: Running BUSCO on polypolish assembly"
            busco -f -i {input.flye_medaka_polypolish} -m genome -l bacteria_odb12 -o {params.busco_outpath}.flye_medaka_polypolish && break || echo "BUSCO medaka attempt $i failed"
            sleep 10
        done

        # Check if BUSCO medaka succeeded
        if [ ! -f {params.busco_outpath}.flye_medaka_polypolish/{params.flye_medaka_polypolish_busco_out} ]; then
            echo "BUSCO medaka failed after 2 attempts" >&2
            exit 1
        fi

        cp {params.busco_outpath}.flye_medaka_polypolish/{params.flye_medaka_polypolish_busco_out} {params.busco_outpath}.flye_medaka_polypolish/busco_flye_medaka_polypolish.txt
        
        """

