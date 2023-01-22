## i wrote this script to facilitate the processing of RNA modification (m6A) meta-analysis. 

rule all:
    input: 
        "pycra/mapped/bed/{n}.plus.sorted.bedgraph",
        "pycra/mapped/bed/{n}.minus.sorted.bedgraph"

rule reverse_complement:
    input: "{n}_1.trimmed.fastq"
    output: "{n}_1.trimmed.RC.fastq"
    shell: "seqkit seq -r -p {input} > {output}"

rule remove_duplicates:
    input:
        "{n}_1.trimmed.RC.fastq",
        "{n}_2.trimmed.fastq"
    output:
        "pycra/{n}_1.fasta",
        "pycra/{n}_2.fasta"
    shell: "pyFastqDuplicateRemover.py -f {input} -o {output}"

rule concatenate_fasta:
    input:
        "pycra/{n}_1.fasta",
        "pycra/{n}_2.fasta"
    output: "pycra/catfasta/{n}_R1R2.fasta"
    shell: "cat {input} > {output}"

rule bwa_mem:
    input: "pycra/catfasta/{n}_R1R2.fasta"
    output: "pycra/mapped/{n}.sam"
    params: 
        REF_PFX="/storehouse/george/m6a/dm6.fasta",
        READ_GRP="{n}",
    shell: "bwa mem -M -R '@RG\tID:{params.READ_GRP}\tPL:ILLUMINA\tSM:{params.READ_GRP}\tDS:pfx={params.REF_PFX}' {params.REF_PFX} {input} > {output}"

rule sam_to_bam:
    input: "pycra/mapped/{n}.sam"
    output: "pycra/mapped/{n}.bam"
    shell: "samtools view -b -S {input} -o {output}"

rule sort_bam:
    input: "pycra/mapped/{n}.bam"
    output: "pycra/mapped/{n}.sorted.bam"
    shell: "samtools sort {input} -o {output}"

rule index_bam:
    input: "pycra/mapped/{n}.sorted.bam"
    shell: "samtools index {input}"

rule create_plus_bam:
    input: "pycra/mapped/{n}.sorted.bam"
    output: "pycra/mapped/bed/{n}.plus.bam"
    shell: "samtools view -F 16 -b -o {output} {input}"

rule create_minus_bam:
    input: "pycra/mapped/{n}.sorted.bam"
    output
rule calculate_total_reads_plus:
    input: "pycra/mapped/bed/{n}.plus.bam"
    output: "pycra/mapped/bed/{n}.plus.total_reads"
    shell: "samtools view -c {input} > {output}"

rule calculate_scale_plus:
    input: "pycra/mapped/bed/{n}.plus.total_reads"
    output: "pycra/mapped/bed/{n}.plus.scale"
    shell: "echo 'scale=6; $(cat {input})/1000000' | bc > {output}"

rule create_plus_bedgraph:
    input:
        "pycra/mapped/bed/{n}.plus.bam",
        "pycra/mapped/bed/{n}.plus.scale"
    output: "pycra/mapped/bed/{n}.plus.bedgraph"
    shell: "bedtools genomecov -ibam {input[0]} -bg -scale $(cat {input[1]}) > {output}"

rule calculate_total_reads_minus:
    input: "pycra/mapped/bed/{n}.minus.bam"
    output: "pycra/mapped/bed/{n}.minus.total_reads"
    shell: "samtools view -c {input} > {output}"

rule calculate_scale_minus:
    input: "pycra/mapped/bed/{n}.minus.total_reads"
    output: "pycra/mapped/bed/{n}.minus.scale"
    shell: "echo 'scale=6; $(cat {input})/1000000' | bc > {output}"

rule create_minus_bedgraph:
    input:
        "pycra/mapped/bed/{n}.minus.bam",
        "pycra/mapped/bed/{n}.minus.scale"
    output: "pycra/mapped/bed/{n}.minus.bedgraph"
    shell: "bedtools genomecov -ibam {input[0]} -bg -scale $(cat {input[1]}) > {output}"

rule sort_plus_bedgraph:
    input: "pycra/mapped/bed/{n}.plus.bedgraph"
    output: "pycra/mapped/bed/{n}.plus.sorted.bedgraph"
    shell: "sort -k1,1 -k2,2n {input} > {output}"

rule sort_minus_bedgraph:
    input: "pycra/mapped/bed/{n}.minus.bedgraph"
    output: "pycra/mapped/bed/{n}.minus.sorted.bedgraph"
    shell: "sort -k1,1 -k2,2n {input} > {output}"

rule convert_minus_bedgraph:
    input: "pycra/mapped/bed/{n}.minus.sorted.bedgraph"
    output: "pycra/mapped/bed/{n}.minus.sorted_converted.bedgraph"
    shell: "awk '{ $4 *= -1 } 1' {input} > {output}"

rule create_plus_bigwig:
    input: 
        "pycra/mapped/bed/{n}.plus.sorted.bedgraph", 
        "/storehouse/george/m6a/dm6.chrom.sizes"
    output: "pycra/mapped/bed/bw/{n}.plus.bigwig"
    shell: "bedGraphToBigWig {input[0]} {input[1]} {output}"

rule create_minus_bigwig:
    input: 
        "pycra/mapped/bed/{n}.minus.sorted_converted.bedgraph",
        "/storehouse/george/m6a/dm6.chrom.sizes"
    output: "pycra/mapped/bed/bw/{n}.minus.bigwig"
    shell: "bedGraphToBigWig {input[0]} {input[1]} {output}"

rule delete_files:
    input:
        "pycra/mapped/bed/{n}.plus.bedgraph",
        "pycra/mapped/bed/{n}.minus.bedgraph",
        "pycra/mapped/bed/{n}.minus.sorted.bedgraph",
        "pycra/mapped/{n}.bam"
    shell: "rm {input}"
