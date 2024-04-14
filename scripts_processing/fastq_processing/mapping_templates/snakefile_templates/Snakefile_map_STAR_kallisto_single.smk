# Jared Bard 230403
# modified on 230503 to deal with different index files better
# takes single read fastqs, maps using kallisto assuming a mean fragment length of 200
FASTQ1=config["fastq1"]
SAMPLE=config["sample"]
STAR_INDEX_DIR=config["STAR_index_dir"]
KALLISTO_INDEX=config["kallisto_index"]
KALLISTO_DIRECTION=config["kallisto_direction"]
GTF=config["gtf"]

rule all:
    input:
        "mapped_reads/"+SAMPLE+"/"+SAMPLE+"_Aligned_Sorted.out.bam.bai",
        "kallisto_quant/"+SAMPLE,

rule trim_galore:
    input:
        read1=FASTQ1,
    output:
        read1=temp("trimmed/"+SAMPLE+"_trimmed.fq.gz"),
    log:
        "logs/trim_galore/"+SAMPLE+"_trim_galore.log"
    conda:
        "envs/star_kallisto_umi.yaml"
    threads: workflow.cores
    params:
        "--gzip"
    shell:
        "trim_galore {params} --fastqc_args '--outdir fastqc/' -j {threads} "
        "-o trimmed --basename "+SAMPLE+" "
        "{input.read1} &> {log}"

rule STAR_map_paired:
    input:
        gtf=GTF,
        idxdone = STAR_INDEX_DIR+"/Genome",
        read1="trimmed/"+SAMPLE+"_trimmed.fq.gz",
    output:
        temp("mapped_reads/"+SAMPLE+"/"+SAMPLE+"_Aligned.out.bam"),
        "mapped_reads/"+SAMPLE+"/"+SAMPLE+"_Log.final.out"
    params:
        out_prefix="mapped_reads/"+SAMPLE+"/"+SAMPLE+"_",
        genome_dir = STAR_INDEX_DIR,
        zip_params = "--readFilesCommand gunzip -c ",
        sjdb_params="--sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon CDS --sjdbGTFtagExonParentGene gene_id"
    log:
        "logs/STAR_map/"+SAMPLE+"_star_map.log"
    threads:
        workflow.cores
    conda:
        "envs/star_kallisto_umi.yaml"
    shell:
        r"""
        STAR --outSAMtype BAM Unsorted \
        {params.zip_params} \
        --sjdbGTFfile {input.gtf} {params.sjdb_params} \
        --runThreadN {threads} --alignMatesGapMax 20000 --limitBAMsortRAM 1445804817 \
        --genomeDir {params.genome_dir} --outFileNamePrefix {params.out_prefix} \
        --readFilesIn {input.read1} >& {log}
        """

rule samtools_sort:
    input:
        "mapped_reads/"+SAMPLE+"/"+SAMPLE+"_Aligned.out.bam"
    output:
        "mapped_reads/"+SAMPLE+"/"+SAMPLE+"_Aligned_Sorted.out.bam"
    log:
        "logs/samtools_sort/"+SAMPLE+"_samtools_sort.log"
    threads:
        workflow.cores
    conda:
        "envs/star_kallisto_umi.yaml"
    shell:
        "samtools sort -@ {threads} {input} -o {output} >& {log}"

rule samtools_index:
    input:
        "mapped_reads/"+SAMPLE+"/"+SAMPLE+"_Aligned_Sorted.out.bam"
    output:
        "mapped_reads/"+SAMPLE+"/"+SAMPLE+"_Aligned_Sorted.out.bam.bai"
    threads:
        workflow.cores
    conda:
        "envs/star_kallisto_umi.yaml"
    shell:
        "samtools index -@ {threads} {input}"

rule kallisto_map:
# takes deduplicated fastqs and uses kallisto to analyze transcript abundance
    input:
        index=KALLISTO_INDEX,
            read1="trimmed/"+SAMPLE+"_trimmed.fq.gz",
    output:
        directory("kallisto_quant/"+SAMPLE)
    log:
        "logs/kallisto_quant/"+SAMPLE+"_kallisto_map.log"
    conda:
        "envs/star_kallisto_umi.yaml"
    params:
        direction=KALLISTO_DIRECTION,
        single="--single -l 200 -s 1"
    shell:
        r"""
        kallisto quant -i {input.index} -o {output} \
        {params.single} {params.direction} --bootstrap-samples=50 -t {threads} \
        {input.read1} &> {log}
        """
        
        
