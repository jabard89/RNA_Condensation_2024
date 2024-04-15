# Jared Bard 230403
# modified on 230503 to deal with different index files better
FASTQ1=config["fastq1"]
FASTQ2=config["fastq2"]
FASTQ3=config["fastq3"]
SAMPLE=config["sample"]
STAR_INDEX_DIR=config["STAR_index_dir"]
KALLISTO_INDEX=config["kallisto_index"]
KALLISTO_DIRECTION=config["kallisto_direction"]

GTF=config["gtf"]

rule all:
    input:
        "mapped_reads/"+SAMPLE+"/"+SAMPLE+"_Aligned_Sorted.out.bam.bai",
        "kallisto_quant/"+SAMPLE,

rule label_fastq:
# takes UMIs generated from NEBNext Multiple Oligos (UMI Adaptor RNA Set 1) and NEBnext Ultra II Directional RNA library Prep
# UMIs are store in fastq R2, which first contains 8 nt of sample barcode, before the 11nt UMI
# uses umi_tools to store the UMIs in the names of the fastq reads
    input:
        read1=FASTQ1,
        read2=FASTQ3,
        umi_read=FASTQ2,
    output:
        read1=temp("labeled_fastq/"+SAMPLE+"/"+SAMPLE+"_R1.umi.fastq"),
        read2=temp("labeled_fastq/"+SAMPLE+"/"+SAMPLE+"_R3.umi.fastq"),
    conda:
        "envs/star_kallisto_umi.yaml"
    shell:
        r"""
        umi_tools extract --bc-pattern=XXXXXXXXNNNNNNNNNNN -I {input.umi_read} \
        --read2-in={input.read1} --read2-out={output.read1} > /dev/null;
        umi_tools extract --bc-pattern=XXXXXXXXNNNNNNNNNNN -I {input.umi_read} \
        --read2-in={input.read2} --read2-out={output.read2} > /dev/null
        """

rule trim_galore:
    input:
        read1="labeled_fastq/"+SAMPLE+"/"+SAMPLE+"_R1.umi.fastq",
        read2="labeled_fastq/"+SAMPLE+"/"+SAMPLE+"_R3.umi.fastq",
    output:
        read1=temp("trimmed/"+SAMPLE+"_val_1.fq.gz"),
        read2=temp("trimmed/"+SAMPLE+"_val_2.fq.gz"),
    log:
        "logs/trim_galore/"+SAMPLE+"_trim_galore.log"
    conda:
        "envs/star_kallisto_umi.yaml"
    threads: workflow.cores
    params:
        "--paired --gzip"
    shell:
        "trim_galore {params} --fastqc_args '--outdir fastqc/' -j {threads} "
        "-o trimmed --basename "+SAMPLE+" "
        "{input.read1} {input.read2} &> {log}"

rule STAR_map_paired:
    input:
        gtf=GTF,
        idxdone = STAR_INDEX_DIR+"/Genome",
        read1="trimmed/"+SAMPLE+"_val_1.fq.gz",
        read2="trimmed/"+SAMPLE+"_val_2.fq.gz",
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
        --readFilesIn {input.read1} {input.read2} >& {log}
        """

rule samtools_sort_and_index:
    input:
        "mapped_reads/"+SAMPLE+"/"+SAMPLE+"_Aligned.out.bam"
    output:
        bam=temp("mapped_reads/"+SAMPLE+"/"+SAMPLE+"_Aligned_Sorted.out.bam"),
        index=temp("mapped_reads/"+SAMPLE+"/"+SAMPLE+"_Aligned_Sorted.out.bam.bai")
    log:
        "logs/samtools_sort/"+SAMPLE+"_samtools_sort.log"
    threads:
        workflow.cores
    conda:
        "envs/star_kallisto_umi.yaml"
    shell:
        "samtools sort -@ {threads} {input} -o {output.bam} >& {log}; "
        "samtools index -@ {threads} {output.bam}"

rule dedup:
# uses umi_tools to deduplicate reads: collapses reads that map to the same place in the genome and have the 
# same 11nt UMI
    input:
        bam="mapped_reads/"+SAMPLE+"/"+SAMPLE+"_Aligned_Sorted.out.bam",
        index="mapped_reads/"+SAMPLE+"/"+SAMPLE+"_Aligned_Sorted.out.bam.bai"
    output:
        "mapped_reads/"+SAMPLE+"/"+SAMPLE+"_Aligned.sortedByCoord.dedup.out.bam"
    log:
        error="logs/dedup/"+SAMPLE+"_dedup.err",
        log="logs/dedup/"+SAMPLE+"_dedup.log"
    threads:
        workflow.cores
    conda:
        "envs/star_kallisto_umi.yaml"
    shell:
        r"""
        umi_tools dedup --stdin={input.bam} \
        --chimeric-pairs=discard --unpaired-reads=discard --spliced-is-unique \
        --paired -S {output} --log={log.log} &> {log.error}
        """

rule dedup_index:
    input:
        "mapped_reads/"+SAMPLE+"/"+SAMPLE+"_Aligned.sortedByCoord.dedup.out.bam"
    output:
        "mapped_reads/"+SAMPLE+"/"+SAMPLE+"_Aligned.sortedByCoord.dedup.out.bam.bai"
    log:
        "logs/dedup_index/"+SAMPLE+".log"
    conda:
        "envs/star_kallisto_umi.yaml"
    shell:
        "samtools index -@ {threads} {input} {output} >& {log}"

rule split_dedup_bam:
# splits deduplicated BAM back out in order to input into kallisto
    input:
        "mapped_reads/"+SAMPLE+"/"+SAMPLE+"_Aligned.sortedByCoord.dedup.out.bam"
    output:
        read1=temp("mapped_reads/"+SAMPLE+"/"+SAMPLE+"_Aligned_dedup_R1.fastq.gz"),
        read2=temp("mapped_reads/"+SAMPLE+"/"+SAMPLE+"_Aligned_dedup_R3.fastq.gz")
    log:
        "logs/split_dedup_bam/"+SAMPLE+"_split_dedup_bam.log"
    threads:
        workflow.cores
    conda:
        "envs/star_kallisto_umi.yaml"
    shell:
        "samtools collate -@ {threads} -u -O {input} | "
        "samtools fastq -@ {threads} -1 {output.read1} -2 {output.read2} "
        "-0 /dev/null -s /dev/null -n >& {log}"

rule kallisto_map:
# takes deduplicated fastqs and uses kallisto to analyze transcript abundance
    input:
        index=KALLISTO_INDEX,
        read1="mapped_reads/"+SAMPLE+"/"+SAMPLE+"_Aligned_dedup_R1.fastq.gz",
        read2="mapped_reads/"+SAMPLE+"/"+SAMPLE+"_Aligned_dedup_R3.fastq.gz",
    output:
        directory("kallisto_quant/"+SAMPLE)
    log:
        "logs/kallisto_quant/"+SAMPLE+"_kallisto_map.log"
    conda:
        "envs/star_kallisto_umi.yaml"
    params:
        KALLISTO_DIRECTION,
    shell:
        r"""
        kallisto quant -i {input.index} -o {output} \
        {params} --bootstrap-samples=50 -t {threads} \
        {input.read1} {input.read2} &> {log}
        """
        
        
