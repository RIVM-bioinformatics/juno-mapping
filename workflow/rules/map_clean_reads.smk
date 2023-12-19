rule copy_ref:
    input:
        config["reference"],
    output:
        OUT + "/reference/reference.fasta",
    message:
        "Copying reference genome to output directory"
    shell:
        """
cp {input} {output}
        """


rule bwa_index_ref:
    input:
        OUT + "/reference/reference.fasta",
    output:
        OUT + "/reference/reference.fasta.sa",
    message:
        "Indexing reference genome using bwa"
    conda:
        "../envs/bwa_samtools.yaml"
    container:
        "docker://staphb/bwa:0.7.17"
    log:
        OUT + "/log/bwa_index_ref.log",
    threads: config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """
bwa index {input} 2>&1>{log}
       """


rule gatk_index_ref:
    input:
        OUT + "/reference/reference.fasta",
    output:
        OUT + "/reference/reference.dict",
    message:
        "Indexing reference genome using GATK"
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/gatk:4.3.0.0"
    log:
        OUT + "/log/gatk_index_ref.log",
    threads: config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """
gatk CreateSequenceDictionary -R {input} 2>&1>{log}
        """


rule samtools_index_ref:
    input:
        OUT + "/reference/reference.fasta",
    output:
        OUT + "/reference/reference.fasta.fai",
    message:
        "Indexing reference genome using samtools"
    conda:
        "../envs/bwa_samtools.yaml"
    container:
        "docker://staphb/samtools:1.17"
    log:
        OUT + "/log/samtools_index_ref.log",
    threads: config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """
samtools faidx {input} 2>&1>{log}
        """


rule minimap2:
    input:
        ref=OUT + "/reference/reference.fasta",
        r1=OUT + "/clean_fastq/{sample}_pR1.fastq.gz",
        r2=OUT + "/clean_fastq/{sample}_pR2.fastq.gz",
    output:
        sam=temp(OUT + "/mapped_reads/raw/{sample}.sam"),
    message:
        "Mapping reads for {wildcards.sample}"
    params:
        bases_per_batch="100000000",
        verbosity="3",
        softclip_supp_aln="-Y",
        rgid="'@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA'",
    conda:
        "../envs/bwa_samtools.yaml"
    container:
        "docker://staphb/minimap2:2.21"
    log:
        OUT + "/log/bwa_mem/{sample}.log",
    threads: config["threads"]["bwa"]
    resources:
        mem_gb=config["mem_gb"]["bwa"],
    shell:
        """
minimap2 \
-x sr \
-a \
-K {params.bases_per_batch} \
-t {threads} \
{params.softclip_supp_aln} \
-R {params.rgid} \
{input.ref} \
{input.r1} {input.r2} 2>{log} 1>{output}
        """


rule sam_to_sorted_bam:
    input:
        sam=OUT + "/mapped_reads/raw/{sample}.sam",
    output:
        bam=temp(OUT + "/mapped_reads/sorted/{sample}.bam"),
    message:
        "Convert sam to sorted bam for {wildcards.sample}"
    conda:
        "../envs/bwa_samtools.yaml"
    container:
        "docker://staphb/samtools:1.17"
    log:
        OUT + "/log/sam_to_sorted_bam/{sample}.log",
    threads: config["threads"]["samtools"]
    resources:
        mem_gb=config["mem_gb"]["samtools"],
    shell:
        """
samtools view -b -@ {threads} {input.sam} 2>{log} | \
samtools sort -@ {threads} - 1> {output.bam} 2>>{log}
        """


rule MarkDuplicates:
    input:
        OUT + "/mapped_reads/sorted/{sample}.bam",
    output:
        bam=OUT + "/mapped_reads/duprem/{sample}.bam",
        metrics=OUT + "/mapped_reads/duprem/{sample}.metrics",
    message:
        "Marking and removing optical duplicates for {wildcards.sample}"
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/picard:2.27.5"
    params:
        use_singularity=config["use_singularity"],
    threads: config["threads"]["picard"]
    resources:
        mem_gb=config["mem_gb"]["picard"],
    log:
        OUT + "/log/MarkDuplicates/{sample}.log",
    shell:
        """
if [ {params.use_singularity} == True ]
then
    EXEC=\"java -jar /usr/picard/picard.jar\"
else
    EXEC=picard
fi

$EXEC \
MarkDuplicates \
INPUT={input} \
OUTPUT={output.bam} \
METRICS_FILE={output.metrics} \
VALIDATION_STRINGENCY=SILENT \
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
ASSUME_SORT_ORDER="queryname" \
CLEAR_DT="false" \
ADD_PG_TAG_TO_READS=false 2>&1>{log}
        """


rule index_bam:
    input:
        bam=OUT + "/mapped_reads/duprem/{sample}.bam",
    output:
        bai=OUT + "/mapped_reads/duprem/{sample}.bam.bai",
    message:
        "Indexing bam file of {wildcards.sample}"
    conda:
        "../envs/bwa_samtools.yaml"
    container:
        "docker://staphb/samtools:1.17"
    threads: config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    log:
        OUT + "/log/index_bam/{sample}.log",
    shell:
        """
samtools index {input.bam} 2>&1>{log}
        """
