rule multiqc:
    input:
        expand(
            OUT + "/qc_raw_fastq/{sample}_{read}_fastqc.zip",
            sample=SAMPLES,
            read="R1 R2".split(),
        ),
        expand(
            OUT + "/qc_clean_fastq/{sample}_p{read}_fastqc.zip",
            sample=SAMPLES,
            read="R1 R2".split(),
        ),
        expand(OUT + "/clean_fastq/{sample}_fastp.json", sample=SAMPLES),
        expand(OUT + "/log/clean_fastq/clean_fastq_{sample}.log", sample=SAMPLES),
        expand(
            OUT + "/qc_mapping/insertsize/{sample}_metrics.txt",
            sample=SAMPLES,
        ),
        expand(
            OUT + "/identify_species/reads/{sample}/{sample}_bracken_species.kreport2",
            sample=SAMPLES,
        ),
        # expand(
        #     OUT + "/qc_mapping/samtools_stats/{sample}.txt",
        #     sample=SAMPLES,
        # ),
        expand(
            OUT + "/qc_mapping/CollectAlignmentSummaryMetrics/{sample}.txt",
            sample=SAMPLES,
        ),
        expand(
            OUT + "/qc_mapping/CollectGcBiasMetrics/{sample}.txt",
            sample=SAMPLES,
        ),
        expand(
            OUT + "/qc_mapping/CollectQualityYieldMetrics/{sample}.txt",
            sample=SAMPLES,
        ),
        expand(
            OUT + "/qc_mapping/CollectWgsMetrics/{sample}.txt",
            sample=SAMPLES,
        ),
        expand(
            OUT + "/qc_variant_calling/bcftools_stats/{sample}.txt",
            sample=SAMPLES,
        ),
        OUT + "/qc_variant_calling/report_filter_status_mqc.tsv"
    output:
        OUT + "/multiqc/multiqc.html",
        json = OUT + "/multiqc/multiqc_data/multiqc_data.json",
        txt = OUT + "/multiqc/multiqc_data/multiqc_fastqc.txt",
    message:
        "Making MultiQC report."
    conda:
        "../envs/multiqc.yaml"
    container:
        "docker://quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"
    threads: config["threads"]["multiqc"]
    resources:
        mem_gb=config["mem_gb"]["multiqc"],
    params:
        config_file = "config/multiqc_config.yaml",
        output_dir = OUT + "/multiqc",
    log:
        OUT + "/log/multiqc/multiqc.log",
    shell:
        """
        multiqc --interactive --force --config {params.config_file} \
            -o {params.output_dir} \
            -n multiqc.html {input} &> {log}
        """