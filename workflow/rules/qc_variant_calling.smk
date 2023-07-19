rule get_filter_status:
    input:
        vcf=OUT + "/variants_raw/FMC_af_depth_masked/{sample}.vcf",
    output:
        tsv=OUT + "/qc_variant_calling/get_filter_status/{sample}.tsv",
    message:
        "Writing filter status of variants to table for {wildcards.sample}"
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/gatk:4.3.0.0"
    log:
        OUT + "/log/get_filter_status/{sample}.log",
    threads: config["threads"]["filter_variants"]
    resources:
        mem_gb=config["mem_gb"]["filter_variants"],
    shell:
        """
gatk VariantsToTable -V {input.vcf} \
-F CHROM \
-F POS \
-F TYPE \
-F REF \
-F ALT \
-F DP \
-F FILTER \
--show-filtered \
-O {output.tsv} 2>&1>{log}
        """


rule combine_filter_status:
    input:
        tsv=expand(
            OUT + "/qc_variant_calling/get_filter_status/{sample}.tsv", sample=SAMPLES
        ),
        header="files/report_filter_status.mqc",
    output:
        OUT + "/qc_variant_calling/report_filter_status_mqc.tsv",
    message:
        "Combining variant QC reports"
    log:
        OUT + "/log/combine_filter_status.log",
    threads: config["threads"]["other"]
    resources:
        mem_gb=config["mem_gb"]["other"],
    shell:
        """
python workflow/scripts/combine_variant_tables.py --input {input.tsv} --output {output} --fields FILTER --mqc {input.header}
        """


rule bcftools_stats:
    input:
        ref=OUT + "/reference/reference.fasta",
        vcf=OUT + "/variants/{sample}.vcf",
    output:
        txt=OUT + "/qc_variant_calling/bcftools_stats/{sample}.txt",
    container:
        "docker://staphb/bcftools:1.16"
    conda:
        "../envs/bcftools.yaml"
    log:
        OUT + "/log/bcftools_stats/{sample}.log",
    threads: config["threads"]["filter_variants"]
    resources:
        mem_gb=config["mem_gb"]["filter_variants"],
    shell:
        """
bcftools stats \
--fasta-ref {input.ref} \
{input.vcf} \
1>{output.txt} \
2>{log}
        """
