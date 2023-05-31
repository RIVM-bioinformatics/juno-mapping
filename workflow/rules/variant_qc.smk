rule get_filter_status:
    input:
        vcf = OUT + "/variants_raw/af_FMC_depth_masked/{sample}.vcf",
    output:
        tsv = OUT + "/variant_qc/get_filter_status/{sample}.tsv"
    message: "Writing filter status of variants to table for {wildcards.sample}"
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/gatk:4.3.0.0"
    log:
        OUT + "/log/get_filter_status/{sample}.log"
    threads:
        config["threads"]["filter_variants"]
    resources:
        mem_gb = config["mem_gb"]["filter_variants"]
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
        expand(OUT + "/variant_qc/get_filter_status/{sample}.tsv", sample = SAMPLES)
    output:
        OUT + "/variant_qc/report_filter_status.tsv"
    message: "Combining variant QC reports"
    log:
        OUT + "/log/combine_filter_status.log"
    threads:
        config["threads"]["other"]
    resources:
        mem_gb = config["mem_gb"]["other"]
    shell:
        """
python workflow/scripts/combine_variant_tables.py --input {input} --output {output} --fields FILTER
        """
