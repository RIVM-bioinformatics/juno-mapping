rule FilterMutectCalls:
    input:
        vcf=OUT + "/variants_raw/raw/{sample}.vcf",
        stats=OUT + "/variants_raw/raw/{sample}.vcf.stats",
        ref=OUT + "/reference/reference.fasta",
    output:
        vcf=OUT + "/variants_raw/FMC/{sample}.vcf",
    message:
        "Marking low confidence variants for {wildcards.sample}, based on FilterMutectCalls microbial mode"
    container:
        "docker://broadinstitute/gatk:4.3.0.0"
    conda:
        "../envs/gatk_picard.yaml"
    params:
        min_reads_per_strand=config["min_reads_per_strand"],
    log:
        OUT + "/log/FilterMutectCalls/{sample}.log",
    threads: config["threads"]["filter_variants"]
    resources:
        mem_gb=config["mem_gb"]["filter_variants"],
    shell:
        """
gatk FilterMutectCalls -V {input.vcf} \
-R {input.ref} \
-O {output.vcf} \
--stats {input.stats} \
--min-reads-per-strand {params.min_reads_per_strand} \
--microbial-mode 2>&1>{log}
        """


rule hard_filter_af:
    input:
        vcf=OUT + "/variants_raw/FMC/{sample}.vcf",
        ref=OUT + "/reference/reference.fasta",
    output:
        vcf=OUT + "/variants_raw/FMC_afhard/{sample}.vcf",
    message:
        "Hard filtering variants with very low allele frequency for {wildcards.sample}"
    container:
        "docker://staphb/bcftools:1.16"
    conda:
        "../envs/bcftools.yaml"
    params:
        min_af=config["hard_filter_minimum_allele_frequency"],
    log:
        OUT + "/log/hard_filter_af/{sample}.log",
    threads: config["threads"]["filter_variants"]
    resources:
        mem_gb=config["mem_gb"]["filter_variants"],
    shell:
        """
bcftools filter \
--exclude \"FORMAT/AF < {params.min_af}\" \
{input.vcf} \
1>{output.vcf} \
2>{log}
        """


rule soft_filter_af:
    input:
        vcf=OUT + "/variants_raw/FMC_afhard/{sample}.vcf",
        ref=OUT + "/reference/reference.fasta",
    output:
        vcf=OUT + "/variants_raw/FMC_afhard_afsoft/{sample}.vcf",
    message:
        "Marking minority variants for {wildcards.sample}"
    container:
        "docker://staphb/bcftools:1.16"
    conda:
        "../envs/bcftools.yaml"
    params:
        min_af=config["soft_filter_minimum_allele_frequency"],
    log:
        OUT + "/log/filter_af/{sample}.log",
    threads: config["threads"]["filter_variants"]
    resources:
        mem_gb=config["mem_gb"]["filter_variants"],
    shell:
        """
bcftools filter \
--exclude \"FORMAT/AF < {params.min_af}\" \
--soft-filter "min_af_{params.min_af}" \
{input.vcf} \
1>{output.vcf} \
2>{log}
        """


rule filter_depth:
    input:
        vcf=OUT + "/variants_raw/FMC_afhard_afsoft/{sample}.vcf",
        ref=OUT + "/reference/reference.fasta",
    output:
        vcf=OUT + "/variants_raw/FMC_afhard_afsoft_depth/{sample}.vcf",
    container:
        "docker://staphb/bcftools:1.16"
    conda:
        "../envs/bcftools.yaml"
    params:
        min_depth=config["minimum_depth"],
    log:
        OUT + "/log/filter_depth/{sample}.log",
    threads: config["threads"]["filter_variants"]
    resources:
        mem_gb=config["mem_gb"]["filter_variants"],
    shell:
        """
bcftools filter \
--exclude \"FORMAT/DP < {params.min_depth}\" \
--soft-filter "min_depth_{params.min_depth}" \
{input.vcf} \
1>{output.vcf} \
2>{log}
        """


if config["disable_mask"] == "True":
    # config value is parsed as str

    rule copy_mask:
        output:
            touch(OUT + "/variants_raw/no_mask.bed"),
        log:
            OUT + "/log/copy_mask.log",
        shell:
            """
    echo "Masking disabled, making empty file {output}" > {log}
            """

    rule filter_mask:
        input:
            vcf=OUT + "/variants_raw/FMC_afhard_afsoft_depth/{sample}.vcf",
            ref=OUT + "/reference/reference.fasta",
            mask=OUT + "/variants_raw/no_mask.bed",
        output:
            vcf=OUT + "/variants/{sample}.vcf",
        log:
            OUT + "/log/filter_depth/{sample}.log",
        threads: config["threads"]["filter_variants"]
        resources:
            mem_gb=config["mem_gb"]["filter_variants"],
        shell:
            """
cp {input.vcf} {output.vcf}
echo "Masking disabled, copying {input.vcf} to {output.vcf}" > {log}
            """

else:

    rule copy_mask:
        input:
            mask=config["mask_bed"],
        output:
            mask=OUT + "/variants_raw/mask.bed",
        log:
            OUT + "/log/copy_mask.log",
        threads: config["threads"]["other"]
        resources:
            mem_gb=config["mem_gb"]["other"],
        shell:
            """
    cp {input.mask} {output.mask}
            """

    rule filter_mask:
        input:
            vcf=OUT + "/variants_raw/FMC_afhard_afsoft_depth/{sample}.vcf",
            ref=OUT + "/reference/reference.fasta",
            mask=OUT + "/variants_raw/mask.bed",
        output:
            vcf=OUT + "/variants/{sample}.vcf",
        container:
            "docker://staphb/bcftools:1.16"
        conda:
            "../envs/bcftools.yaml"
        log:
            OUT + "/log/filter_depth/{sample}.log",
        threads: config["threads"]["filter_variants"]
        resources:
            mem_gb=config["mem_gb"]["filter_variants"],
        shell:
            """
bcftools filter \
--mask-file {input.mask} \
--soft-filter "masked_region" \
{input.vcf} \
1>{output.vcf} \
2>{log}
            """


rule remove_low_confidence_variants:
    input:
        vcf=OUT + "/variants/{sample}.vcf",
        ref=OUT + "/reference/reference.fasta",
    output:
        vcf=OUT + "/variants_clean/{sample}.vcf",
    message:
        "Remove low confidence variants for {wildcards.sample}"
    conda:
        "../envs/gatk_picard.yaml"
    container:
        "docker://broadinstitute/gatk:4.3.0.0"
    log:
        OUT + "/log/remove_low_confidence_variants/{sample}.log",
    threads: config["threads"]["filter_variants"]
    resources:
        mem_gb=config["mem_gb"]["filter_variants"],
    shell:
        """
gatk SelectVariants \
-V {input.vcf} \
-O {output.vcf} \
-R {input.ref} \
--exclude-filtered 2>&1>{log}
        """


rule select_snps:
    input:
        ref=OUT + "/reference/reference.fasta",
        vcf=OUT + "/variants_clean/{sample}.vcf",
    output:
        vcf=OUT + "/variants_snps_only/{sample}.snps.vcf",
    container:
        "docker://broadinstitute/gatk:4.3.0.0"
    conda:
        "../envs/gatk_picard.yaml"
    log:
        OUT + "/log/select_snps/{sample}.log",
    threads: config["threads"]["filter_variants"]
    resources:
        mem_gb=config["mem_gb"]["filter_variants"],
    shell:
        """
gatk SelectVariants \
-R {input.ref} \
-V {input.vcf} \
--select-type-to-include SNP \
--exclude-filtered \
-O {output.vcf} 2>&1>{log}
        """


# borrowed from pathogen-profiler:
# https://github.com/jodyphelan/pathogen-profiler/blob/v4.0.0/pathogenprofiler/bam.py#L64
rule filter_large_deletions:
    input:
        OUT + "/variants_raw/delly_raw/{sample}.bcf",
    output:
        filtered=OUT + "/deletions/{sample}.vcf",
    container:
        "docker://staphb/bcftools:1.16"
    conda:
        "../envs/bcftools.yaml"
    log:
        OUT + "/log/filter_large_deletions/{sample}.log",
    threads: config["threads"]["filter_variants"]
    resources:
        mem_gb=config["mem_gb"]["filter_variants"],
    shell:
        """
bcftools view -c 2 {input} -Oz 2>{log} |\
bcftools view -e '(INFO/END-POS)>=100000' -Ov -o {output.filtered} 2>&1>>{log}
        """
