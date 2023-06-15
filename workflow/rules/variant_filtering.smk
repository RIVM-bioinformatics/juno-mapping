rule filter_af:
    input:
        vcf=OUT + "/variants_raw/raw/{sample}.vcf",
        stats=OUT + "/variants_raw/raw/{sample}.vcf.stats",
        ref=OUT + "/reference/reference.fasta",
    output:
        vcf=OUT + "/variants_raw/af/{sample}.vcf",
    message:
        "Marking minority variants for {wildcards.sample}"
    container:
        "docker://staphb/bcftools:1.16"
    conda:
        "../envs/bcftools.yaml"
    params:
        min_af=config["minimum_allele_frequency"],
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


rule FilterMutectCalls:
    input:
        vcf=OUT + "/variants_raw/af/{sample}.vcf",
        stats=OUT + "/variants_raw/raw/{sample}.vcf.stats",
        ref=OUT + "/reference/reference.fasta",
    output:
        vcf=OUT + "/variants_raw/af_FMC/{sample}.vcf",
    message:
        "Marking low confidence variants for {wildcards.sample}, based on FilterMutectCalls microbial mode"
    container:
        "docker://broadinstitute/gatk:4.3.0.0"
    conda:
        "../envs/gatk_picard.yaml"
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
--microbial-mode 2>&1>{log}
        """


rule filter_depth:
    input:
        vcf=OUT + "/variants_raw/af_FMC/{sample}.vcf",
        ref=OUT + "/reference/reference.fasta",
    output:
        vcf=OUT + "/variants_raw/af_FMC_depth/{sample}.vcf",
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
        vcf=OUT + "/variants_raw/af_FMC_depth/{sample}.vcf",
        ref=OUT + "/reference/reference.fasta",
        mask=OUT + "/variants_raw/mask.bed",
    output:
        vcf=OUT + "/variants_raw/af_FMC_depth_masked/{sample}.vcf",
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
        vcf=OUT + "/variants_raw/af_FMC_depth_masked/{sample}.vcf",
        ref=OUT + "/reference/reference.fasta",
    output:
        vcf=OUT + "/variants/{sample}.vcf",
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
        vcf=OUT + "/variants/{sample}.vcf",
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
