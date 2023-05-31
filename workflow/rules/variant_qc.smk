rule filter_snps:
    input:
        ref = OUT + "/reference/reference.fasta",
        vcf = OUT + "/variants_raw/including_minority_variants/{sample}.vcf",
    output:
        vcf = OUT + "/variants/{sample}.vcf",
    conda:
        "../envs/bcftools.yaml"
    params:
        min_af = config["minimum_allele_frequency"],
        min_depth = config["minimum_depth"],
    log:
        OUT + "/log/filter_snps/{sample}.log"
    threads:
        config["threads"]["filter_variants"]
    resources:
        mem_gb = config["mem_gb"]["filter_variants"]
    shell:
        """
bcftools filter \
--include \"FORMAT/AF >= {params.min_af} && FORMAT/DP >= {params.min_depth}\" \
{input.vcf} \
1>{output.vcf} \
2>{log}
        """

rule select_snps:
    input:
        ref = OUT + "/reference/reference.fasta",
        vcf = OUT + "/variants/{sample}.vcf",
    output:
        vcf = OUT + "/variants_snps_only/{sample}.snps.vcf",
    conda:
        "../envs/gatk_picard.yaml"
    log:
        OUT + "/log/select_snps/{sample}.log"
    threads:
        config["threads"]["filter_variants"]
    resources:
        mem_gb = config["mem_gb"]["filter_variants"]
    shell:
        """
gatk SelectVariants \
-R {input.reference} \
-V {input.vcf} \
--select-type-to-include SNP \
--exclude-filtered \
-O {output.vcf} 2>&1>{log}
        """