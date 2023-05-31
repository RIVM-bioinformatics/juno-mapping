rule Mutect2:
    input:
        bam = OUT + "/mapped_reads/duprem/{sample}.bam",
        bai = OUT + "/mapped_reads/duprem/{sample}.bam.bai",
        ref = OUT + "/reference/reference.fasta",
        ref_gatk_index = OUT + "/reference/reference.dict",
        ref_samtools_index = OUT + "/reference/reference.fasta.fai",
    output:
        vcf = OUT + "/variants_raw/raw/{sample}.vcf",
        stats = OUT + "/variants_raw/raw/{sample}.vcf.stats",
    params:
        dangling_bases = 1,
        annotations = "--annotation StrandBiasBySample --annotation AlleleFraction",
    conda:
        "../envs/gatk_picard.yaml"
    threads: 4
    log:
        OUT + "/log/variant_calling/{sample}.log"
    threads:
        config["threads"]["gatk"]
    resources:
        mem_gb = config["mem_gb"]["gatk"]
    shell:
        """
gatk Mutect2 \
-R {input.ref} \
-I {input.bam} \
-O {output.vcf} \
{params.annotations} \
--num-matching-bases-in-dangling-end-to-recover {params.dangling_bases} \
--max-reads-per-alignment-start 75 2>&1>{log}
        """

rule FilterMutectCalls:
    input:
        vcf = OUT + "/variants_raw/raw/{sample}.vcf",
        stats = OUT + "/variants_raw/raw/{sample}.vcf.stats",
        ref = OUT + "/reference/reference.fasta",
    output:
        vcf = OUT + "/variants_raw/including_minority_variants/{sample}.vcf",
    conda:
        "../envs/gatk_picard.yaml"
    log:
        OUT + "/log/FilterMutectCalls/{sample}.log"
    threads:
        config["threads"]["filter_variants"]
    resources:
        mem_gb = config["mem_gb"]["filter_variants"]
    shell:
        """
gatk FilterMutectCalls -V {input.vcf} \
-R {input.ref} \
-O {output.vcf} \
--stats {input.stats} \
--microbial-mode 2>&1>{log}
        """

