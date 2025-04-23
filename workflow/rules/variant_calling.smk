rule Mutect2:
    input:
        bam=OUT + "/mapped_reads/duprem/{sample}.bam",
        bai=OUT + "/mapped_reads/duprem/{sample}.bam.bai",
        ref=OUT + "/reference/reference.fasta",
        ref_gatk_index=OUT + "/reference/reference.dict",
        ref_samtools_index=OUT + "/reference/reference.fasta.fai",
    output:
        vcf=OUT + "/variants_raw/raw/{sample}.vcf",
        stats=OUT + "/variants_raw/raw/{sample}.vcf.stats",
    params:
        dangling_bases=1,
        annotations="--annotation StrandBiasBySample --annotation AlleleFraction",
    container:
        "docker://broadinstitute/gatk:4.3.0.0"
    conda:
        "../envs/gatk_picard.yaml"
    log:
        OUT + "/log/variant_calling_mutect2/{sample}.log",
    threads: config["threads"]["gatk"]
    resources:
        mem_gb=config["mem_gb"]["gatk"],
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

rule extract_mean_coverage:
    input:
        wgs_metrics=OUT + "/qc_mapping/CollectWgsMetrics/{sample}.txt",
    output:
        mean_coverage=OUT + "/coverage/{sample}_mean_coverage.txt",
    # conda:
    #     "../envs/pandas.yaml"
    resources:
        mem_gb=config["mem_gb"]["calculate_coverage"],
    log:
        OUT + "/log/extract_mean_coverage/{sample}.log",
    script:
        "../scripts/parse_wgs_metrics.py"

checkpoint check_mean_coverage:
    input:
        mean_coverage=OUT + "/coverage/{sample}_mean_coverage.txt",
    output:
        flag=OUT + "/coverage/{sample}_coverage_flag.txt",
    resources:
        mem_gb=config["mem_gb"]["calculate_coverage"],
    run:
        with open(input.mean_coverage) as f:
            mean_coverage = float(f.read().strip())
        if mean_coverage < 1:
            flag = "low"
        else:
            flag = "sufficient"
        with open(output.flag, "w") as f:
            f.write(flag)

rule delly_call:
    input:
        bam=OUT + "/mapped_reads/duprem/{sample}.bam",
        bai=OUT + "/mapped_reads/duprem/{sample}.bam.bai",
        ref=OUT + "/reference/reference.fasta",
        flag=lambda wildcards: checkpoints.check_mean_coverage.get(sample=wildcards.sample).output.flag,
    output:
        bcf=OUT + "/variants_raw/delly_raw/{sample}.bcf",
    container:
        "docker://quay.io/biocontainers/delly:1.2.6--hb7e2ac5_0"
    conda:
        "../envs/delly.yaml"
    log:
        OUT + "/log/variant_calling_delly/{sample}.log",
    threads: config["threads"]["delly"]
    resources:
        mem_gb=config["mem_gb"]["delly"],
    shell:
        """
        if [ $(cat {input.flag}) = "low" ]; then
            echo "large deletions not filtered for the sample" > {log}
            echo "##fileformat=VCFv4.2" > {output.bcf}
            echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> {output.bcf}
        else
            delly call -t DEL -g {input.ref} {input.bam} -o {output.bcf} 2>&1>{log}
        fi
        """

# rule delly_call:
#     input:
#         bam=OUT + "/mapped_reads/duprem/{sample}.bam",
#         bai=OUT + "/mapped_reads/duprem/{sample}.bam.bai",
#         ref=OUT + "/reference/reference.fasta",
#     output:
#         bcf=OUT + "/variants_raw/delly_raw/{sample}.bcf",
#     container:
#         "docker://quay.io/biocontainers/delly:1.2.6--hb7e2ac5_0"
#     conda:
#         "../envs/delly.yaml"
#     log:
#         OUT + "/log/variant_calling_delly/{sample}.log",
#     threads: config["threads"]["delly"]
#     resources:
#         mem_gb=config["mem_gb"]["delly"],
#     shell:
#         """
# delly call \
# -t DEL \
# -g {input.ref} \
# {input.bam} \
# -o {output.bcf} \
# 2>&1>{log}
#         """