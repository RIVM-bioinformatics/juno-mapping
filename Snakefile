import yaml


sample_sheet = config["sample_sheet"]
with open(sample_sheet) as f:
    SAMPLES = yaml.safe_load(f)

for param in ["threads", "mem_gb"]:
    for k in config[param]:
        config[param][k] = int(config[param][k])

# print(SAMPLES)

OUT = config["output_dir"]


localrules:
    all,


include: "workflow/rules/fastqc_raw_data.smk"
include: "workflow/rules/clean_fastq.smk"
include: "workflow/rules/fastqc_clean_data.smk"
include: "workflow/rules/identify_species.smk"
include: "workflow/rules/map_clean_reads.smk"
include: "workflow/rules/qc_mapping.smk"
include: "workflow/rules/variant_calling.smk"
include: "workflow/rules/variant_filtering.smk"
include: "workflow/rules/qc_variant_calling.smk"
include: "workflow/rules/multiqc_report.smk"


rule all:
    input:
        vcf=expand(OUT + "/variants/{sample}.vcf", sample=SAMPLES),
        vcf_snps=expand(OUT + "/variants_snps_only/{sample}.snps.vcf", sample=SAMPLES),
        multiqc=OUT + "/multiqc/multiqc.html",
        filter_report=OUT + "/qc_variant_calling/report_filter_status_mqc.tsv",
