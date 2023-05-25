import yaml


sample_sheet=config["sample_sheet"]
with open(sample_sheet) as f:
    SAMPLES = yaml.safe_load(f)

for param in ["threads", "mem_gb"]:
    for k in config[param]:
        config[param][k] = int(config[param][k])

# print(SAMPLES)

OUT = config["output_dir"]

localrules:
    all,


include: "workflow/rules/rule.smk"


rule all:
    input:
        expand(OUT + "/{sample}_combined.fastq", sample=SAMPLES),
