rule template_rule:
    input:
        lambda wc: SAMPLES[wc.sample]["R1"],
        lambda wc: SAMPLES[wc.sample]["R2"],
    output:
        OUT + "/{sample}_combined.fastq",
    log:
        OUT + "/log/{sample}_template_rule.log"
    message:
        "Merging {input}."
    resources:
        mem_gb=config["mem_gb"]["template_rule"],
    params: script = "workflow/scripts/script.py"
    threads: config["threads"]["template_rule"]
    shell: """
    python {params.script} {input} > {output}
    """
