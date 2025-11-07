# ============================
#      Snakemake Workflow
# ============================

configfile: "config.yaml"

PROJECT = config["project"]
THREADS = config["threads"]
MARKERS = config["markers"]
ANNOTATIONS = config["annotations"]

rule all:
    input:
        "results/motif_analysis.done"

# ---- Step 1 ----
rule preprocess:
    output: "results/preprocess.done"
    shell:
        "Rscript preprocess.R --project {PROJECT} && touch {output}"

# ---- Step 2 ----
rule filter:
    input: "results/preprocess.done"
    output: "results/filter.done"
    shell:
        "Rscript filter.R --project {PROJECT} && touch {output}"

# ---- Step 3 ----
rule addUMAP:
    input: "results/filter.done"
    output: "results/addUMAP.done"
    shell:
        "Rscript addUMAP.R --project {PROJECT} && touch {output}"

# ---- Step 4 ----
rule plotMarkers:
    input: "results/addUMAP.done"
    output: "results/plotMarkers.done"
    shell:
        "Rscript plotMarkers.R --project {PROJECT} --markers {MARKERS} && touch {output}"

# ---- Step 5 ----
rule annotate:
    input: "results/plotMarkers.done"
    output: "results/annotate.done"
    shell:
        "Rscript annotate.R --project {PROJECT} --annotations {ANNOTATIONS} && touch {output}"

# ---- Step 6 ----
rule DGE:
    input: "results/annotate.done"
    output: "results/DGE.done"
    shell:
        "Rscript DGE.R --project {PROJECT} && touch {output}"

# ---- Step 7 ----
rule peak_analysis:
    input: "results/DGE.done"
    output: "results/peak_analysis.done"
    threads: THREADS
    shell:
        "Rscript peak_analysis.R --project {PROJECT} --threads {threads} && touch {output}"

# ---- Step 8 ----
rule motif_analysis:
    input: "results/peak_analysis.done"
    output: "results/motif_analysis.done"
    threads: THREADS
    shell:
        "Rscript motif_analysis.R --project {PROJECT} --threads {threads} && touch {output}"

