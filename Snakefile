
"""
All analyses for phenotime

Notes:
    - Data downloaded from https://osf.io/gqrz9/

Kieran R Campbell, February 2017
"""

R_opts = "--vanilla"


rule all:
    input:
        "data/COAD/sce_coad.rds",
        "data/COAD/sce_coad_clvm.rds",
        "data/COAD/clvm_results.rds",
        "data/OV/sce_ov.rds",
        "data/OV/sce_ov_clvm.rds",
        "data/OV/clvm_results.rds"


## ------ COAD -----

rule construct_scecoad:
    input:
        "data/COAD/TCGA_COAD_tpm.tsv.gz",
        "data/COAD/TCGA_COAD_counts.tsv.gz",
        "data/TCGA_ID_MAP.csv",
        "data/COAD/nm.4191-S3.xlsx",
        "data/COAD/coadread_tcga_clinical_data.tsv"
    output:
        "data/COAD/sce_coad.rds"
    shell:
        "Rscript {R_opts} scripts/COAD/0_coad_to_sceset.R"


rule prepare_coad:
    """
    This rule takes the transcript level quantification in data/COAD/sce_coad.Rdata,
    collapses it to gene level, pulls out interesting genes and turns into
    an SCESet with
    """
    input:
        "data/COAD/sce_coad.rds"
    output:
        "data/COAD/sce_coad_clvm.rds"
    shell:
        "Rscript -e \"rmarkdown::render('analysis/COAD/prepare_for_clvm.Rmd')\""

rule coad_clvm:
    input:
        "data/COAD/sce_coad_clvm.rds"
    output:
        "data/COAD/clvm_results.rds"
    shell:
        "Rscript scripts/run_cavi.R {input} {output}"

## ---- OV ----

rule construct_sceov:
    input:
        "data/OV/TCGA_OV_tpm.tsv.gz",
        "data/OV/TCGA_OV_counts.tsv.gz"
    output:
        "data/OV/sce_ov.rds"
    shell:
        "Rscript {R_opts} scripts/OV/0_ov_to_sceset.R"

rule prepare_ov:
    input:
        "data/OV/sce_ov.rds"
    output:
        "data/OV/sce_ov_clvm.rds"
    shell:
        "Rscript -e \"rmarkdown::render('analysis/OV/prepare_for_clvm.Rmd')\""

rule ov_clvm:
    input:
        "data/OV/sce_ov_clvm.rds"
    output:
        "data/OV/clvm_results.rds"
    shell:
        "Rscript scripts/run_cavi.R {input} {output}"
