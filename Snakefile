
"""
All analyses for phenotime

Notes:
    - Data downloaded from https://osf.io/gqrz9/

Kieran R Campbell, February 2017
"""

R_opts = "--vanilla"


rule all:
    input:
        "data/COAD/sce_coad.Rdata",
        "data/COAD/sce_coad_clvm.Rdata"

rule construct_scecoad:
    input:
        "data/COAD/TCGA_COAD_tpm.tsv.gz",
        "data/COAD/TCGA_COAD_counts.tsv.gz",
        "data/TCGA_ID_MAP.csv",
        "data/COAD/nm.4191-S3.xlsx",
        "data/COAD/coadread_tcga_clinical_data.tsv"
    output:
        "data/COAD/sce_coad.Rdata"
    shell:
        "Rscript {R_opts} scripts/COAD/0_coad_to_sceset.R"


rule prepare_coad:
    """
    This rule takes the transcript level quantification in data/COAD/sce_coad.Rdata,
    collapses it to gene level, pulls out interesting genes and turns into
    an SCESet with
    """
    input:
        "data/COAD/sce_coad.Rdata"
    output:
        "data/COAD/sce_coad_clvm.Rdata"
    shell:
        "Rscript -e \"rmarkdown::render('analysis/COAD/prepare_for_clvm.Rmd')\""
