
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
        "data/OV/clvm_results.rds",
	       "data/BRCA/sce_brca.rds",
        "data/BRCA/sce_brca_clvm.rds", "data/BRCA/sce_brca_gene_level.rds",
        "data/BRCA/clvm_results.rds",
        "data/BRCA/expressed_genes.csv",
        "data/shalek/sce_shalek.rds",
	"data/shalek/sce_shalek_clvm.rds",
	"data/shalek/clvm_results.rds"


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
        "Rscript scripts/run_cavi.R {input} {output} 1"

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
        "Rscript scripts/run_cavi.R {input} {output} 1"

## ---- BRCA ----

rule construct_scebrca:
    input:
        "data/BRCA/TCGA_BRCA_tpm.tsv.gz",
        "data/BRCA/TCGA_BRCA_counts.tsv.gz"
    output:
        "data/BRCA/sce_brca.rds"
    shell:
        "Rscript {R_opts} scripts/BRCA/0_brca_to_sceset.R"

rule prepare_brca:
    input:
        "data/BRCA/sce_brca.rds"
    output:
        "data/BRCA/sce_brca_clvm.rds",
        "data/BRCA/sce_brca_gene_level.rds"
    shell:
        "Rscript -e \"rmarkdown::render('analysis/BRCA/prepare_for_clvm.Rmd')\""


rule brca_clvm:
    input:
        "data/BRCA/sce_brca_clvm.rds"
    output:
        "data/BRCA/clvm_results.rds"
    shell:
        "Rscript scripts/run_cavi.R {input} {output} 1"

rule brca_expressed_genes:
    input:
        "data/BRCA/sce_brca_gene_level.rds"
    output:
        "data/BRCA/expressed_genes.csv"
    shell:
        "Rscript scripts/find_expressed_genes.R {input} {output}"


## ---- Shalek ----

rule shalek_to_sceset:
    input:
        "data/shalek/GSE48968_allgenesTPM_GSM1189042_GSM1190902.txt.gz"
    output:
        "data/shalek/sce_shalek.rds"
    shell:
        "Rscript scripts/shalek/shalek_to_sceset.R"

rule prepare_shalek:
    input:
        "data/shalek/sce_shalek.rds"
    output:
        "data/shalek/sce_shalek_clvm.rds"
    shell:
        "Rscript -e \"rmarkdown::render('analysis/shalek/prepare_for_clvm.Rmd')\""

rule shalek_clvm:
    input:
        "data/shalek/sce_shalek_clvm.rds"
    output:
        "data/shalek/clvm_results.rds"
    shell:
        "Rscript scripts/run_cavi.R {input} {output} 2"
