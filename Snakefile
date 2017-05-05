
"""
All analyses for phenotime

Notes:
    - Data downloaded from https://osf.io/gqrz9/
    - If pandoc gives an error, open Rstudio, type "sys.getenv("RSTUDIO_PANDOC")"
    and append the value to the $PATH environment variable
    e.g.
    export PATH=$PATH:/Applications/RStudio.app/Contents/MacOS/pandoc

Kieran R Campbell, February-May 2017
"""

R_opts = "--vanilla"

plot_names = ['limma_plot', 'goplot', 'gene_plot', 'lplot']
plot_files = [p + ".rds" for p in plot_names]
plot_files_coad = ["figs/coad/" + pf for pf in plot_files]
plot_files_brca = ["figs/brca/" + pf for pf in plot_files]
plot_files_brca.append("figs/brca/cross_plot.rds")
plot_files_brca.append("figs/brca/crossover_thesis.rds")
plot_files_coad.append("figs/coad/tregs.rds")


N = 200
samples = ["%02d" % (x+1) for x in range(N)]
strand = [1,2]

fasta_files = expand("data/simulations/fasta/sample_{samm}_{st}.fasta", samm = samples, st = strand)

fastq_str1 = expand("data/simulations/fastq/sample_{samp}_1.fastq", samp = samples)
fastq_str2 = expand("data/simulations/fastq/sample_{samp}_2.fastq", samp = samples)
kallisto_quants = expand("data/simulations/quant/sample_{sam}/abundance.tsv", sam = samples)

rule all:
    input:
        # "data/COAD/sce_coad.rds",
        # "data/COAD/sce_coad_clvm.rds",
        # "data/COAD/clvm_results.rds",
        # "data/OV/sce_ov.rds",
        # "data/OV/sce_ov_clvm.rds",
        # "data/OV/clvm_results.rds",
        "data/BRCA/sce_brca.rds",
        "data/BRCA/sce_brca_clvm.rds",
        "data/BRCA/sce_brca_gene_level.rds",
        "data/BRCA/clvm_results.rds",
        # "data/BRCA/expressed_genes.csv",
        # "data/shalek/sce_shalek.rds",
	# "data/shalek/sce_shalek_clvm.rds",
        # "data/shalek/clvm_results.rds",
        "data/BRCA/clvm_results_threecov.rds",
        "data/BRCA/clvm_er_pos_results.rds",
        "data/BRCA/clvm_tripleneg_results.rds",
        "analysis/BRCA/clvm_analysis_tripleneg.html",
        # plot_files_coad,
        plot_files_brca,
        "figs/supplementary_crossover.png",
        "figs/supplementary_crossover_2.png",
        "data/BRCA/brca_interactions.csv",
        # "data/COAD/coad_interactions.csv",
        # fasta_files,
        # fastq_str1, fastq_str2,
        # kallisto_quants,
        # "analysis/simulations/simulations.html",
        # "figs/shalek.png",
        # "figs/s_compare_monocle_dpt.png",
        # "figs/shalek.png",
        # "figs/s_compare_monocle_dpt.png",
        "figs/coad_figure.png",
        "figs/brca_figure.png",
        "figs/brca/crossover_thesis.png"


## ------ COAD -----

# rule construct_scecoad:
#     input:
#         "data/COAD/TCGA_COAD_tpm.tsv.gz",
#         "data/COAD/TCGA_COAD_counts.tsv.gz",
#         "data/TCGA_ID_MAP.csv",
#         "data/COAD/nm.4191-S3.xlsx",
#         "data/COAD/coadread_tcga_clinical_data.tsv"
#     output:
#         "data/COAD/sce_coad.rds"
#     shell:
#         "Rscript {R_opts} scripts/COAD/0_coad_to_sceset.R"
#
#
# rule prepare_coad:
#     """
#     This rule takes the transcript level quantification in data/COAD/sce_coad.Rdata,
#     collapses it to gene level, pulls out interesting genes and turns into
#     an SCESet with
#     """
#     input:
#         "data/COAD/sce_coad.rds"
#     output:
#         "data/COAD/sce_coad_clvm.rds"
#     shell:
#         "Rscript -e \"rmarkdown::render('analysis/COAD/prepare_for_clvm.Rmd')\""
#
# rule coad_clvm:
#     input:
#         "data/COAD/sce_coad_clvm.rds"
#     output:
#         "data/COAD/clvm_results.rds"
#     shell:
#         "Rscript scripts/run_cavi.R {input} {output} 1"

rule coad_analysis:
    input:
        "data/COAD/clvm_results.rds",
        "data/COAD/sce_coad_clvm.rds"
    output:
        plot_files_coad,
        "data/COAD/coad_interactions.csv"
    shell:
        "Rscript -e \"rmarkdown::render('analysis/COAD/clvm_analysis.Rmd')\""
#
# ## ---- OV ----
#
# rule construct_sceov:
#     input:
#         "data/OV/TCGA_OV_tpm.tsv.gz",
#         "data/OV/TCGA_OV_counts.tsv.gz"
#     output:
#         "data/OV/sce_ov.rds"
#     shell:
#         "Rscript {R_opts} scripts/OV/0_ov_to_sceset.R"
#
# rule prepare_ov:
#     input:
#         "data/OV/sce_ov.rds"
#     output:
#         "data/OV/sce_ov_clvm.rds"
#     shell:
#         "Rscript -e \"rmarkdown::render('analysis/OV/prepare_for_clvm.Rmd')\""
#
# rule ov_clvm:
#     input:
#         "data/OV/sce_ov_clvm.rds"
#     output:
#         "data/OV/clvm_results.rds"
#     shell:
#         "Rscript scripts/run_cavi.R {input} {output} 1"

# ## ---- BRCA ----

rule construct_scebrca:
    input:
        "data/BRCA/TCGA_BRCA_tpm.tsv.gz",
        "data/BRCA/TCGA_BRCA_counts.tsv.gz",
	"data/BRCA/brca_tcga_clinical_data.tsv",
	"data/BRCA/supplementary.xls"
    output:
        "data/BRCA/sce_brca.rds"
    shell:
        "Rscript {R_opts} scripts/BRCA/0_brca_to_sceset.R"

rule prepare_brca:
    input:
        "data/BRCA/sce_brca.rds"
    output:
        "data/BRCA/sce_brca_clvm.rds",
        "data/BRCA/sce_brca_clvm_er_pos.rds",
        "data/BRCA/sce_brca_gene_level.rds",
        "data/BRCA/brca_pca_plot.rds"
    shell:
        "Rscript -e \"rmarkdown::render('analysis/BRCA/prepare_for_clvm.Rmd')\""
#
#
# rule brca_clvm:
#     input:
#         "data/BRCA/sce_brca_clvm.rds"
#     output:
#         "data/BRCA/clvm_results.rds"
#     shell:
#         "Rscript scripts/run_cavi.R {input} {output} 1"

rule brca_clvm_tripleneg:
    input:
        "data/BRCA/sce_brca_clvm.rds"
    output:
        "data/BRCA/clvm_tripleneg_results.rds"
    shell:
        "Rscript scripts/run_cavi.R {input} {output} 1 is_triple_neg"

rule brca_clvm_er_pos:
    input:
        "data/BRCA/sce_brca_clvm_er_pos.rds"
    output:
        "data/BRCA/clvm_er_pos_results.rds"
    shell:
        "Rscript scripts/run_cavi.R {input} {output} 1"

#
# rule brca_expressed_genes:
#     input:
#         "data/BRCA/sce_brca_gene_level.rds"
#     output:
#         "data/BRCA/expressed_genes.csv"
#     shell:
#         "Rscript scripts/find_expressed_genes.R {input} {output}"
#
rule brca_analysis:
    input:
        "data/BRCA/clvm_results.rds",
        "data/BRCA/sce_brca_clvm.rds"
    output:
        plot_files_brca,
        "figs/supplementary_crossover.png",
        "figs/supplementary_crossover_2.png",
        "data/BRCA/brca_interactions.csv"
    shell:
        "Rscript -e \"rmarkdown::render('analysis/BRCA/clvm_analysis.Rmd')\""

rule brca_tripleneg_analysis:
    input:
        "data/BRCA/clvm_tripleneg_results.rds",
        "data/BRCA/sce_brca_clvm.rds"
    output:
        "analysis/BRCA/clvm_analysis_tripleneg.html",
        "figs/tripleneg/gene_plot.rds"
    shell:
        "Rscript -e \"rmarkdown::render('analysis/BRCA/clvm_analysis_tripleneg.Rmd')\""

rule brca_tripleneg_figure:
    input:
        "figs/tripleneg/gene_plot.rds",
    output:
        "figs/triple_neg_figure.png"
    shell:
        "Rscript scripts/BRCA/2_make_tripleneg_figure.R"

rule brca_threecov_clvm:
    input:
        "data/BRCA/sce_brca_clvm.rds"
    output:
        "data/BRCA/clvm_results_threecov.rds"
    shell:
        "Rscript scripts/BRCA/1_brca_three_covariate_clvm.R {input} {output}"




#
#
# ## ---- Shalek ----
#
# rule shalek_to_sceset:
#     input:
#         "data/shalek/GSE48968_allgenesTPM_GSM1189042_GSM1190902.txt.gz"
#     output:
#         "data/shalek/sce_shalek.rds"
#     shell:
#         "Rscript scripts/shalek/shalek_to_sceset.R"
#
# rule prepare_shalek:
#     input:
#         "data/shalek/sce_shalek.rds"
#     output:
#         "data/shalek/sce_shalek_clvm.rds"
#     shell:
#         "Rscript -e \"rmarkdown::render('analysis/shalek/prepare_for_clvm.Rmd')\""
#
# rule shalek_clvm:
#     input:
#         "data/shalek/sce_shalek_clvm.rds"
#     output:
#         "data/shalek/clvm_results.rds"
#     shell:
#         "Rscript scripts/run_cavi.R {input} {output} 1"



rule prepare_shalek:
    input:
        "data/shalek/sce_shalek.rds"
    output:
        "data/shalek/sce_shalek_clvm.rds"
    shell:
        "Rscript -e \"rmarkdown::render('analysis/shalek/prepare_for_clvm.Rmd')\""

# rule shalek_analysis:
#     input:
#         "data/shalek/clvm_results.rds"
#     output:
#         "figs/shalek.png",
#         "figs/s_compare_monocle_dpt.png"
#     shell:
#         "Rscript -e \"rmarkdown::render('analysis/shalek/shalek_clvm_analysis.Rmd')\""

## ---- Overall cancer figure

rule cancer_figure:
    input:
        plot_files_coad, plot_files_brca
    output:
        "figs/brca_figure.png",
        "figs/coad_figure.png",
        "figs/brca/crossover_thesis.png"
    shell:
        "Rscript scripts/cancer_figure.R"



## ---- Simulations


# fastq_files = ['data/simulations/fasta/sample_' + str(i + 1) + "_" + str(j + 1) + ".fasta" for i in range(N) for j in range(2)]


# rule simulate_fasta:
#     output:
#         "data/simulations/ref/chr22_small.fa",
#         "data/simulations/gene_pars.csv",
#         "data/simulations/pdata.csv",
#         fasta_files
#     shell:
#         "Rscript analysis/simulations/simulate_from_phenopath.R"


# rule fasta_to_fastq_str1:
#     input:
#         "data/simulations/fasta/sample_{sam}_1.fasta"
#     output:
#         "data/simulations/fastq/sample_{sam}_1.fastq"
#     shell:
#         "python3 scripts/fasta_to_fastq.py {input} {output}"

# rule fasta_to_fastq_str2:
#     input:
#         "data/simulations/fasta/sample_{sam}_2.fasta"
#     output:
#         "data/simulations/fastq/sample_{sam}_2.fastq"
#     shell:
#         "python3 scripts/fasta_to_fastq.py {input} {output}"

# rule kallisto_index:
#     input:
#         "data/simulations/ref/chr22_small.fa"
#     output:
#         "data/simulations/ref/chr22_small.idx"
#     shell:
#         "kallisto index -i {output} {input}"

# rule kallisto:
#     input:
#         str1="data/simulations/fastq/sample_{sam}_1.fastq",
#         str2="data/simulations/fastq/sample_{sam}_2.fastq",
#         ref="data/simulations/ref/chr22_small.idx"
#     output:
#         "data/simulations/quant/sample_{sam}/abundance.tsv"
#     shell:
#         "kallisto quant -i data/simulations/ref/chr22_small.idx -o data/simulations/quant/sample_{wildcards.sam} -b 100 {input.str1} {input.str2}"

rule simulations_rmd:
    input:
        kallisto_quants
    output:
        "analysis/simulations/simulations.html"
    shell:
        "Rscript -e \"rmarkdown::render('analysis/simulations/simulations.Rmd')\""
