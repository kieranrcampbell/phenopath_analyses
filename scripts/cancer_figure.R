library(ggplot2)
library(cowplot)
library(splines)


plot_paper_grid <- function(limma_plot, goplot, gene_plot, lplot,
                            A = "A", B = "B", C = "C", D = "D") {
  lsize <- 11
  bottom_left <- plot_grid(plot_grid(NULL, limma_plot, NULL, nrow = 1, rel_widths = c(1,6, 1)), 
                           goplot,
                           ncol = 1, labels = c(B, C),
                           label_size = lsize)
  bottom_grid <- plot_grid(bottom_left, gene_plot, nrow = 1, labels = c("", D),
                           label_size = lsize)
  
  p <- plot_grid(lplot, bottom_grid, ncol = 1, 
                 rel_heights = c(2,1.2), labels = c(A,""),
                 label_size = lsize)
}

coad_limma_plot <- readRDS("figs/coad/limma_plot.rds")
coad_goplot <- readRDS("figs/coad/goplot.rds")
coad_gene_plot <- readRDS("figs/coad/gene_plot.rds")
coad_lplot <- readRDS("figs/coad/lplot.rds")

brca_limma_plot <- readRDS("figs/brca/limma_plot.rds")
brca_goplot <- readRDS("figs/brca/goplot.rds")
brca_gene_plot <- readRDS("figs/brca/gene_plot.rds")
brca_lplot <- readRDS("figs/brca/lplot.rds")

coad_plot <- plot_paper_grid(coad_limma_plot, coad_goplot, coad_gene_plot, coad_lplot)
brca_plot <- plot_paper_grid(brca_limma_plot, brca_goplot, brca_gene_plot, brca_lplot,
                             "E", "F", "G", "H")

coad_plot <- coad_plot + ggtitle("Colorectal cancer") + panel_border()
brca_plot <- brca_plot + ggtitle("Breast cancer") + panel_border()

cancer_plot <- plot_grid(coad_plot, brca_plot, nrow = 1)

ggsave("figs/cancer_figure.png", cancer_plot, width = 16, height = 10)
