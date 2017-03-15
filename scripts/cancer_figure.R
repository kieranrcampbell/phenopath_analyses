library(ggplot2)
library(cowplot)

coad_lplot <- readRDS("figs/coad/lplot.rds")

coad_plot <- readRDS("figs/coad.rds")
brca_plot <- readRDS("figs/brca.rds")


cancer_plot <- plot_grid(coad_plot, brca_plot, nrow = 1)

ggsave("figs/cancer_figure.png", cancer_plot, width = 16, height = 10)
