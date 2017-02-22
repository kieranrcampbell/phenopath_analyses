
# Run CAVI for CLVM and save result

library(clvm)

args <- commandArgs(trailingOnly = TRUE)

input_rds <- args[1]
output_file <- args[2]

sce <- readRDS(input_rds)

y <- scale(t(exprs(sce)))
x <- cbind(sce$x)

pcavi <- clvm(y, x)

saveRDS(pcavi, file = output_file)



