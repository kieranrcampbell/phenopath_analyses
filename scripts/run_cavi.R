
# Run CAVI for CLVM and save result

# Format:
# Rscript clvm [input sceset] [output file] [pc to initialise] [optional covariate name]

library(clvm)
library(scater)

args <- commandArgs(trailingOnly = TRUE)

input_rds <- args[1]
output_file <- args[2]

pc_initialise <- as.numeric(args[3])

cov_name <- NULL

if(length(args) > 3) {
  cov_name <- args[4]
} else {
  cov_name <- 'x'
}

sce <- readRDS(input_rds)

y <- scale(t(exprs(sce)))
x <- cbind(pData(sce)[[ cov_name ]])


pcavi <- clvm(y, x)

saveRDS(pcavi, file = output_file)



