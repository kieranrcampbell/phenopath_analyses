
# Run CAVI for CLVM and save result

library(clvm)
library(scater)

args <- commandArgs(trailingOnly = TRUE)

input_rds <- args[1]
output_file <- args[2]

cov_name <- NULL

if(length(args) > 2) {
  cov_name <- args[3]
} else {
  cov_name <- 'x'
}

sce <- readRDS(input_rds)

y <- scale(t(exprs(sce)))
x <- cbind(pData(sce)[[ cov_name ]])

pcavi <- clvm(y, x)

saveRDS(pcavi, file = output_file)



