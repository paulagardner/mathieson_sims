# generate phenotypes for doing a gwas. make them non-heritable: effect is entirely due to environment
# following the 'noge' version of their phenotypes to make an environmental effect
# https://github.com/Arslan-Zaidi/popstructure/blob/master/code/gwas/grid/tau100/scripts/simphenotype/simphenotype_noge.R

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# hardcoding it in-- this means that the shell script will not need to take in args
popfile <- "test.pop"
outputfile <- "test_phenotype_simulation.txt"

install.packages("data.table") #based on suggestion here: https://stackoverflow.com/questions/52284128/no-fread-function
library("data.table")
pop <- fread(popfile)
colnames(pop) <- c("FID", "IID", "deme")

demes <- unique(pop$deme)
pop$random <- rnorm(nrow(pop), mean = 0, sd = 1)

print(pop)


