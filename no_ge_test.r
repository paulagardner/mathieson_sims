#!/usr/bin/env Rscript
#also, in terminal: chmod +x no_ge_test.r

# generate phenotypes for doing a gwas. make them non-heritable: effect is entirely due to environment
# following the 'noge' version of their phenotypes to make an environmental effect
# https://github.com/Arslan-Zaidi/popstructure/blob/master/code/gwas/grid/tau100/scripts/simphenotype/simphenotype_noge.R

#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #commenting this out for use on the command line
# hardcoding it in-- this means that the shell script will not need to take in args
popfile <- "test.pop"
# outputfile <- "test_phenotype_simulation.txt"

# install.packages("tidyverse") #based on suggestion here: https://stackoverflow.com/questions/52284128/no-fread-function
# install.packages("data.table", INSTALL_opts = c('--no-lock'))
# install.packages("data.table") #https://github.com/Rdatatable/data.table/wiki/

# install.packages("data.table")
# library("data.table")
# pop <- fread(popfile)
# colnames(pop) <- c("FID", "IID", "deme")
# to get this to work, had to first install pkg-config: https://zoomadmin.com/HowToInstall/UbuntuPackage/pkg-config
# then, I was having issues where data.table was not detecting OpenMP (as stated here) https://github.com/Rdatatable/data.table/wiki/Installation
# getting stuck on the openMP step, but this link says what I was thinking-- think I need to redo the r install https://stackoverflow.com/questions/57825428/rstudio-server-error-bin-sh-x86-64-conda-cos6-linux-gnu-cc-command-not-found
# https://github.com/RcppCore/Rcpp/issues/770  ########THIS is where the exact issue I seemed to be having was happening- having r through a conda install
# installing r with conda install gxx_linux-64 unfortunately didn't work for me

# putting a pin in this because I think I can work around it with read.csv, but I will almost certainly need to be able to use the tidyverse at some point
# it seems the advantage of fread instaed of read.csv is a performance one, maybe some options


# the R file on github that I'm following makes a couple of different phenotypes with no genetic basis
# and puts them in different columns. Just doing one, completely random one for now
pop <- read.csv(popfile, sep = "\t", )
pop
colnames(pop) <- c("FID", "IID", "deme")
pop
demes <- unique(pop$deme)
pop$random <- rnorm(nrow(pop), mean = 0, sd = 1)

# pop=pop[,c("FID","IID","random")] #they remove the deme ids but I think I want them to be able to do the GWAS by deme????

print(pop)


# write.csv(pop, file = 'pop.csv')
# https://www.biostars.org/p/491009/
# may instead need to write as a .txt file to match the downstream

pop <- pop[, c("FID", "IID", "random")] # same as above-confused as to how I'll be able to do things w/o the deme ID but phenotype has to be the third column, apparently, so...
pop
write.table(pop, "pop.txt", append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
