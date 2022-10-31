#!/bin/bash
#following the zaidi github, again: https://github.com/Arslan-Zaidi/popstructure/wiki/PRS_calculation

#however skipping filtering SNPs at first- just going to work on getting PGS working with all variants. (they filter to remove SNPs with large P values, filter by genome window)

geno_file_prefix=${1}
phenotype_file=${2}
output_file_prefix=${3}
eigenvec_file=${4}

# plink2 --pfile convert_2_pgen \
# --read-freq convert_2_pgen.snps.frq.afreq \
# --score gwas.phenotype1.glm.linear 3 4 9 'header' \ 
# --out phenotype1 

######the above doesn't work, but this does! 
plink2 --pfile convert_2_pgen --read-freq convert_2_pgen.snps.frq.afreq --score gwas.phenotype1.glm.linear 3 4 9 'header' --pheno pop.txt --out phenotype1  #note that when trying to use --score-col-nums 3,4,9 as a separate flag, it did not work. Also, if you don't tell plink2 to disregard the first line of the glm.linear files with 'header', it still works but you will get an error message that 1 entry was skipped. I think functionally it ends up the same but better to not be ignoring error messages...
# plink2 --pfile convert_2_pgen --read-freq convert_2_pgen.snps.frq.afreq --score gwas.phenotype1.glm.linear 3 4 9 'header' --pheno pop.txt --out phenotype1 #zaidi & mathieson does not include phenotype data, so I'm putting it in with pop.txt, but 

plink2 --pfile convert_2_pgen --read-freq convert_2_pgen.snps.frq.afreq --score gwas.phenotype2.glm.linear 3 4 9 'header' --pheno pop.txt --out phenotype2 

#zaidi et al uses flags to reduce the number of columns: plink2 --pfile convert_2_pgen --read-freq convert_2_pgen.snps.frq.afreq --score gwas.phenotype2.glm.linear 3 4 9 'header' cols=dosagesum,scoresums --out phenotype2