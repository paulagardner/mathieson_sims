#!/bin/bash
conda init bash
conda activate mathieson_sims


python model.py



#convert to PLINK format. skipping splitting data step for now
plink2 --double-id --make-pgen --out convert_2_pgen --vcf output_geno.vcf.gz



# plink2 --vcf output_geno.vcf.gz --set-missing-var-ids @:# --make-pgen --double-id --out convert_2_pgen_test  #the second half shouldn't read as commented out.... but this is what makes the ID column issue you were having diappear
#https://www.cog-genomics.org/plink/2.0/data#set_all_var_ids



#START skip here for just running a GWAS without anything filtered
#just copying their filter rare variants step
# plink2 --mac 2 --max-mac 4 --out rare_variants_filter --pfile convert_2_pgen_test --write-snplist
plink2 --mac 2 --max-mac 4 --out rare_variants_filter --pfile convert_2_pgen --write-snplist allow-dups

#filter common variants
plink2 --maf 0.05 --out common_variants_filter --pfile convert_2_pgen --write-snplist allow-dups


#common variants PCA according to their format:
plink2 --pfile convert_2_pgen --extract common_variants_filter.snplist --pca --out common_PCA

#rare variants PCA according to their format:
plink2 --pfile convert_2_pgen --extract rare_variants_filter.snplist --pca --out rare_PCA 



###########END skip
#their step again: calculate allele frequencies:
#And now it seems they went back and calculated allele frequencies for everything by just making a snplist of EVERYTHING, no variants filtered, then just running freq on it: 
plink2 --pfile convert_2_pgen --mac 1 --write-snplist allow-dups --out convert_2_pgen.snps

#actually running freq: 
plink2 --pfile convert_2_pgen --extract convert_2_pgen.snps.snplist --freq --out convert_2_pgen.snps.frq


# chmod +x no_ge_test.r 
# Rscript no_ge_test.r #note- it's important that you're using the shebang version that specifies path, and not the setwd path that I use in the vscode GUI
#! /home/pdgardne/miniconda3/envs/mathieson_sims/bin/Rscript
##GWAS time 
#generated phenotypes using no_ge_test.r

#"script generates genetic effects and phenotypes with different patterns of stratification"


#not working yet
#bash gwas.sh convert_2_pgen_test pop.txt gwas.test #https://www.biostars.org/p/491009/ --ran into an error before wehere I had pop written as pop.csv, 
    #but the --pheno flag in plink needs it to be a space or comma delimited file

Rscript "PCA_plotting.r"

