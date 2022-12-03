#!/bin/bash
conda init bash
conda activate mathieson_sims


# python model.py

#running the model script should output a .pop and .vcf file.  

#convert the VCF into a format that plink can use (pgen). Skipping a splitting step the mathieson sims uses
#you will get 'no phenotype data present' even though you have it--- it's just in the pop.txt file, instead of the vcf! I *think* it's not necessary just for the PCA. 
plink2 --vcf genos.vcf --double-id --make-pgen --out convert_2_pgen 
#out: a .psam, .pvar, .pgen




#old step used to fill in IDs before I was assigning id #s to the 
# plink2 --vcf output_geno.vcf.gz --set-missing-var-ids @:# --make-pgen --double-id --out convert_2_pgen_test  #the second half shouldn't read as commented out.... but this is what makes the ID column issue you were having disappear
#https://www.cog-genomics.org/plink/2.0/data#set_all_var_ids



#START skipping this for just running a GWAS without anything filtered. However, filtering by rare and common variants means you can do different PCAs with them. 

#just copying their filter rare variants step
# plink2 --mac 2 --max-mac 4 --out rare_variants_filter --pfile convert_2_pgen_test --write-snplist
plink2 --mac 2 --max-mac 4 --out rare_variants_filter --pfile convert_2_pgen --write-snplist allow-dups

#filter common variants
plink2 --maf 0.05 --out common_variants_filter --pfile convert_2_pgen --write-snplist allow-dups


running PCAs
#common variants PCA according to their format:
plink2 --pfile convert_2_pgen --extract common_variants_filter.snplist --pca --out common_PCA

#rare variants PCA according to their format:
plink2 --pfile convert_2_pgen --extract rare_variants_filter.snplist --pca --out rare_PCA 

#overall PCA: 
# plink2 --pfile convert_2_pgen --extract convert_2_pgen.snps.snplist --pca --out overall_PCA

###########END skip




#their step again: calculate allele frequencies:
#And now it seems they went back and calculated allele frequencies for everything by just making a snplist of EVERYTHING, no variants filtered, then just running freq on it: 

#get snps (if not filtering, this should just be the list of mutations in the vcf). 
#--mac 1 was used in the mathieson_sims doc
plink2 --pfile convert_2_pgen --mac 1 --write-snplist allow-dups --out convert_2_pgen.snps
#out: snps.snplist 
# Running into an issue where plink thinks my mutations are all fixed, so I can't trust the betas I get. (my phenotypes are variable, so I think it's associating fixed mutations with indivduals that have variable phenotypes-- which would be like simulating something entirely environmental, since plink thinks everyone is genetically identical). but it doesn't seem to be happening all the time?

#actually running freq: obtain allele frequencies from the list of snps
plink2 --pfile convert_2_pgen --extract convert_2_pgen.snps.snplist --freq --out convert_2_pgen.snps.freq

# plink2 --pfile convert_2_pgen --freq --out convert_2_pgen.snps.freq


# chmod +x no_ge_test.r 
# Rscript no_ge_test.r #note- to use r, it's important that you're using the shebang version that specifies path, and not the setwd path that I use in the vscode GUI
#! /home/pdgardne/miniconda3/envs/mathieson_sims/bin/Rscript
##GWAS time 
#generated phenotypes using no_ge_test.r

#"script generates genetic effects and phenotypes with different patterns of stratification"


#not working yet
# bash gwas.sh convert_2_pgen_test pop.txt gwas.test #https://www.biostars.org/p/491009/ --ran into an error before wehere I had pop written as pop.csv, 
    #but the --pheno flag in plink needs it to be a space or comma delimited file

#Rscript "PCA_plotting.r"


#GWAS time. You will run a glm on each phenotype, and will get a file containing what kind of test (with a linear model, an additive one)
#"beta: Regression coefficient (for A1 if additive test)."




#uncorrected
plink2 --glm  --pfile convert_2_pgen --read-freq convert_2_pgen.snps.freq.afreq  --pheno pop.txt --out gwas  #mathieson et al has --hide-covar, if you remove it





#common variants
plink2 --glm   --pfile convert_2_pgen --read-freq convert_2_pgen.snps.freq.afreq  --pheno pop.txt --covar common_PCA.eigenvec  --out gwas_common_correction
#the output file from this gets pretty tricky

#rare variants
plink2 --glm  --pfile convert_2_pgen --read-freq convert_2_pgen.snps.freq.afreq  --pheno pop.txt --covar rare_PCA.eigenvec  --out gwas_rare_correction



#you could include a covariates file if you'd like-- may try to do that to have the deme id as a covariate. The mathieson paper does not do this, but they do include correction w/ eigenvectors file. (you can get that from doing PCAs)
#tried using --covar pop.txt, and it loads "3 covariates" which I'm assuming end up being deme_id, phenotype1 and phenotype2. I want deme_id to be a covariate, but probably not the phenotypes?

# plink2 --glm --pfile convert_2_pgen --read-freq convert_2_pgen.snps.freq.afreq  --pheno pop.txt --covar pop.txt --covar-col-nums 3 --out gwas 
#specifying covar and taking only the deme id column out (it's the 3rd column, counting from 1, b/c not python) to try to have this as a covariate, but the .glm.linear file you get from that is a bit odd. 

#from plink2 --help --glm "To add covariates which are not constant across all variants, add the'local-covar=', 'local-pvar=', and 'local-psam=' modifiers (note: to the --glm flag), and use full filenames for each.""
# plink2 --glm local-covar= covar.txt --pfile convert_2_pgen --read-freq convert_2_pgen.snps.freq.afreq  --pheno pop.txt --out gwas

#calculating a polygenic score:
#doing it with all variants to start with: 
# plink2 --score convert_2_pgen.snps.freq.afreq list-variants --pfile convert_2_pgen --pheno pop.txt --out testing --score-col-nums 3-12
# plink2 --score gwas.phenotype1.glm.linear list-variants --pfile convert_2_pgen --pheno pop.txt --out testing --score-col-nums 3-12
# plink2 --score genos.vcf list-variants header --pfile convert_2_pgen --pheno pop.txt --out testing --score-col-nums 3-12

#how to get PGS from your GLM: you specify the gwas you want to get PGS for, and three columns: variant ID, allele code (kevin suggested I use the minor allele), and the column that actually contains the numbers to score from. For this, it'll be the betas.   
# plink2 --score gwas.phenotype1.glm.linear 3 5 9 --pfile convert_2_pgen --pheno pop.txt --out PGS

#make effect size vs. p-value plots
#phenotype 1: 

# plink2 --score betas.file list-variants --pfile convert_2_pgen --pheno pop.txt --out PGS