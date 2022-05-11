#!/bin/bash
conda init bash
conda activate mathieson_sims

#simulate the genotypes (later, plan to have this script also generate the phenotypes)
python model.py


#have this also in the gwas.sh file
#had to manually cut the ID column from the frq.afreq file to try to correct for "Error: --read-freq variant ID '.' appears multiple times in main dataset", due to tskit's write_vcf not giving position names, they're all " . "
#first, unzip
gunzip output_geno.vcf.gz 
#now do it
cut --complement -f3 output_geno.vcf #from https://unix.stackexchange.com/questions/222121/how-to-remove-a-column-or-multiple-columns-from-file-using-shell-command. You'd use a flag for space delimited files


#convert to PLINK format. skipping splitting data step for now
plink2 --double-id --make-pgen --out convert_2_pgen_test --vcf output_geno.vcf.gz



#START skip here for just running a GWAS without anything filtered
#just copying their filter rare variants step
# plink2 --mac 2 --max-mac 4 --out rare_variants_filter --pfile convert_2_pgen_test --write-snplist
plink2 --mac 2 --max-mac 4 --out rare_variants_filter --pfile convert_2_pgen_test --write-snplist allow-dups

#filter common variants
plink2 --maf 0.05 --out common_variants_filter --pfile convert_2_pgen_test --write-snplist allow-dups


#common variants PCA according to their format:
plink2 --pfile convert_2_pgen_test --extract common_variants_filter.snplist --pca --out common_PCA

#rare variants PCA according to their format:
plink2 --pfile convert_2_pgen_test --extract rare_variants_filter.snplist --pca --out rare_PCA 



###########END skip
#their step again: calculate allele frequencies:
#And now it seems they went back and calculated allele frequencies for everything by just making a snplist of EVERYTHING, no variants filtered, then just running freq on it: 
plink2 --pfile convert_2_pgen_test --mac 1 --write-snplist allow-dups --out convert_2_pgen_test.snps

#actually running freq: 
plink2 --pfile convert_2_pgen_test --extract convert_2_pgen_test.snps.snplist --freq --out convert_2_pgen_test.snps.frq


##GWAS time 
#generated phenotypes using no_ge_test.r

#"script generates genetic effects and phenotypes with different patterns of stratification"
bash gwas.sh convert_2_pgen_test pop.txt gwas.test #https://www.biostars.org/p/491009/ --ran into an error before wehere I had pop written as pop.csv, 
    #but the --pheno flag in plink needs it to be a space or comma delimited file