#copied from https://github.com/Arslan-Zaidi/popstructure/blob/master/code/gwas/grid/tau-9/scripts/gwas/gwas.sh 
#args for gwas.sh:
#1. genotype file
#2. phenotype file
#3. output file prefix
#4. (if correction required) eigenvector file


#had to manually cut the ID column from the frq.afreq file to try to correct for "Error: --read-freq variant ID '.' appears multiple times in main dataset", due to tskit's write_vcf not giving position names, they're all " . "
#cut --complement --field=2 convert_2_pgen_test.snps.frq.afreq  #from https://unix.stackexchange.com/questions/222121/how-to-remove-a-column-or-multiple-columns-from-file-using-shell-command. You'd use a flag for space delimited files
    
#saving it: 
# echo "$(awk '{cut --complement -f2 convert_2_pgen_test.snps.frq.afreq}' convert_2_pgen_test.snps.frq.afreq)" > convert_2_pgen_test.snps.frq.afreq #https://stackoverflow.com/questions/14660079/how-to-save-the-output-of-this-awk-command-to-file
# plink2 --pfile convert_2_pgen_test --read-freq convert_2_pgen_test.snps.frq.afreq2 --glm hide-covar --pheno pop.txt --out gwas


# cat convert_2_pgen_test.snps.frq.afreq | cut --complement -f2 > list.txt #trying this #this worked (thank Skylar!) but the gwas is not happy with not having an ID column

# plink2 --pfile convert_2_pgen_test --read-freq convert_2_pgen_test.snps.frq.afreq --glm hide-covar --pheno test.txt --out gwas
plink2 --pfile convert_2_pgen_test --read-freq convert_2_pgen_test.snps.frq.afreq --glm hide-covar --pheno pop.txt --out gwas #"Error: --read-freq variant ID '.' appears multiple times in main dataset."
#plink2 --pfile convert_2_pgen_test --read-freq list.txt --glm hide-covar --pheno pop.txt --out gwas #"Error: Missing column(s) in --read-freq file (ID, REF, ALT[1] usually required).""
#didn't work, says I need ID column. Maybe try replacing with 'NA'? 



#exit w/ incorrect arguments if there are less than 3 provided (eigenvector file is optional)
if [ "$#" -lt 3  ]; then
    echo "usage: bash gwas.sh <path to pgen genotype file - prefix> <path to phenotype file> <path to output file> <eigenvector file - optional>"
    exit 1
fi


geno_file_prefix=${1}
phenotype_file=${2}
output_file_prefix=${3}
eigenvec_file=${4}

if [ "$#" == 3 ];
then
  plink2 --pfile ${geno_file_prefix} \
  --read-freq ${geno_file_prefix}.frq.afreq \ 
  --glm hide-covar \
  --pheno ${phenotype_file} \
  --out ${output_file_prefix}
  #--memory 7000 ###########WAS GETTING AN ERROR ON LAPTOP ABOUT MEMORY, TEMP WORKAROUND
  #I'm having issues where my frq.afreq and pfile have different prefixes (convert_2_pgen_test vs convert_2_pgen_test.snps, so tried using: plink2 --pfile convert_2_pgen_test --read-freq convert_2_pgen_test.snps.frq.afreq --glm hide-covar --pheno pop.txt --out gwas )

fi

if [ "$#" == 4 ];
then
  plink2 --pfile ${geno_file_prefix} \
  --read-freq ${geno_file_prefix}.frq.afreq \
  --glm hide-covar \
  --pheno ${phenotype_file} \
  --covar ${eigenvec_file} \
  --covar-col-nums 3-102 \
  --out ${output_file_prefix}
fi
