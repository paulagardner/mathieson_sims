#copied from https://github.com/Arslan-Zaidi/popstructure/blob/master/code/gwas/grid/tau-9/scripts/gwas/gwas.sh 
#args for gwas.sh:
#1. genotype file
#2. phenotype file
#3. output file prefix
#4. (if correction required) eigenvector file

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
  --memory 7000 ###########WAS GETTING AN ERROR ON LAPTOP ABOUT MEMORY, TEMP WORKAROUND
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
