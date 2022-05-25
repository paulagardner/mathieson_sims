import msprime
import io
import demes 
import numpy as np
import gzip
import pandas as pd
import scipy
from scipy import stats
from scipy.stats import norm
from scipy.stats import multivariate_normal

demography = msprime.Demography()
graph = demes.load("no_ancestry.yaml")

mu = 1e-6
rho = 1e-3

def simulate(mu, rho, graph):
    ts = msprime.sim_ancestry(demography = msprime.Demography.from_demes(graph), recombination_rate=rho, sequence_length =1000, samples={"A": 50, "B": 50}, random_seed=1234, discrete_genome = False) #add random seed so it's replicable
    #not including migration info since I'm just doing two demes, so what they did (defining which demes are next to each other) not necessary here.. I think! at least for now

    mts = msprime.sim_mutations(ts, rate = mu, discrete_genome = False) #discrete_genome = false gives infinite sites, basically, so simplifies downstream bc you don't have to account for mult mutations at a site
    return mts


def read_vcf(vcf_path):
    with open(vcf_path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')] 
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
            'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    )#.rename(columns={'#CHROM': 'CHROM'}) #do not want to rename this! Plink2 is very unhappy if you take the #away 


def get_header(vcf_path): #rewrite to try to incorporate into read_vcf, returning both the df and header. ask skylar about it-- should be similar logic
    #to demography parser, sorting through the objects that are returned. This matters because the id definition is counting on one object being returned
    with open(vcf_path, 'r') as f:
        header = [l for l in f if l.startswith('##')]
    return header



print("simulating genotypes under demographic model")
ts = simulate(mu, rho, graph)


print("writing treefile for (potential) downstream analysis")
ts.dump("output.trees")


print("writing vcf")
n_dip_indv = int(ts.num_samples / 2)
indv_names = [f"tsk_{str(i)}indv" for i in range(n_dip_indv)]
vcf_path = "output_geno.vcf"

#convert to a pandas df
with open(vcf_path, "w") as vcf_file:  ts.write_vcf(vcf_file, individual_names=indv_names) 
df = read_vcf(vcf_path) 

#remove the header
header = get_header(vcf_path)
# print(header)

#make an id for all mutations
id  = ["rs" + str(i) for i,_ in enumerate(df.ID)]
#replace the pandas ID column with new IDs
df['ID'] = id
print(df)


print("gzipping file")
with gzip.open(vcf_path+".gz", 'wt') as testfile: 
    for l in header:
        testfile.write(l)

    df.to_csv(testfile, sep = '\t', index = False) 

#make the .txt file that will contain FID and IID    
d = 2 # replace hard-coded with argparser or taking the 'populations' value from the tskit table
#diploid sample size within each deme
ss = 50 #again, hardcoded
deme_id=[[i]*ss for i in range(0,d)] #https://github.com/Arslan-Zaidi/popstructure/blob/master/code/simulating_genotypes/grid/generate_genos_grid.py
#flatten
deme_id=[item for sublist in deme_id for item in sublist] #changes 2 arrays of, say, length 50 into one array of length 100 (for example, will vary depending on deme # and sample sizes))

fid=[f"tsk_{str(i)}indv" for i in range(0,(ss*d))]
iid=[f"tsk_{str(i)}indv" for i in range(0,(ss*d))]

print("writing pop file: FID, IID, deme ids for each individual")
popdf=pd.DataFrame({"FID":fid,
                  "IID":iid,
                  "POP":deme_id})

popdf.to_csv("test"+".pop",sep="\t",header=False,index=False)


#make phenotypes file
print("simulating phenotype")
np.random.seed(10)
random = norm.rvs(0, 1, (len(iid)))
mult_random = multivariate_normal.rvs(0, 1, (len(iid)))
phenotype_number = np.array([random, mult_random])

#random = scipy.stats.norm(0) #start with simulating a random gaussian
#popdf["phenotype"] = random
popdf["phenotype2"] = mult_random #simulating multivariate gaussian
popdf.to_csv("pop"+".txt",sep="\t",header=True,index=False,)







