#taking a LOT from the generate_genos_grid.py
# file, but stripping it back so I can just hardcode to learn w/ 
#the tau = -9/ no shared ancestry model

import msprime

demography = msprime.Demography()
# demography.add_population(name="A", initial_size=100)
# demography.add_population(name="B", initial_size=100)

# mu = 1e-8
#temporarily using a higher value to make sure sims are working correctly
mu = 1e-6
rho = 1e-3

import demes 
graph = demes.load("no_ancestry.yaml")

def simulate(mu, rho, graph):
    ts = msprime.sim_ancestry(demography = msprime.Demography.from_demes(graph), recombination_rate=rho, sequence_length =1000, samples={"A": 50, "B": 50}, random_seed=1234, discrete_genome = False) #add random seed so it's replicable
    #not including migration info since I'm just doing two demes, so what they did (defining which demes are next to each other) not necessary here.. I think! at least for now

    mts = msprime.sim_mutations(ts, rate = mu, discrete_genome = False) #discrete_genome = false gives infinite sites, basically, so simplifies downstream bc you don't have to account for mult mutations at a site
    return mts




print("simulating genotypes under demographic model")

ts = simulate(mu, rho, graph)

print("printing tree sequence")
print(ts)
print("")

import json

print("printing population table")
print(ts.tables.populations)
print("")

print("printing individuals table")
print(ts.tables.individuals)
print("")

print("printing nodes table")
print(ts.tables.nodes)

# population = ts.population()
# print(population)

#provenance = json.loads(ts.provenance(0).record) #for if you wanted to dump to treefile then load again
# provenance = ts.provenance(0).record
# print("PROVENANCE IS")
# print(provenance)

# print("writing treefile")
# import tskit

print("writing treefile")
ts.dump("output.trees")
# new_ts = tskit.load("output.trees")
# print(new_ts) #this print is the same output as print(ts), which makes sense


print("writing to vcf")

#they do omething interesting where it seems like they simulate chromosomes individually- since chromosome number is an argparse argument and vcf files contain only information for whichever chromosome is simulated


    



import numpy #I import it below, but I want to have everything centralized for making the requirements file)
#following the file formats, writing etc made here: https://github.com/Arslan-Zaidi/popstructure/blob/master/code/simulating_genotypes/grid/generate_genos_grid.py

#had this working a while ago, now no data is being written under the header
#I think what happened for a bunch of reps was that things were technically *coded* correctly, but since I had a quite low mutation rate for the genome length parameter I was using, sometimes no mutations would happen and then I believe without that your vcf won't work
n_dip_indv = int(ts.num_samples / 2)
indv_names = [f"tsk_{str(i)}indv" for i in range(n_dip_indv)]
import gzip
import pandas as pd
vcf_path = "output_geno.vcf"
with open(vcf_path, "w") as vcf_file: #normal formatting would be with open("output.vcf", "w") as vcf_file:, here gzipping for when the analysis gets big. ALSO, kevin suggested "wb" as the flag, but tskit docs say "wt" and it didn't give errors like wb did: "TypeError: memoryview: a bytes-like object is required, not 'str'"
    ts.write_vcf(vcf_file, individual_names=indv_names) 

    import io
    import os
    import pandas as pd


    def read_vcf(vcf_path):
        with open(vcf_path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')] #removes the header columns. you will need to figure out how to 
        return pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t'
        ).rename(columns={'#CHROM': 'CHROM'})


    def get_header(vcf_path): #rewrite to try to incorporate into read_vcf, returning both the df and header. ask skylar about it-- should be similar logic
        #to demography parser, sorting through the objects that are returned. This matters because the id definition is counting on one object being returned
        with open(vcf_path, 'r') as f:
            header = [l for l in f if l.startswith('##')]
        return header



    df = read_vcf(vcf_path) 
    header = get_header(vcf_path)
    print(header)
    # print(df)

    id  = ["rs" + str(i) for i,_ in enumerate(df.ID)]
    # print(id)
    # id_df = pd.DataFrame(id)

    # print(len(df))
    # print(len(id_df))

    df['ID'] = id
    print(df)


    with open('a.vcf', 'w') as testfile: #writes to end of file, but logic works
        for l in header:
            testfile.write(l)
    
        df.to_csv(testfile, sep = '\t')




    # with open(vcf_path, 'w') as vcf_file: #writes to end of file, but logic works
    #     for l in header:
    #         vcf_file.write(l)

    # df.to_csv(vcf_file, sep = '\t')



    # with open(vcf_path, "w") as vcf_file:
    #     df.to_csv(vcf_path, sep = '\t')
    # print(vcf_file)



    #tskit docs specified some weird, extraneous-looking stuff because of these docs on tskit for write_vcf having you do this so plink doesn't freak out 
    #but took out ploidy = 2 b/c ValueError: Cannot specify ploidy when individuals present






d = 2 # replace hard-coded with argparser or taking the 'populations' value from the tskit table
#args.deme originally

#diploid sample size within each deme
ss = 50 #again, hardcoded. args.ss originally
deme_id=[[i]*ss for i in range(0,d)]
#flatten
deme_id=[item for sublist in deme_id for item in sublist] #changes 2 arrays of length 50 into one array of length 100 (for example, will vary depending on deme # and sample sizes))
print(deme_id)
#fid and iid
# fid=["msp_"+str(i) for i in range(0,(ss*d))] #mismatch between the above and this was causing my ID's to not work, so I would get a plink error: NO entries in pop.txt correspond to loaded sample IDs.
# iid=["msp_"+str(i) for i in range(0,(ss*d))]
fid=[f"tsk_{str(i)}indv" for i in range(0,(ss*d))]
iid=[f"tsk_{str(i)}indv" for i in range(0,(ss*d))]
# print(fid)
# print(iid)

import pandas as pd

print("checking why valueerror: all arrays must be of the same length")
# print(len(fid))
# print(len(iid))
# print(len(deme_id))

print("writing pop file: FID, IID, deme ids for each individual")
popdf=pd.DataFrame({"FID":fid,
                  "IID":iid,
                  "POP":deme_id})

#print(popdf)

popdf.to_csv("test"+".pop",sep="\t",header=False,index=False)
# popdf.to_csv("test"+".pop",sep="\t")

#
# 
# 
# print("generate phenotypes")

# Goal: generate phenotypes for the GWAS, think I need to 
# make a different file w/ different format for this AND have variance-covariance

#the zaidi paper just makes heritable phenotypes with a .txt file. try that:
#phenotype = scipy.stats.multivariate_normal
#with gzip.open("output_pheno.txt", "wt") as phenotype_file:



print("simulating phenotype")
# print(popdf)

import scipy
from scipy import stats
from scipy.stats import norm
from scipy.stats import multivariate_normal
import numpy as np

np.random.seed(10)
random = norm.rvs(0, 1, (len(iid)))
mult_random = multivariate_normal.rvs(0, 1, (len(iid)))
phenotype_number = np.array([random, mult_random])

#random = scipy.stats.norm(0) #start with simulating a random gaussian
#popdf["phenotype"] = random
popdf["phenotype2"] = mult_random

# print(popdf)
popdf.to_csv("pop"+".txt",sep="\t",header=True,index=False,)







