#taking a LOT from the generate_genos_grid.py
# file, but stripping it back so I can just hardcode to learn w/ 
#the tau = -9/ no shared ancestry model

import msprime

demography = msprime.Demography()
# demography.add_population(name="A", initial_size=100)
# demography.add_population(name="B", initial_size=100)

# mu = 1e-8
#temporarily using a higher value to make sure sims are working correctly
mu = 1e-5
rho = 1e-3

import demes 
graph = demes.load("no_ancestry.yaml")

def simulate(mu, rho, graph):
    ts = msprime.sim_ancestry(demography = msprime.Demography.from_demes(graph), recombination_rate=rho, sequence_length =1000, samples={"A": 50, "B": 50}) 
    #not including migration info since I'm just doing two demes, so what they did (defining which demes are next to each other) not necessary here.. I think! at least for now

    mts = msprime.sim_mutations(ts, rate = mu)
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

# ts.dump("output.trees")
# new_ts = tskit.load("output.trees")
# print(new_ts) #this print is the same output as print(ts), which makes sense


print("writing to vcf")

#they do something interesting where it seems like they simulate chromosomes individually- since chromosome number is an argparse argument and vcf files contain only information for whichever chromosome is simulated




#following the file formats, writing etc made here: https://github.com/Arslan-Zaidi/popstructure/blob/master/code/simulating_genotypes/grid/generate_genos_grid.py

#had this working a while ago, now no data is being written under the header
#I think what happened for a bunch of reps was that things were technically *coded* correctly, but since I had a quite low mutation rate for the genome length parameter I was using, sometimes no mutations would happen and then I believe without that your vcf won't work
n_dip_indv = int(ts.num_samples / 2)
indv_names = [f"tsk_{str(i)}indv" for i in range(n_dip_indv)]
import gzip
with gzip.open("output_geno.vcf.gz", "wt") as vcf_file: #normal formatting would be with open("output.vcf", "w") as vcf_file:, here gzipping for when the analysis gets big. ALSO, kevin suggested "wb" as the flag, but tskit docs say "wt" and it didn't give errors like wb did: "TypeError: memoryview: a bytes-like object is required, not 'str'"
    ts.write_vcf(vcf_file, individual_names=indv_names) 
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
print(len(fid))
print(len(iid))
print(len(deme_id))
popdf=pd.DataFrame({"FID":fid,
                  "IID":iid,
                  "POP":deme_id})

print(popdf)

print("writing pop file- deme ids for each individual")
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





