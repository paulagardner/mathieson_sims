from re import A
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
from numpy.random import default_rng
import tskit

demography = msprime.Demography()
graph = demes.load("no_ancestry.yaml")

mu = ((1e-7)*2)
rho = 1e-7

deme_size = 10
def simulate(mu, rho, graph):
    ts = msprime.sim_ancestry(demography = msprime.Demography.from_demes(graph), recombination_rate=rho, sequence_length = 500, samples={"A": deme_size, "B": deme_size}, random_seed=1234, discrete_genome = False) #add random seed so it's replicable
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
#print(df)


print("gzipping file")
with gzip.open(vcf_path+".gz", 'wt') as testfile: 
    for l in header:
        testfile.write(l)

    df.to_csv(testfile, sep = '\t', index = False) 

#make the .txt file that will contain FID and IID    
demes = 2 # replace hard-coded with argparser or taking the 'populations' value from the tskit table. in this case, you just want the demes of the samples present at the end of the sim. (I think). SO regardless of pop_split.yaml or no_ancestry.yaml, demes A and B split off.
#diploid sample size within each deme
#should be a way to get indices from defining ts above


# ss = 10 #NOTE take this out because I've changed it so it's a little less hardcoding
deme_id=[[i]*deme_size for i in range(0,demes)] #https://github.com/Arslan-Zaidi/popstructure/blob/master/code/simulating_genotypes/grid/generate_genos_grid.py
#flatten
deme_id=[item for sublist in deme_id for item in sublist] #changes 2 arrays of, say, length 50 into one array of length 100 (for example, will vary depending on deme # and sample sizes))

fid=[f"tsk_{str(i)}indv" for i in range(0,(deme_size*demes))]
iid=[f"tsk_{str(i)}indv" for i in range(0,(deme_size*demes))] #number of individuals in the sample

print("writing pop file: FID, IID, deme ids for each individual")
popdf=pd.DataFrame({"FID":fid,
                  "IID":iid,
                  "POP":deme_id})

popdf.to_csv("test"+".pop",sep="\t",header=False,index=False)


#make phenotypes file
# print("simulating phenotypes that are purely environmental")
# np.random.seed(10)
random = norm.rvs(0, 1, (len(iid)))
mult_random = multivariate_normal.rvs(0, 1, (len(iid)))
phenotype_ID = np.array([random, mult_random]) #for the numpy array to work as expected, you'll need to make sure this line is accurate to how many phenotypes you want
num_phenotypes = len(phenotype_ID)
#random = scipy.stats.norm(0) #start with simulating a random gaussian
#popdf["phenotype"] = random
popdf["phenotype2"] = mult_random #simulating multivariate gaussian
popdf.to_csv("pop"+".txt",sep="\t",header=True,index=False,)


phenos = 2
#make the variance-covariance matrix of phenotypic effects: 
var_covar = np.identity(phenos)
print(var_covar)

#make the means matrix for the multivariate generator
means = np.zeros(phenos)

#create an empty numpy array to put the effect sizes in. It will return a table with #phenotypes column and #individuals rows.
effects = np.zeros(n_dip_indv*num_phenotypes).reshape(n_dip_indv,num_phenotypes)
print("empty array for effect sizes:", effects)

 

#make effect size distributtions:
print("simulating phenotypes that are purely environmental")
mult_random = multivariate_normal.rvs(0, 1, (len(iid))) #assigns a random effect size to every individual in one big array. NOTE that this is the same as above, remove this redundancy.
# print(mult_random)

print("number of mutations:", len(ts.tables.mutations))
#make a table with mutations indexed from 0 to m-1
mut_index = [i for i in range(len(ts.tables.mutations))] #goes from 0 to #mutations -1 (again, a consequence of python counting from 0)
print("mutation index array: ", mut_index)
# index = mut_index.index
# print(index)

treesequence = tskit.load("output.trees")
print("sample list:", treesequence.samples()) #20 individuals so 40 samples happen at the end 
sample_list = treesequence.samples()


#these have to be outside the loop below, or else you will create a clear list every loop until the last one.
samples_under_muts = []
muteffects = []

for tree in ts.trees():
    for mutation in tree.mutations(): #this is where kevin suggests to incorporate the array of effect sizes  
        # print(mutation)
        # print(node)
        # for mutation in enumerate(mut_index): 
        #     print("samples with mutation", mutation, ":")
        print("samples with mutation",mutation.id,":") 
        node = mutation.node #ts.node gives you something that looks like this: Node(id=6, flags=1, time=0.0, population=0, individual=3, metadata=b''). notice the individual number is there. ..
        #...this is the individual that contained the node in which the mutation originated
        mutlists = []
        samples_under_muts.append(mutlists)   

        
        rng = default_rng()
        # effects.append(multivariate_normal.rvs([0,0], cov = var_covar))
        muteffects.append(rng.multivariate_normal(mean = means, cov = var_covar).tolist()) #was returning array items before; https://www.journaldev.com/32797/python-convert-numpy-array-to-list
        #also, for above: https://stackoverflow.com/questions/16016959/scipy-stats-seed
        for sample in tree.samples(node): #each sample is a haploid genome, so aren't 1:1 with individuals. Only prints samples with a node that contains a mutation. If you don't specify node, all samples with ANY mutation node will be printed
            #I want to assign a given effect size to all samples that have the same mutation
            #do I need to figure out which node has the mutation and assign that effect to all children? Or can I just use samples to find individuals

            #write something here to replace the entries in blank_array with effect sizes sampled from mut_index # https://stackoverflow.com/questions/2582138/finding-and-replacing-elements-in-a-list ? 
            # print(sample,"node",node) #I've made the sample size really low so that you can see all the samples that contain mutations easily. 
            print(sample)

            mutlists.append(sample) #this makes it so I have a list of mutation lists, which contain the samples that contain those mutations. 
print()
print()
# print("samples with a given mutation:", dict(enumerate(samples_under_muts)))
print("samples with a given mutation:", samples_under_muts)

# print("effect sizes for each mutation:", dict(enumerate(muteffects)))
print("effect sizes for each mutation:", muteffects)

#ALTERNATE WAY: to make a dict
# effects = multivariate_normal.rvs(0, 1, (len(ts.tables.mutations
# effects = dict(enumerate(multivariate_normal.rvs(0, 1, (len(ts.tables.mutations))))) #this way will eventually lead to problems: 
# # Traceback (most recent call last):
# #   File "/home/pdgardne/mathieson_sims/model.py", line 184, in <module>
# #     effects = dict(enumerate(multivariate_normal.rvs(0, 1, (len(ts.tables.mutations)))))
# # TypeError: 'numpy.float64' object is not iterable
# print("effect sizes for each mutation:",effects)


# rng = default_rng()
# effects = rng.multivariate_normal([0,0], cov = var_covar)
# print("effects:",effects)



# def chunker(seq, size): #https://stackoverflow.com/questions/434287/what-is-the-most-pythonic-way-to-iterate-over-a-list-in-chunks/434411#434411
#     return (seq[pos:pos + size] for pos in range(0, len(seq), size))
# samples_by_indv = [] 
# for group in chunker(sample_list, 2): #again, from https://stackoverflow.com/questions/434287/what-is-the-most-pythonic-way-to-iterate-over-a-list-in-chunks/434411#434411
#     samples_by_indv.append(group)
# print(samples_by_indv)



#numpy solutions: 
# samples_by_indv = np.array_split(sample_list, deme_size*2) #https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks

samples_by_indv = []
for int1, int2 in sample_list.reshape(-1,2): #https://stackoverflow.com/questions/434287/what-is-the-most-pythonic-way-to-iterate-over-a-list-in-chunks/434411#434411
    samples_by_indv.append([int1, int2])

print("samples grouped by their individual:", samples_by_indv)
print()

#make a dict indexing individual by the samples they contain
# print(samples_by_indv)
# individuals_dict = dict(enumerate(samples_by_indv))
# print("samples keyed to individuals:",individuals_dict)



""""""

# print(samples_under_muts)
# mutations_dict = dict(enumerate(samples_under_muts))
# print("samples keyed to mutations:",mutations_dict)


# for indv in samples_by_indv:
#     for sample in indv: 


for mutation in samples_under_muts:
    for sample in mutation:
        for effect in muteffects: 
            print(effect)
        




# x = {mut: {} }



# print(multivariate_normal.rvs(0, 1, (len(ts.tables.mutations))))
# print(samples_under_muts)



# print()
# for mut in samples_under_muts:
#     # print(mut)
#     print()
#     for spl in mut: 
#         print(spl)



#get individual id# from samples list
# for i in range(len(samples_by_indv)):
#     print(i)

# for i in enumerate(samples_by_indv):
#     print(i)


#what needs to happen: 
#make effect sizes for each mutation
#be able to access which individual has those mutations, based on the samples
#write the effect sizes to the blank array. 
#make it pleiotropic 

#put effect sizes into nempty numpy array
#make an additive model 

#will have to make a list or dict of effect sizes for each mutation.
#then, will have to iterate over samples_under_muts to add that effect size to every sample with the mutation
#then, maybe do the dictionary grouping thing again, but on the effect sizes to get effects for the individuals?? and then you can access that to write it to the phenotype? 

#custom distribution? for their effect sizes function? https://scicomp.stackexchange.com/questions/1658/define-custom-probability-density-function-in-python
#same as above? https://stackoverflow.com/questions/4265988/generate-random-numbers-with-a-given-numerical-distribution