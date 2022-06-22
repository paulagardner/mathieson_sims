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
from itertools import groupby

demography = msprime.Demography()
graph = demes.load("no_ancestry.yaml")

mu = ((1e-7))
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

###########33may not be necessary, made for working w/ the tuple
def list_duplicates(sequence): #https://stackoverflow.com/questions/9835762/how-do-i-find-the-duplicates-in-a-list-and-create-another-list-with-them/ This only works with this tuple pipeline i'm doing, since making a list of lists below renders it unhashable
    seen = set()
    seen_add = seen.add
    seen_mult = set(x for x in sequence if x in seen or seen_add(x))
    return(list(seen_mult))




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

        # print("samples with mutation",mutation.id,":") 
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

# samples_under_muts = np.array(samples_under_muts, dtype = list)
print("samples under mutations:", samples_under_muts)
samples_under_muts_dict = dict(enumerate(samples_under_muts))
print("samples with a given mutation:", samples_under_muts_dict) 
print()

# muteffects = np.array(muteffects)
print("mutation effects:", muteffects)
mutation_effects_dict = dict(enumerate(muteffects))
print("mutations and effects dictionary:", mutation_effects_dict)
print()

#numpy solutions: 
# samples_by_indv = np.array_split(sample_list, deme_size*2) #https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks

samples_by_indv = []
for int1, int2 in sample_list.reshape(-1,2): #https://stackoverflow.com/questions/434287/what-is-the-most-pythonic-way-to-iterate-over-a-list-in-chunks/434411#434411
    samples_by_indv.append([int1, int2])

#Printing the below will confirm that you can group nodes to individuals in the above manner.
# for individual in ts.individuals():
#     print(individual)

print("samples grouped by their individual:", samples_by_indv)



#make a dict indexing individual by the samples they contain
# print(samples_by_indv)
# print()
individuals_dict = dict(enumerate(samples_by_indv))
print("samples keyed to individuals:",individuals_dict)




# print()

a = [] #https://stackoverflow.com/questions/31250129/python-numpy-array-of-numpy-arrays "Never append to numpy arrays in a loop: it is the one operation that NumPy is very bad at compared with basic Python"
b = []
# muteffects = np.asarray(muteffects)
# print("muteffects as array:",muteffects)


for samples_with_mut, effects_for_mut in zip(samples_under_muts, muteffects):
    # print(samples_with_mut, effects_for_mut)
    for sample in samples_with_mut:
        for k, v in individuals_dict.items():
            if sample in v:
                print("individual",k,"has an allele with the following effect:",effects_for_mut)
                a.append(k)
                b.append(effects_for_mut)
print()
# print()   
# #problem with using a normal dict as below: dictionary keys cannot be repeated, so if an individual has multiple mutations that we want to add together, the dict will only preserve the LAST it encounters. (last: https://www.guru99.com/python-dictionary-append.html). This could be useful for getting the info into                  
# # a = zip(samples_under_muts, muteffects):
# # dic = dict(zip(samples_under_muts, muteffects))
# # print(dic)

l = list(zip(a, b)) #workaround with tuple: https://stackoverflow.com/questions/33593556/zipped-array-returns-as-zip-object-at-0x02b6f198
# print("tuple list:",l)
# print()




dups = (list_duplicates(a))
print("duplicate individuals:",dups)
print()


# #so, this does indeed look a lot like l to begin with, but it's useful in that it's just items where the i entries appear multiple times.
h = []
print("###########################")
for i, j in l:
    if i in dups:
        print(i,j)


# #so, have a function that checks your tuple for duplicates. make your dictionary as you would the list above, and then append the duplicate in (and overwrite it)
# dic = dict(zip(a,b))
# # print(dic)

















# a = []

# #DOING THE ABOVE, BUT MAKES LIST INSTEAD OF THE TUPLE (And makes your dup loop not work:)
# for samples_with_mut, effects_for_mut in zip(samples_under_muts, muteffects):
#     # print(samples_with_mut, effects_for_mut)
#     for sample in samples_with_mut:
#         for k, v in individuals_dict.items():
#             if sample in v:
#                 print("individual",k,"has an allele with the following effect:",effects_for_mut)
#                 a.append([k, effects_for_mut])
#                 # b.append(effects_for_mut)
# print(a)


# way to make lists (instead of tuples) that won't mess up the find_dups loop: 
# z = list(list(a) for a in zip(a,b)) #https://blog.finxter.com/zip-with-list-output-instead-of-tuple-most-pythonic-way/#:~:text=Short%20answer%3A%20Per%20default%2C%20the,a%20new%20nested%20list%20object.
# # z =np.hstack(a,b)
# print("list of indivs and allele effects:", z)

# dups = list_duplicates(a)
# print("duplicate individuals:",dups)

#####make array of just duplicate values to sum them up? and then you could make an array of just the unique stuff below and append the two? but I think it would be easier to just find the duplicates in numpy and deal with the whole array of nonzero effects
# c = []
# d = []
# for i, j in z:
#     if i in dups:
#         c.append(i)
#         d.append(j)
# c = np.asarray(c)
# d = np.asarray(d)
# e = np.column_stack((c,d))
# print(e)
# array_v = np.vstack(z)
# print(array_v)






#works without the chunk above... 
##########Dissue pops up where 'a' values (representing individuals) come up as floats
a = np.asarray(a)
# print(a)
# print(b)
b = np.asarray(b)
# print(a)
# print(b)
# print(b.ndim)

print()
c = np.column_stack((a.astype(int),b))
print("array with individual, effect sizes:",c)

dups = list_duplicates(a) #so this works with the array format--nice
print("duplicate individuals:",dups)

#this causes an error when there's no mutations, so put it in an if/elif 
# d = np.array((c[:,1].sum(), c[:,2].sum())) #example to sum by column-- use this in, presumably, a for loop to match the duplicates 
# print("sum by columns test",d)
# m = np.reshape(A, (-1, num_phenotypes))
# print(m)

q = pd.DataFrame(c)
# print(q)
r = (q.groupby([0]).aggregate(sum)) #https://stackoverflow.com/questions/35961416/how-to-merge-and-sum-duplicates-in-python #note that this doesn't preserve your first column. For my purposes this saves me the step of removing it later, but you could probably always add a column in based on index number. 
# if d[:,0] in dups:
#     print(d)
print("summed effects",r)

df = r.reindex(np.arange(n_dip_indv), fill_value = 0)
print("phenotypes as numpy array",df)

isblank = []
if dups == isblank:
    print("no dups")
else:
    s= np.array(df.iloc[:,[0,1]])
    print(s)


# effects = np.asarray(df)
# print(effects)
# numpy.where()

#https://stackoverflow.com/questions/33929389/python-how-to-find-duplicates-and-sum-their-values
#https://stackoverflow.com/questions/28706209/sum-second-value-in-tuple-for-each-given-first-value-in-tuples-using-python






















# print(effects)




# # blank = np.zeros
# for key, value in samples_under_muts_dict.items():
#     for spl in value:
#         for k, v in mutations_dict.items():
#                 print(v)

#something along the lines of: if i



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