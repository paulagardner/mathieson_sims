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

mu = ((1e-6))
rho = 1e-7

deme_size = 10
def simulate(mu, rho, graph):
    ts = msprime.sim_ancestry(demography = msprime.Demography.from_demes(graph), recombination_rate=rho, sequence_length = 500, samples={"A": deme_size, "B": deme_size}, random_seed=1234, discrete_genome = False) #add random seed so it's replicable
    #not including migration info since I'm just doing two demes, so what they did (defining which demes are next to each other) not necessary here.. I think! at least for now

    mts = msprime.sim_mutations(ts, rate = mu, discrete_genome = False) #discrete_genome = false gives infinite sites, basically, so simplifies downstream bc you don't have to account for mult mutations at a site
    return mts

def make_vcf(vcf_path, indv_names):
    #make a vcf:
    with open(vcf_path, "w") as vcf_file:  ts.write_vcf(vcf_file, individual_names=indv_names) 
    #
    with open(vcf_path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')] 
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
            'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    )#.rename(columns={'#CHROM': 'CHROM'}) #do not want to rename this! Plink2 is very unhappy if you take the #away 

    #isolate the header
    with open(vcf_path, 'r') as f:
        header = [l for l in f if l.startswith('##')]
    return header

    #make an id for all mutations
    id  = ["rs" + str(i) for i,_ in enumerate(df.ID)]
    #replace the pandas ID column with new IDs
    df['ID'] = id

    #gzip vcf file




###########33may not be necessary, made for working w/ the tuple
def list_duplicates(sequence): #https://stackoverflow.com/questions/9835762/how-do-i-find-the-duplicates-in-a-list-and-create-another-list-with-them/ This only works with this tuple pipeline i'm doing, since making a list of lists renders it unhashable
    seen = set()
    seen_add = seen.add
    seen_mult = set(x for x in sequence if x in seen or seen_add(x))
    return(list(seen_mult))




print("simulating genotypes under demographic model")
ts = simulate(mu, rho, graph)

print("writing treefile for downstream analysis")
ts.dump("output.trees")


# print("writing vcf")
vcf_path = "output_geno.vcf"
n_dip_indv = int(ts.num_samples / 2)
indv_names = [f"tsk_{str(i)}indv" for i in range(n_dip_indv)]
make_vcf(vcf_path, indv_names)



demes = 2 # replace hard-coded with argparser or taking the 'populations' value from the tskit table. in this case, you just want the demes of the samples present at the end of the sim. (I think). SO regardless of pop_split.yaml or no_ancestry.yaml, demes A and B split off.
deme_id=[[i]*deme_size for i in range(0,demes)] #https://github.com/Arslan-Zaidi/popstructure/blob/master/code/simulating_genotypes/grid/generate_genos_grid.py
#flatten
deme_id=[item for sublist in deme_id for item in sublist] #changes 2 arrays of, say, length 50 into one array of length 100 (for example, will vary depending on deme # and sample sizes))

fid=[f"tsk_{str(i)}indv" for i in range(0,(deme_size*demes))]
iid=[f"tsk_{str(i)}indv" for i in range(0,(deme_size*demes))] #number of individuals in the sample


phenos = 2 #make this into args...

#make the variance-covariance matrix of phenotypic effects: 
var_covar = np.identity(phenos) #this is essentially no pleitropy- no covariance b/w the two phenotypes means they are independent
var_covar = np.array([1, 0, 0, 1]).reshape(2,2)


print(var_covar)

#make the means matrix for the multivariate generator
means = np.zeros(phenos)

#create an empty numpy array to put the effect sizes in. It will return a table with #phenotypes column and #individuals rows.
effects = np.zeros(n_dip_indv*phenos).reshape(n_dip_indv,phenos) #not sure I actually end up using this....
# print("empty array for effect sizes:", effects)

 

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


samples_by_indv = []
for int1, int2 in sample_list.reshape(-1,2): #https://stackoverflow.com/questions/434287/what-is-the-most-pythonic-way-to-iterate-over-a-list-in-chunks/434411#434411
    samples_by_indv.append([int1, int2])


print("samples grouped by their individual:", samples_by_indv)



#make a dict indexing individual by the samples they contain
individuals_dict = dict(enumerate(samples_by_indv))
print("samples keyed to individuals:",individuals_dict)


indv_with_mutations = [] #https://stackoverflow.com/questions/31250129/python-numpy-array-of-numpy-arrays "Never append to numpy arrays in a loop: it is the one operation that NumPy is very bad at compared with basic Python"
muteffects_in_individuals = []
# muteffects = np.asarray(muteffects)
# print("muteffects as array:",muteffects)


for samples_with_mut, effects_for_mut in zip(samples_under_muts, muteffects):
    # print(samples_with_mut, effects_for_mut)
    for sample in samples_with_mut:
        for key, value in individuals_dict.items():
            if sample in value:
                print("individual",key,"has an allele with the following effect:",effects_for_mut)
                indv_with_mutations.append(key)
                muteffects_in_individuals.append(effects_for_mut)


#pairing the two lists above (individuals and ONE of the mutation effects they posess. might be able to do this within the loop)
paired_tuple = list(zip(indv_with_mutations, muteffects_in_individuals)) #workaround with tuple: https://stackoverflow.com/questions/33593556/zipped-array-returns-as-zip-object-at-0x02b6f198
print("tuple list:",paired_tuple)
# print()

dups = (list_duplicates(indv_with_mutations))
print("duplicate individuals:",dups)
print()

print('###################################################')

# #so, this does indeed look a lot like l to begin with, but it's useful in that it's just items where the i entries appear multiple times.
# h = []
# for i, j in paired_tuple:
#     if i in dups:
#         print(i,j)


##########Dissue pops up where 'a' values (representing individuals) come up as floats
indv_with_mutations = np.asarray(indv_with_mutations)
muteffects_in_individuals = np.asarray(muteffects_in_individuals)

print()
indvs_all_effects_array = np.column_stack((indv_with_mutations.astype(int),muteffects_in_individuals))
print("array with individual, effect sizes:",indvs_all_effects_array)

dups = list_duplicates(indv_with_mutations) 
print("duplicate individuals:",dups)

indvs_all_effects_df = pd.DataFrame(indvs_all_effects_array)
summed_effects = (indvs_all_effects_df.groupby([0]).aggregate(sum)) #https://stackoverflow.com/questions/35961416/how-to-merge-and-sum-duplicates-in-python #note that this doesn't preserve your first (individual number) column. For my purposes this saves me the step of removing it later, but you could probably always add a column in based on index number. 

print("summed effects", summed_effects)

population_effects = summed_effects.reindex(np.arange(n_dip_indv), fill_value = 0)
print("phenotypes as numpy array", population_effects)


isblank = []
if dups == isblank:
    print("no dups")
    population_effects = np.array(population_effects)
else:
    array_for_dups = np.array(population_effects.iloc[:,[0,1]])
    print("phenotypes as numpy array", array_for_dups)

# overall_phenotype = 

env_random_pheno = np.array(norm.rvs(0, 0.25, (2*(len(iid))))).reshape(20,2) #my cheap way of traying to just make totally random phenotypes that are the same length as the two-dimensional np array for mut effects on phenotypes. You will have to figure out the scaling of how to not have the environmental phenos overwhelm your mutational effects
print(env_random_pheno)

#and back to pandas again to get this all written as a txt file with headers... not sure if necessary for downstream . checking now


# def write_txt(pheno_1, pheno_2, input_array, txt_name):
#     pheno_1, pheno_2 = array_name.T #https://stackoverflow.com/questions/30820962/splitting-columns-of-a-numpy-array-easily
#     popdf = np.transpose([fid, iid, deme_id, pheno_1, pheno_2]) 
#     with open(txt_name, 'wb') as f:
#         f.write(b'FID\tIID\tdeme_id\tphenotype1\tphenotype2\n') 
#         np.savetxt(f, popdf, fmt = '%s') #why %s? https://stackoverflow.com/questions/48230230/typeerror-mismatch-between-array-dtype-object-and-format-specifier-18e




if muteffects == isblank:

    # print("no mutations, phenotypes entirely environmental", env_random_pheno)
    # input_array = env_random_pheno
    # txt_name = "pop.txt"
    # write.txt(pheno_1, pheno_2, input_array= input_array, txt_name = txt_name)
    #splitting np array by its columns in order to use np.transpose to write to the .txt\
    txt_name = "pop.txt"
    pheno_1, pheno_2 = env_random_pheno.T #https://stackoverflow.com/questions/30820962/splitting-columns-of-a-numpy-array-easily
    popdf = np.transpose([fid, iid, deme_id, pheno_1, pheno_2]) 
    with open(txt_name, 'wb') as f:
        f.write(b'FID\tIID\tdeme_id\tphenotype1\tphenotype2\n') 
        np.savetxt(f, popdf, fmt = '%s') #why %s? https://stackoverflow.com/questions/48230230/typeerror-mismatch-between-array-dtype-object-and-format-specifier-18e
elif  dups == isblank:

    print("no duplicate mutations", population_effects)
    input_array = population_effects
    txt_name = "pop.txt"
    # write.txt(pheno_1, pheno_2, input_array= input_array, txt_name = txt_name)
    #splitting np array by its columns in order to use np.transpose to write to the .txt
    pheno_1, pheno_2 = population_effects.T #https://stackoverflow.com/questions/30820962/splitting-columns-of-a-numpy-array-easily
    popdf = np.transpose([fid, iid, deme_id, pheno_1, pheno_2]) 
    with open('pop.txt', 'wb') as f:
        f.write(b'FID\tIID\tdeme_id\tphenotype1\tphenotype2\n') 
        np.savetxt(f, popdf, fmt = '%s') #why %s? https://stackoverflow.com/questions/48230230/typeerror-mismatch-between-array-dtype-object-and-format-specifier-18e

else:
    sum_pheno = np.add(array_for_dups, env_random_pheno)

    print("overall phenotypes", sum_pheno)
    input_array = sum_pheno
    txt_name = "pop.txt"
    # write.txt(pheno_1, pheno_2, input_array= input_array, txt_name = txt_name)
    
    #splitting np array by its columns in order to use np.transpose to write to the .txt
    pheno_1, pheno_2 = sum_pheno.T #https://stackoverflow.com/questions/30820962/splitting-columns-of-a-numpy-array-easily
    popdf = np.transpose([fid, iid, deme_id, pheno_1, pheno_2]) 
    with open('pop.txt', 'wb') as f:
        f.write(b'FID\tIID\tdeme_id\tphenotype1\tphenotype2\n') 
        np.savetxt(f, popdf, fmt = '%s') #why %s? https://stackoverflow.com/questions/48230230/typeerror-mismatch-between-array-dtype-object-and-format-specifier-18e



#custom distribution? for their effect sizes function? https://scicomp.stackexchange.com/questions/1658/define-custom-probability-density-function-in-python
#same as above? https://stackoverflow.com/questions/4265988/generate-random-numbers-with-a-given-numerical-distribution #