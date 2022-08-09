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
    #not including migration info since I'm just doing two demes, so what mathieson et al did (defining which demes are next to each other) not necessary here.. at least for now

    mts = msprime.sim_mutations(ts, rate = mu, discrete_genome = False) #discrete_genome = false gives infinite sites, basically, so simplifies downstream bc you don't have to account for mult mutations at a site
    return mts

def make_vcf(vcf_path, indv_names):
    #make a vcf:
    with open(vcf_path, "w") as vcf_file:  ts.write_vcf(vcf_file, individual_names=indv_names) 
    with open(vcf_path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')] 
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
            'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    )#.rename(columns={'#CHROM': 'CHROM'}) #do not want to rename this! Plink2 is very unhappy if you take the # away 

    #isolate the header
    with open(vcf_path, 'r') as f:
        header = [l for l in f if l.startswith('##')]
    return header

    #make an id for all mutations
    id  = ["rs" + str(i) for i,_ in enumerate(df.ID)]
    #replace the pandas ID column with new IDs
    df['ID'] = id

    #gzip vcf file

def make_phenotypes():

    phenos  = 2 #make this into args...
    #make effect size distributions:
    print("simulating phenotypes as if purely environmental")
    mult_random = multivariate_normal.rvs(0, 1, (len(iid))) #assigns a random effect size to every individual in one big array. NOTE that this is the same as above, remove this redundancy.

    print("simulating mutation effects")
    # print("number of mutations:", len(ts.tables.mutations))

    #make a table with mutations indexed from 0 to m - 1
    mut_index = [i for i in range(len(ts.tables.mutations))] #goes from 0 to #mutations - 1 (again, a consequence of python counting from 0)
    # print("mutation index array: ", mut_index)

    treesequence = tskit.load("output.trees")
    # print("sample list:", treesequence.samples()) #if 20 individuals, 40 samples present at end

    #may not be needed, potentially remove (if you change the downstream)
    sample_list = treesequence.samples() 

    #make a random number generator to assign effect sizes to mutations
    rng = default_rng()


    #make the variance-covariance matrix of phenotypic effects: 
    var_covar = np.identity(phenos) #this is essentially no pleitropy- no covariance b/w the two phenotypes means they are independent
    var_covar = np.array([1, 0, 0, 1]).reshape(2,2) #here, I've redefined the variable but it's still an identity matrix. BUT this is how you might make a NON-zero covariance matrix

    #make the means matrix for the multivariate generator
    means = np.zeros(phenos)    

    #create multivariate effects. to do so, you must have mean effect sizes, the covariance matrix (from above). Since it's a matrix of effects, the matrix size must be specified also. 
    muteffects = rng.multivariate_normal(mean = means, cov = var_covar, size = (len(mut_index)))

    # print("effect sizes for each mutation:", "\n", muteffects)
    # print(len(muteffects))  
    # print(phenotypes_array)

    #may not be needed, potentially remove (if you change the downstream)
    individual_id = []

    #print statement to see which nodes go with which individuals
    # for sample, i in enumerate(ts.tables.nodes.individual):
    #     print(sample,i)

    # print(ts.tables.nodes)

    #turn the above into a variable for the for() loop below(????). NOT URGENT
    indv_array = ts.tables.nodes.individual
    # print(indv_array)

    nodes_array = []
    for node, individual in enumerate(ts.tables.nodes.individual):
        nodes_array.append(node)
    # print(nodes_array)

    array = np.column_stack((indv_array, nodes_array))
    # print(array)



    #create an empty numpy array to put the effect sizes in. It will return a table with # of phenotypes column and # of individuals rows
    phenotypes_array = np.zeros(n_dip_indv*phenos).reshape(n_dip_indv,phenos) 
    #print("empty phenotypes array:", phenotypes_array, "\n")


    #make an array to track whether an individual is homozygous or heterozygous for a mutation at a locus

    print("looping over trees to find which individuals have which mutations at which nodes")
    for tree in ts.trees(): #update_samples = True errors out with TypeError

        #loop over all mutations that appear in a tree
        for mutation in tree.mutations():
            # print("newline every time we loop over a new mutation", "\n")
            # print("mutation id", mutation.id)
            nodes = mutation.node
            # print("node the mutation originates from:", nodes)
            #alternate way to get the above:
            # print(mutation.node)

            #make a 1-D numpy array with as many zeros as there are phenotypes
            ncopies = np.zeros(len(phenotypes_array)) #one zero for every individual in the sample at the end of the sim

            #loop over all samples that have a given mutation. If you don't include nodes (mutation.node) here, you loop over ALL samples that contain ANY mutation, so it will be the same full list of samples, iterated over m # of mutations times. Including nodes does it by the sample.
            for sample_node in tree.samples(mutation.node):
                # print("this mutation is in sample node",sample_node)

                #find which 
                if sample_node in array[:,1:]:
                    item = (array[sample_node]) #now I have the individual and node # together for each one that has a mutation
                    individual_with_mut = item[:1]
                    # print(*individual_with_mut) #entechin.com/how-to-print-a-list-without-square-brackets-in-python/#:~:text=use%20asterisk%20'*'%20operator%20to%20print%20a%20list%20without%20square%20brackets
                    ncopies[individual_with_mut] += 1
                    # print("copies of the mut present in each individual:", ncopies)
                    # print("phenotypic value of the indiv preexisting from environmental effect:",phenotypes_array[individual_with_mut])
                    # print("mutation value", muteffects[mutation.id])
                    phenotypes_array[individual_with_mut] += muteffects[mutation.id]


    print("summed phenotypes:", "\n", phenotypes_array)
    return phenotypes_array


def make_popfile(phenotypes_array):

    deme_id=[[i]*deme_size for i in range(0,demes)] #https://github.com/Arslan-Zaidi/popstructure/blob/master/code/simulating_genotypes/grid/generate_genos_grid.py
    #flatten
    deme_id=[item for sublist in deme_id for item in sublist] #changes 2 arrays of, say, length 50 into one array of length 100 (for example, will vary depending on deme # and sample sizes)). Necessary to make the array the correct size for the below 



    phenotypes_array = make_phenotypes()
    txt_name = "pop.txt"

    pheno_1, pheno_2 = phenotypes_array.T #https://stackoverflow.com/questions/30820962/splitting-columns-of-a-numpy-array-easily
    popdf = np.transpose([fid, iid, deme_id, pheno_1, pheno_2]) 
    with open(txt_name, 'wb') as f:
        f.write(b'FID\tIID\tdeme_id\tphenotype1\tphenotype2\n') 
        np.savetxt(f, popdf, fmt = '%s') #why %s? https://stackoverflow.com/questions/48230230/typeerror-mismatch-between-array-dtype-object-and-format-specifier-18e




print("simulating genotypes under demographic model")
ts = simulate(mu, rho, graph)

print("writing treefile for downstream analysis")
ts.dump("output.trees")


print("making vcf for the sim")
vcf_path = "output_geno.vcf"
n_dip_indv = int(ts.num_samples / 2) 
indv_names = [f"tsk_{str(i)}indv" for i in range(n_dip_indv)]
make_vcf(vcf_path, indv_names)


demes = 2 # replace hard-coded with argparser or taking the 'populations' value from the tskit table. in this case, you just want the demes of the samples present at the end of the sim. (I think). SO regardless of pop_split.yaml or no_ancestry.yaml (whether ancestry is shared), demes A and B split off.

fid=[f"tsk_{str(i)}indv" for i in range(0,(deme_size*demes))]
iid=[f"tsk_{str(i)}indv" for i in range(0,(deme_size*demes))] #number of individuals in the sample
print("simulating phenotypes, from environmental and genetic effects")
a = make_phenotypes()

print("making .txt file that contains individuals and phenotypes")
make_popfile(a)

