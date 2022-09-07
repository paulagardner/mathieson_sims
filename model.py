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
import json
import demesdraw
import matplotlib as plt
from dataclasses import dataclass

demography = msprime.Demography()
graph = demes.load("no_ancestry.yaml")

mu = ((1e-7)*2)
rho = 1e-7
bases = 1000

deme_size = 25
def simulate(mu, rho, graph):
    ts = msprime.sim_ancestry(demography = msprime.Demography.from_demes(graph), recombination_rate=rho, sequence_length = bases, samples={"A": deme_size, "B": deme_size}, random_seed=1234, discrete_genome = False) #add random seed so it's replicable
    #not including migration info since I'm just doing two demes, so what mathieson et al did (defining which demes are next to each other) not necessary here.. at least for now


    mts = msprime.sim_mutations(ts, rate = mu, discrete_genome = False) #discrete_genome = false gives infinite sites, basically, so simplifies downstream bc you don't have to account for mult mutations at a site

    return mts


#make a vcf that will have each row signifying a mutation. (thus it's keeping track of the same info as make_phenotypes does: accounting for which individuals have which mutations, and how many copies)
def make_vcf(vcf_path, indv_names):
    with open(vcf_path, "w") as vcf_file:  
        ts.write_vcf(vcf_file, individual_names=indv_names) 

    #isolate each part of the vcf so you can collect all the elements and write a vcf from it, with ID fields
    #get column headers
    with open(vcf_path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')] 
        
        df = pd.read_csv(io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
            'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t')#.rename(columns={'#CHROM': 'CHROM'}) #do not want to rename this! Plink2 is very unhappy if you take the # away. 
    
    
    #isolate header info rows
    with open(vcf_path, 'r') as vcf_file:
        header = [l for l in vcf_file if l.startswith('##')]
    
    header =''.join(header)

    print(header)

    #make an id for all mutations
    id  = ["rs" + str(i) for i in range(len(ts.tables.mutations))] #calling it rs _ because i see it that way on other VCF

    #replace the vcf ID column with new IDs
    df['ID'] = id
    # print(df)

    blankIndex=[''] * len(df)
    df.index=blankIndex
    # print(df)

    with open(vcf_path, 'w') as vcf: #this is almost correct but the headers don't work
        vcf.write(header)#use this instead of df.to_csv because that will sometimes mess up your column headers
    df.to_csv(vcf_path, sep="\t", mode='a', index=False)#if this is indented in with_open() sometimes mess up your column headers




class phenotype_constructor:
    def __init__(self):

        phenos  = 2 #make this into args...
        #make effect size distributions:
        # print("simulating phenotypes as if purely environmental")
        mult_random = multivariate_normal.rvs(0, 1, (len(iid))) #assigns a random effect size to every individual in one big array. NOTE that this is the same as above, remove this redundancy.
        
        # print("simulating mutation effects")
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
        #double check how the above relates to the phenotypes- do the non-diagonal elements behave as you're expecting?dR
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
        # print("empty phenotypes array:", phenotypes_array, "\n")
        


        #make an array to track whether an individual is homozygous or heterozygous for a mutation at a locus
        # print("looping over trees to find which individuals have which mutations at which nodes")
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

                    #find which individuals have 0, 1, or 2 mutations
                    if sample_node in array[:,1:]:
                        item = (array[sample_node]) #now I have the individual and node # together for each one that has a mutation
                        individual_with_mut = item[:1]
                        # print(*individual_with_mut) #entechin.com/how-to-print-a-list-without-square-brackets-in-python/#:~:text=use%20asterisk%20'*'%20operator%20to%20print%20a%20list%20without%20square%20brackets
                        ncopies[individual_with_mut] += 1
                        
                        # print("phenotypic value of the indiv preexisting from environmental effect:",phenotypes_array[individual_with_mut])
                        # print("mutation value", muteffects[mutation.id])
                        phenotypes_array[individual_with_mut] += muteffects[mutation.id]
                # print("copies of the mut present in each individual:", ncopies)
        
        # print("summed phenotypes:", "\n", phenotypes_array)
        # return phenotypes_array, muteffects
            self.phenotypes = phenotypes_array
            self.muts = muteffects

def make_phenotypes():
    return phenotype_constructor()


def make_popfile(phenotypes_array):

    deme_id=[[i]*deme_size for i in range(0,demes)] #https://github.com/Arslan-Zaidi/popstructure/blob/master/code/simulating_genotypes/grid/generate_genos_grid.py
    #flatten
    deme_id=[item for sublist in deme_id for item in sublist] #changes 2 arrays of, say, length 50 into one array of length 100 (for example, will vary depending on deme # and sample sizes)). Necessary to make the array the correct size for the below 


    phenotypes_array = t.phenotypes
    txt_name = "pop.txt"

    pheno_1, pheno_2 = phenotypes_array.T #https://stackoverflow.com/questions/30820962/splitting-columns-of-a-numpy-array-easily
    popdf = np.transpose([fid, iid, deme_id, pheno_1, pheno_2]) 
    with open(txt_name, 'wb') as f:
        f.write(b'FID\tIID\tdeme_id\tphenotype1\tphenotype2\n') 
        np.savetxt(f, popdf, fmt = '%s') #why %s? https://stackoverflow.com/questions/48230230/typeerror-mismatch-between-array-dtype-object-and-format-specifier-18e

    return popdf
    # return phenotypes_array

# def make_covar(): #an attempt to streamline the GWAS process
    # deme_id=[[i]*deme_size for i in range(0,demes)] #https://github.com/Arslan-Zaidi/popstructure/blob/master/code/simulating_genotypes/grid/generate_genos_grid.py
    # #flatten
    # deme_id=[item for sublist in deme_id for item in sublist] #changes 2 arrays of, say, length 50 into one array of length 100 (for example, will vary depending on deme # and sample sizes)). Necessary to make the array the correct size for the below 

    # phenotypes_array = make_phenotypes()[0] #why the indexing? If you return multiple values from the function (such as make_phenotypes), gotta specify which variable you want to access
    # txt_name = "covar.txt"
    
  
    # print(deme_id)
    # covar_df = np.column_stack((fid, iid, deme_id))
    # with open(txt_name, 'wb') as f:
    #     f.write(b'FID\tIID\tdeme_id')
    #     np.savetxt(f, covar_df, fmt = '%s') #why %s? https://stackoverflow.com/questions/48230230/typeerror-mismatch-between-array-dtype-object-and-format-specifier-18e

    # return covar_df

#check that the YAML is doing what you'd expect:
ax = demesdraw.tubes(graph)
ax.figure.savefig("A.svg")


print("simulating genotypes under demographic model")
ts = simulate(mu, rho, graph)

print("writing treefile for downstream analysis")
ts.dump("output.trees")


print("making vcf for the sim")
vcf_path = "genos.vcf"
n_dip_indv = int(ts.num_samples / 2) 
indv_names = [f"tsk_{str(i)}indv" for i in range(n_dip_indv)]
make_vcf(vcf_path, indv_names)


demes = 2 # replace hard-coded with argparser or taking the 'populations' value from the tskit table. in this case, you just want the demes of the samples present at the end of the sim. (I think). SO regardless of pop_split.yaml or no_ancestry.yaml (whether ancestry is shared), demes A and B split off.

fid=[f"tsk_{str(i)}indv" for i in range(0,(deme_size*demes))]
iid=[f"tsk_{str(i)}indv" for i in range(0,(deme_size*demes))] #number of individuals in the sample
print("simulating phenotypes from environmental and genetic effects")


t = make_phenotypes()
print("summed phenotypes", t.phenotypes)

print("mutation effect sizes",  t.muts)

print("making .txt file that contains individuals and phenotypes")
make_popfile(t.phenotypes)




# make_covar()





# print(ts.tables)
# print(ts.tables.mutations.derived_state)




###############SCHEMAS 

# print(ts.table_metadata_schemas)
# print("\n","UNMODIFIED TABLES:", "\n")
# # print(ts.tables) #what this and the above tell me is that there is metadata for the population, but not for the nodes, edges, sites, mutations, migrations, or individuals. So I want to make a schema for the individuals if i'm wanting to invoke phenotypes/add to it

# schema = tskit.MetadataSchema({
#     'codec': 'json',
#     'additionalProperties': False,
#     'properties': {'phenotype1': {'description': 'first phenotype value','type': 'number'}, #why not float? JSON does not like it. https://json-schema.org/understanding-json-schema/reference/numeric.html#number
#                    'phenotype2': {'description': 'second phenotype value','type': 'number'}},
#     'required': ['phenotype1', 'phenotype2'],
#     'type': 'object',
# }) #https://tskit.dev/tutorials/metadata.html#sec-tutorial-metadata. pretty certain this schema is correct







# ###################writing following https://tskit.dev/tutorials/metadata.html#sec-tutorial-metadata. this initiates new tables and writes to the NEW individuals table, and NOT the one you access through the ts variable returned from the actual sim. It's separate.
# tables = tskit.TableCollection(sequence_length=bases)  # make a new, empty set of tables
# tables.individuals.metadata_schema = schema


# # print(tables)
# # print(ts)
# # print("schemas:")
# print(json.dumps(schema.asdict(), indent=4))
# print(ts.table_metadata_schemas)  #this does not show a schema for the metadata
# row_id = tables.individuals.add_row(0, metadata={"phenotype1": 1, "phenotype2": 2,}) #their unspecified argument is for the flags, default value is 0 https://tskit.dev/tutorials/metadata.html#sec-tutorial-metadata
# print(f"Row {row_id} added to the individuals table")
# row_id = tables.individuals.add_row(1, metadata={"phenotype1": 1, "phenotype2": 2,}) #their unspecified argument is for the flags, default value is 0 https://tskit.dev/tutorials/metadata.html#sec-tutorial-metadata
# print(f"Row {row_id} added to the individuals table")


# row_id = tables.mutations.add_row(node=node, site=site, metadata={"phenotype1": 1, "phenotype2": 2,})

# print(tables.individuals)
# print(tables)
# print(tables.individuals.metadata_schema) #we were able to write stuff BUT this cleared all the other tables out, and doesn't work if I wan't to call metadata from ts.individual().metadata.:
# print("Metadata for individual 0:", ts.individual(0).metadata)


# print(ts.tables.individuals) #so we're not actually writing to this! we are writing to a new table collection that is NOT the same as the tables that were saved when we ran the simulate() function. 
# print(ts.tables)
# 



# ################let's see if we can force it to write to individuals that already exists
# print("\n","\n","\n","\n")
# print(json.dumps(schema.asdict(), indent=4))
# ts.tables.individuals.metadata_schema = schema
# row_id = ts.tables.individuals.add_row(0, metadata={"phenotype1": 1, "phenotype2": 2,})




# phenotype_index = [i for i in range(len(phenotypes))]


# for i in phenotype_index: 
#     phenotype1_array = phenotypes[:,0]
#     phenotype2_array = phenotypes[:,1]
#     # print(phenotype1_array[i])
#     # print(phenotype1_array[i].astype)
#     # print(np.asscalar(phenotype1_array[i]))

#     # tables.individuals.add_row(i, metadata= {"phenotype1": phenotype1_array[i], "phenotype2": phenotype2_array[i]})


# print("tskit.IndividualTable:",tskit.IndividualTable)
# print("ts.tables.individuals:",ts.tables.individuals)


# t_m = ts.tables.mutations
# m = len(ts.tables.mutations) #https://tskit.dev/tskit/docs/stable/python-api.html
# site = ts.tables.mutations.site
# print(site)
# node=ts.tables.mutations.node
# derived_state=ts.tables.mutations.derived_state
# time = ts.tables.mutations.time
# d, off = tskit.pack_strings(["1"] * m)

# site = np.arange(m, dtype=np.int32)

# node = np.zeros(m, dtype=np.int32)
# node=ts.tables.nodes
# print(node)



# metadata_schema = tskit.MetadataSchema({
#     'codec': 'json',
#     'additionalProperties': False,
#     'properties': {'phenotype1_effect': {'description': 'first phenotype value','type': 'number'}, #why not float? JSON does not like it. https://json-schema.org/understanding-json-schema/reference/numeric.html#number
#                    'phenotype2_effect': {'description': 'second phenotype value','type': 'number'}},
#     'required': ['phenotype1_effect', 'phenotype2_effect'],
#     'type': 'object',
# }) #https://tskit.dev/tutorials/metadata.html#sec-tutorial-metadata. pretty certain this schema is correct


# ts.tables.mutations.metadata_schema = metadata_schema






# m = len(ts.tables.mutations) #https://tskit.dev/tskit/docs/stable/python-api.html
# d, off = tskit.pack_strings(["1"] * m)
# # row_id = tables.mutations.add_row(site=0,node=0, derived_state=1,parent=-1, metadata={"phenotype1": 1, "phenotype2": 2,})
# # print(tables.mutations)

# print(muteffects)
# muteffects = muteffects.flatten() #flatten, as according to tskit docs, so we can use it in the metadata arg below
# print(muteffects)
# #need to have the list be in bytes format to be able to pack the bytes
# muteffects_bytes = muteffects.tobytes() #this doesn't seem to be formatted correctly for tskit.pack_bytes, as the bytes format you get looks quite different from what's specified here:
# print(muteffects_bytes)
# # a2 = np.frombuffer(muteffects_bytes, dtype=muteffects.dtype)
# # print(a2)
# a = tskit.pack_bytes(muteffects_bytes)



# ts.tables.mutations.set_columns(site=site, node=node, derived_state=derived_state,derived_state_offset=off, time=time)



#put your own data into metadata so that when you make your downstream, you can pull all of the information you need from the tree sequence metadata instead of writing to a million files


# Define a Python class to represent
# your metadata.
@dataclass
class MutationMetadata:
    effect_sizes: list


# Define a method to allow the built-in json
# module to encode your class as json.
# This simply returns the dict representation.python 
# This is a modification from stuff at the json
# module doc site.
class MutationMetadataEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, MutationMetadata):
            return obj.__dict__
        # Base class default() raises TypeError:
        return json.JSONEncoder.default(self, obj)


#make a new table collection to hold everything
tables = tskit.TableCollection(1.0)

tables.nodes.add_row(0, 0)  # flags, time
tables.sites.add_row(1.0, "")  # position, ancestral state
tables.mutations.add_row(0, 0, "")  # site, node, derived state

# Define schema
mutation_metadata_schema = tskit.metadata.MetadataSchema(
    {
        "codec": "json",
        "type": "object",
        "name": "Mutation metadata",
        # See tskit docs for what proprties are allowed
        "properties": {"effect_sizes": {"type": "array"}},
    }
)

# Add schema to table
tables.mutations.metadata_schema = mutation_metadata_schema

# Create a metadata object for the first mutation
md = MutationMetadata([1, 2, 3, 4])

# convert it to json
md_as_json = bytes(json.dumps(md, cls=MutationMetadataEncoder), "utf-8")

# pack it
# NOTE: I put md_as_json into a list here!!!
#       You pass in a list of things to pack
md, offsets = tskit.pack_bytes([md_as_json])

# Efficiently rebuild the mutation table in-place
tables.mutations.set_columns(
    site=tables.mutations.site,
    derived_state=tables.mutations.derived_state,
    derived_state_offset=tables.mutations.derived_state_offset,
    node=tables.mutations.node,  # This is the last required one
    metadata=md,
    metadata_offset=offsets,
)

# Confirm that things work
for m in tables.mutations:
    print(m)

print(tables.mutations)
#print(tables)


# # Method 2
# mutations = tables.mutations.asdict()
# mutations["metadata"] = md
# mutations["metadata_offset"] = offsets
# tables.mutations.set_columns(**mutations)
# # Confirm that things work
# for m in tables.mutations:
#     print(m)
    