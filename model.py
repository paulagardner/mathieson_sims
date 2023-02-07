from re import A
import msprime
import io
import demes 
import numpy as np
from numpy.random import default_rng
rng = default_rng()
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
import matplotlib.pyplot as plt #if it's just matplotlib, your demesdraw subplots part won't work
import seaborn as sns
from dataclasses import dataclass
import random
import argparse



# def make_parser() -> argparse.ArgumentParser:
#     # make an argument parser that can output help.
#     # The __file__ is the full path to this file
#     ADHF = argparse.ArgumentDefaultsHelpFormatter
#     parser = argparse.ArgumentParser(__file__, formatter_class=ADHF)

#     # Add an argument
#     parser.add_argument("--YAML", type=str, default = None help="Yaml file input")
#     parser.add_argument(
#         "--treefile",
#         "-t",
#         type=str,
#         default=None,
#         help="Tree file output name (tskit format)",
#     )


#     parser.add_argument("--N", type=int, help="Number of individuals")

#     parser.add_argument("--MU", "--u", type=float, help="Mutation rate")

#     parser.add_argument(
#         "--POPT", type=float, default=0, help="Population optimum trait value"
#     )

#     parser.add_argument(
#         "--VS",
#         type=float,
#         help="Variance of S, the inverse strength of stabilizing selection",
#     )

#     parser.add_argument(
#         "--E_SD",
#         type=float,
#         help="Environmental effect distribution's standard deviation",
#     )

#     parser.add_argument(
#         "--E_MEAN", type=float, help="Environmental effect distribution's mean"
#     )

#     return parser


# def validate_args(args: argparse.Namespace):
#     if args.treefile is None:
#         raise ValueError(f"treefile to be written cannot be None")

#     # if args.seed < 0:
#     # raise ValueError(f"invalid seed value: {args.seed}")

#     if args.POPT is None:
#         raise ValueError(f"Population optimum trait value cannot be None")

#     if args.N is None:
#         raise ValueError(
#             f"Number of individuals cannot be None"
#         )  # this was giving me trouble. Error was AttributeError: 'Namespace' object has no attribute 'N'

#     if args.MU is None:
#         raise ValueError(f"Mutation rate cannot be None")

#     if args.VS is None:
#         raise ValueError(
#             f"In this simulation using stabilizing selection, VS cannot be None"
#         )

demography = msprime.Demography()
graph = demes.load("pop_split.yaml")
# graph = demes.load("no_ancestry.yaml")
graph2 = demes.load("no_ancestry.yaml")

mu = (1e-4) #normally would have this set to 1e-7 *2, but putting deme size down and #individuals up to be able to visualize tree sequences for working on placing mutations issue.
rho = 1e-7
bases = 1000
#these values of mu and bases don't make a whole lot of sense re: biology, but the math should work out regardless and this will consistently give you small sims with actual mutations showing up.


deme_size = 5 #scaling this to 1000 doesn't affect much 

populations = 2 # replace hard-coded with argparser or taking the 'populations' value from the tskit table. in this case, you just want the demes of the samples present at the end of the sim. (I think). SO regardless of pop_split.yaml or no_ancestry.yaml (whether ancestry is shared), demes A and B split off.
def simulate(mu, rho, graph):
    ts = msprime.sim_ancestry(demography = msprime.Demography.from_demes(graph), recombination_rate=rho, sequence_length = bases, samples={"A": deme_size, "B": deme_size}, random_seed=3456, discrete_genome = False) #add random seed so it's replicable
    #not including migration info since I'm just doing two demes, so what mathieson et al did (defining which demes are next to each other) not necessary here.. at least for now
    return ts

    

    

    # mts = msprime.sim_mutations(ts, rate = mu, discrete_genome = False, ) #discrete_genome = false gives infinite sites, basically, so simplifies downstream bc you don't have to account for mult mutations at a site
    # #random_seed = 1234--- how do I include the seed?
    # return mts 
    
    


def mutation_placer():
    newtables = ts.dump_tables() #modifying a copy of the original tables in order to be able to overwrite 
    #make a new table collection to hold everything
    edge_table = newtables.edges
    mut_table = newtables.mutations
    node_table = newtables.nodes
    site_table = newtables.sites
    for tree in ts.trees():
        branch_lengths = np.zeros(tree.tree_sequence.num_nodes)
        # print(len(branch_lengths))
        dict = {}
        summed_branch_lengths = 0
        nodes_list = np.zeros(tree.tree_sequence.num_nodes)

        #find the interval the tree exists on 
        interval = tree.interval 
        # print("interval of the tree", interval)

        #find the fraction of the total genome this occupies/number of bases
        span = tree.interval.right-tree.interval.left
        # print("span of this tree", span)        

        number_mutations = 0
        for node in tree.nodes(): #type error: tree object is not iterable-- this happens w/ you frequently b/c you're forgetting to include the parenthesis to call the method properly, so the code tries to iterate over the method and not the data structure it's referencing 
            
            #add 
            dict[node] = tree.branch_length(node) #https://stackoverflow.com/questions/43505921/python-trying-to-create-a-dictionary-through-a-for-loop
            
            branch_lengths[node] =+ tree.branch_length(node) 
            nodes_list[node] =+ node
            #kevin's suggestion of how to do it- note that you've largely figured it out already in the line above 
        

            
            child_time = node_table.time[node]
            for index, j in enumerate(edge_table.child): #I think the child node is the appropriate one to pick, since we chose nodes BELOW branches. Mutations are associated with the node below them (according to tskit data model page)
                # print(index, j) #this print statement coupled with the one below lets you check if it's behaving as intended
                if j == node:
                    entry = edge_table[index]
                    parent = entry.parent
                    parent_time = node_table.time[parent]
                    # print("parent time", parent_time, 'node', parent)
                    branch_length = (parent_time - child_time)
                    summed_branch_lengths = summed_branch_lengths + branch_length 

                    

                    ####place mutations on nodes proportional to branch lengths, and proportional to how long of a span of the genome a given tree occupies
                    lineage_mutations = ((rng.poisson(lam=(mu*branch_length*(span/bases)))))  #my idea to populate the sim with the expected number of mutations given the parameters I'm using: get mu * t in here, after the hudson 2015 paper having lineage mu * branch length be the number of mutations occuring along that lineage. Summing them all up to get expectations for each tree. I'm multiplying by the fraction of the genome, since I believe mu* branch length as the expectation for the # of mutations assumes the whole genome


                    #i'm choosing here to use expectations I see in the hudson (2015) and, ultimately, waterson (1975) around mu to be per-base, instead of per-gene (which is what the waterson paper specifies). Come back to this to see if that is justifiable, but I feel that the hudson math should work for per-base instead of per-gene as both should be poisson processes, and whatever way I get to number of mutations, if I have the # of mutations right the SFS should work out
                    number_mutations += lineage_mutations

                    # print("number of mutations expected to go on this node", lineage_mutations)
            

        print("total number of mutations that should happen on this tree", number_mutations)
        print("fraction of the genome this tree occupies", span/bases)
                # time_span = (child_time, parent_time)
            #get branch lengths as a proportion of total
            # child_time = node_table.time
            # print("child time", child_time)
            # entry = edge_table[node]
            # parent_node = entry.parent
           
            # print("parent time", parent)
        print("summed branch lengths:", summed_branch_lengths)
        # total_branch_length = sum(tree.branch_length(u) for u in tree.nodes()) #https://tskit.dev/tutorials/analysing_trees.html note how this compares to the summed branch length
        # print(f"Total branch length: {total_branch_length}")
        proportional_branch_lengths = branch_lengths/summed_branch_lengths 
        # print(proportional_branch_lengths)




        

        


        #randomly sample nodes according to the branch lengths above them F
        #to do that, I'll need to have k (the number of draws numpy.choice/choices is making) specified by a probability distribution that takes in mutation rate
        #mu is from above
        # print("size of distribution", mu*bases)

        p = np.asarray(proportional_branch_lengths).astype('float64')
        p = p/np.sum(p) #not sure if doing this as a workaround to the proportional branch lengths array not summing to 0 is the best idea, but to get it working for now
        # print("probabilities",p, "sums", sum(p))

        #####you need to determine the size parameter using the mutation rate, the span of the tree that the mutations are happening on (and how large of a fraction of the bases that is)
        #notes from meeting: the mu that I'm normally thinking of is a per base pair mean, so you'll want to have a mean number of mutations:

        #mu, in terms of the values you're used to, is a per-base mutation rate. 
        #number of mutataions on any one branch 
        # number_mutations = np.random.poisson(lam=(mu*bases*))
        # print("NUMBER MUTATIONS", number_mutations) #this is a smaller number than what the msprime sim gives you- figure out why 

        #number of mutations: u times the EXPECTED sum of lengths of the branches of the tree-- equation 4 from the hudson paper, aka theta times the sum from i to n-1 of 1/i, i being the number of branches???


        #mutation rate according to span:
        # for 


        nodes_sampled =  rng.choice(nodes_list, p=p, size = int(number_mutations)) #this can sample the same node twice, which IS what you want. however need to figure out a way to determine which parent node was associated with a given branch length, so you can select a site that makes sense given the edges possible for a child:parent node relationship 
        # nodes_sampled =  np.random.choice(nodes_list, p=p, size = np.random.choice(expected mutations))

        nodes_sampled = nodes_sampled.astype(int) ##these need to be integers to be able to use them to access indices/values in the edge table. See if you can make it such when the  variable is created 
        # print("nodes sampled from distribution", nodes_sampled)



        #pick out the spot in the genome that the tree we have iterated over actually occupies (with edges) 

        print("interval this tree occupies", tree.interval) #access this to pick which sites are possible
        # print("left side of interval", interval.left)
        for node in nodes_sampled: #since I do a pretty much identical thing above, there's gotta be a way to combine the two, but I'll leave the stunts for later 
            ######draw a site from the genome that HASN'T ALREADY had a mutation assigned to it. (can probably skip making the sites infinite for now)
            #use these nodes to write mutations to the mutations table

            #kevin suggests I figure out a solution to do the above with sets or dictionaries, which should be much faster than comparing arrays, which will scale up linearly, but sets should not. IN addition, you can't just replace the range argument with interval.left and interval.right, as they're not integer and numpy does not appreciate that
            #try using set differences 
            possible_bases = set(range((int(interval.left)), int(interval.right))) #set them to integers for now?
            bases_already_mutated = set(int(i) for i in site_table.position)
            # print("sites from site table", bases_already_mutated)

            unmutated_sites =(possible_bases.difference(bases_already_mutated)) #get the difference between two sets- it returns only those bases NOT in bases_already_mutated
            # print("set difference test", x)

            site = random.choice(tuple(unmutated_sites))#suggestion from one stackoverflow entry for this(https://stackoverflow.com/questions/15837729/random-choice-from-set-python) links to another to try to make this process constant time???: https://stackoverflow.com/questions/15993447/python-data-structure-for-efficient-add-remove-and-random-choice

            # print(possible_bases)



            #####add binary mutations in- 0/1, whether it's been mutated or not. Or could use just two nucleotides, for instance
            # nucleotides = ["A", "T", "C", "G"] #when I wanted to code them in
            site_table.add_row(position=site, ancestral_state = "0")

            #####get times for the mutations- randomly select them from the interval between parent and child node times 

            #get child node time
            child_time = node_table.time[node]
            # print("node ID", node,  ts.tables.edges.child[node])#https://stackoverflow.com/questions/176918/finding-the-index-of-an-item-in-a-
            
            for index, j in enumerate(edge_table.child): #I think the child node is the appropriate one to pick, since we chose nodes BELOW branches. Mutations are associated with the node below them (according to tskit data model page)
                # print(index, j) #this print statement coupled with the one below lets you check if it's behaving as intended
                if j == node:
                    entry = edge_table[index]
                    # print("node", node, "edges associated with the node selected", entry)
                    current_tree_interval = (interval.left, interval.right)

                    if (entry.left, entry.right) == current_tree_interval: 
                        parent_node = entry.parent
                        # print("tuple matched, parent node", parent_node)

                    elif (entry.left, entry.right) == (0, bases): #this happens when the node is not associated with a recombination event 
                        # print("this node has the same parent regardless of tree")
                        parent_node=entry.parent

            #now, find actual parent time:
            parent_time = node_table.time[parent_node]
            # print("parent time", parent_time)
            time_span = (child_time, parent_time)
            # print("time interval between a node and its parent", time_span)
            time=random.choice(time_span) #you get a suspiciously high amount of '0' times. 

            #####assign all variables to the mutations table 
            mut_table.add_row(site=site, node=node, time=time, derived_state = "1") #this doesn't work as with the above-- perhaps something with the strings. 
            # print(tskit.unpack_strings(site_table.ancestral_state, site_table.ancestral_state_offset))

            # print(np.setdiff1d(nucleotides, site_table.ancestral_state))


        

    print("LOOP BREAK", "\n")
  
    expected_sfs = []
    theta = 4*(populations*deme_size)*(mu)
    for i in range(1,10):
        expected_sfs.append(theta/i)
    print("Expected site frequency spectrum", expected_sfs)
    # simulated_sfs = []
    # a = ts.allele_frequency_spectrum(sample_sets=[ts.samples(population=1), ts.samples(population=2)], windows=None, mode='site', span_normalise=True, polarised=False)
    # print("allele frequency spectrum function", a)


    # ts = newtables.tree_sequence()
    print(newtables)
    return(newtables)




class phenotype_constructor:
    def __init__(self):

        phenos  = 2 #make this into args...
        #make effect size distributions:


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
        # var_covar = np.identity(phenos) #this is essentially no pleitropy- no covariance b/w the two phenotypes means they are independent
        # var_covar = np.array([1, 0, 0, 1]).reshape(2,2) #here, I've redefined the variable but it's still an identity matrix. BUT this is how you might make a NON-zero covariance matrix
        # #double check how the above relates to the phenotypes- do the non-diagonal elements behave as you're expecting?dR
        # #make the means matrix for the multivariate generator

        # #no correlation b/w pleiotropic effects
        # correlation_matrix = np.array([1, 0, 0, 1]).reshape(2,2)  #var-covar = 0.8, 0, 0, 0.8
        # #after the fwdpy11 vignette: you need to convert this correlation matrix into a covariance matrix. To do that, you'll need a vector of standard deviations
        # sd = np.array([0.894427191, 0.894427191]) #the variances in the var-covar matrix come from this vector. since variance is sd^2, you must take the square root of your desired heritability to put in here

        #1/2 correlation:
        correlation_matrix = np.array([1, 0.5, 0.5, 1]).reshape(2,2) #var-covar = 0.8, 0.4, 0.4, 0.8
        sd = np.array([0.894427191, 0.894427191])

       # 1:1 correlation
        # correlation_matrix = np.array([1, 1, 1, 1]).reshape(2,2) #var-covar = 0.8, 0.8, 0.8, 0.8
        # sd = np.array([0.894427191, 0.894427191])







        D = np.identity(2) #argument is however many phenotypes there are 
        np.fill_diagonal(D, sd)
        var_covar = np.matmul(np.matmul(D, correlation_matrix), D)
        print("variance-covariance matrix", var_covar)
        



        means = np.zeros(phenos) 
        environmental_means = np.zeros(phenos) #assign mean environmental effect with argparser 
           

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

        indvs_and_nodes_array = np.column_stack((indv_array, nodes_array))
        # print(indvs_and_nodes_array)



        #create an empty numpy array to put the effect sizes in. It will return a table with # of phenotypes column and # of individuals rows
        genetic_effects_array = np.zeros(n_dip_indv*phenos).reshape(n_dip_indv,phenos) 
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
                ncopies = np.zeros(len(genetic_effects_array)) #one zero for every individual in the sample at the end of the sim

                #loop over all samples that have a given mutation. If you don't include nodes (mutation.node) here, you loop over ALL samples that contain ANY mutation, so it will be the same full list of samples, iterated over m # of mutations times. Including nodes does it by the sample.
                for sample_node in tree.samples(mutation.node):
                    # print("this mutation is in sample node",sample_node)

                    #find which individuals have 0, 1, or 2 mutations
                    if sample_node in indvs_and_nodes_array[:,1:]:
                        item = (indvs_and_nodes_array[sample_node]) #now I have the individual and node # together for each one that has a mutation
                        individual_with_mut = item[:1]
                        # print(*individual_with_mut) #entechin.com/how-to-print-a-list-without-square-brackets-in-python/#:~:text=use%20asterisk%20'*'%20operator%20to%20print%20a%20list%20without%20square%20brackets
                        ncopies[individual_with_mut] += 1
                        
                        # print("phenotypic value of the indiv preexisting from environmental effect:",phenotypes_array[individual_with_mut])
                        # print("mutation value", muteffects[mutation.id])
                        genetic_effects_array[individual_with_mut] += muteffects[mutation.id]
                # print("copies of the mut present in each individual:", ncopies)
        
        


        #make environmental effects:
        #getting a seed for the random number generator??? https://stackoverflow.com/questions/16016959/scipy-stats-seed
        environmental_effects = np.zeros(n_dip_indv*phenos).reshape(n_dip_indv,phenos)
        # env_generator = random.normal(mean=0, seed =1234)
        
        env_generator = rng.normal(loc=environmental_means, scale=1, size = (n_dip_indv,phenos)) #scale is the standard deviation
        # print("environmental effects:",env_generator)
        # print("genetic effects:", genetic_effects_array)
        
        # environmental_effects 
        # environmental_effects = multivariate_normal.rvs(0, 1, (len(iid))) #assigns a random,environmental value to every individual in one big array. 
        #the .rvs method here confuses me and I'd like to be able to set a seed- figure that out at some point 
        # print("ENVIRONMENTAL EFFECTS:", environmental_effects)
        #conceptually, you're going to have to decide how the environmental effects are going to play out. For now, just have the environment have identical effects on both phenotypes for an individual. 


        summed_phenotypes = np.add(genetic_effects_array, env_generator)
        # print(summed_phenotypes)

        # print("summed phenotypes:", "\n", phenotypes_array)
        # return phenotypes_array, muteffects
        self.genetic = genetic_effects_array
        self.environmental = env_generator
        self.phenotypes = summed_phenotypes
        self.muts = muteffects #make sure these are outside of the for() loop!
        self.environmental_means = environmental_means
        # print("PHENOTYPES",phenotypes_array) 


def make_phenotypes():
    return phenotype_constructor()


def add_metadata_to_treefile(newtables):
    # Define a Python class to represent
    # your metadata.
    @dataclass
    class MutationMetadata:
        effect_sizes: list

    @dataclass 
    class IndividualMetadata:
        full_phenotypes: list
        environmental_effects: list
        genetic_effects: list
        # env_effects: list   #might have to do some work to get  


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

    class IndividualMetadataEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, IndividualMetadata):
                return obj.__dict__
            return json.JSONEncoder.default(self, obj)

    #original place for this 
    # #make a new table collection to hold everything
    # newtables = ts.dump_tables() #modifying a copy of the original tables in order to be able to overwrite 



    # tables.nodes.add_row(0, 0)  # flags, time
    # tables.sites.add_row(1.0, "")  # position, ancestral state
    # tables.mutations.add_row(0, 0, "")  # site, node, derived state

    # Define schema. It has to be a dictionary of JSON objects.
    mutation_metadata_schema = tskit.metadata.MetadataSchema(
        {
            "codec": "json",
            "type": "object",
            "name": "Mutation metadata",
            # See tskit docs for what proprties are allowed
            "properties": {"effect_sizes": {"type": "array"}},
        }
    )

    #making a dictionary with two key:value pairs. Each will be a list of values for each phenotype class
    individual_metadata_schema = tskit.metadata.MetadataSchema(
        {
            "codec": "json",
            "type": "object",
            "name": "Individual metadata",
            # See tskit docs for what proprties are allowed
            "properties": {"full_phenotypes": {"type": "array"}, "environmental_effects":{"type":"array"}, "genetic_effects":{"type":"array"}},
        }
    )

    
    
    # "environmental_effects": {"type":"array"}


    # Add schema to table
    newtables.mutations.metadata_schema = mutation_metadata_schema
    newtables.individuals.metadata_schema = individual_metadata_schema



    ###################################################making mutataions into a format that works
    muts = constructor.muts.tolist() #taking just the muts from the constructor class
    metadata_column = []
    # print(muts)
    for index, mut in enumerate(muts):
        # newtables.nodes.add_row(0, 0)  # flags, time
        # newtables.sites.add_row(index, "")  # position, ancestral state
        # newtables.mutations.add_row(index, 0, "")  # site, node, derived state #these are necessary if you're starting from a new table, but since we're modifying the existing table instead, no need
        
        #Create a metadata object for each mutation in the list you're looping over
        md = MutationMetadata(mut)
        # print(md)
        #convert it to json
        md_as_json = bytes(json.dumps(md, cls=MutationMetadataEncoder), "utf-8")
        # md_as_json = json.dumps(md, cls=MutationMetadataEncoder)
        # print(md_as_json)
        metadata_column.append(md_as_json)
        
        # # pack it
        # # NOTE: I put md_as_json into a list here!!!
        # #       You pass in a list of things to pack
        # md, md_offset = tskit.pack_bytes([md_as_json]) #this doesn't work inside the loop, as each time the tables get rese
        




        # print("LOOP BREAK", "\n")

    ##########################validate_and_encode step that appears to be unnecessary? 
    # encoded_metadata_column = [
    #     newtables.individuals.metadata_schema.validate_and_encode_row(r) for r in metadata_column
    # ]
    # md, md_offset = tskit.pack_bytes(encoded_metadata_column)
    #this is the same as the step inside the loop..
    md, md_offset = tskit.pack_bytes(metadata_column)#looks like it doesn't matter which one you use, so this may circumvent the issue you're having with encoded_metadata_col;umn using the individuals table and issuing an object problem??????
    #this works, in a way that the encode_and_validate step seems to not once you include individual metadata (something to do with how the mutation metadata is referencing the individual metadata). 
    


    # Efficiently rebuild the mutation table in-place (method 1)
    newtables.mutations.set_columns(
        site=ts.tables.mutations.site,
        derived_state=ts.tables.mutations.derived_state,
        derived_state_offset=ts.tables.mutations.derived_state_offset,
        node=ts.tables.mutations.node,  # This is the last required one
        time= ts.tables.mutations.time,
        metadata=md,
        metadata_offset=md_offset,
    ) #rebuild using original tables for most arguments





    ###################################################################set it up for environmental effects 
    
    phenotypes = constructor.phenotypes.tolist()
    genetic_effects = constructor.genetic.tolist()
    environmental_effects = constructor.environmental.tolist()


    #print(phenotypes)
    #print("get environmental effects from constructor function", environmental_effects)
    metadata_column = []
    # print(muts)
    # for index, phenotype in enumerate(phenotypes):
    #     md = IndividualMetadata(phenotype,environmental_effects)#this almost works, but what it appends is the correct phenotype, then ALL environmental effects for all individuals.
    #     md_as_json = bytes(json.dumps(md, cls=IndividualMetadataEncoder), "utf-8")
    #     metadata_column.append(md_as_json)
        

    # for index, phenotype in enumerate(phenotypes):
    #     print(phenotypes)


    for index, phenotype in enumerate(phenotypes):
        md = IndividualMetadata(phenotype,environmental_effects[index], genetic_effects[index])#this almost works, but what it appends is the correct phenotype, then ALL environmental effects for all individuals.
        md_as_json = bytes(json.dumps(md, cls=IndividualMetadataEncoder), "utf-8")
        metadata_column.append(md_as_json)

    # for index, environmental_effect in enumerate(environmental_effects):
    #     md = IndividualMetadata(environmental_effect)#this almost works, but what it appends is the correct phenotype, then ALL environmental effects for all individuals.
    #     md_as_json = bytes(json.dumps(md, cls=IndividualMetadataEncoder), "utf-8")
    #     metadata_column.append(md_as_json)

 

    md, md_offset = tskit.pack_bytes(metadata_column)#looks like it doesn't matter which one you use, so this may circumvent the issue you're having with encoded_metadata_col;umn using the individuals table and issuing an object problem??????
    # print(ts.tables.individuals.flags)
    newtables.individuals.set_columns(
        flags=ts.tables.individuals.flags,
        location=None,
        location_offset=None,
        parents=None,
        parents_offset=None,
        metadata=md,
        metadata_offset=md_offset)




    # print(newtables)

    # print(newtables)
    # print(newtables)
    newtables.dump("output.trees")




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
ax = demesdraw.tubes(graph2)
ax.title.set_text('no ancestry')
ax.figure.savefig("no_ancestry.svg")

ax = demesdraw.tubes(graph)
ax.title.set_text('population split')
ax.figure.savefig("pop_split.svg")

#attempt to plot them side-by side
# fig, (ax1, ax2) = plt.subplots(1,2)
# ax1 = demesdraw.tubes(graph)
# ax2 = demesdraw.tubes(graph2)
# plt.savefig("A.svg")

# fig = plt.figure()
# ax1 = demesdraw.tubes(graph)
# ax2 = demesdraw.tubes(graph2)

# fig = plt.figure()
# ax1 = fig.add_subplot(demesdraw.tubes(graph))
# plt.savefig("A.svg")




print("simulating genotypes under demographic model")

ts = simulate(mu, rho, graph) #NOTE: not the ts WITHIN def(simulate), or mts. new variable for the tree sequence containing mutations. The function returns mts, which you've renamed back to ts for simplicity (it's the tree sequence to be used in the rest of the file)
# print(ts) #there's no sites on the original ts

newtables = mutation_placer()
# print(newtables)
# newtables.sort() #https://tskit.dev/tutorials/tables_and_editing.html
# tskit.TableCollection.sort(self=newtables)

new_ts = newtables.tree_sequence()
print(new_ts)

print("writing treefile for downstream analysis")
new_ts.dump("output.trees")




# print("making vcf for the sim") #hopefully not necessary 
# vcf_path = "genos.vcf"
n_dip_indv = int(ts.num_samples / 2) 
# indv_names = [f"tsk_{str(i)}indv" for i in range(n_dip_indv)]
# make_vcf(vcf_path, indv_names)



fid=[f"tsk_{str(i)}indv" for i in range(0,(deme_size*populations))]
iid=[f"tsk_{str(i)}indv" for i in range(0,(deme_size*populations))] #number of individuals in the sample
print("simulating phenotypes from environmental and genetic effects")


constructor = make_phenotypes()
# print("summed phenotypes", constructor.phenotypes)
print("environmental means:", constructor.environmental_means)

# print("mutation effect sizes",  constructor.muts)

# print(constructor.phenotypes)
# print("making .txt file that contains individuals and phenotypes")
# make_popfile(constructor.phenotypes)

print("dumping treefile with metadata")
add_metadata_to_treefile(newtables)



ts.draw_svg("treesequence_visualization.svg")
print("tree sequence data:","\n", new_ts)
#view trees in-terminal 
# tskit.TableCollection.sort(edge_start=0, self=newtables, site_start=0, mutation_start=0)
# ts = newtables.tree_sequence()
# for t in new_ts.trees():
#   print(t.draw(format='unicode'))


