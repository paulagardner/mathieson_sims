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
import matplotlib.pyplot as plt #if it's just matplotlib, your demesdraw subplots part won't work
import seaborn as sns
from dataclasses import dataclass
import random

demography = msprime.Demography()
graph = demes.load("pop_split.yaml")
# graph = demes.load("no_ancestry.yaml")
graph2 = demes.load("no_ancestry.yaml")

mu = ((1e-7)*20) #normally would have this set to 1e-7 *2, but putting deme size down and #individuals up to be able to visualize tree sequences for working on placing mutations issue 
rho = 1e-7
bases = 1000


deme_size = 5 #scaling this to 1000 doesn't affect much 
def simulate(mu, rho, graph):
    ts = msprime.sim_ancestry(demography = msprime.Demography.from_demes(graph), recombination_rate=rho, sequence_length = bases, samples={"A": deme_size, "B": deme_size}, random_seed=1234, discrete_genome = False) #add random seed so it's replicable
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
        nodes_list = np.zeros(tree.tree_sequence.num_nodes)
        for node in tree.nodes(): #type error: tree object is not iterable-- this happens w/ you frequently b/c you're forgetting to include the parenthesis to call the method properly, so the code tries to iterate over the method and not the data structure it's referencing 
            
            #add 
            dict[node] = tree.branch_length(node) #https://stackoverflow.com/questions/43505921/python-trying-to-create-a-dictionary-through-a-for-loop
            
            branch_lengths[node] =+ tree.branch_length(node) 
            nodes_list[node] =+ node
            #kevin's suggestion of how to do it- note that you've largely figured it out already in the line above 
        
            
            #get branch lengths as a proportion of total
            proportional_branch_lengths = branch_lengths/tree.total_branch_length


        #randomly sample nodes according to the branch lengths above them 
        #to do that, I'll need to have k (the number of draws numpy.choice/choices is making) specified by a probability distribution that takes in mutation rate
        #mu is from above
        # print("size of distribution", mu*bases)
        #conceptually: the probability of a mutation happening, look at coalescent theory book
        nodes_sampled =  np.random.choice(nodes_list, p=proportional_branch_lengths, size = 3) #this can sample the same node twice, which IS what you want. however need to figure out a way to determine which parent node was associated with a given branch length, so you can select a site that makes sense given the edges possible for a child:parent node relationship 
        
        nodes_sampled = nodes_sampled.astype(int) ##these need to be integers to be able to use them to access indices/values in the edge table. See if you can make it such when the  variable is created 
        # print("nodes sampled from distribution", nodes_sampled)

        #pick out the spot in the genome that the tree we have iterated over actually occupies (with edges) 

        print("interval this tree occupies", tree.interval)
        interval = tree.interval #access this to pick which sites are possible
        # leftint= tree.interval(left)
        # rightint= tree.interval(right)
        # print(ts.tables.edges)
        print(interval)
        # print("left side of interval", interval.left)
        for node in nodes_sampled:
            ######draw a site from the genome that HASN'T ALREADY had a mutation assigned to it. (can probably skip making the sites infinite for now)
            #use these nodes to write mutations to the mutations table
            # print(ts.tables.mutations)
            
            # site =np.random.choice(np.setdiff1d(range(0, bases), site_table.position)) #resource for what I'm wanting to do to simulate infinite sites: sample bases without replacement, but ACROSS loops (bases has a replacement parameter, that only works within calls of the function) #https://stackoverflow.com/questions/44507803/in-numpy-how-can-i-randomly-choose-a-number-from-a-range-that-excludes-a-random #need to test to see if it's working properly. range(1, bases is so you're actually getting values up to the integer value you've set as bases (print bases above to check)). 
            # # print("node:", node, "site chosen from bases:", site)
            # #add the site to the sites table

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
            #get the entries from the edges tables that correspond to the nodes you've chosen 

            #get child node time
            child_time = node_table.time[node]
            #get parent node time
            #make the span of that another list you can randomly select from 

            # print("node ID", node,  ts.tables.edges.child[node])#https://stackoverflow.com/questions/176918/finding-the-index-of-an-item-in-a-
            
            for index, j in enumerate(edge_table.child): #I think the child node is the appropriate one to pick, since we chose nodes BELOW branches. Mutations are associated with the node below them (according to tskit data model page)
                # print(index, j) #this print statement coupled with the one below lets you check if it's behaving as intended
                if j == node:
                    entry = edge_table[index]
                    print("node", node, "edges associated with the node selected", entry)
                    current_tree_interval = (interval.left, interval.right)

                    if (entry.left, entry.right) == current_tree_interval: 
                        parent_node = entry.parent
                        print("tuple matched, parent node", parent_node)

                    elif (entry.left, entry.right) == (0, bases): #this happens when the node is not associated with a recombination event 
                        print("this node has the same parent regardless of tree"
                        )
                        parent_node=entry.parent

            #now, find actual parent time:
            parent_time = node_table.time[parent_node]
            print("parent time", parent_time)
            time_span = (child_time, parent_time)
            print("time interval between a node and its parent", time_span)
            time=random.choice(time_span) #you get a suspiciously high amount of '0' times. 

            #####assign all variables to the mutations table 
            mut_table.add_row(site=site, node=node, time=time, derived_state = "1") #this doesn't work as with the above-- perhaps something with the strings. 
            # print(tskit.unpack_strings(site_table.ancestral_state, site_table.ancestral_state_offset))

            # print(np.setdiff1d(nucleotides, site_table.ancestral_state))

            
            print("LOOP BREAK", "\n")
  

        # print(ts.tables.edges)
        # print(node_table)
        # print(mut_table)
        # print(site_table)
        # print(edge_table)
    # print(site_table)
    # print(mut_table)
    print(newtables)
    return(newtables)
                
            #check where that site should go in terms of edges (eg which tree it'll end up falling on)


           

            #conceptually: a way to add rows
            #mut_table.add_row(site=(randomly select a site from the span of this tree's 'genome without replacement), node=node, derived_state=(randomly select what the derived state is), parent=NULL, metadata=None, time=(select a time interval on the branch length)) 
            # site = #randomly sample between left and right edges 

            #may also need to figure out (from theory) if branch lengths should be affecting mutation effect sizes 
            




         
        #write that to the mutations table to actually have the mutations
        
        # print("normalized branch lengths", proportional_branch_lengths)
        
        # print(branch_lengths)


    # for t in ts.trees():
    #     print(t)
    #     print(type(t))
    #     # t.branch_length(x) 


#make a vcf that will have each row signifying a mutation. (thus it's keeping track of the same info as make_phenotypes does: accounting for which individuals have which mutations, and how many copies)
# def make_vcf(vcf_path, indv_names):
#     with open(vcf_path, "w") as vcf_file:  
#         ts.write_vcf(vcf_file, individual_names=indv_names) 

#     #isolate each part of the vcf so you can collect all the elements and write a vcf from it, with ID fields
#     #get column headers
#     with open(vcf_path, 'r') as f:
#         lines = [l for l in f if not l.startswith('##')] 
        
#         df = pd.read_csv(io.StringIO(''.join(lines)),
#         dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
#             'QUAL': str, 'FILTER': str, 'INFO': str},
#             sep='\t')#.rename(columns={'#CHROM': 'CHROM'}) #do not want to rename this! Plink2 is very unhappy if you take the # away. 
    
    
#     #isolate header info rows
#     with open(vcf_path, 'r') as vcf_file:
#         header = [l for l in vcf_file if l.startswith('##')]
    
#     header =''.join(header)

#     # print(header)

#     #make an id for all mutations
#     id  = ["rs" + str(i) for i in range(len(ts.tables.mutations))] #calling it rs _ because i see it that way on other VCF

#     #replace the vcf ID column with new IDs
#     df['ID'] = id
#     # print(df)

#     blankIndex=[''] * len(df)
#     df.index=blankIndex
#     # print(df)

#     with open(vcf_path, 'w') as vcf: #this is almost correct but the headers don't work
#         vcf.write(header)#use this instead of df.to_csv because that will sometimes mess up your column headers
#     df.to_csv(vcf_path, sep="\t", mode='a', index=False)#if this is indented in with_open() sometimes mess up your column headers



    


        # for ts.nodes in tree:
        #     print(nodes)
    # for tree in ts.trees():
    #     print(ts.nodes)
        # print(tree.nodes)
        # nodes = nodes.id
        # # for node in tree:
        # #     print(node)
        # print(nodes)



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


# def make_popfile(phenotypes_array):

#     deme_id=[[i]*deme_size for i in range(0,demes)] #https://github.com/Arslan-Zaidi/popstructure/blob/master/code/simulating_genotypes/grid/generate_genos_grid.py
#     #flatten
#     deme_id=[item for sublist in deme_id for item in sublist] #changes 2 arrays of, say, length 50 into one array of length 100 (for example, will vary depending on deme # and sample sizes)). Necessary to make the array the correct size for the below 

#     # phenotypes_array = constructor.phenotypes
#     txt_name = "pop.txt"

#     pheno_1, pheno_2 = phenotypes_array.T #https://stackoverflow.com/questions/30820962/splitting-columns-of-a-numpy-array-easily
#     popdf = np.transpose([fid, iid, deme_id, pheno_1, pheno_2]) 
#     # print(popdf)
#     with open(txt_name, 'wb') as f:
#         f.write(b'FID\tIID\tdeme_id\tphenotype1\tphenotype2\n') 
#         np.savetxt(f, popdf, fmt = '%s') #why %s? https://stackoverflow.com/questions/48230230/typeerror-mismatch-between-array-dtype-object-and-format-specifier-18e

#     return popdf
#     return phenotypes_array


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
ax.figure.savefig("A.svg")

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
# print(ts.tables)

print("writing treefile for downstream analysis")
ts.dump("output.trees")


newtables = mutation_placer()

# print("making vcf for the sim") #hopefully not necessary 
# vcf_path = "genos.vcf"
n_dip_indv = int(ts.num_samples / 2) 
# indv_names = [f"tsk_{str(i)}indv" for i in range(n_dip_indv)]
# make_vcf(vcf_path, indv_names)


demes = 2 # replace hard-coded with argparser or taking the 'populations' value from the tskit table. in this case, you just want the demes of the samples present at the end of the sim. (I think). SO regardless of pop_split.yaml or no_ancestry.yaml (whether ancestry is shared), demes A and B split off.

fid=[f"tsk_{str(i)}indv" for i in range(0,(deme_size*demes))]
iid=[f"tsk_{str(i)}indv" for i in range(0,(deme_size*demes))] #number of individuals in the sample
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
print(ts)
#view tree sequences in-terminal 
# tskit.TableCollection.sort(edge_start=0, self=newtables, site_start=0, mutation_start=0)
# ts = newtables.tree_sequence()
# for t in ts.trees():
#   print(t.draw(format='unicode'))


