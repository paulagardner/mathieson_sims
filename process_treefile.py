import tskit
import pandas as pd
import io
import numpy as np

 

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

    # print(header)

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




# def reconstruct_phenotypes(): #this should be rendered obsolete by actually putting the phenotypes in the metadata
#     #get effect sizes for each mutation from metadata: 
#     muteffects = []
#     for index, m in enumerate(ts.mutations()):
#         # print(index,m.metadata)
#         for value in m.metadata.values():
#             muteffects.append(value)
#             phenos = len(value) #get number of phenotypes in each list so you can pass it down to the loop that needs it to calculate an array size 

#     # print(muteffects)
#     # print(phenos)





#     #make phenotypes
#     indv_array = ts.tables.nodes.individual
#             # print(indv_array)

#     nodes_array = []
#     for node, individual in enumerate(ts.tables.nodes.individual):
#         nodes_array.append(node)
#     # print(nodes_array)

#     array = np.column_stack((indv_array, nodes_array))
#     # print(array)
#     # phenos = 
#     #create an empty numpy array to put the effect sizes in. It will return a table with # of phenotypes column and # of individuals rows
#     phenotypes_array = np.zeros(n_dip_indv*phenos).reshape(n_dip_indv,phenos) 

#     # print("empty phenotypes array:", phenotypes_array, "\n")



    # #make an array to track whether an individual is homozygous or heterozygous for a mutation at a locus
    # for tree in ts.trees():
    # #     print("a")
    # #     #loop over all mutations that appear in a tree
    #     for mutation in tree.mutations():
    #         nodes = mutation.node
    #         #make a 1-D numpy array with as many zeros as there are phenotypes
    #         ncopies = np.zeros(len(phenotypes_array)) #one zero for every individual in the sample at the end of the sim  
    #         #loop over all samples that have a given mutation. If you don't include nodes (mutation.node) here, you loop over ALL samples that contain ANY mutation, so it will be the same full list of samples, iterated over m # of mutations times. Including nodes does it by the sample.
    #         for sample_node in tree.samples(mutation.node):
    #             #find which individuals have 0, 1, or 2 mutations
    #             if sample_node in array[:,1:]:
    #                 item = (array[sample_node]) #now I have the individual and node # together for each one that has a mutation
    #                 individual_with_mut = item[:1]
    #                 # print(*individual_with_mut) #entechin.com/how-to-print-a-list-without-square-brackets-in-python/#:~:text=use%20asterisk%20'*'%20operator%20to%20print%20a%20list%20without%20square%20brackets
    #                 ncopies[individual_with_mut] += 1
                    
    #                 # print("phenotypic value of the indiv preexisting from environmental effect:",phenotypes_array[individual_with_mut])
    #                 # print("mutation value", muteffects[mutation.id])
    #                 phenotypes_array[individual_with_mut] += muteffects[mutation.id]

    
#     return phenotypes_array

class phenotypes:
    def __init__(self): #I think by writing the phenotypes into the sim, you're making the phenotype reconstructor obsolete????????
        muteffects = []
        for index, m in enumerate(ts.mutations()):
            # print(index,m.metadata)
            for value in m.metadata.values():
                muteffects.append(value)
                phenos = len(value) #get number of phenotypes in each list so you can pass it down to the loop that needs it to calculate an array size 



        phenotypes_array = np.zeros(n_dip_indv*phenos).reshape(n_dip_indv,phenos) 
        print("empty phenotypes array",phenotypes_array)
        # for tree in ts.trees():
        #     for individual in tree.individuals():
        #         print(metadata)

        # print("table", ts.tables.individuals)
        # print(ts.tables.individuals.metadata)
        # for i in ts.tables.individuals:
        #     print(i.metadata)
        environmental_effects = np.zeros(n_dip_indv*phenos).reshape(n_dip_indv,phenos) 
        genetic_effects = np.zeros(n_dip_indv*phenos).reshape(n_dip_indv,phenos)

        for index, i in enumerate(ts.tables.individuals):
            # print(i.metadata)
            # print("loop")
            # print((i.metadata)[value] for value in i.metadata)
            # print("loop")
            # print(*(i.metadata).values()) #you run into issues if you don't include the parenthesis around i.metadata. HOWEVER doing it like this only works with the print statement 
            # print((i.metadata)['full_phenotypes']) #get the values for the key you want    
            phenotypes_array[index] += (i.metadata)['full_phenotypes']
            environmental_effects[index] += (i.metadata)['environmental_effects']
            genetic_effects[index]+= (i.metadata)['genetic_effects']
        
        # print(phenotypes_array)
        # print(environmental_effects)

        
        # for index, i in enumerate(ts.individuals()):
        #     print(i.metadata)
            # for i.metadata():
            #     print('a')

        # array = []
        # for i in ts.individuals():
            # array.append(i.metadata)
            # print(i.metadata)
            # print(i.metadata)
            # for key in i.metadata:
            #     print(key)
            # for key in ts.individuals[i]:
            #     print(key)
        # print(array)
        # print(ts.individuals.metadata)


        self.full_phenotypes = phenotypes_array
        self.environmental_effects = environmental_effects
        self.genetic_effects = genetic_effects


def get_phenotypes():
    return phenotypes()





# def make_popfile(phenotypes_array): #not sure why this works, I would've thought I would have to call self.full_phenotypes here???
def make_popfile(full_phenotypes,environmental_effects,genetic_effects): #not sure why this works, I would've thought I would have to call self.full_phenotypes here???
    deme_id=[[i]*deme_size for i in range(0,populations_at_end)] #https://github.com/Arslan-Zaidi/popstructure/blob/master/code/simulating_genotypes/grid/generate_genos_grid.py
    #flatten
    deme_id=[item for sublist in deme_id for item in sublist] #changes 2 arrays of, say, length 50 into one array of length 100 (for example, will vary depending on deme # and sample sizes)). Necessary to make the array the correct size for the below 

    
    
    # phenotypes_array = phenotypes_array
    txt_name = "pop.txt"

    
    # full_pheno_1, full_pheno_2 = phenotypes_array.T #https://stackoverflow.com/questions/30820962/splitting-columns-of-a-numpy-array-easily
    # popdf = np.transpose([fid, iid, deme_id, full_pheno_1, full_pheno_2]) 
    
    full_pheno_1, full_pheno_2 = full_phenotypes.T #https://stackoverflow.com/questions/30820962/splitting-columns-of-a-numpy-array-easily
    environmental_effect_1, environmental_effect_2 = environmental_effects.T
    genetic_effect_1, genetic_effect_2 = genetic_effects.T 

    popdf = np.transpose([fid, iid, deme_id, full_pheno_1, full_pheno_2, environmental_effect_1, environmental_effect_2, genetic_effect_1, genetic_effect_2]) 
    # print(popdf)

    



    with open(txt_name, 'wb') as f:
        # f.write(b'FID\tIID\tdeme_id\tfull_phenotype1\tfull_phenotype2\n') 
        f.write(b'FID\tIID\tdeme_id\tfull_phenotype1\tfull_phenotype2\tenvironmental_effect1\tenvironmental_effect2\tgenetic_effect1\tgenetic_effect2\n') 
        np.savetxt(f, popdf, fmt = '%s') #why %s? https://stackoverflow.com/questions/48230230/typeerror-mismatch-between-array-dtype-object-and-format-specifier-18e


# def make_datafiles(phenotypes_array, environmental_effects):
#     deme_id=[[i]*deme_size for i in range(0,populations_at_end)
#     deme_id=[item for sublist in deme_id for item in sublist

ts = tskit.load("output.trees") # how did this get deleted???


# print(ts.num_samples)
print("ts.num_samples:", ts.num_samples)


print("making vcf for the sim") 
vcf_path = "genos.vcf"
n_dip_indv = int(ts.num_samples / 2) 
indv_names = [f"tsk_{str(i)}indv" for i in range(n_dip_indv)]
make_vcf(vcf_path, indv_names)

# demes = len(ts.tables.populations) #this will mess up your popfile
# print("NUMBER OF DEMES",demes)


print("making .txt file that contains individuals and phenotypes")
#find how many demes are present at the end of the sim
print(type(ts.tables.nodes))
populations = ts.tables.nodes.population
times = ts.tables.nodes.time
array = np.column_stack((populations, times))
present_nodes = array[array[:,1] == 0] #sort to get only those nodes present at the end of the sim
populations_at_end = len(set(present_nodes[:,0])) #length of the SET of the populations column in filtered set. This gives you number of demes/populations at end, regardless of ancestral nodes

#make a list of individuals in each deme. The demes must be equal in size for this to work
# individuals_at_end = len(present_nodes)
# print(individuals_at_end)
# individuals_at_end = (len(present_nodes)/2)
individuals_at_end = len(present_nodes) #this isn't actually the # of individuals, since this is actually the length of nodes.... but does it mess up the code?
print("individuals_at_end:",individuals_at_end)
deme_size = int(individuals_at_end / (2*populations_at_end)) #divide by 2 because nodes represent haploid genomes in diploid organisms. 


#use deme sizes and number of demes to create IDs that plink will understand 
fid=[f"tsk_{str(i)}indv" for i in range(0,(deme_size*populations_at_end))]
iid=[f"tsk_{str(i)}indv" for i in range(0,(deme_size*populations_at_end))] #number of individuals in the sample

constructor = get_phenotypes() #will have to rewrite this if you're returning multiple values

# phenotypes_array = reconstruct_phenotypes() #obsolete????

print("simulating phenotypes from environmental and genetic effects")
# make_popfile(constructor.full_phenotypes)
make_popfile(constructor.full_phenotypes, constructor.environmental_effects, constructor.genetic_effects)

# print(ts.tables.individuals.metadata)#these will print as the non-decoded versions
# print(ts.tables.mutations.metadata)
# for index, m in enumerate(ts.mutations()):
#     print(index,m.metadata)
# for index, i in enumerate(ts.individuals()):
#     print(i.metadata)