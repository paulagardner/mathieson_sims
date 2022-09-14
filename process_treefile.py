import tskit
import pandas as pd
import io
import numpy as np

ts = tskit.load("output.trees")

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




def reconstruct_phenotypes():
    #get effect sizes for each mutation from metadata: 
    muteffects = []
    for index, m in enumerate(ts.mutations()):
        # print(index,m.metadata)
        for value in m.metadata.values():
            muteffects.append(value)
            phenos = len(value) #get number of phenotypes in each list so you can pass it down to the loop that needs it to calculate an array size 

    print(muteffects)
    print(phenos)





    #make phenotypes
    indv_array = ts.tables.nodes.individual
            # print(indv_array)

    nodes_array = []
    for node, individual in enumerate(ts.tables.nodes.individual):
        nodes_array.append(node)
    # print(nodes_array)

    array = np.column_stack((indv_array, nodes_array))
    # print(array)
    # phenos = 
    #create an empty numpy array to put the effect sizes in. It will return a table with # of phenotypes column and # of individuals rows
    phenotypes_array = np.zeros(n_dip_indv*phenos).reshape(n_dip_indv,phenos) 

    # print("empty phenotypes array:", phenotypes_array, "\n")



    #make an array to track whether an individual is homozygous or heterozygous for a mutation at a locus
    for tree in ts.trees():
    #     print("a")
    #     #loop over all mutations that appear in a tree
        for mutation in tree.mutations():
            nodes = mutation.node
            #make a 1-D numpy array with as many zeros as there are phenotypes
            ncopies = np.zeros(len(phenotypes_array)) #one zero for every individual in the sample at the end of the sim  
            #loop over all samples that have a given mutation. If you don't include nodes (mutation.node) here, you loop over ALL samples that contain ANY mutation, so it will be the same full list of samples, iterated over m # of mutations times. Including nodes does it by the sample.
            for sample_node in tree.samples(mutation.node):
                #find which individuals have 0, 1, or 2 mutations
                if sample_node in array[:,1:]:
                    item = (array[sample_node]) #now I have the individual and node # together for each one that has a mutation
                    individual_with_mut = item[:1]
                    # print(*individual_with_mut) #entechin.com/how-to-print-a-list-without-square-brackets-in-python/#:~:text=use%20asterisk%20'*'%20operator%20to%20print%20a%20list%20without%20square%20brackets
                    ncopies[individual_with_mut] += 1
                    
                    # print("phenotypic value of the indiv preexisting from environmental effect:",phenotypes_array[individual_with_mut])
                    # print("mutation value", muteffects[mutation.id])
                    phenotypes_array[individual_with_mut] += muteffects[mutation.id]
    
    return phenotypes_array








def make_popfile(phenotypes_array):

    deme_id=[[i]*deme_size for i in range(0,demes)] #https://github.com/Arslan-Zaidi/popstructure/blob/master/code/simulating_genotypes/grid/generate_genos_grid.py
    #flatten
    deme_id=[item for sublist in deme_id for item in sublist] #changes 2 arrays of, say, length 50 into one array of length 100 (for example, will vary depending on deme # and sample sizes)). Necessary to make the array the correct size for the below 

    
    
    # phenotypes_array = phenotypes_array
    txt_name = "pop.txt"
    pheno_1, pheno_2 = phenotypes_array.T #https://stackoverflow.com/questions/30820962/splitting-columns-of-a-numpy-array-easily
    # print(phenotypes_array)

    # pheno_1 = phenotypes_array[:,0]
    # pheno_2 = phenotypes_array[:,1]

    # print(pheno_1, pheno_2)
    # popdf = np.hstack([fid, iid, deme_id, pheno_1, pheno_2])
    # print(popdf)
    print(len(fid), len(iid), len(deme_id), len(pheno_1), len(pheno_2))
    popdf = np.vstack([fid, iid, deme_id, pheno_1, pheno_2]) 
    print(popdf)
    
    # popdf = np.vstack([fid, iid, deme_id, pheno_1, pheno_2])
    # with open(txt_name, 'wb') as f:
    #     f.write(b'FID\tIID\tdeme_id\tphenotype1\tphenotype2\n') 
    #     np.savetxt(f, popdf, fmt = '%s') #why %s? https://stackoverflow.com/questions/48230230/typeerror-mismatch-between-array-dtype-object-and-format-specifier-18e


# def make_popfile(phenotypes_array):

#     deme_id=[[i]*deme_size for i in range(0,demes)] #https://github.com/Arslan-Zaidi/popstructure/blob/master/code/simulating_genotypes/grid/generate_genos_grid.py
#     #flatten
#     deme_id=[item for sublist in deme_id for item in sublist] #changes 2 arrays of, say, length 50 into one array of length 100 (for example, will vary depending on deme # and sample sizes)). Necessary to make the array the correct size for the below 


#     # phenotypes_array = constructor.phenotypes
#     txt_name = "pop.txt"

#     pheno_1, pheno_2 = phenotypes_array.T #https://stackoverflow.com/questions/30820962/splitting-columns-of-a-numpy-array-easily
#     popdf = np.transpose([fid, iid, deme_id, pheno_1, pheno_2]) 
#     with open(txt_name, 'wb') as f:
#         f.write(b'FID\tIID\tdeme_id\tphenotype1\tphenotype2\n') 
#         np.savetxt(f, popdf, fmt = '%s') #why %s? https://stackoverflow.com/questions/48230230/typeerror-mismatch-between-array-dtype-object-and-format-specifier-18e

#     return popdf
#     # return phenoty

    






# print(phenotypes_array)



# for m in ts.mutations():
#     print(m.metadata)

# print(ts.tables.mutations.metadata)
# print(ts.table_metadata_schemas.mutations)
# print("Metadata for individual 0:", ts.individual(0).metadata)
# print("Metadata for individual 0:", ts.tables.individuals[0].metadata)  # Table access

# for m in ts.tables.mutations.metadata: #this prints the array used to construct the metadata, so not what you want
#     print(m)













print(ts.tables)
print(ts.num_samples)


print("making vcf for the sim") 
vcf_path = "genos.vcf"
n_dip_indv = int(ts.num_samples / 2) 
indv_names = [f"tsk_{str(i)}indv" for i in range(n_dip_indv)]
make_vcf(vcf_path, indv_names)

# demes = len(ts.tables.populations) #this will mess up your popfile
# print("NUMBER OF DEMES",demes)


print("making .txt file that contains individuals and phenotypes")
deme_size  = 25 #harcoding first to check if it works
demes = 2 ###########need to figure out how to not have these hardcoded
fid=[f"tsk_{str(i)}indv" for i in range(0,(deme_size*demes))]
iid=[f"tsk_{str(i)}indv" for i in range(0,(deme_size*demes))] #number of individuals in the sample
phenotypes_array = reconstruct_phenotypes()
print(phenotypes_array)


print("simulating phenotypes from environmental and genetic effects")
make_popfile(phenotypes_array)
