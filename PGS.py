import tskit
import numpy as np
import pandas as pd
#get treefile from directory
ts = tskit.load("output.trees")

 
# using loadtxt()
# arr = np.loadtxt("gwas.full_phenotype1.glm.linear",
                #  delimiter="\t", dtype=str)
# print(arr)
arr = np.loadtxt("gwas.full_phenotype1.glm.linear",
                 delimiter="\t", dtype=str)
print(arr)









#get a column (the betas) from a 3d numpy array 
# https://stackoverflow.com/questions/40742046/extract-certain-columns-of-3d-numpy-array
# print("betas from the numpy array from our GWAS", arr[:, 8:9])

# betas = arr[:, (2,8)] #version where you're pulling out ID column as well 
#need to manually change the array value to a string so that numpy can add the effect sizes in your loop below, since above you've read things in as strings (since a bunch of them are!)
betas = arr[:, 8]
betas = betas.flatten()
betas = betas.astype(float) #https://stackoverflow.com/questions/3877209/how-to-convert-an-array-of-strings-to-an-array-of-floats-in-numpy
# print("column of beta values for the given phenotype",betas)


#do the loop from model.py where you iterate over the tree to add effect sizes
        #make an array to track whether an individual is homozygous or heterozygous for a mutation at a locus
        # print("looping over trees to find which individuals have which mutations at which nodes")


#############sample of how muteffects was constructed: presumably you can just feed in a numpy array of betas from the GWAS 
# muteffects = rng.multivariate_normal(mean = means, cov = var_covar, size = (len(mut_index)))

# for tree in ts.trees(): #update_samples = True errors out with TypeError
#     #loop over all mutations that appear in a tree
#     for mutation in tree.mutations():
#         # print("newline every time we loop over a new mutation", "\n")
#         # print("mutation id", mutation.id)
#         nodes = mutation.node
#         # print("node the mutation originates from:", nodes)
#         #alternate way to get the above:
#         # print(mutation.node)

#         #make a 1-D numpy array with as many zeros as there are phenotypes
#         ncopies = np.zeros(len(genetic_effects_array)) #one zero for every individual in the sample at the end of the sim

#         #loop over all samples that have a given mutation. If you don't include nodes (mutation.node) here, you loop over ALL samples that contain ANY mutation, so it will be the same full list of samples, iterated over m # of mutations times. Including nodes does it by the sample.
#         for sample_node in tree.samples(mutation.node):
#             # print("this mutation is in sample node",sample_node)

#             #find which individuals have 0, 1, or 2 mutations
#             if sample_node in indvs_and_nodes_array[:,1:]:
#                 item = (indvs_and_nodes_array[sample_node]) #now I have the individual and node # together for each one that has a mutation
#                 individual_with_mut = item[:1]
#                 # print(*individual_with_mut) #entechin.com/how-to-print-a-list-without-square-brackets-in-python/#:~:text=use%20asterisk%20'*'%20operator%20to%20print%20a%20list%20without%20square%20brackets
#                 ncopies[individual_with_mut] += 1
                
#                 # print("phenotypic value of the indiv preexisting from environmental effect:",phenotypes_array[individual_with_mut])
#                 # print("mutation value", muteffects[mutation.id])
#                 genetic_effects_array[individual_with_mut] += muteffects[mutation.id]
#         print("copies of the mut present in each individual:", ncopies)
        


indv_array = ts.tables.nodes.individual
nodes_array = []
for node, individual in enumerate(ts.tables.nodes.individual):
    nodes_array.append(node)
        # print(nodes_array)



#get true mutation effect sizes so you can calculate residuals. Do this before everything else so you can access the number of phenotypes from the mutation
muteffects = []
for index, m in enumerate(ts.mutations()):
    # print(index,m.metadata)
    for value in m.metadata.values():
        muteffects.append(value)
        # phenos = len(value) #get number
#convert muteffects array so we can plot them (get the first phenotype value only for the test case)
muteffects = np.array(muteffects)
# phenos = muteffects.shape[1]
# print(phenos)
true_muteffects = muteffects[:,0] 
# print("true mutation effects", true_muteffects)
# print("estimated mutation effects", betas)


phenos = 2 #DO NOT HARDCODE THIS LATER
n_dip_indv = int(ts.num_samples / 2) 


# # true_genetic_effects_array = np.zeros((phenos*n_dip_indv,phenos)) #Multiple phenotypes version
# # true_genetic_effects_array = np.zeros((phenos,n_dip_indv)) 
# true_genetic_effects_array = np.zeros((n_dip_indv*phenos)).reshape(n_dip_indv*phenos)
# print("true genetic effects array:", true_genetic_effects_array)
# print(len(true_genetic_effects_array))
# # true_genetic_effects_array = np.zeros(n_dip_indv).reshape(n_dip_indv)
# # true_genetic_effects_array = np.zeros(n_dip_indv*phenos).reshape(n_dip_indv*phenos)
# pgs_array = np.zeros(n_dip_indv*phenos).reshape(n_dip_indv*phenos)
# # pgs_array = np.zeros((phenos,n_dip_indv))
# indvs_and_nodes_array = np.column_stack((indv_array, nodes_array))
# # print(indvs_and_nodes_array)
# # print(true_genetic_effects_array)
# # print(pgs_array)




true_genetic_effects_array = np.zeros((n_dip_indv))
# print("true genetic effects array:", true_genetic_effects_array)
# print(len(true_genetic_effects_array))
pgs_array = np.zeros(n_dip_indv)
indvs_and_nodes_array = np.column_stack((indv_array, nodes_array))


##this is my attempt to be having the two phenotypes, but it's not important for committee meeting, where I"ll be using a symmetric correlation matrix-- so the GWAS on phenotype1 should be more or less the same as the GWAS you'd get on the second phenotype
# true_genetic_effects_array = np.zeros((n_dip_indv,phenos)) #Multiple phenotypes version
# print(true_genetic_effects_array)
# print(len(true_genetic_effects_array))
# pgs_array = np.zeros((n_dip_indv, phenos))
# indvs_and_nodes_array = np.column_stack((indv_array, nodes_array))




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
        ncopies = np.zeros(len(indv_array)) #one zero for every individual in the sample at the end of the sim

        #loop over all samples that have a given mutation. If you don't include nodes (mutation.node) here, you loop over ALL samples that contain ANY mutation, so it will be the same full list of samples, iterated over m # of mutations times. Including nodes does it by the sample.
        for sample_node in tree.samples(mutation.node):
            # print("this mutation is in sample node",sample_node)

            #find which individuals have 0, 1, or 2 mutations
            if sample_node in indvs_and_nodes_array[:,1:]: #looking through the second column of the 3D matrix
                item = (indvs_and_nodes_array[sample_node]) #now I have the individual and node # together for each one that has a mutation
                
                individual_with_mut = item[:1]
                # print(individual_with_mut)
                # print(*individual_with_mut) #entechin.com/how-to-print-a-list-without-square-brackets-in-python/#:~:text=use%20asterisk%20'*'%20operator%20to%20print%20a%20list%20without%20square%20brackets
                
                ncopies[individual_with_mut] += 1
                # print(ncopies[individual_with_mut])
                
                # print("phenotypic value of the indiv preexisting from environmental effect:",phenotypes_array[individual_with_mut])
                # print("mutation value", muteffects[mutation.id])

                # print(genetic_effects_array[individual_with_mut])
                # genetic_effects_array[individual_with_mut] += betas[]
                # for idx, x in enumerate(betas):
                #     # print(idx)
                #     genetic_effects_array[individual_with_mut] += betas[idx]
                
                # true_genetic_effects_array[individual_with_mut] += muteffects[mutation.id]


                #get genetic values calculated off of the predicted effect sizes (the betas)
                pgs_array[individual_with_mut] += betas[mutation.id] #this is just going to be your polygenic score as the mathieson paper's done it, because they calculate the sums of (estimated effect size)*(# of DERIVED alleles. However, since you're centering your phetnoype at 0, any allele with an effect is a derived one, so adding up mutations per individual is just making a polygenic score)
                # pgs_array[individual_with_mut,0] += betas[mutation.id]
                true_genetic_effects_array[individual_with_mut] += true_muteffects[mutation.id] #this wasn't working with the normal array of muteffects, for some reason. mysterious


                #version of this for when you're dealing w/ 2 phenotype-GWAS (for when it matters, when the matrix is not symmetrical)
                # pheno1_pgs_array[individual_with_mut,0] += betas[mutation.id] 
                # # pgs_array[individual_with_mut,0] += betas[mutation.id]
                # pheno1_true_genetic_effects_array[individual_with_mut,0] += true_muteffects[mutation.id] 




# print("polygenic scores calculated from genetic_effects_array")
# print(ts.tables.mutations)
# print(ts.tables.individuals) #compare this with your estimated_genetic_effects_array w/ chris tomorrow



import seaborn as sns
import matplotlib.pyplot as plt

# concatenated = pd.concat([betas.assign(dataset='betas'), phenotype1_muteffects.assign(dataset='phenotype1_muteffects')]) #https://stackoverflow.com/questions/51732867/seaborn-plot-two-data-sets-on-the-same-scatter-plot
# concatenated = pd.concat([betas, phenotype1_muteffects]) #https://stackoverflow.com/questions/64145095/how-to-plot-data-from-multiple-dataframes-with-seaborn-relplot


betas = pd.DataFrame(betas)
true_muteffects = pd.DataFrame(true_muteffects)
# df = pd.concat([pd.DataFrame[betas], pd.DataFrame[phenotype1_muteffects]], axis = 1) #doesn't work although the above does????
df = pd.concat([betas, true_muteffects], axis = 1, keys = ['betas', 'true_muteffects'])
# print(true_muteffects)
# print(betas)
print(df)
plot = sns.scatterplot(data=df)
plt.show()
plt.savefig('output.png')


# fig, ax = plt.subplots() #https://stackoverflow.com/questions/59243174/creating-a-seaborn-factor-plot-using-two-data-sets
# sns.scatterplot(data=betas, ax=ax) # first dataset
# sns.scatterplot(data=phenotype1_muteffects, ax=ax) # second dataset


# get a matrix of predicted genetic effects
# print("phenotype1 polygenic score array:", pgs_array)
# print("true genetic effects array:",true_genetic_effects_array)

# subtract genetic effects FROM polygenic scores 
residuals_array = (pgs_array) - (true_genetic_effects_array)
# print("array of residuals", residuals_array)

#make a .txt file with the residuals so you can plot them later 

# np.savetxt("no_correlation_no_ancestry.txt", residuals_array)
# np.savetxt("no_correlation_including_ancestry.txt", residuals_array)

# np.savetxt("half_correlation_no_ancestry.txt", residuals_array)
np.savetxt("half_correlation_including_ancestry.txt", residuals_array)

# np.savetxt("full_correlation_no_ancestry.txt", residuals_array)
# np.savetxt("full_correlation_including_ancestry.txt", residuals_array)







