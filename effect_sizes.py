
#load treesequence to get 
import tskit 
ts = tskit.load("output.trees")
#print(ts)


#effect sizes sampling file: https://github.com/Arslan-Zaidi/popstructure/blob/master/code/shared_scripts/simgeffects.R

# for ts in ts.trees():
#     for m in ts.mutations():
#         for sample in ts.samples():
#             for individual in ts.nodes():

#print(ts.trees)
#print(ts.mutations)
#print(ts.samples)
#print(ts.nodes)


#load in frequency file
freq_file = "convert_2_pgen_test.snps.frq.afreq"

# import numpy as np
# np.loadtxt(freq_file)

import pandas as pd
dataframe = pd.read_csv(freq_file, sep = "\t")
# print(dataframe)

list = []
for i in ts.individuals():
    list.append(i.id)
I = len(list)
print(I)

# print(ts.tables.mutations)
mutlist = []
for m in ts.mutations():
    mutlist.append(m.id)
M = len(mutlist)
print(M)



#making phenotypic values




P = 2 #when specifying phenotype in the future, it should be easy to take in the number of phenotypes invoked if using argparser. for now, hardcoding it
print(P)

import numpy as np




#make the empty array of numpy objects 
#will result in an array that is I lists containing P values in each one
array = np.zeros(I*P).reshape(I,P)
# print(array)

# I = 2
# P = 5
# array2 = np.zeros(I*P) #this just makes one list of length I*P
# print(array2)
# array2 = np.zeros(I*P).reshape(I,P)
# print(array2)



#βk∼N(0,σl^2[pk(1−pk)]^α) #their effect size distribution
#they calculate sigma^2 aka overall genetic variance/heritability at 0.8, might as well do the same

#directly transcribing into python the R and comments from https://github.com/Arslan-Zaidi/popstructure/blob/master/code/shared_scripts/simgeffects.R
sigma^2_l = 0.8 / sum( sapply( causal.variants$ALT_FREQS,function(x){
  beta= ( 2*x*(1-x)) ^ (1-0.4)
  return(beta)
}))

sigma^2_l = 0.8 / sum(map())
}))