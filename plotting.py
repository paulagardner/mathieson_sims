import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

#using seaborn b/c R and I are not on speaking terms

# common_pca = open("common_PCA.eigenvec", "r")
# print(common_pca.readline())

# common_pca = pd.read_csv("common_PCA.eigenvec", sep = '\t')


# phenotype1_glm = pd.read_csv("gwas.phenotype1.glm.linear", sep = '\t')
# beta1 = phenotype1_glm["BETA"]
# abs_beta1 = abs(beta1)
# phenotype1_glm['ABS BETA'] = abs_beta1
# # popfile = np.loadtxt("pop.txt")
# # popfile = pd.read_csv("pop.txt")
# # popfile = pd.read_csv('pop.txt', sep=" ", header=None)
# # data.columns = ["FID", "IID", "deme_id", "phenotype_1", 'phenotype_2']
# # print(popfile)
# # popfile = pd.read_fwf('pop.txt')
# popfile = pd.read_csv('pop.txt', sep=" ", header=None, names=["FID", "IID", "deme_id", "phenotype_1", "phenotype_2"], skiprows = 1)
# # print(popfile)

# # common_pca = open("common_PCA.eigenvec", "r")
# # print(common_pca.readline())
# common_pca = pd.read_csv("common_PCA.eigenvec", sep = '\t')
# # print(common_pca)



# #GWAS plots
# phenotype1_glm = pd.read_csv("gwas.phenotype1.glm.linear", sep = '\t')
# beta1 = phenotype1_glm["BETA"]
# abs_beta1 = abs(beta1)
# phenotype1_glm['ABS BETA'] = abs_beta1

# # sns.scatterplot(data= phenotype1_glm, x="P", y=abs_beta1)

# phenotype2_glm = pd.read_csv("gwas.phenotype2.glm.linear", sep = '\t')
# beta2 = phenotype1_glm["BETA"]
# abs_beta2 = abs(beta2)
# phenotype2_glm['ABS BETA'] = abs_beta2

# #https://stackoverflow.com/questions/69352701/seaborn-implot-combine-datasets-into-one-plot

# concatenated = pd.concat([phenotype1_glm.assign(dataset='phenotype1_glm'), phenotype2_glm.assign(dataset='phenotype2_glm')])
# # print(concatenated)
# sns.scatterplot(x='P', y='ABS BETA', data=concatenated, hue="dataset")

# # plt.show()
# plt.savefig("beta_vs_pvalue.png")



#plot cases against each other for each demography
# no_correlation_no_ancestry = np.loadtxt("no_correlation_no_ancestry.txt", dtype=float) #probably good to make your script name the items whatever you have as the file name
# no_correlation_including_ancestry = np.loadtxt("no_correlation_including_ancestry.txt", dtype=float)
# print(no_correlation_no_ancestry)

# colnames = ['with_ancestry']
no_correlation_no_ancestry = pd.read_csv("no_correlation_no_ancestry.txt", names=['no correlation, no ancestry'] ,dtype=float) 
#probably good to make your script name the items whatever you have as the file name
print(no_correlation_no_ancestry)
no_correlation_including_ancestry = pd.read_csv("no_correlation_including_ancestry.txt", names=['no correlation, with ancestry'],dtype=float)
print(no_correlation_including_ancestry)



# no_correlation_no_ancestry_common_correction = pd.read_csv("no_correlation_no_ancestry_common_correction.txt", names=['with ancestry'],dtype=float)
# no_correlation_including_ancestry_common_correction = pd.read_csv("no_correlation_including_ancestry_common_correction.txt", names=['with ancestry'],dtype=float)
# print(no_correlation_including_ancestry_common_correction)

half_correlation_no_ancestry = pd.read_csv("half_correlation_no_ancestry.txt", names=['.5 correlation, no ancestry'],dtype=float)
half_correlation_including_ancestry = pd.read_csv("half_correlation_including_ancestry.txt", names=['.5 correlation, with ancestry'],dtype=float)

full_correlation_no_ancestry = pd.read_csv("full_correlation_no_ancestry.txt", names=['full correlation, no ancestry'],dtype=float)
full_correlation_including_ancestry = pd.read_csv("full_correlation_including_ancestry.txt", names=['full correlation, with ancestry'],dtype=float)






# concatenated = pd.concat([no_correlation_no_ancestry.assign(dataset='no_correlation_no_ancestry'), no_correlation_including_ancestry.assign(dataset='no_correlation_including_ancestry')])
# print(concatenated)
# sns.scatterplot(x='P', y='ABS BETA', data=concatenated, hue="dataset")

fig, ax = plt.subplots(1,1, sharey=True)
sns.scatterplot(data=no_correlation_no_ancestry, palette = ['blue'])
sns.scatterplot(data=no_correlation_including_ancestry, palette = ['teal'])
# sns.scatterplot(data=no_correlation_no_ancestry_common_correction, palette = ['orange'])
# sns.scatterplot(data=no_correlation_including_ancestry_common_correction, palette = ['red'])
sns.scatterplot(data=half_correlation_no_ancestry, palette = ['yellow'])
sns.scatterplot(data=half_correlation_including_ancestry, palette = ['orange'])
sns.scatterplot(data=full_correlation_no_ancestry, palette = ['red'])
sns.scatterplot(data=full_correlation_including_ancestry, palette = ['purple'])
ax.set_xlabel('Individual')

# plt.show()
plt.savefig("z.png")



fig, ax = plt.subplots(1,1, sharey=True)
sns.scatterplot(data=no_correlation_no_ancestry, palette = ['blue'])
sns.scatterplot(data=no_correlation_including_ancestry, palette = ['teal'])
# ax.set_xlabel('Individual')
plt.savefig("z1.png")


fig, ax = plt.subplots(1,1, sharey=True)
sns.scatterplot(data=half_correlation_no_ancestry, palette = ['yellow'])
sns.scatterplot(data=half_correlation_including_ancestry, palette = ['orange'])
# ax.set_xlabel('Individual')
plt.savefig("z2.png")

fig, ax = plt.subplots(1,1, sharey=True)
sns.scatterplot(data=full_correlation_no_ancestry, palette = ['red'])
sns.scatterplot(data=full_correlation_including_ancestry, palette = ['purple'])
# ax.set_xlabel('Individual')
plt.savefig("z3.png")

#problems with ipykernel, prob bc I've pip installed some things and conda installed others #https://github.com/microsoft/vscode-jupyter/wiki/Failure-to-start-Kernel-due-to-Modules-not-installed , but it still just isn't working
# #%%
# msg='hello world'
# print(msg)




