import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

#using seaborn b/c R and I are not on speaking terms

# common_pca = open("common_PCA.eigenvec", "r")
# print(common_pca.readline())
common_pca = pd.read_csv("common_PCA.eigenvec", sep = '\t')


phenotype1_glm = pd.read_csv("gwas.phenotype1.glm.linear", sep = '\t')
beta1 = phenotype1_glm["BETA"]
abs_beta1 = abs(beta1)
phenotype1_glm['ABS BETA'] = abs_beta1

# sns.scatterplot(data= phenotype1_glm, x="P", y=abs_beta1)

phenotype2_glm = pd.read_csv("gwas.phenotype2.glm.linear", sep = '\t')
beta2 = phenotype1_glm["BETA"]
abs_beta2 = abs(beta2)
phenotype2_glm['ABS BETA'] = abs_beta2

#https://stackoverflow.com/questions/69352701/seaborn-implot-combine-datasets-into-one-plot

concatenated = pd.concat([phenotype1_glm.assign(dataset='phenotype1_glm'), phenotype2_glm.assign(dataset='phenotype2_glm')])
print(concatenated)
sns.scatterplot(x='P', y='ABS BETA', data=concatenated, hue="dataset")

# plt.show()
plt.savefig("beta_vs_pvalue.png")


