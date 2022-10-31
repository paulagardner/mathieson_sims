
# %%
%pip install seaborn #had issues where seaborn wasn't working. This is from: https://stackoverflow.com/questions/65280038/modulenotfounderror-no-module-named-seaborn-in-jupyter-notebook. Not sure if this will permanently fix it in this env or what  

# import seaborn as sns 
# import pandas as pd 
# #trying to plot in seaborn


# phenotype1 = pd.read_csv('phenotype1.sscore', sep='\t')
# print(phenotype1)
# phenotype2 = pd.read_csv('phenotype2.sscore', sep='\t')


# sns.set_theme()
# g = sns.relplot(data=phenotype1, x = "phenotype1", y = "phenotype2")


# %%
import seaborn as sns 
import pandas as pd 
#trying to plot in seaborn


phenotype1 = pd.read_csv('phenotype1.sscore', sep='\t')
phenotype2 = pd.read_csv('phenotype2.sscore', sep='\t')


sns.set_theme()
pheno_comp = sns.relplot(data=phenotype1, x = "phenotype1", y = "phenotype2")


# %%
