import msprime

demography = msprime.Demography()
mu = 1e-6
rho = 1e-3

import demes 
graph = demes.load("no_ancestry.yaml")


def simulate(mu, rho, graph):
    ts = msprime.sim_ancestry(demography = msprime.Demography.from_demes(graph), recombination_rate=rho, sequence_length =1000, samples={"A": 50, "B": 50}, random_seed=1234, discrete_genome = False) #add random seed so it's replicable
    #not including migration info since I'm just doing two demes, so what they did (defining which demes are next to each other) not necessary here.. I think! at least for now

    mts = msprime.sim_mutations(ts, rate = mu, discrete_genome = False) #discrete_genome = false gives infinite sites, basically, so simplifies downstream bc you don't have to account for mult mutations at a site
    return mts




print("simulating genotypes under demographic model")

ts = simulate(mu, rho, graph)





n_dip_indv = int(ts.num_samples / 2)
indv_names = [f"tsk_{str(i)}indv" for i in range(n_dip_indv)]

vcf_path = "output_geno.vcf"

with open(vcf_path, "w") as vcf_file: #normal formatting would be with open("output.vcf", "w") as vcf_file:, here gzipping for when the analysis gets big. ALSO, kevin suggested "wb" as the flag, but tskit docs say "wt" and it didn't give errors like wb did: "TypeError: memoryview: a bytes-like object is required, not 'str'"
    ts.write_vcf(vcf_file, individual_names=indv_names) 

    import io
    import os
    import pandas as pd


    def read_vcf(vcf_path):
        with open(vcf_path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')] #removes the header columns. you will need to figure out how to 
        return pd.read_csv(
            io.StringIO(''.join(lines)),
            dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                'QUAL': str, 'FILTER': str, 'INFO': str},
            sep='\t'
        ).rename(columns={'#CHROM': 'CHROM'})


    def get_header(vcf_path): #kind of a silly way to do it when read_vcf could potentially be the whole function, returning both the df and header. ask skylar about it-- should be similar logic
        #to demography parser, sorting through the objects that are returned. This matters because the id definition is counting on one object being returned
        with open(vcf_path, 'r') as f:
            header = [l for l in f if l.startswith('##')]
        return header

    df = read_vcf(vcf_path)
    header = get_header(vcf_path)
    print(header)
    # print(df)

    id  = ["rs" + str(i) for i,_ in enumerate(df.ID)]
    print(id)
    # id_df = pd.DataFrame(id)

    # print(len(df))
    # print(len(id_df))

    df['ID'] = id
    print(df)
    df.to_csv('a.vcf', sep = '\t')


    # n = df.columns[4]
    # df.drop(n, axis = 1, inplace = True)
    # df[4] = id
    # print(df)
    # df['ID'] = id_df[0]
    # print(df)
    

    # df.to_csv('a.vcf', sep = '\t')

    # with open(vcf_path, "w") as vcf_file:
    #     df.to_csv(vcf_path, sep = '\t')
    # print(vcf_file)

        













# def read_vcf(vcf_file):
#     with open("output_geno.vcf.gz", "w") as vcf_file:
#         lines = [l for l in vcf_file if not vcf_file.startswith('##')]
#     return pd.read_csv("output_geno.vcf.gz",
#         dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
#             'QUAL': str, 'FILTER': str, 'INFO': str},
#         sep='\t'
#         ).rename(columns={'#CHROM': 'CHROM'})

# print(vcf_file)



