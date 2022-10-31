import msprime
import tskit 

ts = tskit.load("output.trees")
print(ts.node(6))

print(ts.samples())