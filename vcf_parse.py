import sys
import allel
import pandas as pd
import os

deletions = 0
insertion = 0
inversion = 0
translocations = 0
duplications = 0

df = allel.vcf_to_dataframe(sys.argv[1])
list_var = []

ID = list((df.ID))
for element in ID:
    variant = element.split(".")
    if variant[1] == 'DEL':
        deletions += 1
    if variant[1] == 'INS':
        insertion += 1
    if variant[1] == 'INV':
        inversion += 1
    if variant[1] == 'BND':
        translocations += 1
    if variant[1] == 'DUP':
        duplications += 1

list_var = [deletions, insertion, inversion, duplications, translocations]
index_data = ['deletions', 'insertions', 'inversions', 'duplications', 'translocations']
data = pd.DataFrame(list_var, index = index_data)
data.to_csv(f"{sys.argv[2]}")




#df = df[df.CHROM == "VMNF01000006.1"]
#print(df.head())

#df.to_csv("UK_092_vcf_VMNF06.csv")
