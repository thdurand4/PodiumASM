## Arg : fasta input file, new name of the new file

from Bio import SeqIO
import sys

ancien_names = []
new_names = []
dict = {}
data = SeqIO.parse(sys.argv[1], "fasta")

for seq in data:
    dict[len(str(seq.seq))] = seq
dict = (sorted(dict.items(), key = lambda t:t[0], reverse=True))

list_records = []
for k in range(len(dict)):
    list_records.append(dict[k][1])

x=[]
for k in range(len(dict)):
    x.append("contig_"+str(k+1))

for k in range(len(list_records)):
    list_records[k].id = x[k]
    list_records[k].description = ""

SeqIO.write(list_records, sys.argv[2], 'fasta')
