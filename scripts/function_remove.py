from Bio import SeqIO
import pandas as pd
import sys

dict = {}

for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
    seq = str(seq_record.seq)
    dict[seq_record.id] = (seq.count("N") / len(seq) * 100)

L1 = list((dict.keys()))
L2 = list(dict.values())
data = pd.DataFrame(L2, index=L1)
data.columns = ["Npercent"]
data = data[data.Npercent > float(sys.argv[4])] # Fait un dataframe avec les listes et enlève les contigs avec %N<70%
L = list(data.index) ## L est la liste des contigs qu'on enlève

dictseq = {}

sequences = SeqIO.parse(sys.argv[2], "fasta")

for i in sequences:
    dictseq[i.id] = i  ## dictionnaire dictseq avec en clé les id et en valeurs les seq records
Ldel = []  # Ldel même chose que L

for j in dictseq.keys():
    if j in L:
        Ldel.append(j)
for count in range(len(L)):
    del dictseq[L[count]]

SeqIO.write(dictseq.values(), sys.argv[3], 'fasta')