import sys
from collections import defaultdict
from Bio import SeqIO
import os
import pandas as pd

idtocontig = defaultdict(str)
contigtoN = defaultdict(int)
contigto_percent = defaultdict(int)
contigs = []
taille = []
n_percent = []
GC_contig = defaultdict(str)
GC = []
depth = []
start_telo = []
end_telo = []

for record in SeqIO.parse(sys.argv[1], "fasta"):
    idtocontig[record.id] = record.seq



for cle in idtocontig:
    contigtoN[cle]=idtocontig[cle].count("N")
    contigto_percent[cle]=(contigtoN[cle]/len(idtocontig[cle])*100)
    contigs.append(cle)
    taille.append(len(idtocontig[cle]))
    n_percent.append(contigtoN[cle]/len(idtocontig[cle])*100)

for sequences in SeqIO.parse(sys.argv[2], "fasta"):
    GC_contig[sequences.id] = sequences.seq

for cle in GC_contig:
    G = GC_contig[cle].count("G")
    C = GC_contig[cle].count("C")
    GC.append((G+C)/len(GC_contig[cle])*100)

with open(sys.argv[3],"r") as f1:
    for lignes in f1:
        test = lignes.split("\t")
        depth.append(test[6])
depth.pop(0)


with open(sys.argv[4],"r") as f1:
    for lignes in f1:
        test = lignes.split("\t")
        start_telo.append(test[4])
        end_telo.append(test[5])

start_telo.pop(0)
end_telo.pop(0)

df = pd.DataFrame({"Chromosome":contigs, "Taille_en_bp": taille, "Pourcentage_de_N": n_percent,"GC":GC, "Meand_Depth":depth, "Start_Telomere": start_telo, "End_Telomere": end_telo})

fasta = sys.argv[1]
strain = fasta.split(".")

df.to_csv(sys.argv[5],index=False)

#print(df)

