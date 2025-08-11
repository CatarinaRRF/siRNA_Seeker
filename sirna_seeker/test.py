import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction


sequencia = Seq("CCGCCACCATGACCACCTCCATC")
comp = sequencia.complement()
trasc = Bio.Seq.transcribe(sequencia)
cg = round(gc_fraction(sequencia)*100, 2)
#print(len(trasc))
print(trasc)
print(cg)
#print(comp)