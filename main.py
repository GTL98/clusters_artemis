from Bio import SeqIO
from json import load
from obter_clusters import obter_clusters

with open('genes_bifido_pullorum.json', 'r') as doc:
    dados = load(doc)

registro = SeqIO.read('Bifidobacterium_pullorum.gb', 'genbank')
x = obter_clusters(registro, 'Beta-galactosidase', 'ESN35_10250', dados)
print(x)
