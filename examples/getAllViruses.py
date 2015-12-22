#There is a lineage paramater that may be helpful
from Bio import Entrez, SeqIO
import os #TEMP
os.chdir('..')
import pyviko.core as core

def extractJoin(location):
	return ((int(location[location.index('[')+1:location.index(':')]),int(location[location.index(':')+1:location.index(']')])),(int(location[location.rindex('[')+1:location.rindex(':')]),int(location[location.rindex(':')+1:location.rindex(']')])))
	
def findGenes(geneList, positive):
	over = []
	toKO = []
	if positive:
		for gene in geneList:
			for otherGene in geneList:
				if gene <> otherGene:
					if gene[0] > otherGene[0] and gene[0] < otherGene[1]:
						over.append(otherGene)
						toKO.append(gene)
	else:
		for gene in geneList:
			for otherGene in geneList:
				if gene <> otherGene:
					if gene[0] < otherGene[0] and gene[0] > otherGene[1]:
						over.append(otherGene)
						toKO.append(gene)
	return over, toKO

Entrez.email = ""
#handle = Entrez.esearch(db="nuccore", term='"complete genome"[All Fields] AND viruses[filter] NOT segment[All Fields] NOT partial[All Fields]', retmax=48770)
handle = Entrez.esearch(db="nuccore", term='"complete genome"[All Fields] AND viruses[filter] NOT segment[All Fields] NOT partial[All Fields]', retmax=31)
record = Entrez.read(handle)
idlist = record["IdList"]
handle.close()
addedCounter = 0
handle = Entrez.efetch(db="nuccore", id=idlist, rettype="gb", retmode="text")
results = SeqIO.parse(handle, "gb")
for seq_record in results:
	plus_genes = []
	minus_genes = []
	ccc = 0
	for i in seq_record.features:
		if i.type == "CDS":
			ccc+=1
			if "join" in str(i.location):
				fancyJoin = extractJoin(str(i.location))
				if i.location.strand == 1:
					plus_genes += [(fancyJoin[0][0], fancyJoin[1][1], seq_record.seq[fancyJoin[0][0]:]+seq_record.seq[:fancyJoin[1][1]])]
				elif i.location.strand == -1:
					minus_genes += [(fancyJoin[1][1], fancyJoin[0][0], core.reverseComplement(seq_record.seq[fancyJoin[0][0]]+seq_record.seq[fancyJoin[1][1]]))]
			else:
				if i.location.strand == 1:
					plus_genes += [(int(i.location.start), int(i.location.end), str(seq_record.seq[i.location.start:i.location.end]))]
				elif i.location.strand == -1:
					minus_genes += [(int(i.location.end), int(i.location.start), str(core.reverseComplement(seq_record.seq[i.location.start:i.location.end])))]

	overPlus, toKoPlus = findGenes(plus_genes, True)
	overMinus, toKoMinus = findGenes(minus_genes,False)
	over = overPlus + overMinus
	toKO = toKoPlus + toKoMinus
	#print seq_record.id, len(over), len(toKO), ccc
handle.close()

'''records = Entrez.read(handle)
for i in range(0,10):
	print records[1]['GBSeq_feature-table'][i]
	print
handle.close()
'''