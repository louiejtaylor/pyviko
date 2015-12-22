#There is a lineage paramater that may be helpful
from Bio import Entrez, SeqIO
import os #TEMP
os.chdir('..')
import pyviko.core as core
import time

def extractJoin(location):
	return ((int(location[location.index('[')+1:location.index(':')]),int(location[location.index(':')+1:location.index(']')])),(int(location[location.rindex('[')+1:location.rindex(':')]),int(location[location.rindex(':')+1:location.rindex(']')])))
	
def findGenes(geneList, positive):
	over = []
	toKO = []
	goodGenes = []
	for gene in geneList:
		trigger = False
		for notNt in ['R', 'Y', 'W', 'S', 'M', 'K', 'B', 'D', 'H', 'V', 'N', 'F']:
			if notNt in gene[2]:
				print "Bad nt caught: " + notNt				
				trigger = True
		if trigger == False:
			goodGenes.append(gene)
	geneList = [g for g in goodGenes]
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
handle = Entrez.esearch(db="nuccore", term='"complete genome"[All Fields] AND viruses[filter] NOT segment[All Fields] NOT partial[All Fields]', retmax=48770)
#handle = Entrez.esearch(db="nuccore", term='"complete genome"[All Fields] AND viruses[filter] NOT segment[All Fields] NOT partial[All Fields]', retmax=131)
record = Entrez.read(handle)
idlist = record["IdList"]
handle.close()
addedCounter = 0
finum = 0
fiKO = open('test/dem/ko/'+str(finum)+'.fasta', 'w')
fiOver = open('test/dem/over/'+str(finum)+'.fasta', 'w')
for jjj in range(1,int(len(idlist)/100)+2):
	print 
	print "Page: " + str(jjj) + " of " +str(int(len(idlist)/100)+1)
	time.sleep(2)
	handle = Entrez.efetch(db="nuccore", id=idlist[jjj*100-100:jjj*100], rettype="gb", retmode="text")
	results = SeqIO.parse(handle, "gb")
	for seq_record in results:
		plus_genes = []
		minus_genes = []
		for i in seq_record.features:
			try:
				if i.type == "CDS":
					if "join" in str(i.location):
						fancyJoin = extractJoin(str(i.location))
						if i.location.strand == 1:
							plus_genes += [(fancyJoin[0][0], fancyJoin[1][1], str(seq_record.seq[fancyJoin[0][0]:]+seq_record.seq[:fancyJoin[1][1]]))]
						elif i.location.strand == -1:
							minus_genes += [(fancyJoin[1][1], fancyJoin[0][0], core.reverseComplement(seq_record.seq[fancyJoin[0][0]]+seq_record.seq[fancyJoin[1][1]]))]
					else:
						if i.location.strand == 1:
							plus_genes += [(int(i.location.start), int(i.location.end), str(seq_record.seq[i.location.start:i.location.end]))]
						elif i.location.strand == -1:
							minus_genes += [(int(i.location.end), int(i.location.start), str(core.reverseComplement(seq_record.seq[i.location.start:i.location.end])))]
			except:
				print "Error with: "
				print i
		overPlus, toKoPlus = findGenes(plus_genes, True)
		overMinus, toKoMinus = findGenes(minus_genes,False)
		over = overPlus + overMinus
		toKO = toKoPlus + toKoMinus
		
		if len(over) <> len(toKO):
			print "ERROR: LENGTHS NOT EVEN"
			print len(over), len(toKO), addedCounter
		else:
			for i in range(len(over)):
				addedCounter += 1
				if addedCounter % 100 == 0:
					fiKO.close()
					fiOver.close()
					finum += 1
					fiKO = open('test/dem/ko/'+str(finum)+'.fasta', 'w')
					fiOver = open('test/dem/over/'+str(finum)+'.fasta', 'w')
				fiKO.write('>'+str(seq_record.id)+' '+str(toKO[i][0])+':'+str(toKO[i][1]) +'\n'+toKO[i][2]+'\n')
				fiOver.write('>'+str(seq_record.id)+' '+str(over[i][0])+':'+str(over[i][1])+'\n'+over[i][2]+'\n')
	print "Added: " + str(addedCounter)
	handle.close()
	
fiKO.close()
fiOver.close()
'''
			try:
				print core.findOverlap(toKO[i][2],over[i][2])
			except core.SequenceError:
				print "Ohno", toKO[i], over[i]
			except IndexError:
				print "we dont care"
'''

#print seq_record.id, len(over), len(toKO), ccc
'''
try:
	print core.findOverlap(toKO[0][-1],over[0][-1])
except core.SequenceError:
	print "Ohno"
except IndexError:
	next
'''

'''records = Entrez.read(handle)
for i in range(0,10):
	print records[1]['GBSeq_feature-table'][i]
	print
handle.close()
'''