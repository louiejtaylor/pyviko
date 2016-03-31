'''
Script to extract overprinted gene pairs from NCBI Nucleotide database.
'''
from Bio import Entrez, SeqIO
import os
os.chdir('..')
import pyviko.core as core
import time
	
def findGenes(geneList, positive):
	'''
	Extracts gene pairs from Entrez records.
	'''
	over = []
	toKO = []
	goodGenes = []
	for gene in geneList:
		trigger = False
		for notNt in ['R', 'Y', 'W', 'S', 'M', 'K', 'B', 'D', 'H', 'V', 'N', 'F']:
			if notNt in gene[2]:
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

def extractJoin(location):
	'''
	Extracts the locations of genes and their 
	overprinted counterparts from Entrez record.
	'''
	return ((int(location[location.index('[')+1:location.index(':')]),int(location[location.index(':')+1:location.index(']')])),(int(location[location.rindex('[')+1:location.rindex(':')]),int(location[location.rindex(':')+1:location.rindex(']')])))

#setup
Entrez.email = "your_email_here@university.edu"
#query
handle = Entrez.esearch(db="nuccore", term='"complete genome"[All Fields] AND viruses[filter] NOT segment[All Fields] NOT partial[All Fields]', retmax=48770)
record = Entrez.read(handle)
idlist = record["IdList"]
handle.close()
#setting up files to store collected genes
addedCounter = 0
finum = 0
fiKO = open('test/dem/ko/'+str(finum)+'.fasta', 'w')
fiOver = open('test/dem/over/'+str(finum)+'.fasta', 'w')
for jjj in range(1,int(len(idlist)/100)+2):
	print 
	print "Page: " + str(jjj) + " of " +str(int(len(idlist)/100)+1)
	time.sleep(1)
	#grab genes
	handle = Entrez.efetch(db="nuccore", id=idlist[jjj*100-100:jjj*100], rettype="gb", retmode="text")
	results = SeqIO.parse(handle, "gb")
	for seq_record in results:
		plus_genes = []
		minus_genes = []
		for i in seq_record.features:
			try:
				#find coding sequences
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
				#rare errors, but just in case
				print "Error with: "
				print i
		overPlus, toKoPlus = findGenes(plus_genes, True)
		overMinus, toKoMinus = findGenes(minus_genes,False)
		#store total genes
		over = overPlus + overMinus
		toKO = toKoPlus + toKoMinus
		
		if len(over) <> len(toKO):
			#if gene pair lists are not of even size (for whatever reason), we don't want to store this round
			print "Error: file lengths uneven"
			print len(over), len(toKO), addedCounter
		else:
			#add genes to files
			for i in range(len(over)):
				addedCounter += 1
				if addedCounter % 100 == 0:
					#generate multiple files so that we don't have a single giant FASTA
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
