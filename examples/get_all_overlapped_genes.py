'''
Script to extract overprinted gene pairs from NCBI Nucleotide database.
'''
from Bio import Entrez, SeqIO
import os
import pyviko.core as core
import time
	
def find_genes(gene_list, positive):
	'''
	Extracts gene pairs from Entrez records.
	'''
	over = []
	to_ko = []
	good_genes = []
	for gene in gene_list:
		trigger = False
		for not_nt in ['R', 'Y', 'W', 'S', 'M', 'K', 'B', 'D', 'H', 'V', 'N', 'F']:
			if not_nt in gene[2]:
				trigger = True
		if trigger == False:
			good_genes.append(gene)
	gene_list = [g for g in good_genes]
	if positive:
		for gene in gene_list:
			for other_gene in gene_list:
				if gene <> other_gene:
					if gene[0] > other_gene[0] and gene[0] < other_gene[1]:
						over.append(other_gene)
						toKO.append(gene)
	else:
		for gene in gene_list:
			for other_gene in gene_list:
				if gene <> other_gene:
					if gene[0] < other_gene[0] and gene[0] > other_gene[1]:
						over.append(other_gene)
						to_ko.append(gene)
	return over, to_ko

def extract_join(location):
	'''
	Extracts the locations of genes and their 
	overprinted counterparts from Entrez record.
	'''
	return ((int(location[location.index('[')+1:location.index(':')]),int(location[location.index(':')+1:location.index(']')])),(int(location[location.rindex('[')+1:location.rindex(':')]),int(location[location.rindex(':')+1:location.rindex(']')])))

# setup
Entrez.email = "your_email_here@university.edu"
# query
handle = Entrez.esearch(db="nuccore", term='"complete genome"[All Fields] AND viruses[filter] NOT segment[All Fields] NOT partial[All Fields]', retmax=48770)
record = Entrez.read(handle)
idlist = record["IdList"]
handle.close()
# setting up files to store collected genes
addedCounter = 0
finum = 0
fiKO = open('test/dem/ko/'+str(finum)+'.fasta', 'w')
fiOver = open('test/dem/over/'+str(finum)+'.fasta', 'w')
for jjj in range(1,int(len(idlist)/100)+2):
	print("Page: " + str(jjj) + " of " +str(int(len(idlist)/100)+1))
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
						fancy_join = extract_join(str(i.location))
						if i.location.strand == 1:
							plus_genes += [(fancy_join[0][0], fancy_join[1][1], str(seq_record.seq[fancy_join[0][0]:]+seq_record.seq[:fancy_join[1][1]]))]
						elif i.location.strand == -1:
							minus_genes += [(fancy_join[1][1], fancy_join[0][0], core.reverse_complement(seq_record.seq[fancy_join[0][0]]+seq_record.seq[fancy_join[1][1]]))]
					else:
						if i.location.strand == 1:
							plus_genes += [(int(i.location.start), int(i.location.end), str(seq_record.seq[i.location.start:i.location.end]))]
						elif i.location.strand == -1:
							minus_genes += [(int(i.location.end), int(i.location.start), str(core.reverse_complement(seq_record.seq[i.location.start:i.location.end])))]
			except:
				#rare errors, but just in case
				print("Error with:", i)
		over_plus, to_ko_plus = find_genes(plus_genes, True)
		over_minus, to_ko_minus = find_genes(minus_genes,False)
		#store total genes
		over = over_plus + over_minus
		to_ko = to_ko_plus + to_ko_minus
		
		if len(over) <> len(to_ko):
			# if gene pair lists are not of even size (for whatever reason), we don't want to store this round
			print("Error: file lengths uneven")
			print(len(over), len(to_ko), added_counter)
		else:
			#add genes to files
			for i in range(len(over)):
				added_counter += 1
				if added_counter % 100 == 0:
					#generate multiple files so that we don't have a single giant FASTA
					fi_ko.close()
					fi_over.close()
					finum += 1
					fi_ko = open('test/dem/ko/'+str(finum)+'.fasta', 'w') #better system for this
					fi_over = open('test/dem/over/'+str(finum)+'.fasta', 'w')
				fi_ko.write('>'+str(seq_record.id)+' '+str(to_ko[i][0])+':'+str(to_ko[i][1]) +'\n'+to_ko[i][2]+'\n')
				fi_over.write('>'+str(seq_record.id)+' '+str(over[i][0])+':'+str(over[i][1])+'\n'+over[i][2]+'\n')
	print("Added: " + str(added_counter))
	handle.close()
	
fi_ko.close()
fi_over.close()
