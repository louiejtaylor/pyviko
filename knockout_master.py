# -*- coding: utf-8 -*-

'''
Proof-of-concept module for pyViKO. [Examples and source](https://github.com/louiejtaylor/pyViKO)
'''

from pyviko.core import codonify, seqify, insertMutation, translate, stopCodons
from pyviko.restriction import findNcutters

# workhorse: find stop codons that can be generated given codons, mutation
def findPossibleStopCodons(codons, n):
	
	if codons[-1] in stopCodons:
		codons = codons[:-1] # remove c-terminal stop codon
	
	almostStopCodons = {}
	# build dict of codons that can be mutated to a stop codon 
	
	for c in stopCodons:
		for i in range(0,3):
			for nt in 'ACTG':
				if c[:i]+nt+c[i+1:] not in stopCodons:
					try:
						almostStopCodons[c[:i]+nt+c[i+1:]].append(c)
					except KeyError:
						almostStopCodons[c[:i]+nt+c[i+1:]] = [c]
	if n == 2:
		for c in almostStopCodons.keys():
			for i in range(0,3):
				for nt in 'ACTG':
					if c[:i]+nt+c[i+1:] not in stopCodons:
						try:
							almostStopCodons[c[:i]+nt+c[i+1:]] += almostStopCodons[c]
						except KeyError:
							 almostStopCodons[c[:i]+nt+c[i+1:]] = almostStopCodons[c]
		for i in almostStopCodons.keys():
			almostStopCodons[i] = list(set(almostStopCodons[i]))
	
	
	# creates a list of tuples of the form (index, ['list', 'of', 'stop', 'codons']) to be further pre-processed	
	preMatches = [(i, almostStopCodons[codons[i]]) for i in range(0,len(codons)) if codons[i] in almostStopCodons.keys()]
	
	matches = []  
	
	# further processing to create actual tuples (index, 'codon')	
	for m in preMatches:
		for codon in m[1]:
			matches.append((m[0],codon))
			
	print(len(matches))
	
	return matches

# finds overprinted gene, given input sequence, frameshift, and bool startsBefore 
# TODO: change startsBefore to an **index** for gene	
def findOverprintedGene(seq, frame, startsBefore):	
	codons = codonify(seq[frame - 1:])[:-1] # remove last (incomplete) codon
	for i in range(0,len(codons)):
		if codons[i] in stopCodons:
			codons = codons[:i]
			break
	
	if not startsBefore:
		x = 1
		# TODO: complete
		
	return codons

def findNonHarmfulMutations(seq, frame, startsBefore, numMutations):
	codons = codonify(seq)
	stops = findPossibleStopCodons(codons, numMutations)
	overAAs = translate(findOverprintedGene(seq, frame, startsBefore))
	winners = []
	for poss in stops:
		nCodons = [codon for codon in codons]
		newCodons = insertMutation(nCodons, poss)
		newOverAAs = translate(findOverprintedGene(seqify(newCodons), frame, startsBefore))
		if newOverAAs == overAAs:
			winners.append(poss)
	return winners

def findRestrictionSiteChanges(seq, frame, startsBefore, numMutations):
	safeMutations = findNonHarmfulMutations(seq, frame, startsBefore, numMutations)
	rSites = findNcutters(seq,6)
	newSites = []
	for i in safeMutations:
		newSites.append((i, findNcutters(seqify(insertMutation(codonify(seq),i)), 6)))
	winners = []
	print len(newSites)
	for j in newSites:
		if j[1] <> rSites:
			re = [r for r in rSites]
			for site in j[1]:		  
				try:
					re.remove(site)
				except ValueError:
					re.append((site,'+++'))
			k = (j[0], re)
			winners.append(k)
	for w in winners:
		print w
	return winners
	
sequence = '''ATGGAACAGGCACCAGAAGATCAAGGACCACAGAGGGAGCCATACAACGAATGGGCTTTAGAATTGTTGGAAGACCTAAAGAATGAGGCTCTGCGCCACTTTCCTCGGCCTTGGCTACATGGACTAGGGCAATACTTCTATAATACATATGGAGATACCTGGGAGGGAGTAGAGGCCATCATTAGGACACTACAACAACTGTTGTTTATACATTATAGGATTGGCTGTCAACATAGCAGGATAGGAATCACTCCTCAAAGGAGAAGGAATGGAGCCAGTAGATCCTGA'''

if __name__ == "__main__":
	mutFile = open('test/mutations.fasta',"w")
	mutFile.write("> Original sequence\n")
	mutFile.write(sequence+"\n")
	for f in findRestrictionSiteChanges(sequence, 2, True, 1):
		mutFile.write('> Mutant at codon ' + str(f[0][0]+1) +' '+ str(f[1])+'\n')
		#print f #for debugging, prints output to screen
		mutFile.write(seqify(insertMutation(codonify(sequence),f[0]))+'\n')
	mutFile.close()
