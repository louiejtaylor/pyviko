# -*- coding: utf-8 -*-

'''
Proof-of-concept module for pyViKO. [Examples and source](https://github.com/louiejtaylor/pyViKO)
'''

from pyviko.core import codonify, seqify, insert_mutation, translate, stop_codons
from pyviko.restriction import find_n_cutters

# workhorse: find stop codons that can be generated given codons, mutation
def find_possible_stop_codons(codons, n):
	
	if codons[-1] in stop_codons:
		codons = codons[:-1] # remove c-terminal stop codon
	
	almost_stop_codons = {}
	# build dict of codons that can be mutated to a stop codon 
	
	for c in stop_codons:
		for i in range(0,3):
			for nt in 'ACTG':
				if c[:i]+nt+c[i+1:] not in stop_codons:
					try:
						almost_stop_codons[c[:i]+nt+c[i+1:]].append(c)
					except KeyError:
						almost_stop_codons[c[:i]+nt+c[i+1:]] = [c]
	if n == 2:
		for c in almost_stop_codons.keys():
			for i in range(0,3):
				for nt in 'ACTG':
					if c[:i]+nt+c[i+1:] not in stop_codons:
						try:
							almost_stop_codons[c[:i]+nt+c[i+1:]] += almost_stop_codons[c]
						except KeyError:
							 almostStopCodons[c[:i]+nt+c[i+1:]] = almost_stop_codons[c]
		for i in almost_stop_codons.keys():
			almost_stop_codons[i] = list(set(almost_stop_codons[i]))
	
	
	# creates a list of tuples of the form (index, ['list', 'of', 'stop', 'codons']) to be further pre-processed	
	pre_matches = [(i, almost_stop_codons[codons[i]]) for i in range(0,len(codons)) if codons[i] in almost_stop_codons.keys()]
	
	matches = []  
	
	# further processing to create actual tuples (index, 'codon')	
	for m in pre_matches:
		for codon in m[1]:
			matches.append((m[0],codon))
			
	print(len(matches))
	
	return matches

# finds overprinted gene, given input sequence, frameshift, and bool startsBefore 
# TODO: change startsBefore to an **index** for gene	
def find_overprinted_gene(seq, frame, starts_before):	
	codons = codonify(seq[frame - 1:])[:-1] # remove last (incomplete) codon
	for i in range(0,len(codons)):
		if codons[i] in stop_codons:
			codons = codons[:i]
			break
	
	if not starts_before:
		x = 1
		# TODO: complete whatever this is
		
	return codons

def find_non_harmful_mutations(seq, frame, starts_before, num_mutations):
	codons = codonify(seq)
	stops = find_possible_stop_codons(codons, num_mutations)
	overAAs = translate(find_overprinted_gene(seq, frame, starts_before))
	winners = []
	for poss in stops:
		n_codons = [codon for codon in codons]
		new_codons = insert_mutation(n_codons, poss)
		new_over_AAs = translate(find_overprinted_gene(seqify(new_codons), frame, starts_before))
		if new_over_AAs == over_AAs:
			winners.append(poss)
	return winners

def find_restriction_site_changes(seq, frame, starts_before, num_mutations):
	safe_mutations = find_non_harmful_mutations(seq, frame, starts_before, num_mutations)
	r_sites = find_n_cutters(seq,6)
	new_sites = []
	for i in safe_mutations:
		new_sites.append((i, find_n_cutters(seqify(insert_mutation(codonify(seq),i)), 6)))
	winners = []
	print(len(new_sites))
	for j in new_sites:
		if j[1] <> r_sites:
			res = [r for r in r_sites]
			for site in j[1]:		  
				try:
					res.remove(site)
				except ValueError:
					res.append((site,'+++'))
			k = (j[0], res)
			winners.append(k)
	for w in winners:
		print(w)
	return winners
	
sequence = '''ATGGAACAGGCACCAGAAGATCAAGGACCACAGAGGGAGCCATACAACGAATGGGCTTTAGAATTGTTGGAAGACCTAAAGAATGAGGCTCTGCGCCACTTTCCTCGGCCTTGGCTACATGGACTAGGGCAATACTTCTATAATACATATGGAGATACCTGGGAGGGAGTAGAGGCCATCATTAGGACACTACAACAACTGTTGTTTATACATTATAGGATTGGCTGTCAACATAGCAGGATAGGAATCACTCCTCAAAGGAGAAGGAATGGAGCCAGTAGATCCTGA'''

if __name__ == "__main__":
	mut_file = open('test/mutations.fasta',"w")
	mut_file.write("> Original sequence\n")
	mut_file.write(sequence+"\n")
	for f in find_restriction_site_changes(sequence, 2, True, 1):
		mut_file.write('> Mutant at codon ' + str(f[0][0]+1) +' '+ str(f[1])+'\n')
		#print f #for debugging, prints output to screen
		mut_file.write(seqify(insert_mutation(codonify(sequence),f[0]))+'\n')
	mut_file.close()
