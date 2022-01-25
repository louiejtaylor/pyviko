from pyviko import core, restriction

class OverGene:
	'''
	Class representing the overprinted gene.
	'''
	frame = 1
	start_nucleotide_index = -1
	pre_sequence = '' # includes 1-2nt removed by core.find_over_gene
	gene_sequence = ''
	post_sequence = '' # includes 1-2nt removed by core.find_over_gene
	over_aas = ''
	
	def __init__(self, over_seq, start_nt_index, seq, frame_over = 1):
		if over_seq != '':
			ol = core.find_overlap(seq, over_seq)
			if ol[0] == 0: #overprinted gene starts before
				start_nt_index = -1
				frame_over = 4-(ol[1]%3)
				self.pre_sequence = over_seq[ol[1]-(-frame_over+4):ol[1]] + seq[:3-(-frame_over+4)]
			else: #overprinted gene starts after
				start_nt_index = ol[0]
		self.gene_sequence = over_seq
		self.start_nucleotide_index = start_nt_index
		self.comb_sequence = seq
		self.frame = frame_over
		self.over_aas = core.translate(self.pre_sequence + core.seqify(core.find_overprinted_gene(seq, start_nt_index, frame_over)) + self.post_sequence)

class Mutant:
	'''
	Class to mutate target gene.
	'''
	seq = ''
	codons = []
	n_mut = 1
	over_gene = False
	mutants = []
	
	def __init__(self, sequence, num_mutations = 1, regex = False):
		self.seq = sequence
		self.n_mut = num_mutations
		self.codons = core.codonify(sequence)
		self.regex = regex
		
	def set_over_gene(self, over_seq = '', start_nt_index = -1, over_frame = 1):
		'''
		Adds the overprinted gene to the current `Mutant` object.
		'''
		if over_seq == '' and start_nt_index == -1 and over_frame == 1:
			raise core.SequenceError("You must provide either the sequence of an overprinted gene, or its start position/frame in the knockout sequence.")
		else:
			self.over_gene = OverGene(over_seq, start_nt_index, self.seq, over_frame)
		
	def vector(self, sequence):
		'''
		Adds the vector sequence to the current `Mutant` object 
		(primarily for making primers of early knockouts).
		'''
		if self.seq in sequence:
			self.vector_seq = sequence
		else:
			raise core.SequenceError("Could not find target sequence in vector sequence.")
			
	def find_mutants(self, ignore_rx_sites = True, r_site_length = 6, r_sites = restriction.default_enzymes()):
		'''
		Returns a list of mutants that add a premature stop codon 
		(or mutate the start codon) without changing the overprinted 
		gene, and which add or remove a restriction site.
		'''

		stops = find_possible_stop_codons(self.codons, self.nMut)

		if self.over_gene:
			if len(self.over_gene.gene_sequence) > 0:
				stops = mutate_start_codon(self.codons, self.n_mut) + stops
			safe_mutations = []
			new_pre_sequence = '';
			for poss in stops:
				n_codons = [codon for codon in self.codons]
				new_codons = core.insert_mutation(nCodons, poss)
				if self.over_gene.gene_sequence != '':
					new_pre_sequence = self.over_gene.pre_sequence[:4-self.over_gene.frame] + new_codons[0][:self.overGene.frame - 1]
				new_over_aas = core.translate(new_pre_sequence + core.seqify(core.find_overprinted_gene(
                                    core.seqify(new_codons), self.over_gene.start_nucleotide_index, self.over_gene.frame)))
				if new_over_aas == self.overGene.over_aas:
					safe_mutations.append(poss)
		else:
			safe_mutations = stops
		final_winners = [s for s in safe_mutations]
		if not ignore_rx_sites:
			### Two approaches: regex and non-regex.
			restriction_site_lengths = list(set([len(k) for k in r_sites.keys()]))
			
			if r_site_length == 'all':
				temp_restriction_sites = r_sites
			elif r_site_length >= min(restriction_site_lengths) and r_site_length <= max(restriction_site_lengths):
				r_keys = [k for k in r_sites.keys() if len(k) == r_site_length]
				temp_restriction_sites = {} # Reduce size of dict. searched
				for site in r_keys:
					temp_restriction_sites[site] = r_sites[site]
				restriction_site_lengths = [r_site_length]
			else:
				raise core.SequenceError("Invalid restriction site length.")		
			
			new_sites = [] # list of lists

			# should do all one way--why is this inconsistent?
			### Regex:
			if self.regex:
				base_sites = restriction.re_find_enzymes(self.seq)
				
				for mut in safe_mutations:
					new_sites.append([])
					new_sites[-1] += restriction.re_find_enzymes(core.seqify(core.insert_mutation(self.codons, mut)))
					
			### Non-regex:
			else:
				
				base_sites = []
				for length in restriction_site_lengths:
					base_sites += restriction.find_n_cutters(self.seq, length)
					
				for mut in safe_mutations:
					new_sites.append([])
					for length in restriction_site_lengths:
						new_sites[-1] += restriction.find_n_cutters(core.seqify(core.insert_mutation(self.codons, mut)), length)

			winners = {}
			for l in new_sites:
				if l != base_sites: #this is why I should use sets
					temp_sites = [c for c in base_sites]
					temp_added_sites = []
					for site in l: #basically, removing everything in the new list from the old list to get the differences
						try: #set logic would improve this
							temp_sites.remove(site)
						except ValueError:
							temp_added_sites.append(site)

					for i in range(0,len(temp_sites)):
						temp_sites[i] = (temp_sites[i][0], temp_sites[i][1], '-')
						
					for i in range (0, len(temp_added_sites)):
						temp_added_sites[i] = (temp_added_sites[i][0], temp_added_sites[i][1], '+')
					
					winners[safe_mutations[new_sites.index(l)]] = temp_sites + temp_added_sites
			
			final_winners = [(x,winners[x]) for x in sorted(winners.keys(), key=lambda x: x[0])]		
		# Rx can be either ignored or not
		return final_winners

def find_stop_codon_mutants(codons, n):
        '''
        Given a list `codons`, finds individual codons that
        can be mutated to a stop codon given `n` mutations. Returns
        a list of tuples of the form `(index, 'codon')` where `'codon'`
        is the mutated codon.
        '''
        #TODO: can optimize this as in js
        if codons[-1] in core.stop_codons:
                codons = codons[:-1] #remove c-terminal stop codon

        almost_stop_codons = {}
        #build dict of codons that can be mutated to a stop codon

        for c in core.stop_codons:
                for i in list(range(0,3)):
                        for nt in 'ACTG':
                                if c[:i]+nt+c[i+1:] not in core.stop_codons:
                                        try:
                                                almost_stop_codons[c[:i]+nt+c[i+1:]].append(c)
                                        except KeyError:
                                                almost_stop_codons[c[:i]+nt+c[i+1:]] = [c]

        if n == 2:
                for c in almost_stop_codons.keys():
                        for i in list(range(0,3)):
                                for nt in 'ACTG':
                                        if c[:i]+nt+c[i+1:] not in core.stop_codons:
                                                try:
                                                        almost_stop_codons[c[:i]+nt+c[i+1:]] += almost_stop_codons[c]
                                                except KeyError:
                                                         almost_stop_codons[c[:i]+nt+c[i+1:]] = almost_stop_codons[c]
                for i in almost_stop_codons.keys():
                        almost_stop_codons[i] = list(set(almost_stop_codons[i]))

        # creates a list of tuples of the form (index, ['list', 'of', 'stop', 'codons']) to be further pre-processed
        pre_matches = [(i, almost_stop_codons[codons[i]]) for i in range(0,len(codons)) if codons[i] in almost_stop_codons.keys()]
        matches = []

        # further processing to create actual tuples (index, 'codon')
        for m in pre_matches:
                for codon in m[1]:
                        matches.append((m[0],codon))

        return matches

def find_start_codon_mutants(codons, n):
	'''
	Given a list `codons`, makes up to `n` mutations (n<2)
	to destroy the start codon. Returns a formatted list
	of tuples in the form `(index, 'mutated codon')`.
	'''
	start = codons[0]
	muts = []
	for i in list(range(0,3)):
		for nt in 'ACTG':
			mut_codon = start[:i]+nt+start[i+1:]
			if mut_codon != start and mut_codon != 'ATG':
				muts.append(mut_codon)
	new_muts = [z for z in muts]
	if n == 2:
		for e in muts:
			for i in list(range(0,3)):
				for nt in 'ACTG':
					mut_codon = e[:i]+nt+e[i+1:]
					if mut_codon != start and mut_codon != 'ATG' and mut_codon not in new_muts:
						new_muts.append(mut_codon)
					
	return [(0,m) for m in new_muts]
