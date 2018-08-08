from pyviko import core, restriction

class OverGene:
	'''
	Class for the overprinted gene.
	'''
	frame = 1
	startNucleotideIndex = -1
	preSequence = '' # includes 1-2nt removed by core.findOverGene
	geneSequence = ''
	postSequence = '' # includes 1-2nt removed by core.findOverGene
	overAAs = ''
	
	def __init__(self, overSeq, startNtIndex, seq, frameOver = 1):
		if overSeq != '':
			ol = core.findOverlap(seq, overSeq)
			if ol[0] == 0: #overprinted gene starts before
				startNtIndex = -1
				frameOver = 4-(ol[1]%3)
				self.preSequence = overSeq[ol[1]-(-frameOver+4):ol[1]] + seq[:3-(-frameOver+4)]
			else: #overprinted gene starts after
				startNtIndex = ol[0]
		self.geneSequence = overSeq
		self.startNucleotideIndex = startNtIndex
		self.combSequence = seq
		self.frame = frameOver
		self.overAAs = core.translate(self.preSequence + core.seqify(core.findOverprintedGene(seq, startNtIndex, frameOver)) + self.postSequence)

class Mutant:
	'''
	Class for target gene mutagenesis.
	'''
	seq = ''
	codons = []
	nMut = 1
	overGene = False
	mutants = []
	
	def __init__(self, sequence, numMutations = 1, regEx = False):
		self.seq = sequence
		self.nMut = numMutations
		self.codons = core.codonify(sequence)
		self.regex = regEx
		
	def setOverGene(self, overSeq = '', startNtIndex = -1, overFrame = 1):
		'''
		Adds the overprinted gene to the current `Mutant` object.
		'''
		if overSeq == '' and startNtIndex == -1 and overFrame == 1:
			raise core.SequenceError("You must provide either the sequence of an overprinted gene, or its start position/frame in the knockout sequence.")
		else:
			self.overGene = OverGene(overSeq, startNtIndex, self.seq, overFrame)
		
	def vector(self, sequence):
		'''
		Adds the vector sequence to the current `Mutant` object 
		(primarily for making primers of early knockouts).
		'''
		if self.seq in sequence:
			self.vectorSeq = sequence
		else:
			raise core.SequenceError("Could not find target sequence in vector sequence.")
			
	def findMutants(self, ignoreRxSites = True, rSiteLength = 6, rSites = restriction.defaultEnzymes()):
		'''
		Returns a list of mutants that add a premature stop codon 
		(or mutate the start codon) without changing the overprinted 
		gene, and which add or remove a restriction site.
		'''
		stops = findPossibleStopCodons(self.codons, self.nMut)

		if self.overGene != False:
			if len(self.overGene.geneSequence) > 0:
				stops = mutateStartCodon(self.codons, self.nMut) + stops
			safeMutations = []
			newPreSequence = '';
			for poss in stops:
				nCodons = [codon for codon in self.codons]
				newCodons = core.insertMutation(nCodons, poss)
				if self.overGene.geneSequence != '':
					newPreSequence = self.overGene.preSequence[:4-self.overGene.frame] + newCodons[0][:self.overGene.frame - 1]
				newOverAAs = core.translate(newPreSequence + core.seqify(core.findOverprintedGene(core.seqify(newCodons), self.overGene.startNucleotideIndex, self.overGene.frame)))
				if newOverAAs == self.overGene.overAAs:
					safeMutations.append(poss)
		else:
			safeMutations = stops
		finalWinners = [s for s in safeMutations]
		if not ignoreRxSites:
			### Two approaches: regex and non-regex.
			restrictionSiteLengths = list(set([len(k) for k in rSites.keys()]))
			
			if rSiteLength == 'all':
				tempRestrictionSites = rSites
			elif rSiteLength >= min(restrictionSiteLengths) and rSiteLength <= max(restrictionSiteLengths):
				rKeys = [k for k in rSites.keys() if len(k) == rSiteLength]
				tempRestrictionSites = {} # Reduce size of dict. searched
				for site in rKeys:
					tempRestrictionSites[site] = rSites[site]
				restrictionSiteLengths = [rSiteLength]
			else:
				raise core.SequenceError("Invalid restriction site length.")		
			
			newSites = [] # list of lists
			
			### Regex:
			if self.regex:
				baseSites = restriction.reFindEnzymes(self.seq)
				
				for mut in safeMutations:
					newSites.append([])
					newSites[-1] += restriction.reFindEnzymes(core.seqify(core.insertMutation(self.codons, mut)))
					
			### Non-regex:
			else:
				
				baseSites = []
				for length in restrictionSiteLengths:
					baseSites += restriction.findNcutters(self.seq, length)
					
				for mut in safeMutations:
					newSites.append([])
					for length in restrictionSiteLengths:
						newSites[-1] += restriction.findNcutters(core.seqify(core.insertMutation(self.codons, mut)), length)
			
			winners = {}
			for l in newSites:
				if l != baseSites: #this is why I should use sets
					tempSites = [c for c in baseSites]
					tempAddedSites = []
					for site in l: #basically, removing everything in the new list from the old list to get the differences
						try:
							tempSites.remove(site)
						except ValueError:
							tempAddedSites.append(site)
					
					for i in range(0,len(tempSites)):
						tempSites[i] = (tempSites[i][0], tempSites[i][1], '-')
						
					for i in range (0, len(tempAddedSites)):
						tempAddedSites[i] = (tempAddedSites[i][0], tempAddedSites[i][1], '+')
					
					winners[safeMutations[newSites.index(l)]] = tempSites + tempAddedSites
			
			finalWinners = [(x,winners[x]) for x in sorted(winners.keys(), key=lambda x: x[0])]		
			
		return finalWinners
	
def findPossibleStopCodons(codons, n):
	'''
	Given a list `codons`, finds individual codons that
	can be mutated to a stop codon given `n` mutations. Returns 
	a list of tuples of the form `(index, 'codon')` where `'codon'`
	is the mutated codon. 
	'''
	#TODO: can optimize this as in js
	if codons[-1] in core.stopCodons: 
		codons = codons[:-1] #remove c-terminal stop codon
	
	almostStopCodons = {}
	#build dict of codons that can be mutated to a stop codon 
	
	for c in core.stopCodons:
		for i in list(range(0,3)):
			for nt in 'ACTG':
				if c[:i]+nt+c[i+1:] not in core.stopCodons:
					try:
						almostStopCodons[c[:i]+nt+c[i+1:]].append(c)
					except KeyError:
						almostStopCodons[c[:i]+nt+c[i+1:]] = [c]
	
	if n == 2:
		for c in almostStopCodons.keys():
			for i in list(range(0,3)):
				for nt in 'ACTG':
					if c[:i]+nt+c[i+1:] not in core.stopCodons:
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
			
	return matches

def mutateStartCodon(codons, n):
	'''
	Given a list `codons`, makes up to `n` mutations (n<2)
	to destroy the start codon. Returns a formatted list
	of tuples in the form `(index, 'mutated codon')`.
	'''
	start = codons[0]
	muts = []
	for i in list(range(0,3)):
		for nt in 'ACTG':
			mutCodon = start[:i]+nt+start[i+1:]
			if mutCodon != start and mutCodon != 'ATG':
				muts.append(mutCodon)
	newMuts = [z for z in muts]
	if n == 2:
		for e in muts:
			for i in list(range(0,3)):
				for nt in 'ACTG':
					mutCodon = e[:i]+nt+e[i+1:]
					if mutCodon != start and mutCodon != 'ATG' and mutCodon not in newMuts:
						newMuts.append(mutCodon)
					
	return [(0,m) for m in newMuts]
