from pyviko import core, restriction

# TODO - scrambling, excision KO?
# TODO - allow for input of full overprinted gene (DONE: Overlap finder, now just need to handle input)

class OverGene:
	
	frame = 1
	startNucleotide = 0
	sequence = ''
	overAAs = ''
	
	def __init__(self, startNt, seq, frameOver = 1):
		self.startNucleotide = startNt
		self.sequence = seq
		self.frame = frameOver
		self.overAAs = core.translate(core.findOverprintedGene(seq, startNt-1, frameOver))

class Mutant:
	
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
		
	def setOverGene(self, startNt, overFrame = 1):
		'''
		Adds the overprinted gene to the current `Mutant` object.
		'''
		self.overGene = OverGene(startNt, self.seq, overFrame)
		
	def vector(self, sequence):
		'''
		Adds the vector sequence to the current `Mutant` object 
		(primarily for making primers of early knockouts).
		'''
		if self.seq in sequence:
			self.vectorSeq = sequence
		else:
			raise core.SequenceError("Could not find target sequence in vector sequence.")
			
	def addMutant(self, mut):
		'''
		Add individual point mutants by hand, mut tuple
		'''
		#Should change all places I use `mut` to either ntMut or codMut for less ambiguity
		
	def findMutants(self, rSiteLength = 6, rSites = restriction.defaultEnzymes()):
		'''
		Returns a list of mutants that add a premature stop codon 
		(or mutate the start codon) without changing the overprinted 
		gene, and which add or remove a restriction site.
		'''
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
			
		stops = findPossibleStopCodons(self.codons, self.nMut)
			
		if self.overGene:
			safeMutations = []
			for poss in stops:
				nCodons = [codon for codon in self.codons]
				newCodons = core.insertMutation(nCodons, poss)
				newOverAAs = core.translate(core.findOverprintedGene(core.seqify(newCodons), self.overGene.startNucleotide-1, self.overGene.frame))
				if newOverAAs == self.overGene.overAAs:
					safeMutations.append(poss)
		else:
			safeMutations = stops
		
		###### Two approaches: regex and non-regex.
		newSites = [] # list of lists
		
		### Regex:
		
		if self.regex:
			baseSites = restriction.reFindEnzymes(self.seq)
			
			for mut in safeMutations:
				newSites.append([])
				for length in restrictionSiteLengths:
					newSites[-1] += restriction.reFindEnzymes(core.seqify(core.insertMutation(self.codons, mut)))
				
		### Non-regex:
		else:
			
			baseSites = []
			for length in restrictionSiteLengths:
				baseSites += restriction.findNcutters(self.seq, length)
				
			# baseSites = set(baseSites) #for some reason this will not work... Collection? Should use sets though
			for mut in safeMutations:
				newSites.append([])
				for length in restrictionSiteLengths:
					newSites[-1] += restriction.findNcutters(core.seqify(core.insertMutation(self.codons, mut)), length)
					
					
		winners = {}
		for l in newSites:
			if l <> baseSites: #this is why I should use sets
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
		
		finalWinners = winners
		
		return finalWinners
		
		# also look for start codon KOs
		# next step: post process `winners` to make the tuples look pretty (- and +), store in Mutant object
	
def findPossibleStopCodons(codons, n):

	'''
	Given a list of `codons`, finds individual codons that
	can be mutated to a stop codon given `n` mutations. Returns 
	a list of tuples of the form `(index, 'codon')` where 
	'''
	if codons[-1] in core.stopCodons: 
		codons = codons[:-1] #remove c-terminal stop codon
	
	almostStopCodons = {}
	#build dict of codons that can be mutated to a stop codon 
	
	for c in core.stopCodons:
		for i in range(0,3):
			for nt in 'ACTG':
				if c[:i]+nt+c[i+1:] not in core.stopCodons:
					try:
						almostStopCodons[c[:i]+nt+c[i+1:]].append(c)
					except KeyError:
						almostStopCodons[c[:i]+nt+c[i+1:]] = [c]
	
	if n == 2:
		for c in almostStopCodons.keys():
			for i in range(0,3):
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

print restriction.defaultEnzymes()