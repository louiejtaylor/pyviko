from pyviko import core, restriction

class OverGene:
	# dummy variables
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
	
	# dummy variables
	seq = ''
	nMut = 1
	isOverGene = False
	mutants = [] #will hold tuple?!
	
	def __init__(self, sequence, numMutations = 1):
		self.seq = sequence
		self.nMut = numMutations
		
	def overGene(self, startNt, overFrame = 1):
	'''
	Adds the overprinted gene to the current sequence object.
	'''
		self.overGene = OverGene(startNt, self.seq, overFrame)
		self.isOverGene = True
		
	def findMutants(self, rSiteLength = 6, rSites = restriction.restrictionSites):
	'''
	Returns a list of mutants that add a premature stop codon 
	(or mutate the start codon) without changing the overprinted 
	gene, and which add or remove a restriction site.
	'''
		if rSiteLength == 'all':
			x = 7
		elif 5 <= rSiteLength and rSiteLength <= 10:
			x = 8
		else:
			#throw error
			x = 9
	
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

# pyviko.mutation.Mutant('ATGGCGCGGCTAAGGGCCTAA')
