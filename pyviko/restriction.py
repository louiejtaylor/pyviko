import warnings
from pyviko import core

#from future import print_function

try:
	import regex as re
except ImportError:
	warnings.warn("To support overlapping restriction sites, please update to the new regex module.")
	import re

def findNonRegexEnzymeSite(site):
	'''
	Builds a list of sequences that correspond to a given 
	restriction enzyme recognition site (tree walking).
	'''
	possible_seqs = ['']
	for nt in site:
		if nt in 'ACGT':
			for i in xrange(len(possible_seqs)):
				possible_seqs[i] = possible_seqs[i] + nt
		else:
			placeholder = []
			for seq in possible_seqs:
				for possibility in nucleotideMatrix[nt]:
					placeholder.append(seq+possibility)
			possible_seqs = [ss for ss in placeholder]
	return possible_seqs
	
def findEnzymeSiteRegex(site):
	'''
	Returns a naive regular expression for a given
	restriction site.
	'''
	r_site = ''
	for nt in site:
		if nt in 'ACGT':
			r_site += nt
		else:
			r_site += '['
			try:
				for m in nucleotideMatrix[nt]:
					r_site += m
			except KeyError:
				print("Unknown nucleotide '" + nt  + "' encountered.")
			r_site += ']'
	return r_site

def generateEnzymeDict(enzDict):
	'''
	Function to help pyviko recognize both a restriction
	enzyme site and its reverse complement. Takes as input 
	a dictionary of restriction enzyme sites in the pyviko format
	`{'ntseq':['list','of','cutting,'enzymes'],...}` and returns
	a dictionary in the same format including the reverse compelement 
	of all input sites.
	'''	
	
	newDict = {}
	added = []
	for k in enzDict.keys():
		current = enzDict[k]
		for j in [k,core.reverseComplement(k)]:
			if j not in added:
				newDict[j]=current
				added.append(j)
			else:
				processed = newDict[j]
				for enzyme in current:
					if enzyme not in processed:
						newDict[j].append(enzyme)
	
	return newDict

def findNcutters(seq, n, rSites = None):
	'''
	Find restriction sites of length `n` in a sequence `seq`
	in `O(n*m)` where `n` is the sequence length and `m` is the
	number of restriction enzymes. Returns a list of tuples
	of the form `(site index, 'enzyme name')`.
	'''	
	
	if rSites == None:
		rSites = defaultEnzymes()
		
	tempSites = {}

	for si in rSites.keys():
		for se in findNonRegexEnzymeSite(si):
			tempSites[se] = rSites[si]

	recKeys = tempSites.keys()
	actualSites = []
	for i in list(range(0, len(seq) - (n-1))):
		if seq[i:i+n] in recKeys: 
			actualSites.append((i,seq[i:i+n]))		
			
	return actualSites
	
nucleotideMatrix = {	'R':['A','G'],
					'Y':['C','T'],
					'W':['A','T'],
					'S':['C','G'],
					'M':['A','C'],
					'K':['G','T'],
					'B':['C','G','T'],
					'D':['A','G','T'],			  
					'H':['A','C','T'],
					'V':['A','C','G'],
					'N':['A','C','G','T']}

def reFindEnzymes(seq, rSites=None):
	'''
	Find restriction sites in a sequence `seq`
	using regular expressions. Returns a list of tuples
	of the form `(site index, 'enzyme name')`.
	'''
	if rSites == None:
		rSites = defaultEnzymes()
	actualSites = []
	for site in rSites.keys():
		regexSite = findEnzymeSiteRegex(site)
		try:
			matches = re.finditer(regexSite, seq, overlapped=True)
		except TypeError: #if no new regex module
			matches = re.finditer(regexSite, seq)	
		for result in matches:
			actualSites.append((result.start(),rSites[site]))
			
	return actualSites
	
def defaultEnzymes():
	'''
	Returns the default enzyme set (New England BioLabs)
	as a dictionary.
	'''
	return {'GCGATCGC': ['AsiSI'], 'ACATGT': ['PciI'], 'TCCGGA': ['BspEI'], 'GGCGCGCC': ['AscI'], 'AGTACT': ['ScaI'], 'ACGCGT': ['MluI'], 'CTCAAG': ['BpuEI'], 'GCCGAG': ['NmeAIII'], 'GGTACC': ['KpnI', 'Acc65I'], 'GGCGGA': ['EciI'], 'CAGCTG': ['PvuII'], 'TTCGAA': ['BstBI'], 'CCCAGT': ['BmrI'], 'GATATC': ['EcoRV'], 'CTCCTC': ['BseRI'], 'CTTGAG': ['BpuEI'], 'CTGCAG': ['PstI'], 'CTGAAG': ['AcuI'], 'CTGCAC': ['BsgI'], 'CATCGC': ['BtgZI'], 'GTATCC': ['BciVI'], 'CCGCGG': ['SacII'], 'GATCC': ['AlwI', 'Nt.AlwI'], 'ATTTAAAT': ['SwaI'], 'TCACC': ['HphI'], 'CTCGAG': ['XhoI'], 'GAGTC': ['Nt.BstNBI', 'MlyI', 'PleI'], 'CGGCCG': ['EagI'], 'CTGGAG': ['BpmI'], 'ACGGC': ['BceAI'], 'ACTGG': ['BsrI'], 'GGTCTC': ['BsaI'], 'GATGG': ['BccI'], 'GACTC': ['Nt.BstNBI', 'MlyI', 'PleI'], 'ACTAGT': ['SpeI'], 'TCGCGA': ['NruI'], 'CTTAAG': ['AflII'], 'GCGTC': ['HgaI'], 'CCTAGG': ['AvrII'], 'GCAGGT': ['BspMI'], 'ATTAAT': ['AseI'], 'CCTTC': ['HpyAV'], 'GCTGGG': ['BseYI'], 'CAATTG': ['MfeI'], 'AGATCT': ['BglII'], 'ACTGGG': ['BmrI'], 'AGCGCT': ['AfeI'], 'CTCGGC': ['NmeAIII'], 'GCTAGC': ['BmtI', 'NheI'], 'CATATG': ['NdeI'], 'TCATGA': ['BspHI'], 'ACCTGC': ['BspMI'], 'GCTCTTC': ['Nt.BspQI', 'SapI'], 'CATTGC': ['BsrDI', 'Nb.BsrDI'], 'GTGCAC': ['ApaLI'], 'GTGCAG': ['BsgI'], 'CACGAG': ['BssSI'], 'CTGAG': ['BspCNI'], 'GAAGG': ['HpyAV'], 'GAAGA': ['MboII'], 'ACCGGT': ['AgeI'], 'GGCGCC': ['KasI', 'SfoI', 'NarI', 'PluTI'], 'GTATAC': ['BstZ17I'], 'GTTAAC': ['HpaI'], 'GACGC': ['HgaI'], 'ATCGAT': ['ClaI'], 'GAAGAG': ['EarI'], 'GGATCC': ['BamHI'], 'GAAGAC': ['BbsI'], 'CGTCTC': ['BsmBI'], 'GAGACC': ['BsaI'], 'CTCTTC': ['EarI'], 'TACGTA': ['SnaBI'], 'CCTCAGC': ['Nb.BbvCI', 'Nt.BbvCI', 'BbvCI'], 'CACTGC': ['BtsI', 'Bts-alpha-I'], 'GTCCC': ['BsmFI'], 'GCATC': ['SfaNI'], 'AAGCTT': ['HindIII'], 'TCCGCC': ['EciI'], 'GCGGG': ['FauI'], 'CTCGTG': ['BssSI'], 'GCCGT': ['BceAI'], 'CCCGGG': ['SmaI', 'XmaI'], 'GAGGAG': ['BseRI'], 'TGATCA': ['BclI'], 'GACGTC': ['AatII', 'ZraI'], 'GCAGC': ['BbvI'], 'TTTAAA': ['DraI'], 'GACGTG': ['BmgBI'], 'CATCC': ['BtsCI', 'FokI'], 'CTCCAG': ['BpmI'], 'AACGTT': ['AclI'], 'CCCGC': ['FauI'], 'GTTTAAAC': ['PmeI'], 'CCCAGC': ['BseYI'], 'TCTTC': ['MboII'], 'CCATGG': ['NcoI'], 'GAGAC': ['Nt.BsmAI', 'BsmAI'], 'CCATC': ['BccI'], 'CCGCTC': ['BsrBI'], 'GAAGAGC': ['Nt.BspQI', 'SapI'], 'CTGCTG': ['EcoP15I'], 'GGTGA': ['HphI'], 'GCGCGC': ['BssHII'], 'GGGAC': ['BsmFI'], 'GGATAC': ['BciVI'], 'GCATGC': ['SphI'], 'GTCGAC': ['SalI'], 'GCAGTG': ['BtsI', 'Bts-alpha-I'], 'TCTAGA': ['XbaI'], 'GAGCTC': ['SacI', 'Eco53kI'], 'GCTGAGG': ['Nb.BbvCI', 'Nt.BbvCI', 'BbvCI'], 'CGTACG': ['BsiWI'], 'GCATTC': ['Nb.BsmI', 'BsmI'], 'GAATGC': ['Nb.BsmI', 'BsmI'], 'CAGCAG': ['EcoP15I'], 'GGATG': ['BtsCI', 'FokI'], 'GGATC': ['AlwI', 'Nt.AlwI'], 'TGCGCA': ['FspI'], 'GCGGCCGC': ['NotI'], 'TGTACA': ['BsrGI'], 'GCCGGC': ['NaeI', 'NgoMIV'], 'TGGCCA': ['MscI'], 'CTTCAG': ['AcuI'], 'GAGCGG': ['BsrBI'], 'TTAATTAA': ['PacI'], 'GCTGC': ['BbvI'], 'ATGCAT': ['NsiI'], 'GATGC': ['SfaNI'], 'CTCAG': ['BspCNI'], 'AATATT': ['SspI'], 'CACGTG': ['PmlI'], 'CACGTC': ['BmgBI'], 'GAGACG': ['BsmBI'], 'CGATCG': ['PvuI'], 'GGGCCC': ['PspOMI', 'ApaI'], 'CCAGT': ['BsrI'], 'GGCCGGCC': ['FseI'], 'GCAATG': ['BsrDI', 'Nb.BsrDI'], 'CCTGCAGG': ['SbfI'], 'GTCTTC': ['BbsI'], 'GAATTC': ['EcoRI'], 'TTATAA': ['PsiI'], 'GTCTC': ['Nt.BsmAI', 'BsmAI'], 'GCGATG': ['BtgZI'], 'AGGCCT': ['StuI'], 'GTTGGA': ['MmeI'], 'GTCGGA': ['MmeI'], 'TCCAAC': ['MmeI'], 'TCCGAC': ['MmeI']}
