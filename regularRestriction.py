# -*- coding: utf-8 -*-
# take advantage of the degenerate amino acid codon coding
'''
Testing restriction site finder methods. KEY:
	A
	C
	G
	T
	R = A/G (purine)
	Y = C/T (pyrimidine)
	W = A/T (weak)
	S = C/G (strong)
	M = A/C (amino)
	K = G/T (keto)
	B = C/G/T
	D = A/G/T
	H = A/C/T
	V = A/C/G
	N = A/C/G/T
'''

try:
	import regex as re
except ImportError:
	import re
	print "reeee"
from pyviko.restriction import restrictionSites, findNonRegexEnzymeSite, findEnzymeSiteRegex, findNcutters


###
# Unified site testing
sites = ['AAAGGG','AAASSS','SSS','SSSWWW','SWSR', 'ACGTRYSWMKBDHVN']
sequence = 'AAAGGGCCCTTTAGCTAGAGAGACAGACAACGTACGTATCGTAA'

###################################
# Testing enzyme regex search
# This will be important to actually implement as a function!
# findRestrictionSites or something
s = []
for si in sites:
	s.append(findEnzymeSiteRegex(si))
	
for ww in s:
	try:
		q = re.finditer(ww, sequence, overlapped=True)
	except TypeError: #no new regex module
		q = re.finditer(ww, sequence)		
		
	for i in q:
		print (i.start(), i.string[i.start():i.end()]),
	print ''
# We 
# LIMITATION: Does not search overlapping sequences.
# For palindromes, this is not an issue
###################################

###################################
# weak func testing cutters with a list of restriction sites as input
def testNcutters(seq, n, rlist):
	n_mers = []
	for i in range(0, len(seq) - (n-1)):
		n_mers.append((i, seq[i:i+n]))
	actualSites = []
	for s in n_mers:
		if s[1] in rlist:
			actualSites.append((s[0],s[1]))
	return actualSites

# Testing non-regex search. 
# Ideally, should have the same results for regex and non-regex searches
s = []
for sis in sites:
	s.append(findNonRegexEnzymeSite(sis))
	
for ww in s:
	print testNcutters(sequence, len(ww[0]), ww)
	
# looks good! 
# This search is complete (including overlapping sequences)
###################################
