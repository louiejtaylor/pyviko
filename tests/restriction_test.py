try:
	import regex as re
except ImportError:
	import re
	print("Module 'regex' not found")

from pyviko.restriction import restrictionSites, findNonRegexEnzymeSite, findEnzymeSiteRegex, findNcutters

# Unified site testing
sites = ['AAAGGG','AAASSS','SSS','SSSWWW','SWSR', 'ACGTRYSWMKBDHVN']
# sites = {'AAAGGG':1,'AAASSS':2,'SSS':3,'SSSWWW':4,'SWSR':5, 'ACGTRYSWMKBDHVN':6} #for true testing

sequence = 'AAAGGGCCCTTTAGCTAGAGAGACAGACAACGTACGTATCGTAA'

# Find sites w/ regex

s = []
for si in sites:
	s.append(find_enzyme_site_regex(si))

for ww in s:
	try:
		q = re.finditer(ww, sequence, overlapped=True)
	except TypeError: #no new regex module
		q = re.finditer(ww, sequence)
	for i in q:
		print(i.start(), i.string[i.start():i.end()],"\n")

s = []
for sis in sites:
	s.append(expand_ambiguous_sequence(sis))

for ww in s:
	print test_n_cutters(sequence, len(ww[0]), ww)

