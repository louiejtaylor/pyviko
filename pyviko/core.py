stop_codons = ['TAG', 'TAA', 'TGA']
translation = {'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'ACG': 'T', 'CAA': 'Q', 'AGT': 'S', 'CAG': 'Q', 'CCG': 'P', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CCA': 'P', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S', 'TTA': 'L', 'GAC': 'D', 'TCC': 'S', 'GAA': 'E', 'TCA': 'S', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T', 'TTG': 'L', 'CGT': 'R', 'CGC': 'R'}

from Bio.Seq import Seq
from Bio import SeqIO
import warnings, os

class SequenceError(Exception):
	pass

def codonify(sequence):
	'''
	Converts an input DNA sequence (str) to a list of codons.
	'''
	if type(sequence) == list:
		warnings.warn("passed list to codonify, returning unchanged")
		return sequence
	return [sequence[i:i+3] for i in list(range(0,len(sequence),3))]

def seqify(codons):
	'''
	Converts an input list of codons into a DNA sequence (str).
	'''
	if type(codons) == str:
		warnings.warn("passed string to seqify, returning unchanged")
		return codons
	return ''.join(codons)

def translate(seq):
	'''
	Translates a list of DNA codons into the corresponding amino
	acids, stopping translation if a stop codon is encountered.
	'''
	# this should become obsolete and just use the biopython/Seq builtin
	translation = str(Seq(seqify(seq)).translate())
	if "*" in translation: # match existing functionality--truncates at stop codon
		translation = translation[translation.index('*')]
	return translation

def mutate(seq, mut):
	'''
	Takes as input a sequence or list `seq` to mutate
	and a tuple `mut` in the form (index, 'mutated element')
	eg. `(3, 'A')` or `(3, 'AUG')`.
	Returns a modified nucleotide sequence or codon list.
	'''
	return seq[:mut[0]] + mut[1] + seq[mut[0]+1:]

def find_overprinted_gene(seq, startIndex, frame=1):
	'''
	Given a sequence `seq` and the `startIndex` of
	an overprinted gene, returns a list of codons that
	correspond to the overprinted gene. The `frame`
	argument is only necessary if the overprinted
	gene's start codon is before the input sequence, in
	which case `startIndex` must be -1. (0-indexed)
	'''

	if startIndex != -1:
		frame = 1   # In case `frame` argument provided erroneously
		codons = codonify(seq[startIndex:])[:-1] # Remove last (incomplete) codon
	else:
		if frame == 1:
			raise SequenceError("The overprinted sequence is in the same frame as the main coding sequence. Please provide a frame argument.")
		codons = codonify(seq[frame - 1:])[:-1] # Remove last (incomplete) codon
	for i in list(range(1,len(codons))):
		if codons[i] in stopCodons:
			codons = codons[:i]
			break

	if codons[0] != 'ATG' and startIndex != -1:
		warnings.warn("The first codon of your sequence is not a start codon.")
		# no need to error, fine if not, some viral transcripts don't start with ATG

	return codons

def find_overlap(seq1, seq2, min_overlap = 10):
	"""
	Find overlaps between two sequences. Returns index
	of the nucleotide in `seq1` where the overlap starts.
	Index is negative if the overlap is in the reverse frame.
	"""
	# same direction: excluding 5' parts
	for i in range(len(seq2)-min_overlap):
		if seq2[i:] in seq1:
			return seq1.index(seq2[i:])

	# same direction: excluding 3' parts
	for i in range(len(seq2)-min_overlap):
		if seq2[:-i] in seq1:
			return seq1.index(seq2[:-i])

	rc_seq2 = reverse_complement(seq2)
	# opposite direction: 5' exclusion
	for i in range(len(rc_seq2)-min_overlap):
		if rc_seq2[i:] in seq1:
			return -seq1.index(rc_seq2[i:])

	# opposite direction: 3' exclusion
	for i in range(len(rc_seq2)-min_overlap):
		if rc_seq2[:-i] in seq1:
			return -seq1.index(rc_seq2[:-i])

	# no return value here rather than error makes sense

def reverse_complement(input_seq):
	return str(Seq(input_seq).reverse_complement())

def read_fasta(loc):
	'''
	Reads in a FASTA file, returns tuples in the form
	`('> identifying information', 'sequence')`.
	'''
	# TODO: other code handle output of SeqIO.parse to not store in memory
	return [(">"+r.id, str(r.seq)) for r in SeqIO.parse(loc, "fasta")]

# should implement this as a wrapper around Bio.Seq--has a list of mutations and a helper func to generate a list of SeqRecords with appropriate names
#class MutantSeq() here or in mutation.py?

def write_fasta(fname, mutlist, seq, has_rx_sites = False, rloc = "", floc = ""): # wrapper around Bio.SeqIO?
	'''
	Given a filename `fname`, list of mutations `mutlist` input sequence `seq`
	and an optional file location  `rloc` (relative location) or `floc` (absolute location), 
	generates a FASTA file with all mutants in the sequence.
	Accepts input in two formats: if hasRxSites=False, assumes the mutations are of the form
	`(mutant codon index, 'stop codon')`.
	'''
	fname = fname.replace('|', '.')[:30]
	for character in fname:
		if character not in '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz.,\'"()-_':
			fname = fname.replace(character, '_')
	fname = ' '.join([i for i in fname.split(' ') if i != '']) + '.fasta'
	if floc == "":
		dirs = os.listdir(os.getcwd())
	else:
		dirs = os.listdir(floc)
	if fname in dirs:
		i = 1
		while fname[:-6] + '(' + str(i) + ').fasta' in dirs:
			i += 1
		fname = fname[:-6] + '(' + str(i) + ').fasta'
	del i
	fasta = open(fname, 'w')
	for m in mutlist:
		fasta.write(">mutant at codon " + str(m[0][0]+1) +': \n' )
		codons = codonify(seq)
		mut_seq = seqify(codons[:m[0][0]]+[m[0][1]]+codons[m[0][0]+1:])
		fasta.write('\n'.join([mut_seq[i:i+100] for i in range(0,len(mut_seq),100)]) +'\n')
	fasta.close()
	return True
