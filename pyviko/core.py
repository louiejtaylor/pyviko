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

def translate_bio(seq):
	return str(Seq(seqify(seq)).translate().seq)

def translate(codons):
	'''
	Translates a list of DNA codons into the corresponding amino
	acids, stopping translation if a stop codon is encountered.
	'''
	codons = codonify(codons)
	aa = ''
	for i in list(range(0,len(codons))):
		if codons[i] in stop_codons or len(codons[i]) != 3:
			if codons[i] in stop_codons:
				aa = aa + '*'
			break
		try:
			aa = aa + translation[codons[i]]
		except KeyError as e:
			raise SequenceError("Invalid codon: " + e.message)
	return aa

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
		#NOTE: Not all viral genes are initiated with ATG.
		warnings.warn("The first codon of your sequence is not a start codon.")

	return codons

def find_overlap_new(seq1, seq2, min_overlap = 10):
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

def reverse_complement_bio(input_seq):
	return str(Seq(input_seq).reverse_complement().seq)

def reverse_complement(seq):
	'''
	Given a sequence `seq`, returns the reverse complement.
	'''
	seq = seqify(seq)
	pairs = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
	rev = ""
	#Here should add reverse complements for regex sites? i.e. Y -> R
	try:
		for nt in seq[::-1]:
			rev += pairs[nt]
	except KeyError:
		print("Unknown nucleotide '" + nt + "' encountered.")
		return "False"

	return rev

def find_overlap(seq1, seq2, minimum=10):
	'''
	Given two sequences, returns a tuple `(i1, i2)` where `i1`
	is the index in `seq1` where the overlap with `seq2` begins and
	`i2` is the corresponding index in `seq2`. `minimum` is the minimum
	overlap length considered (default 10).
	'''
	i1 = 0
	i2 = 0
	if seq1 in seq2:
		i2 = seq2.index(seq1)
	elif seq2 in seq1:
		i1 = seq1.index(seq2)
	else:
		max12 = 0
		max21 = 0
		l1 = len(seq1)
		l2 = len(seq2)
		overall = min(l1,l2) + 1
		for i in list(range(1, overall)):
			if seq1[:i] == seq2[l2-i:]:
				max21 = i
			if seq2[:i] == seq1[l1-i:]:
				max12 = i
		if max12 > max21:
			if max12 > minimum:
				i1 = l1 - max12
			else:
				raise SequenceError("No overlap detected between input sequences")
		else: #without this small overlaps equal on both sides are mishandled
			if max21 > minimum:
				i2 = l2 - max21
			else:
				raise SequenceError("No overlap detected between input sequences")

	return (i1, i2)

def read_fasta_bio(loc):
	return [(r.accession, str(r.seq)) for r in SeqIO.parse(loc)]

def read_fasta(loc): # Bio.SeqIO
	'''
	Reads in a FASTA file, returns tuples in the form
	`('> identifying information', 'sequence')`.
	'''
	#uses a lot of memory. use Bio's generators
	f = open(loc, 'r')
	seqs = []
	iden = ''
	seq = ''
	for line in f.readlines():
		if iden == '':
			try:
				if line.lstrip()[0] != '>':
					raise SequenceError("Invalid file format: id line doesn't begin with '>'")
				iden = line.strip()
			except IndexError: #blank line
				next
		else:
			try:
				if line.lstrip()[0] == '>':
					seqs.append((iden, seq))
					iden = line.strip()
					seq = ''
				else:
					seq += line.strip().upper().replace(' ', '')
			except IndexError: #blank line
				next
	seqs.append((iden, seq))
	f.close()
	return seqs

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
