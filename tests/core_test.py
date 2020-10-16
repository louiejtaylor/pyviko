import sys
sys.path.append('..')

import unittest
from pyviko.core import *

class coreTestCase(unittest.TestCase):
	"""Tests core function."""

	def setup(self):
		self.seq = 'ATGGCTAAATGACT'
		self.cod = ['ATG', 'GCT', 'AAA', 'TGA', 'CT']
		self.trans = 'MAK'
	
		self.cod_no_start = ['GGG', 'GCT', 'AAA', 'TGA', 'CT']
		self.trans_no_start = 'GAK'
		
		self.cod_no_stop = ['ATG', 'GCT', 'AAA']
		self.trans_no_stop = 'MAK'

		self.cod_invalid = ['ATG', 'GCT', 'XXX', 'TGA', 'CT']
		self.invalid_input = 42

	def teardown(self):
		pass

	# TESTS

	# codonify
	def test_seq_to_codon(self):
		"""Tests codonify takes seq -> cod"""
		self.assertEqual(codonify(self.seq), self.cod)
	def test_codon_to_codon(self):
		"""Tests codonify takes cod -> cod"""
		self.assertEqual(codonify(self.cod), self.cod)
	def test_codonify_invalid(self):
		"""Tests codonify's handling of invalid input types"""
		self.assertRaises(TypeError, codonify(self.invalid_input))

	# seqify 
	def test_codon_to_seq(self):
		"""Tests sequify takes cod -> seq"""
		self.assertEqual(seqify(self.cod), self.seq)
	def test_seq_to_seq(self):
		"""Tests sequify takes seq -> seq"""
		self.assertEqual(seqify(self.seq), self.seq)
	def test_seqify_invalid(self):
		"""Tests sequify's handling of invalid input types"""
		self.assertRaises(TypeError, seqify(self.invalid_input))
	
	# translate
	def test_translate_codon(self):
		"""Tests translate takes cod -> aa"""
		self.assertEqual(translate(self.cod), self.trans)
	def test_invalid_codon(self):
		"""Tests translate rasies SequenceError for invalid codon"""
		self.assertRaises(SequenceError, translate(self.cod_invalid))
	def test_start_codon_warning(self):
		"""Tests translate warns user if no start codon is found"""
	        self.assertWarns(Warning, translate(self.cod_no_start))
                self.assertEqual(translate(self.cod_no_start), self.trans_no_start)
	def test_stop_codon_warning(self):
		"""Tests translate warns user if no stop codon is found"""
                self.assertWarns(Warning, translate(self.cod_no_stop))
                self.assertEqual(translate(self.cod_no_stop), self.trans_no_stop)

	# insertMutation
	def test_successful_insertion(self):
		"""Tests mutation (str) sucessfully inserted to codon list"""
		self.assertEqual(insertMutation(self.cod, (2, 'GCC')), ['ATG', 'GCT', 'GCC', 'TGA', 'CT'])
	def test_incorrect_index(self):
		"""Error when user insertion index exceeds codon list length"""
	        self.assertRaises(IndexError, insertMutation(self.cod), 21)

	def test_invalid_mutation(self): 
		"""Tests mutation's handling of invalid input types"""
		self.assertRaises(TypeError, seqify(self.invalid_input))

	def test_overprinted_gene_id(self):
		"""
		Overprinted gene completely contained within sequence 
		""" 
		self.assertEqual(findOverprintedGene('AAGTTTCGCTTAAC', startIndex=1, frame=1), ['AGT', 'TTC', 'GCT'])
		self.assertEqual(findOverprintedGene('AATGTTCGCTTAA', startIndex=1, frame=1), ['ATG', 'TTC', 'GCT'])
		

if __name__ == '__main__':
	unittest.main()

