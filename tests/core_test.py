import sys
sys.path.append('..')

import unittest
from pyviko.core import *

class coreTestCase(unittest.TestCase):
	"""Tests core function."""

	def setUp(self):
		self.seq = 'ATGGCTAAATGACT'
		self.cod = ['ATG', 'GCT', 'AAA', 'TGA', 'CT']
		self.trans = 'MAK'
	
		self.cod_no_start = ['GGG', 'GCT', 'AAA', 'TGA', 'CT']
		self.trans_no_start = 'GAK'
		
		self.cod_no_stop = ['ATG', 'GCT', 'AAA']
		self.trans_no_stop = 'MAK'

		self.cod_invalid = ['ATG', 'GCT', 'XXX', 'TGA', 'CT']
		self.invalid_input = 42

	def tearDown(self):
		pass

	## TESTS

	## Codonify ##
	def test_codonify_1(self):
		"""Tests codonify takes seq -> cod"""
		self.assertEqual(codonify(self.seq), self.cod)
	def test_codonify_2(self):
		"""Tests codonify takes cod -> cod"""
		self.assertEqual(codonify(self.cod), self.cod)
	def test_codonify_3(self):
		"""Tests codonify's handling of invalid input types"""
		self.assertRaises(TypeError, lambda: codonify(self.invalid_input))
	##############

	## Seqify ##
	def test_seqify_1(self):
		"""Tests sequify takes cod -> seq"""
		self.assertEqual(seqify(self.cod), self.seq)
	def test_seqify_2(self):
		"""Tests sequify takes seq -> seq"""
		self.assertEqual(seqify(self.seq), self.seq)
	def test_seqify_3(self):
		"""Tests sequify's handling of invalid input types"""
		self.assertRaises(TypeError, lambda: seqify(self.invalid_input))
	############
	
	## Translate ##
	def test_translate_1(self):
		"""Tests translate takes cod -> aa"""
		self.assertEqual(translate(self.cod), self.trans)
	def test_translate_2(self):
		"""Tests translate rasies TranslationError for invalid codon"""
		self.assertRaises(TranslationError, lambda: translate(self.cod_invalid))
	#def test_translate_3(self):
	#	"""Tests translate warns user if no start codon is found"""
	#
	#def test_translate_4(self):
	#	"""Tests translate warns user if no stop codon is found"""
	###############

	## insertMutation ##
	def test_insertMutation_1(self):
		"""Tests mutation (str) sucessfully inserted to codon list"""
		self.assertEqual(insertMutation(self.cod, (2, 'GCC')), ['TTC', 'GCT', 'GCC', 'TGA', 'CTA'])
	#def test_insertMutation_2(self):
	#	"""Tests mutation warns user when insertion index exceeds codon list length"""
	#
	#def test_insertMutation_3(self):
	#	"""Tests mutation warns user when the codon inserted may give translation error"""
	#
	def test_insertMutation_4(self): 
		"""Tests mutation's handling of invalid input types"""
		self.assertRaises(TypeError, lambda: seqify(self.invalid_input))
	####################


#	def test_findOverprintedGene(self):
#		"""
#		Various test for the findOverprintedGene function:
#		1. Overprinted gene completely contain in sequence
#		2. 
#		""" 
#		#self.assertEqual(findOverprintedGene('AAGTTTCGCTTAAC', startIndex=1, frame=1), ['AGT', 'TTC', 'GCT'])
#		self.assertEqual(findOverprintedGene('AATGTTCGCTTAA', startIndex=1, frame=1), ['ATG', 'TTC', 'GCT'])
		

if __name__ == '__main__':
	unittest.main()

