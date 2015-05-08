import unittest
from core import codonify

class codonifyTestCase(unittest.TestCase):
    """Tests for `codonify` function."""

    def test_seq2cod(self):
        """Test sequence (str) sucessfully converted to codon (list)"""
        self.assertEqual(codonify('ATGCCGA'), ['ATG', 'CCG', 'A'])

if __name__ == '__main__':
    unittest.main()
