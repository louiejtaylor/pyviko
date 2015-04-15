# -*- coding: utf-8 -*-

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

nucleotideMatrix = {'R':['A','G'],
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

from pyviko.restriction import restrictionSites
from pyviko.base import findOverprintedGene

def findNonRegexEnzymeSites(site):
    '''
    Builds a list of sequences that correspond to a given 
    restriction enzyme recognition site.
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

def findNcutters(seq, n):
    '''
    Find restriction sites of length `n` in a sequence `seq`
    in `O(n*m)` where n is the sequence length and m is the
    number of restriction enzymes.
    '''
    n_mers = []
    for i in range(0, len(seq) - (n-1)):
        n_mers.append((i, seq[i:i+n]))
    actualSites = []
    for s in n_mers:
        if s[1] in restrictionSites.keys():
            actualSites.append((s[0],restrictionSites[s[1]]))
            
    # Limitation: assumes dict restrictionSites is in non-regex form (could be very large)
    return actualSites

'''
# testing regex counter.
testSites = ['AA', 'TATGCG', 'ARA', 'RRR', 'NN', 'ABC', 'SS', 'NAV']
for t in testSites:
    print findNonRegexEnzymeSites(t)
# looks good!
'''
sequence = 'ATGGAACAAGCCCCGGAAGACCAAGGGCCACAAAGAGAGCCATACAATGAATGGACACTAGAATTATTAGATGAACTCAAACAGGAAGCAGTAAGACATTTTCCTAGACAGTGGCTTCATGATTTAGGACAGCACATTTATAACACATATGGAGACACTTGGGCGGGGGTTGAGGCTATCATAAGGATCCTGCAACAATTGCTGTTTATTCATTACAGAATTGGCTGCCAACATAGCAGAATAGGCATTCTGCCACAAGGAAGAAGGAGAAATGGATCCAATAGATCCTAA'

'''
print findNcutters(sequence, 6)
# initial test good
'''
'''
print findOverprintedGene(sequence, -1, 3)
# looks good!
'''