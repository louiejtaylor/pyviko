# -*- coding: utf-8 -*-

'''
Proof-of-concept module for pyViKO. [Examples and source](https://github.com/louiejtaylor/pyViKO)
'''

#from collections import Counter

# import os

try:
    from pyviko.base import codonify, seqify, insertMutation, translate, stopCodons, translation, restrictionSites
except ImportError as e:
    print "ImportError!",
    print e

#workhorse: find stop codons that can be generated given codons, mutation
def findPossibleStopCodons(codons, n):
    
    if codons[-1] in stopCodons:
        codons = codons[:-1] #remove c-terminal stop codon
    
    almostStopCodons = {}
    #build dict of codons that can be mutated to a stop codon 
    
    for c in stopCodons:
        for i in range(0,3):
            for nt in 'ACTG':
                if c[:i]+nt+c[i+1:] not in stopCodons:
                    try:
                        almostStopCodons[c[:i]+nt+c[i+1:]].append(c)
                    except KeyError:
                        almostStopCodons[c[:i]+nt+c[i+1:]] = [c]
                    # does this include duplicates??? Yes!
    if n == 2:
        #TESTING: 2 MUTATIONS
        # issue resolved, just need to take sets instead of using all possible mutants!
        for c in almostStopCodons.keys():
            for i in range(0,3):
                for nt in 'ACTG':
                    if c[:i]+nt+c[i+1:] not in stopCodons:
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
            
    print len(matches)
    #print matches
    return matches

#finds overprinted gene, given input sequence, frameshift, and bool startsBefore 
# TODO: change startsBefore to an **index** for gene    
def findOverprintedGene(seq, frame, startsBefore):    
    codons = codonify(seq[frame - 1:])[:-1] #remove last (incomplete) codon
    for i in range(0,len(codons)):
        if codons[i] in stopCodons:
            codons = codons[:i]
            break
    
    if not startsBefore:
        x = 1
        #
        # what
        #
        #To fill in, although I don't need this case yet
        #
        #
        #
        
    return codons


def findNonHarmfulMutations(seq, frame, startsBefore, numMutations):
    codons = codonify(seq)
    stops = findPossibleStopCodons(codons, numMutations)
    overAAs = translate(findOverprintedGene(seq, frame, startsBefore))
    winners = []
    for poss in stops:
        nCodons = [codon for codon in codons]
        newCodons = insertMutation(nCodons, poss)
        newOverAAs = translate(findOverprintedGene(seqify(newCodons), frame, startsBefore))
        if newOverAAs == overAAs:
            winners.append(poss)
    return winners

def findRestrictionSites(seq):
    all6mers = []
    for i in range(0, len(seq) - 5):
        all6mers.append((i, seq[i:i+6]))
    actualSites = []
    for s in all6mers:
        if s[1] in restrictionSites.keys():
            actualSites.append((s[0],restrictionSites[s[1]]))
    return actualSites

def findRestrictionSiteChanges(seq, frame, startsBefore, numMutations):
    safeMutations = findNonHarmfulMutations(seq, frame, startsBefore, numMutations)
    restrictionSites = findRestrictionSites(seq)
    newSites = []
    for i in safeMutations:
        newSites.append((i, findRestrictionSites(seqify(insertMutation(codonify(seq),i)))))
    winners = []
    for j in newSites:
        if j[1] <> restrictionSites:
            re = [r for r in restrictionSites]
            for site in j[1]:          
                try:
                    re.remove(site)
                except ValueError:
                    re.append((site,'+++'))
            k = (j[0], re)
            winners.append(k)
    
    return winners
        
    # print all6mers[0], all6mers[1], all6mers[-2], all6mers[-1]
    
sequence = '''ATGGAACAAGCCCCGGAAGACCAAGGGCCACAAAGAGAGCCATACAATGAATGGACACTA
GAATTATTAGATGAACTCAAACAGGAAGCAGTAAGACATTTTCCTAGACAGTGGCTTCAT
GATTTAGGACAGCACATTTATAACACATATGGAGACACTTGGGCGGGGGTTGAGGCTATC
ATAAGGATCCTGCAACAATTGCTGTTTATTCATTACAGAATTGGCTGCCAACATAGCAGA
ATAGGCATTCTGCCACAAGGAAGAAGGAGAAATGGATCCAATAGATCCTAA
'''.replace('\n','')

if __name__ == "__main__":
    mutFile = open('test/mutations.fasta',"w")
    #print findNonHarmfulMutations(sequence, 3, True, 1)
    mutFile.write("> Original sequence\n")
    mutFile.write(sequence+"\n")
    for f in findRestrictionSiteChanges(sequence, 3, True, 1):
        mutFile.write('> Mutant at codon ' + str(f[0][0]+1) +' '+ str(f[1])+'\n')
        #print f #for debugging, prints output to screen
        mutFile.write(seqify(insertMutation(codonify(sequence),f[0]))+'\n')
        #break
    mutFile.close()