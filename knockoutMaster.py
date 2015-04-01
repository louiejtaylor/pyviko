# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 17:44:11 2015

@author: louis
"""

from collections import Counter

stopCodons = ['TAG', 'TAA', 'TGA']

translation = {'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'ACG': 'T', 'CAA': 'Q', 'AGT': 'S', 'CAG': 'Q', 'CCG': 'P', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CCA': 'P', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S', 'TTA': 'L', 'GAC': 'D', 'TCC': 'S', 'GAA': 'E', 'TCA': 'S', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T', 'TTG': 'L', 'CGT': 'R', 'CGC': 'R'}

#restrictionSites = {'CCCGGG':['XmaI'], 'TCCAAC':['MmeI'], 'TCCGAC':['MmeI'], 'GTTGGA':['MmeI']}

restrictionSites = {'GCGATCGC': ['AsiSI'], 'ACATGT': ['PciI'], 'TCCGGA': ['BspEI'], 'GGCGCGCC': ['AscI'], 'AGTACT': ['ScaI'], 'ACGCGT': ['MluI'], 'GCCGAG': ['NmeAIII'], 'GGTACC': ['KpnI', 'Acc65I'], 'GGCGGA': ['EciI'], 'CAGCTG': ['PvuII'], 'AAGCTT': ['HindIII'], 'GATATC': ['EcoRV'], 'CTTGAG': ['BpuEI'], 'CTGCAG': ['PstI'], 'GCTCTTC': ['SapI', 'Nt.BspQI'], 'GTATCC': ['BciVI'], 'CCGCGG': ['SacII'], 'ATTTAAAT': ['SwaI'], 'GCATC': ['SfaNI'], 'GAGTC': ['PleI', 'Nt.BstNBI', 'MlyI'], 'CGGCCG': ['EagI'], 'CTGGAG': ['BpmI'], 'ACGGC': ['BceAI'], 'ACTGG': ['BsrI'], 'GGTCTC': ['BsaI'], 'ACTAGT': ['SpeI'], 'TCGCGA': ['NruI'], 'CTTAAG': ['AflII'], 'CCTAGG': ['AvrII'], 'ATTAAT': ['AseI'], 'CCTTC': ['HpyAV'], 'CAATTG': ['MfeI'], 'AGATCT': ['BglII'], 'ACTGGG': ['BmrI'], 'AGCGCT': ['AfeI'], 'GCTAGC': ['NheI', 'BmtI'], 'CATATG': ['NdeI'], 'TCATGA': ['BspHI'], 'ACCTGC': ['BspMI'], 'CTGAAG': ['AcuI'], 'GTGCAC': ['ApaLI'], 'GTGCAG': ['BsgI'], 'CACGAG': ['BssSI'], 'GAAGA': ['MboII'], 'ACCGGT': ['AgeI'], 'CCCGGG': ['SmaI', 'XmaI'], 'GTATAC': ['BstZ17I'], 'GTTAAC': ['HpaI'], 'GACGC': ['HgaI'], 'GAGGAG': ['BseRI'], 'GGATCC': ['BamHI'], 'GAAGAC': ['BbsI'], 'CGTCTC': ['BsmBI'], 'CTCTTC': ['EarI'], 'TACGTA': ['SnaBI'], 'CCTCAGC': ['BbvCI', 'Nb.BbvCI', 'Nt.BbvCI'], 'CTCGAG': ['XhoI'], 'TTCGAA': ['BstBI'], 'GGCGCC': ['NarI', 'KasI', 'SfoI', 'PluTI'], 'ATCGAT': ['ClaI'], 'TGATCA': ['BclI'], 'GACGTC': ['AatII', 'ZraI'], 'GCAGC': ['BbvI'], 'TTTAAA': ['DraI'], 'AACGTT': ['AclI'], 'CCCGC': ['FauI'], 'GTTTAAAC': ['PmeI'], 'CCCAGC': ['BseYI'], 'CCATGG': ['NcoI'], 'CCATC': ['BccI'], 'CCGCTC': ['BsrBI'], 'GGTGA': ['HphI'], 'GCGCGC': ['BssHII'], 'GGGAC': ['BsmFI'], 'GCATGC': ['SphI'], 'GTCGAC': ['SalI'], 'GCAGTG': ['BtsI', 'Bts\xce\xb1I'], 'TCTAGA': ['XbaI'], 'GAGCTC': ['Eco53kI', 'SacI'], 'CGTACG': ['BsiWI'], 'GAATGC': ['BsmI', 'Nb.BsmI'], 'CAGCAG': ['EcoP15I'], 'GGATG': ['FokI', 'BtsCI'], 'GGATC': ['AlwI', 'Nt.AlwI'], 'TGCGCA': ['FspI'], 'GCGGCCGC': ['NotI'], 'TGTACA': ['BsrGI'], 'GCCGGC': ['NaeI', 'NgoMIV'], 'TGGCCA': ['MscI'], 'TTAATTAA': ['PacI'], 'ATGCAT': ['NsiI'], 'CTCAG': ['BspCNI'], 'AATATT': ['SspI'], 'CACGTG': ['PmlI'], 'CACGTC': ['BmgBI'], 'CGATCG': ['PvuI'], 'GGGCCC': ['ApaI', 'PspOMI'], 'GGCCGGCC': ['FseI'], 'GCAATG': ['BsrDI', 'Nb.BsrDI'], 'CCTGCAGG': ['SbfI'], 'GAATTC': ['EcoRI'], 'TTATAA': ['PsiI'], 'GTCTC': ['Nt.BsmAI', 'BsmAI'], 'GCGATG': ['BtgZI'], 'AGGCCT': ['StuI']}

def codonify(seq):
    return [seq[i:i+3] for i in range(0,len(seq),3)]
    
def seqify(cod):
    sequence = ""
    for codon in cod:
        sequence += codon
    return sequence
    
def translate(codons):
    aa = ''
    for c in codons:
        aa = aa + translation[c]
    return aa
    
def findPossibleStopCodons(codons, n):

    
    if codons[-1] in stopCodons:
        codons = codons[:-1] #remove c-terminal stop codon
    
    almostStopCodons = {}
    for c in stopCodons:
        for i in range(0,3):
            for nt in 'ACTG':
                if c[:i]+nt+c[i+1:] not in stopCodons:
                    almostStopCodons[c[:i]+nt+c[i+1:]] = c
                    # does this include duplicates???
                
    #almostStopCodons = [c for c in list(set(almostStopCodons)) if c[0] not in stopCodons]
    
    #matches = [(i, codons[i]) for i in range(0,len(codons)) if codons[i] in almostStopCodons]
    
    matches = [(i, almostStopCodons[codons[i]]) for i in range(0,len(codons)) if codons[i] in almostStopCodons.keys()]
    
    return matches
    
def findOverprintedGene(seq, frame, startsBefore):
    
    codons = codonify(seq[frame - 1:])[:-1] #remove last (incomplete) codon
    for i in range(0,len(codons)):
        if codons[i] in stopCodons:
            codons = codons[:i]
            break
    
    if not startsBefore:
        x = 1
        #
        #
        #
        #To fill in, although I don't need this case yet
        #
        #
        #
        
    return codons

def insertMutation(codons, mut):
    codons.pop(mut[0])
    newCodons = codons[:mut[0]] + [mut[1]] + codons[mut[0]:]
    return newCodons

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
    

sequence = 'ATGGCCTCAGGAGGCTGGCTGCCACCAACAGGAGGGGACCCGCCACAGGACCCTCCAAAAAACCCCAGAGAAGAGATCCCGGGGTGGTTGGAAACATGGGATTTACCAAGAGAGCCTTTCGACGAATGGCTCAGGGACATGTTACAGGATCTCAATGCAGAGGCCCAGTGCCACTTCCCAAGGAATCTCCTTTTCCGGCTTTGGTGGAACATTGTAGAGGAGCCAGCTATTGACAGGGGACAACCCAGACTAGAGGGATGGTATAAATATTATATCATAGTTCAGAGAGCTCTGTTTGTGCATATGAAAGGCAGGTGCTGTAAGCCCAAGACACATCCCGCATATGGCCCAGGAAGAGGGCCTCCAGGTCTGGGAGGAGCTCCAGGAGGAGCTGCAGCGGCCCCTCCAGGCCTGTAA'

#print findNonHarmfulMutations(sequence, 3, True, 1)
for f in findRestrictionSiteChanges(sequence, 3, True, 1):
    print f
    print seqify(insertMutation(codonify(sequence),f[0]))
    print ''