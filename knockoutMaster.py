# -*- coding: utf-8 -*-

#from collections import Counter

import os

stopCodons = ['TAG', 'TAA', 'TGA']

translation = {'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'ACA': 'T', 'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 'ACG': 'T', 'CAA': 'Q', 'AGT': 'S', 'CAG': 'Q', 'CCG': 'P', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R', 'CCA': 'P', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 'GGG': 'G', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S', 'TTA': 'L', 'GAC': 'D', 'TCC': 'S', 'GAA': 'E', 'TCA': 'S', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'TTC': 'F', 'GTT': 'V', 'GCT': 'A', 'ACC': 'T', 'TTG': 'L', 'CGT': 'R', 'CGC': 'R'}

restrictionSites = {'GCGATCGC': ['AsiSI'], 'ACATGT': ['PciI'], 'TCCGGA': ['BspEI'], 'GGCGCGCC': ['AscI'], 'AGTACT': ['ScaI'], 'ACGCGT': ['MluI'], 'CTCAAG': ['BpuEI'], 'GCCGAG': ['NmeAIII'], 'GGTACC': ['KpnI', 'Acc65I'], 'GGCGGA': ['EciI'], 'CAGCTG': ['PvuII'], 'TTCGAA': ['BstBI'], 'CCCAGT': ['BmrI'], 'GATATC': ['EcoRV'], 'CTCCTC': ['BseRI'], 'CTTGAG': ['BpuEI'], 'CTGCAG': ['PstI'], 'CTGAAG': ['AcuI'], 'CTGCAC': ['BsgI'], 'CATCGC': ['BtgZI'], 'GTATCC': ['BciVI'], 'CCGCGG': ['SacII'], 'GATCC': ['AlwI', 'Nt.AlwI'], 'ATTTAAAT': ['SwaI'], 'TCACC': ['HphI'], 'CTCGAG': ['XhoI'], 'GAGTC': ['Nt.BstNBI', 'MlyI', 'PleI'], 'CGGCCG': ['EagI'], 'CTGGAG': ['BpmI'], 'ACGGC': ['BceAI'], 'ACTGG': ['BsrI'], 'GGTCTC': ['BsaI'], 'GATGG': ['BccI'], 'GACTC': ['Nt.BstNBI', 'MlyI', 'PleI'], 'ACTAGT': ['SpeI'], 'TCGCGA': ['NruI'], 'CTTAAG': ['AflII'], 'GCGTC': ['HgaI'], 'CCTAGG': ['AvrII'], 'GCAGGT': ['BspMI'], 'ATTAAT': ['AseI'], 'CCTTC': ['HpyAV'], 'GCTGGG': ['BseYI'], 'CAATTG': ['MfeI'], 'AGATCT': ['BglII'], 'ACTGGG': ['BmrI'], 'AGCGCT': ['AfeI'], 'CTCGGC': ['NmeAIII'], 'GCTAGC': ['BmtI', 'NheI'], 'CATATG': ['NdeI'], 'TCATGA': ['BspHI'], 'ACCTGC': ['BspMI'], 'GCTCTTC': ['Nt.BspQI', 'SapI'], 'CATTGC': ['BsrDI', 'Nb.BsrDI'], 'GTGCAC': ['ApaLI'], 'GTGCAG': ['BsgI'], 'CACGAG': ['BssSI'], 'CTGAG': ['BspCNI'], 'GAAGG': ['HpyAV'], 'GAAGA': ['MboII'], 'ACCGGT': ['AgeI'], 'GGCGCC': ['KasI', 'SfoI', 'NarI', 'PluTI'], 'GTATAC': ['BstZ17I'], 'GTTAAC': ['HpaI'], 'GACGC': ['HgaI'], 'ATCGAT': ['ClaI'], 'GAAGAG': ['EarI'], 'GGATCC': ['BamHI'], 'GAAGAC': ['BbsI'], 'CGTCTC': ['BsmBI'], 'GAGACC': ['BsaI'], 'CTCTTC': ['EarI'], 'TACGTA': ['SnaBI'], 'CCTCAGC': ['Nb.BbvCI', 'Nt.BbvCI', 'BbvCI'], 'CACTGC': ['BtsI', 'Bts-alpha-I'], 'GTCCC': ['BsmFI'], 'GCATC': ['SfaNI'], 'AAGCTT': ['HindIII'], 'TCCGCC': ['EciI'], 'GCGGG': ['FauI'], 'CTCGTG': ['BssSI'], 'GCCGT': ['BceAI'], 'CCCGGG': ['SmaI', 'XmaI'], 'GAGGAG': ['BseRI'], 'TGATCA': ['BclI'], 'GACGTC': ['AatII', 'ZraI'], 'GCAGC': ['BbvI'], 'TTTAAA': ['DraI'], 'GACGTG': ['BmgBI'], 'CATCC': ['BtsCI', 'FokI'], 'CTCCAG': ['BpmI'], 'AACGTT': ['AclI'], 'CCCGC': ['FauI'], 'GTTTAAAC': ['PmeI'], 'CCCAGC': ['BseYI'], 'TCTTC': ['MboII'], 'CCATGG': ['NcoI'], 'GAGAC': ['Nt.BsmAI', 'BsmAI'], 'CCATC': ['BccI'], 'CCGCTC': ['BsrBI'], 'GAAGAGC': ['Nt.BspQI', 'SapI'], 'CTGCTG': ['EcoP15I'], 'GGTGA': ['HphI'], 'GCGCGC': ['BssHII'], 'GGGAC': ['BsmFI'], 'GGATAC': ['BciVI'], 'GCATGC': ['SphI'], 'GTCGAC': ['SalI'], 'GCAGTG': ['BtsI', 'Bts-alpha-I'], 'TCTAGA': ['XbaI'], 'GAGCTC': ['SacI', 'Eco53kI'], 'GCTGAGG': ['Nb.BbvCI', 'Nt.BbvCI', 'BbvCI'], 'CGTACG': ['BsiWI'], 'GCATTC': ['Nb.BsmI', 'BsmI'], 'GAATGC': ['Nb.BsmI', 'BsmI'], 'CAGCAG': ['EcoP15I'], 'GGATG': ['BtsCI', 'FokI'], 'GGATC': ['AlwI', 'Nt.AlwI'], 'TGCGCA': ['FspI'], 'GCGGCCGC': ['NotI'], 'TGTACA': ['BsrGI'], 'GCCGGC': ['NaeI', 'NgoMIV'], 'TGGCCA': ['MscI'], 'CTTCAG': ['AcuI'], 'GAGCGG': ['BsrBI'], 'TTAATTAA': ['PacI'], 'GCTGC': ['BbvI'], 'ATGCAT': ['NsiI'], 'GATGC': ['SfaNI'], 'CTCAG': ['BspCNI'], 'AATATT': ['SspI'], 'CACGTG': ['PmlI'], 'CACGTC': ['BmgBI'], 'GAGACG': ['BsmBI'], 'CGATCG': ['PvuI'], 'GGGCCC': ['PspOMI', 'ApaI'], 'CCAGT': ['BsrI'], 'GGCCGGCC': ['FseI'], 'GCAATG': ['BsrDI', 'Nb.BsrDI'], 'CCTGCAGG': ['SbfI'], 'GTCTTC': ['BbsI'], 'GAATTC': ['EcoRI'], 'TTATAA': ['PsiI'], 'GTCTC': ['Nt.BsmAI', 'BsmAI'], 'GCGATG': ['BtgZI'], 'AGGCCT': ['StuI'], 'GTTGGA': ['MmeI'], 'GTCGGA': ['MmeI'], 'TCCAAC': ['MmeI'], 'TCCGAC': ['MmeI']}

def codonify(sequence):
    '''
    Converts an input DNA sequence (str) to a list of codons.
    '''
    return [sequence[i:i+3] for i in range(0,len(sequence),3)]
    
def seqify(cod):
    '''
    Converts an input list of codons into a DNA sequence (str).
    '''
    sequence = ""
    for codon in cod:
        sequence += codon
    return sequence
    
def translate(codons):
    '''
    Translates a list of DNA codons into the corresponding amino acids (str).
    '''
    aa = ''
    for c in codons:
        aa = aa + translation[c]
    return aa
    
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
    '''
    #TESTING: 2 MUTATIONS
    # Major issue--horribly large amount of possibilities, need a different approach here
    for c in almostStopCodons.keys():
        for i in range(0,3):
            for nt in 'ACTG':
                if c[:i]+nt+c[i+1:] not in stopCodons:
                    try:
                        almostStopCodons[c[:i]+nt+c[i+1:]] += almostStopCodons[c]
                    except KeyError:
                         almostStopCodons[c[:i]+nt+c[i+1:]] = almostStopCodons[c]
    '''
    
    # creates a list of tuples of the form (index, ['list', 'of', 'stop', 'codons']) to be further pre-processed    
    preMatches = [(i, almostStopCodons[codons[i]]) for i in range(0,len(codons)) if codons[i] in almostStopCodons.keys()]
    
    matches = []  
    
    # further processing to create actual tuples (index, 'codon')    
    for m in preMatches:
        for codon in m[1]:
            matches.append((m[0],codon))
            
    #print len(matches)
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
    
sequence = '''ATGGAACAAGCCCCGGAAGACCAAGGGCCACAAAGAGAGCCATACAATGAATGGACACTA
GAATTATTAGATGAACTCAAACAGGAAGCAGTAAGACATTTTCCTAGACAGTGGCTTCAT
GATTTAGGACAGCACATTTATAACACATATGGAGACACTTGGGCGGGGGTTGAGGCTATC
ATAAGGATCCTGCAACAATTGCTGTTTATTCATTACAGAATTGGCTGCCAACATAGCAGA
ATAGGCATTCTGCCACAAGGAAGAAGGAGAAATGGATCCAATAGATCCTAA
'''.replace('\n','')

if __name__ == "__main__":
    mutFile = open('mutations.fasta',"w")
    #print findNonHarmfulMutations(sequence, 3, True, 1)
    mutFile.write("> Original sequence\n")
    mutFile.write(sequence+"\n")
    for f in findRestrictionSiteChanges(sequence, 3, True, 1):
        mutFile.write('> Mutant at codon ' + str(f[0][0]+1) +' '+ str(f[1])+'\n')
        #print f
        mutFile.write(seqify(insertMutation(codonify(sequence),f[0]))+'\n')
        #break
    mutFile.close()