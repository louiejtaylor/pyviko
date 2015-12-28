'''
Data collection from viral corpus:
-First nt mutated
 -including start codon mutants
 -ignoring start codon mutants
 -including start codon mutants and restriction sites
 -excluding start codon mutants and including restriction sites
-collect all as individual position and ratio of position/length (measure of closeness to beginnning/knockout efficacy)
'''

#####
##temporary dev hack
import os
os.chdir('/home/louis/pyViKO')
#####

from pyviko import core, mutation, restriction

allFiles = sorted(os.listdir('test/dem/ko'))
firstSiteNoRx = {}
firstSiteRx = {}
firstSiteNoRxNoZero = {}
firstSiteRxNoZero = {}

firstSiteNoRxNONE = 0
firstSiteRxNONE = 0
firstSiteNoRxNoZeroNONE = 0
firstSiteRxNoZeroNONE = 0

total = 0
oldTotal = 0
err = 0

for f in allFiles:
	fFirstSiteNoRx = open('test/dem/NoRx/fFirstSiteNoRx'+str(total)+'.txt', 'w')
	fFirstSiteRx = open('test/dem/Rx/fFirstSiteRx'+str(total)+'.txt', 'w')
	fFirstSiteNoRxNoZero = open('test/dem/NoRxNoZero/fFirstSiteNoRxNoZero'+str(total)+'.txt', 'w')
	fFirstSiteRxNoZero = open('test/dem/RxNoZero/fFirstSiteRxNoZero'+str(total)+'.txt', 'w')
	ovr = core.readFasta('test/dem/over/'+f)
	toKO = core.readFasta('test/dem/ko/'+f)
	print '.',
	if total > oldTotal + 1000:
		print total,
		oldTotal = total
	
	for i in range(len(toKO)):
		
		if len(toKO[i][1]) > 7000:
			print ",",
		else:
			try:
				m = mutation.Mutant(toKO[i][1],numMutations=1,regEx=True)
				m.setOverGene(overSeq = ovr[i][1])
				geneLength = len(m.codons)
			except IndexError as e:
				#print e
				#print toKO[i]
				#print ovr[i]
				err += 1
			try:
				mutRx = m.findMutants(ignoreRxSites = False, rSiteLength='all')
				mutNoRx = m.findMutants()
				total += 1
				
				try:
					firstRx = mutRx[0][0][0]
					try:
						firstSiteRx[firstRx] += 1
					except KeyError:
						firstSiteRx[firstRx] = 1
					fFirstSiteRx.write(str(firstRx) + ' ' + str(geneLength) + '\n')
				except IndexError:
					firstSiteRxNONE += 1
					
					
				try:
					firstNoRx = mutNoRx[0][0]
					try:
						firstSiteNoRx[firstNoRx] += 1
					except KeyError:
						firstSiteNoRx[firstNoRx] = 1
					fFirstSiteNoRx.write(str(firstNoRx) + ' ' + str(geneLength) + '\n')
				except IndexError:
					firstSiteNoRxNONE += 1
				
				try:
					firstRxNoZero = False
					for m in mutRx:
						if m[0][0] <> 0:
							firstRxNoZero = m[0][0]
							break
						
					#If false gets added to the dict, that means the only possible mutants were start codon mutants
						
					
					try:
						firstSiteRxNoZero[firstRxNoZero] += 1
					except KeyError:
						firstSiteRxNoZero[firstRxNoZero] = 1
					fFirstSiteRxNoZero.write(str(firstRxNoZero) + ' ' + str(geneLength) + '\n')
				except IndexError:
					firstSiteRxNoZeroNONE += 1
					
				try:
					firstNoRxNoZero = False
					for m in mutNoRx:
						if m[0] <> 0:
							firstNoRxNoZero = m[0]
							break
						
					#If false gets added to the dict, that means the only possible mutants were start codon mutants
						
					try:
						firstSiteNoRxNoZero[firstNoRxNoZero] += 1
					except KeyError:
						firstSiteNoRxNoZero[firstNoRxNoZero] = 1
					fFirstSiteNoRxNoZero.write(str(firstNoRxNoZero) + ' ' + str(geneLength) + '\n')
				except IndexError:
					firstSiteNoRxNoZeroNONE += 1
					
			except:
				print "Unknown error with: "
				print toKO[i]
				print ovr[i]
			del m
			del mutRx
			del mutNoRx
	
	fFirstSiteNoRx.close()
	fFirstSiteRx.close()
	fFirstSiteNoRxNoZero.close()
	fFirstSiteRxNoZero.close()
	del fFirstSiteNoRx
	del fFirstSiteRx
	del fFirstSiteNoRxNoZero
	del fFirstSiteRxNoZero
	del toKO
	del ovr

print
print "firstSiteNoRxNONE: " + str(firstSiteNoRxNONE)
print firstSiteNoRx
print
print "firstSiteRxNONE: " + str(firstSiteRxNONE)
print firstSiteRx
print
print "firstSiteNoRxNoZeroNONE: " + str(firstSiteNoRxNoZeroNONE)
print firstSiteNoRxNoZero
print
print "firstSiteRxNoZeroNONE: " + str(firstSiteRxNoZeroNONE)
print firstSiteRxNoZero
print 
print "Total analyzed: " +str(total)
print "Errors: " +str(err)