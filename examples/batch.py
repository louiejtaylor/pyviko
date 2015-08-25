if __name__ == '__main__':
	
	#####
	##temporary dev hack
	import os,time
	os.chdir('..')
	#####
	
	from pyviko import core, mutation, restriction, bio

	##	
	#testing
#	
#	y = mutation.Mutant('ATGGCCCGGGACGCGCGCTTAGTTAGTTTCTCGAGATAG')
#	y.setOverGene(overSeq = 'ATGGATGGCCCGGGACGCGCGCTTAGTTAG')
#	print y.findMutants(rSiteLength = 'all')
#	print
#
#	x = mutation.Mutant('ATGGCCCGGGACGCGCGCTTAGTTAGTTTCTCGAGATAG')
#	#                 '''  M  A  R  D  A  R  L  V  S  F  S  R  -'''
#	#           '''  ATGGATGGCCCGGGACGCGCGCTTAGTTAG'''
#	#           '''    M  D  G  P  G  R  A  L  S  -'''
#
#	x.setOverGene(startNtIndex = -1, overFrame = 3)	
#	print x.findMutants(rSiteLength = 'all')
#	print 
	

	#print core.findOverlap('ATGGCCCGGGACGCGCGCTTAGTTAGTTTCTCGAGATAG','ATGGATGGCCCGGGACGCGCGCTTAGTTAG')
	
	#print mutation.mutateStartCodon(['ATG','GGC'], 1)	
	
	##

	ovr = core.readFasta('examples/over.fasta')
	toKO = core.readFasta('examples/ko.fasta')
	t0 = time.time()
	for i in range(len(toKO)):
		m = mutation.Mutant(toKO[i][1],numMutations=1,regEx=False)
		m.setOverGene(overSeq = ovr[i][1])
		print toKO[i][1]
		print m.findMutants(rSiteLength='all')[:5]
		print "done "+str(i)		
		print
	print time.time()-t0
	
	#overlaps = [core.findOverlap(toKO[i][1],ovr[i][1]) for i in range(len(toKO))]
	#print overlaps
	