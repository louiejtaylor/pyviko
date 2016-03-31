if __name__ == '__main__':
	
	#####
	##temporary dev hack
	import os,time
	os.chdir('..')
	#####
	
	from pyviko import core, mutation, restriction
	
	#ovr = core.readFasta('test/dem/over/2000.fasta')
	#toKO = core.readFasta('test/dem/ko/2000.fasta')
	'''#True batch script
	t1 = time.time()
	for i in range(len(toKO)):
		m = mutation.Mutant(toKO[i][1],numMutations=1,regEx=True)
		m.setOverGene(overSeq = ovr[i][1])
		#print toKO[i][1]
		print m.findMutants(ignoreRxSites = False, rSiteLength='all')[:5]
		print
		print m.findMutants()[:5]
		print "done "+str(i)		
		print
	print time.time()-t1
	'''
	#testing RC search
	ovr=['ATGATTACCCGGGTTTCCCAAAGGGTTTCATCCTAA']
	z='''     TTACCCGGGTTTCCCAAAGGGTTTCAT'''
	toKO  = ['ATGAAACCCTTTGGGAAACCCGGGTAA']
	t1 = time.time()
	for i in range(len(toKO))[:1]:
		m = mutation.Mutant(toKO[i][1],numMutations=1,regEx=True)
		m.setOverGene(overSeq = ovr[i][1])
		#print toKO[i][1]
		print m.findMutants(ignoreRxSites = False, rSiteLength='all')[:5]
		print
		print m.findMutants()[:5]
		print "done "+str(i)		
		print
	print time.time()-t1
	#overlaps = [core.findOverlap(toKO[i][1],ovr[i][1]) for i in range(len(toKO))]
	#print overlaps
	