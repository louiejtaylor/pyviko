if __name__ == '__main__':
	
	from pyviko import core, mutation, restriction
	
	# testing RC search
	ovr=['ATGATTACCCGGGTTTCCCAAAGGGTTTCATCCTAA']
	z='''     TTACCCGGGTTTCCCAAAGGGTTTCAT'''
	toKO  = ['ATGAAACCCTTTGGGAAACCCGGGTAA']
        # benchmarking overlapping in the reverse frame
	for i in range(len(toKO))[:1]:
		m = mutation.Mutant(toKO[i][1],numMutations=1,regEx=True)
		m.setOverGene(overSeq = ovr[i][1])
		print(m.findMutants(ignoreRxSites = False, rSiteLength='all')[:5])
		print(m.findMutants()[:5])
		print("done "+str(i))
	
