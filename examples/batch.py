if __name__ == '__main__':
	
	#####
	##temporary dev hack
	import os
	os.chdir('..')
	#####
	
	from pyviko import core, mutation, restriction, bio

	##	
	#testing
	x = mutation.Mutant('ATGGAACAGGCACCAGAAGATCAAGGACCACAGAGGGAGCCATACAACGAATGGGCTTTAGAATTGTTGGAAGACCTAAAGAATGAGGCTCTGCGCCACTTTCCTCGGCCTTGGCTACATGGACTAGGGCAATACTTCTATAATACATATGGAGATACCTGGGAGGGAGTAGAGGCCATCATTAGGACACTACAACAACTGTTGTTTATACATTATAGGATTGGCTGTCAACATAGCAGGATAGGAATCACTCCTCAAAGGAGAAGGAATGGAGCCAGTAGATCCTGA')
	x.setOverGene(0, 2)
	print x.findMutants(rSiteLength = 'all')
	##
	ovr = core.readFasta('examples/over.fasta')
	toKO = core.readFasta('examples/ko.fasta')
	overlaps = [core.findOverlap(toKO[i][1],ovr[i][1]) for i in range(len(toKO))]
	print overlaps