if __name__ == '__main__':
	from pyviko import core, mutation, restriction	
	# testing RC search
	ovr=['ATGATTACCCGGGTTTCCCAAAGGGTTTCATCCTAA']
	#         TTACCCGGGTTTCCCAAAGGGTTTCAT bounds of gene in the minus dir`
	to_ko = ['ATGAAACCCTTTGGGAAACCCGGGTAA']
        # benchmarking overlapping in the reverse frame
	for i in range(len(to_ko))[:1]:
		m = mutation.Mutant(to_ko[i][1],num_mutations=1,regex=True)
		m.set_over_gene(over_seq = ovr[i][1])
		print(m.find_mutants(ignore_rx_sites = False, r_site_length='all')[:5])
		print(m.find_mutants()[:5])
		print("done "+str(i))
