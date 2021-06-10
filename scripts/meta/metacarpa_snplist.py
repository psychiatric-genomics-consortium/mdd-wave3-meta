from bitarray import bitarray

# remove SNPs with rare patterns of missingness across cohorts

# example error
## ERROR:  Not enough variants with mask ❌ ❌ ❌ ✅ ❌ ❌ ❌ ❌ ✅ ❌ ❌ ❌ ❌ ❌ ❌ ❌ . Increase number of variants used to calculate the matrix, or check for allele issues (error occurred at position  1:906122 ).
# # Not enough variants with mask ❌ ❌ ❌ ❌ ❌ ❌ ✅ ✅ ❌ ❌ ❌ ❌ ❌ ❌ ✅ ❌ . Increase number of variants used to calculate the matrix, or check for allele issues (error occurred at position  4:122799607_A_C ). [4:122799607_C_A]

input = snakemake.input

# number of cohort sumstats
n = len(input)

# minimum count for each unique mask
mask_min = snakemake.params.mask_min

# store as text
# snp_default = '?' * n
# store as bitarray
snp_default = bitarray('0' * n)

# keep track of which SNPs appear in each study using a dictionary
# rsid is key and value is bitarray of whether the SNP was in the ith study
snps_tracker = {}

for i in range(n):
	print(input[i])
	f = open(input[i])
	sumstats_snps = f.read().splitlines()
	f.close()
	print(len(sumstats_snps))
	for snp in sumstats_snps:
		snp_update = snps_tracker.setdefault(snp, snp_default.copy())
		snp_update[i] = True

# tally up how many of each pattern of missingness we've seen
print("Tallying patterns")
patterns_tracker = {}

for bitpattern in snps_tracker.values():
	pattern01 = bitpattern.to01()
	pattern_count = patterns_tracker.setdefault(pattern01, 0)
	patterns_tracker[pattern01] = pattern_count + 1
	

# keep pattern that have been seen `mask_min` or more times and
# markers that appear in two or more cohorts
print("Filtering patterns")
keep_patterns = [pattern[0] for pattern in patterns_tracker.items() if pattern[1] >= mask_min and pattern[0].count('1') > 1]

# find SNPs that have one of these patterns
print("Filtering SNPs")
keep_snps = [snp[0] for snp in snps_tracker.items() if snp[1].to01() in keep_patterns]

print("Keeping")
print(len(keep_snps))

with open(snakemake.output[0], 'w') as out:
	out.writelines(map(lambda snp: snp+'\n', list(keep_snps)))