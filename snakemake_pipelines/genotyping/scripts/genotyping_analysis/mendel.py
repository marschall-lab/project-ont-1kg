import sys
import cyvcf2
import resource
import pickle

vcf = sys.argv[1]
output_path = sys.argv[2]
vcf = cyvcf2.VCF(vcf)
vcf_iter = iter(vcf)
sample_list = vcf.samples

done_looping = False

print("INFO: Reading variants and calculating value...", file = sys.stderr)
multiallelic = 0
multiallelic_not_genotyped = 0
multiallelic_not_genotyped_pos = []
not_genotyped = 0
not_genotyped_pos = []
het = []
aaf = []
tot = []
c = 0
while not done_looping:
    try:
        c += 1
        n_aa = 0
        n_het = 0
        n_tot = 0
        is_multiallelic = False
        x = 0
        position = None
        var = next(vcf_iter)
        genotype = var.genotypes
        for n,s in enumerate(sample_list):
            if n == 0:
                position = str(c-1)+"_"+str(var.CHROM)+"_"+str(var.POS)
            gen = genotype[n][:2]
            if len(var.ALT) > 1:
                is_multiallelic = True
                if gen[0] == -1:
                    x += 1
                continue
            if gen[0] == -1:
                continue
            n_tot += 1
            assert (gen[0] in [0,1]) and (gen[1] in [0,1])
            if gen[0] != gen[1]:
                n_het += 1
            if gen[0] == 1:
                n_aa += 1
            if gen[1] == 1:
                n_aa += 1
        if c%100000 == 0:
            print("INFO: Processed %d variants..."%(c), file = sys.stderr)
        if is_multiallelic:
            multiallelic += 1
            if x == len(sample_list):
                multiallelic_not_genotyped += 1
                multiallelic_not_genotyped_pos.append(position)
            continue
        try:
            het.append(n_het/n_tot)
            aaf.append(n_aa/(2*n_tot))
            tot.append(n_tot)
        except ZeroDivisionError:
            not_genotyped += 1
            not_genotyped_pos.append(position)

    except StopIteration:
        done_looping = True

d = {"het": het,
    "aaf": aaf,
    "tot": tot,
    "multiallelic_skipped": multiallelic,
    "multiallelic_not_genotyped": multiallelic_not_genotyped,
    "not_genotyped": not_genotyped,
    "multiallelic_not_genotyped_pos": multiallelic_not_genotyped_pos,
    "not_genotyped_pos": not_genotyped_pos
    }

with open(output_path, 'wb') as handle:
    pickle.dump(d, handle, protocol=pickle.HIGHEST_PROTOCOL)

print("\nINFO: Result Summary", file=sys.stderr)
print("INFO: Skipped %d multiallelic variant positions"%(multiallelic), file=sys.stderr)
print("INFO: Skipped %d variant positions which could not be geneotyped for any sample"%(not_genotyped), file=sys.stderr)
print("INFO: Maximum Memory Usage: %.3fGB" %(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6), file=sys.stderr)
