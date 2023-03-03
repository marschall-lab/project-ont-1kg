import sys
import cyvcf2
import pickle
import resource
import pandas as pd

vcf = sys.argv[1]
path_to_groups = sys.argv[2]
output_path = sys.argv[3]
vcf = cyvcf2.VCF(vcf)
sample_list = vcf.samples
print("INFO: Reading sample population origin...", file = sys.stderr)
df = pd.read_csv(path_to_groups, sep="\t")
df = df[["SAMPLE","GROUP"]].to_numpy()
group={}
for sample,grp in df:
    try:
        group[grp].append(sample_list.index(sample))
    except:
        group[grp] = [sample_list.index(sample)]

vcf_iter = iter(vcf)

done_looping = False

print("INFO: Reading variants and calculating value...", file = sys.stderr)
het = [0]*len(sample_list)
hom_ref = [0]*len(sample_list)
hom_alt = [0]*len(sample_list)
skipped = [0]*len(sample_list)
c = 0
while not done_looping:
    try:
        c += 1
        var = next(vcf_iter)
        genotype = var.genotypes
        for n,s in enumerate(sample_list):
            gen = genotype[n][:2]
            if gen[0] == -1:
                skipped[n] += 1
                continue
            if gen[0] == gen[1]:
                if gen[0] == 0:
                    hom_ref[n] += 1
                else:
                    hom_alt[n] += 1
            else:
                het[n] += 1
        if c%100000 == 0:
            print("INFO: Processed %d variants..."%(c), file = sys.stderr)

    except StopIteration:
        done_looping = True

vcf.close()

d = {"het": het,
    "hom_alt": hom_alt,
    "hom_ref": hom_ref,
    "skipped": skipped,
    "group": group,
    "sample list": sample_list}

with open(output_path, 'wb') as handle:
    pickle.dump(d, handle, protocol=pickle.HIGHEST_PROTOCOL)

print("INFO: Maximum Memory Usage: %.3fGB" %(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6), file=sys.stderr)
