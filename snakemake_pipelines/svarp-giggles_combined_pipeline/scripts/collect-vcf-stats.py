import argparse
from cyvcf2 import VCF
import sys
from collections import namedtuple
from collections import defaultdict

AlleleStats = namedtuple('AlleleStats','af ac an untyped')
GenotypeStats = namedtuple('GenotypeStats', 'heterozygosity het hom_ref hom_alt total')

def compute_allele_statistics(record):
    """
    Compute allele related statistics.
    """
    an = 0
    ac = 0
    unknown = 0
    for genotype in record.genotypes:
        alleles = genotype[:-1]
        assert 1 <= len(alleles) <= 2
        if len(alleles) == 1:
            # haploid genotype
            alleles.append(alleles[0])
        for a in alleles:
            if a == -1:
                unknown += 1
                continue
            assert a in [0,1]
            an += 1
            ac += a
    if an < 1:
        assert ac < 1
    af = ac / max(1.0, float(an))
    return AlleleStats(str(af), str(ac), str(an), str(unknown))

def compute_genotype_statistics(record, qualities=None):
    """
    Compute genotype related statistics.
    """
    counts = defaultdict(int)
    het_genotypes = 0
    hom_ref_genotypes = 0
    hom_alt_genotypes = 0
    total_genotypes = 0
    gqs = record.format('GQ') if qualities is not None else [None]*len(record.genotypes)
    for genotype, quality in zip(record.genotypes, gqs):
        alleles = genotype[:-1]
        assert 1 <= len(alleles) <= 2
        if len(alleles) == 1:
            # haploid genotype
            alleles.append(alleles[0])
        if not -1 in alleles:
            total_genotypes += 1
            if sum(alleles) == 0:
                assert alleles == [0,0]
                hom_ref_genotypes += 1
            elif sum(alleles) == 1:
                assert 0 in alleles
                assert 1 in alleles
                het_genotypes += 1
            elif sum(alleles) == 2:
                assert alleles == [1,1]
                hom_alt_genotypes += 1
            else:
                sys.stderr.write("Inconsistent allele detected for %s"%(record.INFO['ID']))
            # read GQ
            if qualities is not None:
                for q in qualities:
                    if int(quality) >= q:
                        counts[q] += 1
    genotype_stats = GenotypeStats( str(het_genotypes / max(1.0, float(total_genotypes))), str(het_genotypes), str(hom_ref_genotypes), str(hom_alt_genotypes), str(total_genotypes))
    return genotype_stats, counts

parser = argparse.ArgumentParser(prog='collect-vcf-stats.py', description="Collects the stats from the giggles callset vcf")
parser.add_argument('-panel', metavar='panel', help='Biallelic panel VCF.')
parser.add_argument('-callset', metavar='callset', help='Giggles multisample biallelic VCF')
args = parser.parse_args()

panel_reader = VCF(args.panel)
panel_stats = {}
callset_reader = VCF(args.callset)
callset_stats = {}

for variant in panel_reader:
    # require bi-allelic vcf with IDs
    assert len(variant.ALT) == 1
    var_id = variant.INFO['ID']
    allele_stats = compute_allele_statistics(variant)
    panel_stats[var_id] = allele_stats

sys.stderr.write("Completed generating panel statistics.\n")
sys.stderr.write("Found %d variant IDs.\n"%(len(panel_stats)))

quals = [0,50,100,200]
for variant in callset_reader:
    assert len(variant.ALT) == 1
    var_id = variant.INFO['ID']
    allele_stats = compute_allele_statistics(variant)
    genotype_stats, counts = compute_genotype_statistics(variant, quals)
    callset_stats[var_id] = [allele_stats, genotype_stats, counts]

sys.stderr.write("Completed generating callset statistics.\n")
sys.stderr.write("Qualities used for stat generation: %s.\n"%(','.join([str(q) for q in quals])))
sys.stderr.write("Found %d variant IDs.\n"%(len(callset_stats)))

# print stats for all IDs in genotypes VCF
header = [ 	'variant_id',
        'panel_allele_freq',
        'panel_alternative_alleles',
        'panel_total_alleles',
        'panel_unknown_alleles',

        'giggles_allele_freq',
        'giggles_alternative_alleles',
        'giggles_total_alleles',
        'giggles_unknown_alleles',
        'giggles_heterozygosity',
        'giggles_heterozygous_genotypes',
        'giggles_homozygous_reference_genotypes',
        'giggles_homozygous_alternate_genotypes',
        'giggles_total_genotypes',
    ]

for q in quals:
    header.append('giggles_GQ>=' + str(q))

print('\t'.join(header))

assert len(panel_stats) == len(callset_stats)

for var_id in callset_stats:
    if not var_id in panel_stats:
        continue
        
    line = [var_id,
            panel_stats[var_id].af,
            panel_stats[var_id].ac,
            panel_stats[var_id].an,
            panel_stats[var_id].untyped,

            callset_stats[var_id][0].af,
            callset_stats[var_id][0].ac,
            callset_stats[var_id][0].an,
            callset_stats[var_id][0].untyped,
            callset_stats[var_id][1].heterozygosity,
            callset_stats[var_id][1].het,
            callset_stats[var_id][1].hom_ref,
            callset_stats[var_id][1].hom_alt,
            callset_stats[var_id][1].total,
        ]
        
    # add counts for GQs
    for q in quals:
        line.append(str(callset_stats[var_id][2][q]))
    assert len(line) == len(header)
    print('\t'.join(line))