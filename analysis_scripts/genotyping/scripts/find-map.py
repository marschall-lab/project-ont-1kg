import argparse

parser = argparse.ArgumentParser(prog='find-map.py', description="Find maps from biallelic and multiallelic bubbles")
parser.add_argument('-vcf', metavar='vcf', help='Biallelic VCF')
parser.add_argument('-map', metavar='map', help='Bubble ID to Allele ID map')
parser.add_argument('-outdir', metavar='outdir', help='Output directory')
args = parser.parse_args()

def read_map(map):
    is_bub_biallelic={}
    allele_id_to_bub_id_map = {}
    bub_id_to_allele_id_map={}
    map_reader = open(map, 'r')
    for line in map_reader:
        if line[0] == '#':
            continue
        bub_id, allele_ids = line.strip().split('\t')
        allele_ids = allele_ids.split(',')
        bub_id_to_allele_id_map[bub_id] = allele_ids
        for a in allele_ids:
            allele_id_to_bub_id_map[a] = bub_id
        if len(allele_ids) == 1:
            is_bub_biallelic[bub_id] = True
        else:
            assert len(allele_ids) > 1
            is_bub_biallelic[bub_id] = False
    map_reader.close()

    return is_bub_biallelic, allele_id_to_bub_id_map, bub_id_to_allele_id_map

def get_sv_ids(vcf):
    sv_ids = []
    vcf_reader = open(vcf, 'r')
    for line in vcf_reader:
        if line[0] == '#':
            continue
        sv_id = line.split('\t')[2]
        sv_ids.append(sv_id) 
    vcf_reader.close()

    return sv_ids

is_bub_biallelic, allele_id_to_bub_id_map, bub_id_to_allele_id_map = read_map(args.map)
sv_ids = get_sv_ids(args.vcf)

# creating the stats
biallelic_sv_count=0
multiallelic_sv_count=0

# writing stats
stat_writer = open(args.outdir+'map-stats.txt', 'w')
# writing bi-allelic SVs
bi_writer = open(args.outdir+'biallelic-sv-ids.txt', 'w')
# writing multiallelic SVs
multi_writer = open(args.outdir+'multiallelic-sv-ids.txt', 'w')

for sv in sv_ids:
    if is_bub_biallelic[allele_id_to_bub_id_map[sv]]:
        biallelic_sv_count+=1
        print(sv, file=bi_writer)
    else:
        multiallelic_sv_count+=1
        print(sv, file=multi_writer)

print(f'Biallelic VNTR Contractions: {biallelic_sv_count}', file=stat_writer)
print(f'Multiallelic VNTR Contractions: {multiallelic_sv_count}', file=stat_writer)

stat_writer.close()
bi_writer.close()
multi_writer.close()