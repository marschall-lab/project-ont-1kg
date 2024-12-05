import argparse

parser = argparse.ArgumentParser(prog='find-vntr-map.py', description="Find VNTRs from biallelic and multiallelic bubbles")
parser.add_argument('-con', metavar='con', help='VNTR contraction vcf')
parser.add_argument('-exp', metavar='exp', help='VNTR expansion vcf')
parser.add_argument('-map', metavar='map', help='Bubble ID to Allele ID map')
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
deletion_sv_ids = get_sv_ids(args.con)
insertion_sv_ids = get_sv_ids(args.exp)

all_sv_ids = deletion_sv_ids+insertion_sv_ids

# creating the stats
biallelic_vntr_count=0
multiallelic_vntr_count=0

for sv in deletion_sv_ids:
    if is_bub_biallelic[allele_id_to_bub_id_map[sv]]:
        biallelic_vntr_count+=1
    else:
        multiallelic_vntr_count+=1

print(f'Biallelic VNTR Contractions: {biallelic_vntr_count}')
print(f'Multiallelic VNTR Contractions: {multiallelic_vntr_count}')

biallelic_vntr_count=0
multiallelic_vntr_count=0

for sv in insertion_sv_ids:
    if is_bub_biallelic[allele_id_to_bub_id_map[sv]]:
        biallelic_vntr_count+=1
    else:
        multiallelic_vntr_count+=1

print(f'Biallelic VNTR Expansions: {biallelic_vntr_count}')
print(f'Multiallelic VNTR Expansions: {multiallelic_vntr_count}')