import argparse
from collections import Counter

def vcfformat(tag):
    if tag == 'i':
        return int
    elif tag == 'Z':
        return str
    else:
        exit('Unknown VCF Info format')

def run(gfa = None):
    
    gfa_reader = open(gfa, 'r')
    scaffold_node_dict = {}
    for line in gfa_reader:
        if line[0] != 'S':
            continue
        node_id = line.strip().split('\t')[1]
        fields = line.strip().split('\t')[3:]
        tags = {x.split(':')[0]:vcfformat(x.split(':')[1])(x.split(':')[2]) for x in fields}
        if tags['NO'] == 0:
            scaffold_node_dict[node_id] = tags
    gfa_reader.close()
    sorted(scaffold_node_dict.items(), key=lambda item: item[1]['BO'])
    node_list=list(scaffold_node_dict.items())
    print('#chrom\tstart_pos\tend_pos\tbubble_id')
    for index, values in enumerate(node_list[:-1]):
        s_node = node_list[index][0]
        e_node = node_list[index+1][0]
        if scaffold_node_dict[s_node]['SN'] != scaffold_node_dict[e_node]['SN']:
            continue
        chr = scaffold_node_dict[s_node]['SN']
        start_coord = scaffold_node_dict[s_node]['SO'] + scaffold_node_dict[s_node]['LN']
        end_coord = scaffold_node_dict[e_node]['SO']
        print(f'{chr}\t{start_coord}\t{end_coord}\t>{s_node}>{e_node}')

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='extract-bubbles.py', description="Extract bubbles from a tagged rGFA and outputs a BED file")
    parser.add_argument("-gfa", required=True, help="GFA file")
    
    options = parser.parse_args()

    run(**vars(options))