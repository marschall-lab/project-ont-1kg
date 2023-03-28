'''
Compare the list of bubbles produced by gfatools bubble and the scaffold node method Tobias proposed. Ideally create a Venn diagram to compare.
'''

'''
What I found?
- The bubbles found by the scaffold node definition and this minigraph function is almost the same.
- Regarding the inverse motifs, the minigraph method does not consider it as a scaffold node and hence the allele traversal between it is not defined. That seems to be the major difference.
- Our method also is not able to find the first bubble at the start and end of chromosomes because of the scaffold node definition.
'''

import argparse
import gzip
from collections import defaultdict, namedtuple
import sys
from matplotlib_venn import venn2
import matplotlib.pyplot as plt

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", action="store_true", default=False, help="Print debug messages")
    parser.add_argument("--gfa", required=True, help="GFA file with the sort keys (BO and NO tagged)")
    parser.add_argument("--bed", required=True, help="BED file produced by gfatools bubbles. Column 4 and 5 should have the flanking node ids for the bubble.")
    parser.add_argument("--print-list", help="Enter 'gfatools' to get list bubbles exculively in BED file. Enter 'scaffold' to get list of bubbles exclusively present in our method. Enter 'common' to get the common bubbles.")
    parser.add_argument("--venn", action="store_true", default=False, help="Draw Venn Diagram.")
    options = parser.parse_args()
    
    scaffold_bubbles = read_gfa(options.gfa)
    gfatools_bubbles = read_bed(options.bed)

    if options.venn:
        venn2([set(scaffold_bubbles), set(gfatools_bubbles)], set_labels=['Our Method', 'GFAtools'])
        plt.show()
    
    if options.print_list == 'gfatools':
        bubbles = list(set(gfatools_bubbles).difference(set(scaffold_bubbles)))
        for b in bubbles:
            print(b)
    elif options.print_list == 'scaffold':
        bubbles = list(set(scaffold_bubbles).difference(set(gfatools_bubbles)))
        for b in bubbles:
            print(b)
    elif options.print_list == 'common':
        bubbles = list(set(gfatools_bubbles).intersection(set(scaffold_bubbles)))
        for b in bubbles:
            print(b)
    

def read_bed(bed):
    
    print("Reading BED file.", file=sys.stderr)
    if is_file_gzipped(bed):
        reader = gzip.open(bed, 'rt')
    else:
        reader = open(bed, 'r')
    bubbles = []
    while True:
        line = reader.readline()
        if not line:
            break
        fields = line.rstrip().split('\t')
        bubble = fields[3]+fields[4]
        bubbles.append(bubble)
    reader.close()
    print("\tNumber of Bubbles found = %d"%(len(bubbles)), file=sys.stderr)
    return bubbles

def read_gfa(gfa):
    
    print("Reading GFA file.", file=sys.stderr)
    Node = namedtuple("Node", ["SN", "SO", "SR", "LN", "BO", "NO"], defaults=[-1, -1, None, -1, -1, -1])
    node = defaultdict(lambda: Node())
    if is_file_gzipped(gfa):
        reader = gzip.open(gfa, 'rt')
    else:
        reader = open(gfa, 'r')
    total_nodes = 0
    while True:
        line = reader.readline()
        if not line:
            break
        if line[0] != 'S':
            continue
        total_nodes += 1
        fields = line.split("\t")
        LN = -1
        SN = None
        SO = -1
        SR = -1
        BO = -1
        NO = -1
        for f in fields:
            if f.startswith("LN:i:"):
                LN = int(f[5:])
            if f.startswith("SN:Z:"):
                SN = f[5:]
            if f.startswith("SO:i:"):
                SO = int(f[5:])
            if f.startswith("SR:i:"):
                SR = int(f[5:])
            if f.startswith("BO:i:"):
                BO = int(f[5:])
            if f.startswith("NO:i:"):
                NO = int(f[5:])
        #Get the length from the sequence if the LN tag is not given
        node[fields[1]] = node[fields[1]]._replace(BO=BO, NO=NO, SO=SO, SR=SR, SN=SN, LN=LN)
        if node[fields[1]].LN == -1:
            node[fields[1]] = node[fields[1]]._replace(LN = len(fields[2]))
    reader.close()
    scaffold_nodes = [node_id for node_id, value in sorted(node.items(), key=lambda item: item[1].BO) if value.NO==0]
    bubbles = []
    for i in range(len(scaffold_nodes) - 1):
        # Take into consideration the chromosome switch
        if node[scaffold_nodes[i]].SN != node[scaffold_nodes[i+1]].SN:
            continue
        bubble = ">%s>%s"%(scaffold_nodes[i], scaffold_nodes[i+1])
        bubbles.append(bubble)
    
    print("\tNumber of Bubbles found = %d"%(len(bubbles)), file=sys.stderr)
    
    return bubbles

def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b'\x1f\x8b'

if __name__ == "__main__":
    main()