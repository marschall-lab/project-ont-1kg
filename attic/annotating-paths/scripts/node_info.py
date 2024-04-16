import argparse
import sys
import logging
import resource
from collections import defaultdict, namedtuple, Counter
import gzip

logger = logging.getLogger(__name__)
VariantRecord = namedtuple('VariantRecord', ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--gfa", required=True, help="GFA file")
    parser.add_argument("--node", required=True, help="Node to search. Can provide multiple nodes separated by comma.")
    
    options = parser.parse_args()
    Node = namedtuple("Node", ["SN", "SO", "SR", "LN", "BO", "NO"], defaults=[-1, -1, None, -1, -1, -1])
    nodes = defaultdict(lambda: Node())
    read_gfa(options.gfa, nodes)

    node_list = options.node.split(',')
    for node in node_list:
        print("Node %s: %s"%(node, nodes[node]))

def read_gfa(gfa, node):
    
    if is_file_gzipped(gfa):
        reader = gzip.open(gfa, 'rt')
    else:
        reader = open(gfa, 'r')
    total_nodes = 0
    tagged_nodes = 0
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
                tagged_nodes += 1
            if f.startswith("NO:i:"):
                NO = int(f[5:])
        #Get the length from the sequence if the LN tag is not given
        node[fields[1]] = node[fields[1]]._replace(BO=BO, NO=NO, SO=SO, SR=SR, SN=SN, LN=LN)
        if node[fields[1]].LN == -1:
            node[fields[1]] = node[fields[1]]._replace(LN = len(fields[2]))
            
    reader.close()

def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b'\x1f\x8b'


if __name__ == "__main__":
    main()