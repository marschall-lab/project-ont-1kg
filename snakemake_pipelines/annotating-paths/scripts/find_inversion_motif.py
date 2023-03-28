
import argparse
import sys
import logging
import resource
from collections import defaultdict, namedtuple, Counter
import gzip
import re
import pysam

logger = logging.getLogger(__name__)

class Edge:
    def __init__(self, id1, id2, orient1, orient2) -> None:
        self.id1 = id1
        self.id2 = id2
        self.orient1 = orient1
        self.orient2 = orient2

class Node:
    def __init__(self, id, chrom, SO, LN) -> None:
        self._id = id
        self._chrom = chrom
        self._SO = SO
        self._LN = LN
        self.edges = []
        self.forward_edges = []
        self.back_edges = []
    
    def add_edge(self, id, o1, o2, forward):
        if forward:
            self.forward_edges.append(Edge(self._id, id, o1, o2))
        else:
            self.back_edges.append(Edge(self._id, id, o1, o2))
    
    def check_inversion_motif(self):
        forward_nodes = []
        back_nodes = []
        for f in self.forward_edges:
            forward_nodes.append(f.id2)
        for f in self.back_edges:
            back_nodes.append(f.id2)
        
        common_nodes = set(forward_nodes).intersection(back_nodes)

        return common_nodes    



def setup_logging(debug):
    handler = logging.StreamHandler()
    root = logging.getLogger()
    root.addHandler(handler)
    root.setLevel(logging.DEBUG if debug else logging.INFO)

def validate_arguments(options):
    pass

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", action="store_true", default=False, help="Print debug messages")
    parser.add_argument("--gfa", required=True, help="GFA file with the sort keys (BO and NO tagged)")
    
    options = parser.parse_args()
    setup_logging(options.debug)
    
    validate_arguments(options)

    nodes = defaultdict(lambda: None)
    read_gfa(options.gfa, nodes)
    total_len = 0
    for keys, values in nodes.items():
        common_nodes = values.check_inversion_motif()
        if len(common_nodes) > 0:
            print("Node: %s\tSN: %s\tNumber of nodes with loop: %d\t Length of Node: %d"%(keys, nodes[keys]._chrom, len(common_nodes), nodes[keys]._LN))
            total_len += nodes[keys]._LN
    print("Total Length of Problematic Nodes: %d"%(total_len))
    memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    logger.info("\nMemory Information")
    logger.info("  Maximum memory usage: %.3f GB", memory_kb / 1e6)

def read_gfa(gfa, node):
    
    logger.info("\nINFO: Parsing GFA file and reading sort key information")
    logger.info("INFO: Requires the following tags: SN, SO, SR, BO, and NO. Either the sequence or the LN tag has to be provided also.")
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
        if line[0] == 'S':
            total_nodes += 1
            fields = line.split("\t")
            SN = None
            SO = -1
            BO = -1
            NO = -1
            LN = -1
            for f in fields:
                if f.startswith("SN:Z:"):
                    SN = f[5:]
                if f.startswith("SO:i:"):
                    SO = int(f[5:])
                if f.startswith("LN:i:"):
                    LN = int(f[5:])
                if f.startswith("BO:i:"):
                    BO = int(f[5:])
                    tagged_nodes += 1
                if f.startswith("NO:i:"):
                    NO = int(f[5:])
            node[fields[1]] = Node(id=fields[1], chrom=SN, SO=SO, LN=LN)
        elif line[0] == 'L':
            fields = line.split("\t")
            id1 = fields[1]
            o1 = fields[2]
            id2 = fields[3]
            o2 = fields[4]
            if o1 == "+":
                node[id1].add_edge(id2, o1, o2, True)
            elif o1 == "-":
                node[id1].add_edge(id2, o1, o2, False)
            if o2 == "+":
                node[id2].add_edge(id1, o2, o1, False)
            elif o2 == "-":
                node[id2].add_edge(id1, o2, o1, True)
            
    logger.info("INFO: %d total Nodes Processed"%(total_nodes))
    logger.info("INFO: %d nodes with tags"%(tagged_nodes))
    reader.close()

def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b'\x1f\x8b'

if __name__ == "__main__":
    main()