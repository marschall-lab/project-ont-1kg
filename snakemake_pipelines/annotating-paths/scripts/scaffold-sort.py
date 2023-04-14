'''
Scaffold Sort:

The basis of the script is the BO and NO tags defined by Tobias Marschall in his order_gfa.py script. The script defines scaffold nodes as nodes whose removal increases total number of connected components.
Using this definition of scaffold nodes, the sorting has been done.

The main difference between this sorting and the bubble-sort.py script is that here the the entire alignment is processed and used to determine whether the alignment is reveresed.
Here inversions are not taken into account.
'''

import sys
import argparse
import logging
from collections import defaultdict, namedtuple
import functools
import gzip
import re
from pysam import libcbgzf
import resource
from whatshap.timer import StageTimer

logger = logging.getLogger(__name__)
timers = StageTimer()
def setup_logging(debug):
    handler = logging.StreamHandler()
    root = logging.getLogger()
    root.addHandler(handler)
    root.setLevel(logging.DEBUG if debug else logging.INFO)

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", action="store_true", default=False, help="Print debug messages")
    parser.add_argument("--gfa", required=True, help="GFA file with the sort keys (BO and NO tagged)")
    parser.add_argument("--gaf", required=True, help="Input GAF File")
    parser.add_argument("--output", default=None, help="Output GAF File path (Default: sys.stdout)")
    parser.add_argument("--bgzip", action='store_true', help="Flag to bgzip the output. Can only be given with --output.")
    
    options = parser.parse_args()
    setup_logging(options.debug)
    
    validate_arguments(options)
    
    bubble_sort(options)

    memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    logger.info("\nMemory Information")
    logger.info("  Maximum memory usage: %.3f GB", memory_kb / 1e6)
    logger.info("\nTime Summary:")
    logger.info("  Time to parse GFA file: %.3f"%(timers.elapsed("read_gfa")))
    logger.info("  Total time to sort GAF file: %.3f"%(timers.elapsed("total_sort")))
    logger.info("    Time to parse GAF file: %.3f"%(timers.elapsed("read_gaf")))
    logger.info("    Time to sort GAF file: %.3f"%(timers.elapsed("sort_gaf")))
    logger.info("    Time to write GAF file: %.3f"%(timers.elapsed("write_gaf")))
    

def bubble_sort(options):
    
    if options.output == None:
        writer = sys.stdout
    else:
        if options.bgzip:
            writer = libcbgzf.BGZFile(options.output, 'wb')
        else:
            writer = open(options.output, "w")
    nodes = defaultdict(lambda: [-1,-1])
    with timers("read_gfa"):
        read_gfa(options.gfa, nodes)
    with timers("total_sort"):
        sort(options.gaf, nodes, writer)
    writer.close()
    

def sort(gaf, nodes, writer):
    
    logger.info("\n##### Parsing GAF file and sorting it #####")
    if is_file_gzipped(gaf):
        reader = gzip.open(gaf, 'rt')
    else:
        reader = open(gaf, 'r')
    
    Alignment = namedtuple('Alignment', ['offset', 'BO', 'NO', 'start', 'confused'])
    gaf_alignments = []
    count_inverse = 0
    # First pass: Store all the alignment lines as minimally. Just storing line offset and alignment string.
    with timers("read_gaf"):
        logger.debug("Processing and Storing GAF Alignments...")
        while True:
            offset = reader.tell()
            line = reader.readline()
            if not line:
                break
            line = line.split('\t')
            bo, no, start, confused = process_alignment(line, nodes, offset)
            if confused:
                count_inverse += 1
            gaf_alignments.append(Alignment(offset=offset, BO=bo, NO=no, start=start, confused=confused))
    
    logger.info("Number of alignments with scaffold nodes in both orientation: %d"%(count_inverse))
    # Sorting the alignments based on BO and NO tag
    with timers("sort_gaf"):
        logger.debug("Sorting the alignments...")
        gaf_alignments.sort(key=functools.cmp_to_key(compare_gaf))
    
    # Writing the sorted file
    with timers("write_gaf"):
        logger.debug("Writing Output File...")
        for alignment in gaf_alignments:
            off = alignment.offset
            if alignment.confused:
                logger.debug("[DEBUG|INV] Offset of read in sorted file: %d"%(writer.tell()))
            reader.seek(off)
            line = reader.readline()
            write_to_file(line, writer)
    
    reader.close()

def read_gfa(gfa, node):
    
    logger.info("\n##### Parsing GFA file and reading sort key information #####")
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
        for f in fields:
            if f.startswith("BO:i:"):
                node[fields[1]][0] = int(f[5:])
                tagged_nodes += 1
            if f.startswith("NO:i:"):
                node[fields[1]][1] = int(f[5:])
            
    logger.info("Total Nodes Processed: %d"%(total_nodes))
    logger.info("Nodes with tags: %d"%(tagged_nodes))
    reader.close()


def process_alignment(line, nodes, offset):
    path = list(filter(None, re.split('(>)|(<)', line[5])))
    orient = None
    # If there is no scaffold node present in the alignment, then assuming that it is in the correct orientation.
    # TODO: Need to find a way to deal with such alignments.
    rv = False
    bo = None
    no = None
    start = None
    orient_list = []
    confused = False
    for n in path:
        if n in [">", "<"]:
            orient = n
            continue
        if nodes[n][0] == -1 or nodes[n][1] == -1:
            logger.debug("[ERR]\tOF:i:%d\tND:Z:%s"%(offset,n))
            continue
        if nodes[n][1] != 0:
            continue
        orient_list.append(orient)
    if orient_list.count('>') !=0 and orient_list.count('<') != 0:
        confused = True
    # TODO: What to do when the start node that we have is an untagged node?
    # One argument is to just sort them to the end and completely ignore them since they dont belong to any vcf bubble.
    # That works for my tools but not ideal for general purpose.
    if orient_list.count('>') < orient_list.count('<'):
        l = int(line[6])
        e = int(line[8])
        start = l-e
        n = path[-1]
        bo = nodes[n][0]
        no = nodes[n][1]
    else:
        start = int(line[7])
        n = path[1]
        bo = nodes[n][0]
        no = nodes[n][1]
    return bo, no, start, confused


def compare_gaf(al1, al2):
    # Since we are sorting in ascending order, al1 is above al2 if it comes before al2.
    # Have to consider the case when the start node is untagged and has BO and NO has -1. These should be sorted to the end and not the start
    if al1.BO == -1 and al2.BO == -1:
        if al1.offset < al2.offset:
            return -1
        else:
            return 1
    elif al1.BO == -1:
        return 1
    elif al2.BO == -1:
        return 1

    # Comparing BO tags
    if al1.BO < al2.BO:
        return -1
    if al1.BO > al2.BO:
        return 1
    
    # Comparing NO tags
    if al1.NO < al2.NO:
        return -1
    if al1.BO > al2.BO:
        return 1
    
    # Comparing start position in node
    if al1.start < al2.start:
        return -1
    if al1.start > al2.start:
        return 1
    
    # This will be execulted only when two reads start at the exact same position at the same node. Then we use offset values. So the read which comes first in the GAF file will come higher up.
    if al1.offset < al2.offset:
        return -1
    if al1.offset > al2.offset:
        return 1


def write_to_file(line, writer):
    try:
        writer.write(line)
    except TypeError:
        writer.write(str.encode(line))


def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b'\x1f\x8b'


def validate_arguments(options):
    if options.bgzip and not options.output:
        raise RuntimeError("--bgzip flag has been specified but not output path has been defined. Please define the output path.")


if __name__ == "__main__":
    main()
