#!/usr/bin/env python
import sys
from collections import namedtuple, defaultdict
from argparse import ArgumentParser
import re
from pywfa.align import WavefrontAligner


complement = str.maketrans('ACGT', 'TGCA')

class Node:
    def __init__(self, name, tags, sequence=None):
        self.name = name
        self.tags = tags
        self.sequence = sequence

class Edge:
    def __init__(self, from_node, from_dir, to_node, to_dir, overlap, tags):
        self.from_node = from_node
        self.from_dir = from_dir
        self.to_node = to_node
        self.to_dir = to_dir
        self.overlap = overlap
        self.tags = tags

def parse_tag(s):
    name, type_id, value = s.split(':')
    assert len(name) == 2
    if type_id == 'i':
        return name, int(value)
    elif type_id == 'Z':
        return name, value
    else:
        assert False

def parse_gfa(gfa_filename, with_sequence=False):
    nodes = {}
    edges = defaultdict(list)

    for nr, line in enumerate(open(gfa_filename)):
        fields = line.split('\t')
        if fields[0] == 'S':
            name = fields[1]
            tags = dict(parse_tag(s) for s in fields[3:])
            sequence = None
            if with_sequence and (fields[2] != '*'):
                sequence = fields[2]
            nodes[name] = Node(name,tags,sequence)
        elif fields[0] == 'L':
            from_node = fields[1]
            from_dir = fields[2]
            to_node = fields[3]
            to_dir = fields[4]
            overlap = fields[5]
            tags = dict(parse_tag(s) for s in fields[6:])
            e = Edge(from_node,from_dir,to_node,to_dir,overlap, tags)
            edges[(from_node,to_node)].append(e)

    return nodes, edges

GafLine = namedtuple("GafLine", "query_name query_length query_start query_end strand path path_length path_start path_end residue_matches alignment_block_length mapping_quality")

def parse_gaf(filename):
    for line in open(filename):
        fields = line.split('\t')
        yield GafLine(
            query_name = fields[0],
            query_length = int(fields[1]),
            query_start = int(fields[2]),
            query_end = int(fields[3]),
            strand = fields[4],
            path = fields[5],
            path_length = int(fields[6]),
            path_start = int(fields[7]),
            path_end = int(fields[8]),
            residue_matches = int(fields[9]),
            alignment_block_length = int(fields[10]),
            mapping_quality = int(fields[11]),
        )
    return

FastqEntry = namedtuple("FastqEntry", "name sequence")
def parse_fastq(filename):
    for n, line in enumerate(open(filename)):
        if n % 4 == 0:
            assert line.startswith('@')
            name = line[1:].strip()
        elif n % 4 == 1:
            yield FastqEntry(name, line.strip())

def joint_iterator(gaf_iterator, fastq_iterator):
    gaf_lines = []
    while True:
        try:
            fastq_entry = fastq_iterator.__next__()
            for gaf_line in gaf_lines:
                assert gaf_line.query_name == fastq_entry.name, 'GAF and FASTQ out of sync'
        except StopIteration:
            return
        try:
            while True:
                gaf_line = gaf_iterator.__next__()
                if gaf_line.query_name == fastq_entry.name:
                    gaf_lines.append(gaf_line)
                else:
                    yield fastq_entry, gaf_lines
                    gaf_lines = [gaf_line]
                    break
        except StopIteration:
            yield fastq_entry, gaf_lines
            return

def get_path(nodes, path):
    l = []
    for s in re.findall('[><][^><]+', path):
        node = nodes[s[1:]]
        print(node.name, len(node.sequence))
        if s[0] == '>':
            l.append(node.sequence)
        elif s[0] == '<':
            l.append(node.sequence[::-1].translate(complement))
        else:
            assert False
    return ''.join(l)

def add_arguments(parser):
    arg = parser.add_argument
    arg('graph', metavar='GRAPH', help='Input GFA file')
    arg('gaf', metavar='GAF', help='Input GAF file')
    arg('fastq', metavar='FASTQ', help='Input FASTQ file, order of reads have be match GAF')


if __name__ == "__main__":
    parser = ArgumentParser()
    add_arguments(parser)
    args = parser.parse_args()

    gfa_filename = args.graph

    print('Reading', gfa_filename)
    nodes, edges = parse_gfa(gfa_filename, with_sequence=True)

    for fastq_entry, gaf_lines in joint_iterator(parse_gaf(args.gaf), parse_fastq(args.fastq)):
        print('Processing read', fastq_entry.name)
        for gaf_line in gaf_lines:
            print('RESULT GAF', gaf_line.query_name, gaf_line.query_start, gaf_line.query_end, gaf_line.path, gaf_line.strand)
            path_sequence = get_path(nodes, gaf_line.path)
            #print(path_sequence)
            ref = path_sequence[gaf_line.path_start:gaf_line.path_end]
            query = fastq_entry.sequence[gaf_line.query_start:gaf_line.query_end]
            print('ref:  ', ref[:10], '...', ref[-10:])
            print('query:', query[:10], '...', query[-10:])
            aligner = WavefrontAligner(query)
            score = aligner.wavefront_align(ref)
            aligner.cigar_print_pretty()
