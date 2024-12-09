'''
Compare two GAF files (made from same FASTA and GFA) made with minigraph with and without the --vc flag to see which has correct CIGAR.
'''

'''
What I found?
- Vertex Coordinate GAF is correct
'''

import argparse
import re
import gzip
import logging
from collections import defaultdict, namedtuple, Counter
import sys
from pysam import libcbgzf, FastaFile

handler = logging.StreamHandler()
logger = logging.getLogger()
logger.addHandler(handler)
logger.setLevel(logging.INFO)

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--stable", required=True, help="Stable Coordinate GAF File")
    parser.add_argument("--vertex", required=True, help="Vertex Coordinate GAF File")
    parser.add_argument("--rgfa", required=True, help="rGFA File")
    parser.add_argument("--fasta", required=True, help="Read FASATA File")
    options = parser.parse_args()
    
    logger.info('INFO: Finding the alignment lines with reversed CIGARs')
    pos = find_reversed_cigars(options)

    evaluate_cigar(options, pos)
    #print(pos)


def evaluate_cigar(options, pos):
    logger.info('INFO: Loading FASTA.')
    fasta_file = FastaFile(options.fasta)
    logger.info('INFO: Loading rGFA.')
    gfa = parse_gfa_file(options.rgfa)

    if is_file_gzipped(options.vertex):
        reader = libcbgzf.BGZFile(options.vertex, 'rb')
    else:
        reader = open(options.vertex, 'r')

    vertex_correct_list = []
    stable_correct_list = []
    for p in pos:
        reader.seek(p)
        line = reader.readline()
        try:
            line = line.rstrip().split("\t")
        except TypeError:
            line = line.decode('utf8').rstrip().split("\t")

        query = fasta_file.fetch(line[0])[int(line[2]):int(line[3])]
        reference = get_path_sequence(line[5], gfa)[int(line[7]):int(line[8])]
        for f in line[12:]:
            if f[:5] == "cg:Z:":
                cg = f
                break
        correct = check_cigar_correctness(cg, reference, query)
        vertex_correct_list.append(correct)
        if not correct:
            correct = check_cigar_correctness(reverse_cigar(cg), reference, query)
            stable_correct_list.append(correct)

        
    
    print("Counter for the correctness of Vertex Coordinate GAF CIGARs: ", Counter(vertex_correct_list))
    print("The False instances from the above Counter is then further checked if the Stable GAF CIGAR is the correct CIGAR")
    print("Counter for the correctness of Stable Coordinate GAF CIGARs: ", Counter(stable_correct_list))


        

    reader.close()


def find_reversed_cigars(options):
    if is_file_gzipped(options.stable):
        stable_reader = libcbgzf.BGZFile(options.stable, 'rb')
    else:
        stable_reader = open(options.stable, 'r')
    if is_file_gzipped(options.vertex):
        vertex_reader = libcbgzf.BGZFile(options.vertex, 'rb')
    else:
        vertex_reader = open(options.vertex, 'r')
    
    total_alignments = 0
    wrong_alignments = 0
    reversed_cigar_pos = []
    while True:
        v_pos = vertex_reader.tell()
        s_line = stable_reader.readline()
        v_line = vertex_reader.readline()
        if not s_line:
            assert (not v_line)
            break
        total_alignments += 1

        try:
            s_line = s_line.rstrip().split("\t")
        except TypeError:
            s_line = s_line.decode('utf8').rstrip().split("\t")
        try:
            v_line = v_line.rstrip().split("\t")
        except TypeError:
            v_line = v_line.decode('utf8').rstrip().split("\t")

        assert (all([a == b for a,b in zip(s_line[0:4],v_line[0:4])]))

        s_cg = None
        v_cg = None
        for f in s_line[12:]:
            if f[:5] == "cg:Z:":
                s_cg = f
                break
        
        for f in v_line[12:]:
            if f[:5] == "cg:Z:":
                v_cg = f
                break
        
        if s_line[4] == '-':
            s_cg = reverse_cigar(s_cg)
        
        if v_cg != s_cg:
            assert s_cg == reverse_cigar(v_cg), "CIGAR is different"
            wrong_alignments += 1
            reversed_cigar_pos.append(v_pos)
    
    stable_reader.close()
    vertex_reader.close()

    print("Total Number of Alignments: ", total_alignments)
    print("Alignments with Reversed CIGAR: ", wrong_alignments)

    return reversed_cigar_pos

def get_path_sequence(path, gfa):
    seq = ''
    path = list(filter(None, re.split('(>)|(<)', path)))
    orient = None
    for nd in path:
        if nd in ['<','>']:
            orient = nd
            continue
        if orient == '>':
            seq += gfa[nd].sequence
        else:
            assert orient == '<'
            seq += reverse_complement(gfa[nd].sequence)
    return seq


def check_cigar_correctness(cigar, ref, query):
    cigar = cigar[5:]
    cg_tuples = []
    cg_letter_to_op = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, 'X': 7, '=': 8}
    cg = list(filter(None, re.split("([MIDNSHP=X])", cigar)))
    for i in range(0,len(cg),2):
        l = int(cg[i])
        op = cg_letter_to_op[cg[i+1]]
        cg_tuples.append((op,l))
    correct = True
    ref_pos = 0
    query_pos = 0
    for op,l in cg_tuples:
        if op == 0:
            RuntimeError('Cannot process CIGAR correctness with CIGAR operation M. Need = and X.')
        elif op == 1:
            query_pos += l
        elif op == 2:
            ref_pos += l
        elif op == 3:
            ref_pos += l
        elif op == 4:
            query_pos += l
        elif op == 5 or op == 6:
            pass
        elif op == 7:
            query_pos += l
            ref_pos += l
        elif op == 8:
            ref_in_match = ref[ref_pos:ref_pos+l]
            query_in_match = query[query_pos:query_pos+l]
            if ref_in_match != query_in_match:
                correct = False
                break
            query_pos += l
            ref_pos += l
        else:
            RuntimeError("Unknown CIGAR Operation")
    
    return correct
        
        


def parse_gfa_file(path):
    if is_file_gzipped(path):
        file = libcbgzf.BGZFile(path, 'rb')
    else:
        file = open(path, 'r')
    Node = namedtuple("Node", ['sequence', 'start', 'contig', 'tags'])
    node_dict = {}
    for line in file:
        if line[0] != "S":
            continue
        line = line.rstrip().split("\t")
        node_id = line[1]
        node_seq = line[2]
        tags = {}
        for i in line[3:]:
            i = i.split(":")
            if len(i) != 3:
                continue
            if i[1] == "i":
                tags[i[0]] = int(i[2])
            else:
                tags[i[0]] = i[2]
        # Assuming that these three tags are present
        node_contig = tags.pop("SN")
        node_start = tags.pop("SO")
        node_dict[node_id] = Node(node_seq, node_start, node_contig, tags)
    
    return node_dict

def reverse_complement(seq):
    seq = seq.replace("A", "t").replace(
        "C", "g").replace("T", "a").replace("G", "c")
    seq = seq.upper()
    
    seq = seq[::-1]
    return seq

def reverse_cigar(cg):
    import itertools
    cg = cg[5:]
    all_cigars = ["".join(x) for _, x in itertools.groupby(cg, key=str.isdigit)]
    new_cigar = "cg:Z:"
    for i in range(len(all_cigars), 0, -2):
        new_cigar += str(all_cigars[i-2]) + str(all_cigars[i-1])
    return new_cigar

def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b'\x1f\x8b'

if __name__ == "__main__":
    main()
