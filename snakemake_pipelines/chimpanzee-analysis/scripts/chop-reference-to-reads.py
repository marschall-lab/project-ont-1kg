import argparse
import sys

def chop_sequence(seq, s, o):
    
    start = 0
    reads = []
    while True:
        try:
            read = seq[start:start+s]
        except IndexError:
            read = seq[start:]
            reads.append(read)
            break
        reads.append(read)
        start = start + s - o
    return reads

def run(size = None, overlap = None, ref = None):
    
    references = []
    sequence = ''
    fasta_reader = open(ref, 'r')
    for line in fasta_reader:
        line = line.rstrip()
        if line.startswith('>'):
            references.append(sequence)
            sequence = ''
            continue
        sequence += line
    references = references[1:]
    count = 1
    for ref in references:
        start = 0
        while True:
            try:
                read = ref[start:start+size]
            except IndexError:
                read = ref[start:]
                print(f">read_{count}")
                print(read)
                count += 1    
                break
            start = start + size - overlap        
            print(f">read_{count}")
            print(read)
            count += 1
        

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='chop-reference-to-reads.py', description="Converts the reference into synthetic reads")
    parser.add_argument("-size", required=True, type=int, help="Size of synthetic reads")
    parser.add_argument("-overlap", required=True, type=int, help="Base overlap between sythetic reads")
    parser.add_argument("-ref", required=True, help="Reference FASTA")

    options = parser.parse_args()

    run(**vars(options))