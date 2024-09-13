import argparse
import sys

def run(size = None, overlap = None, ref = None):
    
    print(f"Size of the synthetic reads = {size}", file=sys.stderr)
    print(f"Length of overlap between the synthetic reads = {overlap}", file=sys.stderr)
    size_list = [int(x) for x in size.split(',')]
    overlap_list = [int(x) for x in overlap.split(',')]
    assert len(size_list) ==  len(overlap_list)
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
    for ind, ref in enumerate(references):
        print(f"\nWorking on reference {ind}", file=sys.stderr)
        print(f"Length of reference {ind} = {len(ref)}", file=sys.stderr)
        for s, o in zip(size_list, overlap_list):
            start = 0
            while True:
                if start+s < len(ref):
                    read = ref[start:start+s].upper()
                else:
                    read = ref[start:].upper()
                    print(f">read_{count}")
                    print(read)
                    count += 1
                    break
                start = start + s - o        
                print(f">read_{count}")
                print(read)
                count += 1
        

if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='chop-reference-to-reads.py', description="Converts the reference into synthetic reads")
    parser.add_argument("-size", required=True, help="Comma-separated list of sizes of synthetic reads")
    parser.add_argument("-overlap", required=True, help="Comma-separated list of base overlap between sythetic reads for each lenght")
    parser.add_argument("-ref", required=True, help="Reference FASTA")

    options = parser.parse_args()

    run(**vars(options))