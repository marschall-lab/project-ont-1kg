import argparse
import matplotlib.pyplot as plt
from collections import defaultdict
from pywfa import WavefrontAligner

def run(miller = None, vienna = None, sample = None, output = None):
    
    miller_reader = open(miller, 'r')
    vienna_reader = open(vienna, 'r')
    miller_seq = defaultdict(lambda: [])
    vienna_seq = defaultdict(lambda: [])
    while True:
        line = miller_reader.readline()
        if not line:
            break
        chrom, pos, seq = line.rstrip().split('\t')
        miller_seq[(chrom, int(pos))].append(seq)    
    while True:
        line = vienna_reader.readline()
        if not line:
            break
        chrom, pos, seq = line.rstrip().split('\t')
        vienna_seq[(chrom, int(pos))].append(seq)
    
    scores = []
    aligner = WavefrontAligner(scope='score', distance='affine2p', span='end-to-end', heuristic='adaptive')
    for info in vienna_seq.keys():
        v_seqs = vienna_seq[info]
        m_seqs = miller_seq[info]
        if len(v_seqs) != len(m_seqs):
            continue
        if len(v_seqs) != 2:
            continue
        score_matrix = [[None, None], [None, None]]
        for i, m in enumerate(m_seqs):
            for j, v in enumerate(v_seqs):
                res = aligner(m, v)
                assert aligner.status == 0
                score_matrix[i][j] = res.score
        pair1 =  score_matrix[0][0] + score_matrix[1][1]
        pair2 =  score_matrix[0][1] + score_matrix[1][0]
        if pair1 > pair2:
            scores.append([score_matrix[0][0], score_matrix[1][1]])
        else:
            scores.append([score_matrix[0][1], score_matrix[1][0]])
    
    total_scores_1000 = [x[0]+x[1] for x in scores if x[0]+x[1] > -1000]
    total_scores = [x[0]+x[1] for x in scores]
    
    plt.figure(figsize=(10,10), dpi=200)
    plt.hist(total_scores_1000, bins=100)
    plt.xlabel("WFA Score")
    plt.title(f"{sample}-subset1000")
    plt.yscale('log')
    plt.savefig(output+'-subset1000.png')
    
    plt.figure(figsize=(10,10), dpi=200)
    plt.hist(total_scores, bins=100)
    plt.xlabel("WFA Score")
    plt.title(f"{sample}-full")
    plt.yscale('log')
    plt.savefig(output+'-full.png')
    
    
if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='plot-allele-distance-histogram.py', description="Creates a scatter plot of the repeating unit counts comparing the vamos run on miller vcf and vienna vcf")
    parser.add_argument("-miller", required=True, help="Miller Vamos VCF file")
    parser.add_argument("-vienna", required=True, help="Vienna Vamos VCF file")
    parser.add_argument("-sample", required=True, help="Sample label")
    parser.add_argument("-output", required=True, help="Output histogram prefix.")

    options = parser.parse_args()

    run(**vars(options))