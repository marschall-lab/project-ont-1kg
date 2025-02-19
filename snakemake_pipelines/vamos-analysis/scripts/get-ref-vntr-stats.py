import argparse
import matplotlib.pyplot as plt

def run(vcf = None, plot = None):
    
    # parse reference VNTRs
    len_diff = []
    num_ref_len_mismatch = 0
    num_ref_len_match = 0
    with open(vcf, 'r') as vcffile:
        for line in vcffile:
            if line.startswith('#'):
                continue
            line = line.strip().split('\t')
            chr = line[0]
            start = int(line[1])
            end = None
            rus = None
            altanno = None
            info = line[7]
            assert line[9] == '1/1'
            for field in info.split(';'):
                if field.startswith('END'):
                    end = int(field.split('=')[1])
                if field.startswith('RU'):
                    rus = tuple(field.split('=')[1].split(','))
                if field.startswith('ALTANNO_H1'):
                    altanno = tuple(field.split('=')[1].split(','))
            altanno_len = 0
            for allele in altanno:
                altanno_len += len(rus[int(allele)])
            if altanno_len == (end-start):
                num_ref_len_match += 1
            else:
                num_ref_len_mismatch += 1
                len_diff.append(abs(altanno_len - (end-start)))

    print(f"Number of reference VNTRs with matching length: {num_ref_len_match}")
    print(f"Number of reference VNTRs with mismatching length: {num_ref_len_mismatch}")

    plt.hist(len_diff, bins=range(0, 100, 1))
    plt.xlabel('Length difference')
    plt.ylabel('Number of VNTRs')
    plt.title('Length difference between reference VNTR (defined in TSV) and vamos VNTR (defined in VCF)')
    plt.savefig(plot)


            
    
if __name__=='__main__':
    
    parser = argparse.ArgumentParser(prog='get-ref-vntr-stats.py', description="Checks compatibility between VNTR coordiante in TSV and VCF")
    parser.add_argument("-vcf", required=True, help="Reference VNTRs (made from contigs)")
    parser.add_argument("-plot", required=True, help="Path to the histogram plot")
    
    options = parser.parse_args()

    run(**vars(options))