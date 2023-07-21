'''
Python script to read the assembly-to-graph alignmnets and recognize alleles from that.
Output a GFA file with P and W lines.
P lines - Reference Path
W lines - Haplotype Path
'''

import argparse
import os
import logging
import resource
from collections import defaultdict, namedtuple, abc
import gzip
import re
from whatshap.timer import StageTimer
import pysam

logger = logging.getLogger(__name__)
timers = StageTimer()
VariantRecord = namedtuple('VariantRecord', ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])


class VCFWriter:
    def __init__(self, path, bgzip):
        self.bgzipped = bgzip
        if bgzip:
            self.writer = pysam.libcbgzf.BGZFile(path, 'wb')
        else:
            self.writer = open(path, "w")
    
    def write(self, line):
        line += '\n'
        if not self.bgzipped:
            self.writer.write(line)
        else:
            self.writer.write(str.encode(line))
    
    def close(self):
        self.writer.close()


class Genotype:
    def __init__(self):
        self.maternal = None
        self.paternal = None
    
    def to_string(self):
        return '%s|%s'%(self.paternal, self.maternal)


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
    parser.add_argument("--assembly-list", required=True, help="Text file with the list of assembly-to-graph GAF files. The GAF files should be sorted.")
    parser.add_argument("--output", required=True, help="Output VCF File path")
    parser.add_argument("--bgzip", action='store_true', help="Flag to bgzip the outputs (VCF and GFA). Can only be given with --output.")
    
    options = parser.parse_args()
    setup_logging(options.debug)
    
    validate_arguments(options)

    annotate(options)

    memory_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    logger.info("\nMemory Information")
    logger.info("\tMaximum memory usage:       %.3f GB", memory_kb / 1e6)
    logger.info("\nTime Summary (in seconds)")
    logger.info("\tTime to parse GFA file:     %.3f"%(timers.elapsed("read_gfa")))
    logger.info("\tTime to read assemblies:    %.3f"%(timers.elapsed("read_assemblies")))
    logger.info("\tTime to write VCF file:     %.3f"%(timers.elapsed("write_vcf")))
    


def annotate(options):
    '''
    Options:
        1. gfa: Input GFA that needs to be updated.
        2. assembly_list: a file containing the list of file paths to assembly-to-graph GAFs.
        3. outgfa: path to output GFA file. Will contain the original GFA info with P and W lines.
        4. outvcf: path to output VCF file. This will contain the allele traversal paths.
        5. bgzip: flag to bgzip the output files.
    
    What to keep in the VCF file?
        1. ID will be the start and end node with the orientation
        2. INFO will contain Allele Traversals (AT), Allele Frequencies (AF), Number of Samples (NS). (Anything else?)
        3. FORMAT will contain the GT. (Anything else?)
        4. FILTER can have ns80 (present in 80% of samples). (Anything else?)
        5. contigs have to be added later after the gfa has been processed.
    '''
    
    Node = namedtuple("Node", ["SN", "SO", "SR", "LN", "BO", "NO", "Seq"], defaults=[-1, -1, None, -1, -1, -1, None])
    Edge = namedtuple("Edge", ["next", "prev"], defaults=[None, None])
    nodes = defaultdict(lambda: Node())
    edges = defaultdict(lambda: Edge())
    with timers("read_gfa"):
        read_gfa(options.gfa, nodes, edges)
    
    scaffold_nodes = [node_id for node_id, value in sorted(nodes.items(), key=lambda item: item[1].BO) if value.NO==0]
    with timers("read_assemblies"):
        variants, haplotype_list, contigs = read_assemblies(options.assembly_list, scaffold_nodes, nodes)
    ref_alleles = get_reference_alleles(edges, list(variants.keys()))
    
    writer = VCFWriter(options.output, options.bgzip)
    with timers("write_vcf"):
        write_vcf(writer, variants, ref_alleles, haplotype_list, contigs, options, nodes)
    writer.close()
    

def read_gfa(gfa, node, edges):
    
    logger.info("\nParsing GFA file and reading sort key information")
    logger.info("\tRequires the following tags: SN, SO, SR, BO, and NO. Either the sequence or the LN tag has to be provided also.")
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
        node[fields[1]] = node[fields[1]]._replace(BO=BO, NO=NO, SO=SO, SR=SR, SN=SN, LN=LN, Seq=fields[2])
        if node[fields[1]].LN == -1:
            node[fields[1]] = node[fields[1]]._replace(LN = len(fields[2]))
    logger.info("\t%d total Nodes Processed"%(total_nodes))
    logger.info("\t%d nodes with tags"%(tagged_nodes))
    
    reader.seek(0)
    while True:
        line = reader.readline()
        if not line:
            break
        if line[0] != 'L':
            continue
        fields = line.split('\t')
        n1 = fields[1]
        o1 = fields[2]
        n2 = fields[3]
        o2 = fields[4]
        if o1 != "+" or o2 != "+" or node[n1].SR != 0 or node[n2].SR != 0:
            continue
        if node[n2].SO != (node[n1].SO + node[n1].LN):
            continue
        new_next = n2
        new_previous = edges[n1].prev
        edges[n1] = edges[n1]._replace(next=new_next, prev=new_previous)
        new_next = edges[n2].next
        new_previous = n1
        edges[n2] = edges[n2]._replace(next=new_next, prev=new_previous)


    reader.close()

def read_assemblies(assembly_list, scaffold_nodes, nodes):
    variants = {}
    bubbles = []
    haplotype_list = []
    contigs = set()
    # Generate the variant INFO (the starting and end node of a bubble)
    for i in range(len(scaffold_nodes) - 1):
        # Take into consideration the chromosome switch
        if nodes[scaffold_nodes[i]].SN != nodes[scaffold_nodes[i+1]].SN:
            continue
        contigs.add(nodes[scaffold_nodes[i]].SN)
        bubble = ">%s>%s"%(scaffold_nodes[i], scaffold_nodes[i+1])
        bubbles.append(bubble)
        variants[bubble] = {}
    logger.info("\nNumber of possible bubbles = %d"%(len(bubbles)))
    gaf_files = []
    with open(assembly_list, 'r') as f:
        for line in f:
            gaf_files.append(line.rstrip())
    for gaf in gaf_files:
        logger.info("\nProcessing Assembly File %s"%(gaf))
        haplotype = gaf.split("/")[-1][:-4]
        haplotype_list.append(haplotype)
        find_variant_alleles(gaf, variants, nodes, haplotype)

    haplotype_list.sort()
    contigs = list(contigs)
    contigs.sort()    
            
    return variants, haplotype_list, contigs

def find_variant_alleles(gaf, variants, nodes, haplotype):
    if is_file_gzipped(gaf):
        reader = gzip.open(gaf, 'rt')
    else:
        reader = open(gaf, 'r')
    counts = defaultdict(lambda: 0)
    while True:
        offset = reader.tell()
        alignment = reader.readline()
        if not alignment:
            break
        #if offset != 30432:
        #    continue
        #if (counts['total']+1) % 1 == 0:
        #    logger.info('\tProcessed %d alignments'%(counts['total']))
        alignment = alignment.split('\t')
        logger.debug("\tProcessing contig %s"%(alignment[0]))
        path = list(filter(None, re.split('(>)|(<)', alignment[5])))
        rv = False
        confused = False
        orient = None
        orient_list = []
        scaffold_index = []
        scaffold_list = []
        counts['Total Number of Alignments'] += 1
        # Reading the entire alignment once to check if the alignment is in the reverse direction.
        # Also storing the scaffold node indices to quickls access the variant bubbles.
        for n, p in enumerate(path):
            if p in ['>', '<']:
                orient = p
                continue
            if nodes[p].NO != 0:
                continue
            scaffold_index.append(n)
            scaffold_list.append(p)
            orient_list.append(orient)
        if len(orient_list) == 0:
            counts['Number of Alignments with no Scaffold Nodes'] += 1
            logger.debug("\t\tSkipping contig since no scaffold nodes were found.")
            continue
        
        # Assign reversed if there are more scaffold nodes with < direction
        # Overall orientation has been reversed
        if orient_list.count('>') < orient_list.count('<'):
            rv = True
            logger.debug('\t\tContig alignment identified as reversed.')
            counts['reversed'] += 1
            scaffold_index.reverse()
            new_scaffold_index = []
            for i in scaffold_index:
                new_scaffold_index.append(len(path) - i)
            scaffold_index = new_scaffold_index
            orient_list.reverse()
            new_orient_list = []
            for n,o in enumerate(orient_list):
                if o == ">":
                    new_orient_list.append("<")
                elif o == "<":
                    new_orient_list.append(">")
            orient_list = new_orient_list
            path = reverse_path(path)
        # Assign confused if there are scaffold nodes in both directions.
        breakpoint = []
        if orient_list.count('>') != 0 and orient_list.count('<') != 0:
            # Break up the path into forward and reverse segments.
            logger.debug('\t\tContig alignment contains inverted scaffold nodes.')
            counts['Number of Alignment with Scaffold Nodes in both orientations'] += 1
            orient = None
            for n,o in enumerate(orient_list):
                if orient == None:
                    orient = o
                    breakpoint.append(orient)
                    continue
                if o != orient:
                    breakpoint.append(n)
                    breakpoint.append(o)
                    orient = o
            breakpoint.append(len(orient_list))
        else:
            breakpoint = [">", len(orient_list)]

        start_path = scaffold_index[0]
        start_index = 0
        for n_bp, b in enumerate(breakpoint):
            if b in ['<', '>']:
                orient  = b
                continue
            end_path = scaffold_index[b-1]
            subpath = path[start_path-1:end_path+1]
            subindex = scaffold_index[start_index:b]
            # If there is only one scaffold node, then I can conisder that to be an inversion event with the inverted allele 
            # Cannot be supported since the inverted node is a scaffold node and hence variant cannot be defined between two nodes that are not consecutive
            # >s1<s2>s3: In this case, since all at scaffold nodes, the variants defined are >s1>s2 and >s2>s3. >s1>s3 is not allowed.
            # Instead now I just ignore these short inversions.
            if len(subindex) == 1:
                logger.debug('\t\tContig alignment contains single scaffold node breakpoint.')
            else:
                diff = subindex[0]
                subindex = [s-diff+1 for s in subindex]
                if orient == '<':
                    subpath = reverse_path(subpath)
                    subindex.reverse()
                    new_scaffold_index = []
                    for i in subindex:
                        new_scaffold_index.append(len(subpath) - i)
                    subindex = new_scaffold_index
                

            s = 0
            e = 1
            while e < len(subindex):
                sn = subpath[subindex[s]]
                en = subpath[subindex[e]]
                bub = ">%s>%s"%(sn,en)
                assert(nodes[sn].NO == 0)
                assert(nodes[en].NO == 0)
                var = subpath[subindex[s]-1:subindex[e]+1]
                try:
                    variants[bub][haplotype].add(''.join(var))
                except KeyError:
                    variants[bub][haplotype] = set()
                    variants[bub][haplotype].add(''.join(var))
                s=e
                e=s+1

            # Updating the breakpoint start variables
            try:
                start_path = scaffold_index[b]
                start_index=b
            except IndexError:
                assert b==len(scaffold_index)
    
    for key, value in variants.items():
        try:
            var = value[haplotype]
            if len(var) > 1:
                counts['Number of Variant Bubbles with multiple Allele Traversals found'] += 1
        except KeyError:
            counts['Number of Variant Bubbles not found in the alignments'] += 1
    for key, value in counts.items():
        logger.info("\t%s: %d"%(key, value))
    reader.close()


def reverse_path(path):
    new_path = []
    for p in path:
        if p == '<':
            orient = ">"
            continue
        elif p == '>':
            orient = "<"
            continue
        new_path.insert(0,p)
        new_path.insert(0,orient)
    return new_path


def write_vcf(writer, variants, ref_alleles, haplotypes, contigs, options, nodes):
    write_header(writer, contigs, haplotypes, options)
    write_records(writer, variants, ref_alleles, haplotypes, nodes)
    pass

def write_header(writer, contigs, haplotypes, options):
    writer.write("##fileformat=VCFv4.2")
    writer.write("##reference=%s"%(os.path.abspath(options.gfa)))
    writer.write("##script=%s"%(os.path.abspath(__file__)))
    writer.write("##phased=True")
    
    #Writing INFO Header
    writer.write('##INFO=<ID=CONFLICT,Number=.,Type=String,Description="Assembly names for which there are multiple conflicting allele traversals">')
    writer.write('##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">')
    writer.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">')
    writer.write('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">')
    writer.write('##INFO=<ID=AT,Number=R,Type=String,Description="Allele Traversal as path in graph">')
    
    #Writing FILTER Header
    writer.write('##FILTER=<ID=PASS,Description="All filters passed">')
    writer.write('##FILTER=<ID=NS80,Description="Number of Samples with data less than 80%">')
    
    #Writing FOMRAT Header
    writer.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')

    #Writing Contig Header
    for contig in contigs:
        writer.write(contig_header_to_string(contig))
    
    samples = list(set([x.split(".")[0] for x in haplotypes]))
    writer.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s"%('\t'.join(samples)))


def write_records(writer, variants, ref_alleles, haplotypes, nodes):

    samples = list(set([x.split(".")[0] for x in haplotypes])) 
    for bub in variants.keys():
        variant = variants[bub]
        ref_path = ref_alleles[bub]
        alleles = {ref_path: 0}
        ns = 0
        count = 1
        allele_count = {0: 0}
        conflict = []
        genotypes = []
        for sample in samples:
            #TODO: Hardcoded the identifier for the maternal and paternal haplotypes
            mat = sample+'.2'
            pat = sample+'.1'
            gen = Genotype()
            try:
                mat_allele = variant[mat]
                if len(mat_allele) > 1:
                    mat_allele = ['.']
                    conflict.append(mat)
            except KeyError:
                mat_allele = ['.']
            try:
                pat_allele = variant[pat]
                if len(pat_allele) > 1:
                    pat_allele = ['.']
                    conflict.append(pat)
            except KeyError:
                pat_allele = ['.']
            assert len(mat_allele) == 1 and len(pat_allele) == 1
            mat_allele = list(mat_allele)[0]
            pat_allele = list(pat_allele)[0]
            pat_anum = '.'
            if pat_allele != '.':
                try:
                    pat_anum = alleles[pat_allele]
                    allele_count[pat_anum] += 1
                except KeyError:
                    pat_anum = count
                    alleles[pat_allele] = count
                    count += 1
                    allele_count[pat_anum] = 1
            gen.paternal = pat_anum
            mat_anum = '.'
            if mat_allele != '.':
                try:
                    mat_anum = alleles[mat_allele]
                    allele_count[mat_anum] += 1
                except KeyError:
                    mat_anum = count
                    alleles[mat_allele] = count
                    count += 1
                    allele_count[mat_anum] = 1
            gen.maternal = mat_anum
            genotypes.append(gen.to_string())
            #Determining number of available samples and conflicting samples
            if pat_anum != '.' and mat_anum != '.':
                ns += 1
        #Determing filter
        if ns/len(samples) < 0.8:
            filter = ['NS80']
        else:
            filter = ['PASS']
        #Determine alternate allele count
        ac = {}
        for i in range(len(allele_count)):
            ac[i] = 0
        for gen in genotypes:
            g = gen.split('|')
            try:
                ac[int(g[0])] += 1
            except ValueError:
                pass
            try:
                ac[int(g[1])] += 1
            except ValueError:
                pass
        ac.pop(0)
        ac = list(ac.values)
        #Determine allele traversals
        at = []
        for key,_ in alleles.items():
            at.append(key)
        seq = get_sequences(at, nodes)
        #Determine allele frequencies
        #TODO: What should the denominator be? The number of available alleles or 88?
        af = []
        tot = 0
        for key, value in allele_count.items():
            if key == 0:
                tot += value
                continue
            af.append(value)
            tot += value
        af = [a/tot for a in af]
        nd = nodes[bub.split(">")[1]]
        chr = nd.SN
        pos = nd.SO + nd.LN
        id = bub
        ref = seq[0]
        alt = seq[1:]
        qual = 60
        info={"CONFLICT": conflict, "AC": ac, "AF": af, "NS": ns, "AT": at}
        writer.write(variant_record_to_string(chr, pos, id, ref, alt, qual, filter, info, genotypes))
        

def get_reference_alleles(edges, bubbles):
    ref_alleles = {}
    allele = ''
    for bubble in bubbles:
        b = bubble.split(">")
        n1 = b[1]
        n2 = b[-1]
        allele = ''
        allele += '>'+n1
        n = n1
        while n != n2:
            n = edges[n].next
            allele += '>'+n
        ref_alleles[bubble] = allele
    
    return ref_alleles


def get_sequences(at, nodes):
    alleles = []
    for traversal in at:
        path = list(filter(None, re.split('(>)|(<)', traversal)))
        seq = ''
        orient = None
        for nd in path[2:-2]:
            if nd in ['>', '<']:
                orient = nd
                continue
            if orient == '>':
                seq += nodes[nd].Seq
            elif orient == '<':
                seq += reverse_complement(nodes[nd].Seq)
        if seq == '':
            seq = '*'
        alleles.append(seq)
    return alleles

def reverse_complement(seq):
    seq = seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
    seq = seq.upper()
    # reverse strand
    seq = seq[::-1]
    return seq

def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b'\x1f\x8b'

def write_to_file(line, writer):
    try:
        writer.write(line)
    except TypeError:
        writer.write(str.encode(line))

def format_header_to_string(id, number, type, description):
    return '##FORMAT=<ID=%s,Number=%s,Type=%s,Description="%s">'%(id, number, type, description)

def info_header_to_string(id, number, type, description):
    return '##INFO=<ID=%s,Number=%s,Type=%s,Description="%s">'%(id, number, type, description)

def filter_header_to_string(id, description):
    return '##FILTER=<ID=%s,Description="%s">'%(id, description)

def contig_header_to_string(id):
    return '##contig=<ID=%s>'%(id)

def variant_record_to_string(chr, pos, id, ref, alt, qual, filter, info, genotypes):
    for key, value in info.items():
        if isinstance(value, abc.Iterable):
            info[key] = ','.join([str(a) for a in value])
        else:
            info[key] = str(value)
    return '%s\t%d\t%s\t%s\t%s\t%d\t%s\t%s\tGT\t%s'%(chr, pos, id, ref, ','.join(alt), qual, ','.join(filter), ';'.join(['%s=%s'%(key,values) for key, values in info.items()]), '\t'.join(genotypes))


if __name__ == "__main__":
    main()
