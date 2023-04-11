import gzip
import argparse
import pysam

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", action="store_true", default=False, help="Print debug messages")
    parser.add_argument("--tag-list", required=True, help="Provide read tag list given as whatshap haplotag output.")
    parser.add_argument("--bam", required=True, help="Provide the tagged BAM/CRAM given as whatshap haplotag output/ BAM file used for tagging.")
    parser.add_argument("--output", help="Output path to store the output tagged list with the non-PAR region reads tagged. Default: sys.stdout")

    options = parser.parse_args()
    
    bam = options.bam
    tag_list = options.tag_list
    read_tags = {}
    read_tag_list(tag_list, read_tags)
    
    reader = None
    bam_reader = pysam.AlignmentFile(bam, "rb")
    nonPAR_iter = bam_reader.fetch("chrX", 2781479, 155701383)

    out = options.output
    writer = None
    
    for read in nonPAR_iter:
        start = read.reference_start
        end = read.reference_end
        if start < 2781479:
            continue
        elif end > 155701383:
            continue
        readname = read.query_name
        print(read_tags[readname])

        exit()
        pass

def read_tag_list(file, read):
    
    reader = open(file, 'r')
    while True:
        line = reader.readline()
        if not line:
            break
        if line[0] == "#":
            continue
        readname,haplotype,ps,chr = line.rstrip().split("\t")
        read[readname] = [haplotype, ps, chr]
    reader.close()

def write_to_file(line, writer):
    try:
        writer.write(line)
    except TypeError:
        writer.write(str.encode(line))

def is_file_gzipped(src):
    with open(src, "rb") as inp:
        return inp.read(2) == b'\x1f\x8b'

if __name__ == "__main__":
    main()