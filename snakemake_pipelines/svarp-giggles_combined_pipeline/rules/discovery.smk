# run svarp
rule svarp:
    input:
        reads='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz',
        reads_gzi='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz.gzi',
        reads_fai='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/data/fasta/{sample}.fasta.gz.fai',
        haplotag='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/haplotagging-results/GRCh38/{sample}/{sample}.tsv',
        alignment='/gpfs/project/projects/medbioinf/data/share/globus/1000g-ont/gaf/{sample}.gaf.gz',
        gfa='/gpfs/project/projects/medbioinf/users/spani/files/gfa/1000GP/chm13-90c.r518_tagged.gfa'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/{sample}_svtigs_H1.fa',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/{sample}_svtigs_H2.fa',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/{sample}_svtigs_untagged.fa'
    params:
        outdir='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}'
    log:
        stdout='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/log/{sample}.stdout',
        stderr='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/log/{sample}.stderr'
    conda:
        "../envs/svarp.yml"
    resources:
        mem_total_mb=1024*400,
        runtime_hrs=24*5,
        runtime_min=1
    priority: 2
    shell:
        '''
        module load Python/3.11.3
        module load Minimap2/2.17
        module load SamTools/1.6
        module load gcc/10.2.0
        export PATH=$PATH:/gpfs/project/projects/medbioinf/users/spani/packages/wtdbg2
        /gpfs/project/projects/medbioinf/projects/1000g-ont/svarp/build/svarp -a {input.alignment} -g {input.gfa} --fasta {input.reads} -i {wildcards.sample} --phase {input.haplotag} -o {params.outdir} 2> {log.stderr} 1> {log.stdout}
        module unload gcc/10.2.0
        module unload Minimap2/2.17
        module unload Python/3.11.3
        module unload SamTools/1.6
        '''

# map svtigs with minimap2
rule svtigs_minimap2:
    input:
        ref='/gpfs/project/projects/medbioinf/users/ebler/long-read-1kg/data/1KG_ONT_VIENNA_hg38.fa',
        fasta='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/{sample}_svtigs_{haplotype}.fa'
    output:
        temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/{sample}_{haplotype}.sam')
    wildcard_constraints:
        haplotype='H1|H2'
    conda:
        '../envs/svarp_processing.yml'
    resources:
        mem_total_mb=30000
    threads: 2
    priority:3
    shell:
        'minimap2 -a -x asm5 --cs -r2k -t {threads} {input.ref} {input.fasta}  > {output}'

# sort aligned SAM file
rule sort_svtig_alignment:
    input:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/{sample}_{haplotype}.sam'
    output:
        temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/{sample}_{haplotype}.sorted.bam')
    wildcard_constraints:
        haplotype='H1|H2'
    conda:
        '../envs/svarp_processing.yml'
    threads: 2
    priority:3
    shell:
        'samtools sort -m4G -@{threads} -o {output} {input}'

# index sorted BAM
rule index_sorted_svtig_alignment:
    input:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/{sample}_{haplotype}.sorted.bam'
    output:
        temp('/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/{sample}_{haplotype}.sorted.bam.bai')
    wildcard_constraints:
        haplotype='H1|H2'
    conda:
        '../envs/svarp_processing.yml'
    shell:
        'samtools index {input}'

# run svim-asm to get the vcf
rule run_svim_asm:
    input:
        ref='/gpfs/project/projects/medbioinf/users/ebler/long-read-1kg/data/1KG_ONT_VIENNA_hg38.fa',
        bam1='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/{sample}_H1.sorted.bam',
        bam1_index='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/{sample}_H1.sorted.bam.bai',
        bam2='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/{sample}_H2.sorted.bam',
        bam2_index='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/{sample}_H2.sorted.bam.bai'
    output:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/variants.vcf',
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/sv-lengths.png'
    params:
        outdir='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/svimasm/'
    conda:
        '../envs/svarp_processing.yml'
    priority: 3
    shell:
        'svim-asm diploid --sample {wildcards.sample} {params.outdir} {input.bam1} {input.bam2} {input.ref}'


###### PAV rules from Arda ######

rule prepare_pav_assembly_table:
    input:
        svtigs_h1='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/{sample}_svtigs_H1.fa',
        svtigs_h2='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/{sample}_svtigs_H2.fa',
    output:
        svtigs_h1= '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/pav/svtigs_h1.asm.fa',
        svtigs_h2= '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/pav/svtigs_h2.asm.fa',
        assm_tsv= '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/pav/assemblies.tsv'
    run:
        import pathlib as pl
        import shutil
        # I don't trust Snakemake keeping a sort-order intact
        input_assemblies = [input.svtigs_h1, input.svtigs_h2
                        ]
        output_assemblies = [output.svtigs_h1, output.svtigs_h2
                        ]
        # relative path to working directory for PAV assemblies.tsv
        pav_input_assemblies = [fp.replace('results/', '', 1) for fp in output_assemblies]
        names = [pl.Path(fp).name.split('_')[0] for fp in output_assemblies]
        with open(output.assm_tsv, 'w') as table:
            _ = table.write('\t'.join(['NAME', 'HAP1', 'HAP2']) + '\n')
            _ = table.write('\t'.join([names[0], pav_input_assemblies[0], pav_input_assemblies[1]])  + '\n')
        for infile, outfile in zip(input_assemblies, output_assemblies):
            shutil.copy(infile, outfile)

rule prepare_pav_config:
    input:
        assm_tsv='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/pav/assemblies.tsv',
        ref='/gpfs/project/projects/medbioinf/users/ebler/long-read-1kg/data/1KG_ONT_VIENNA_hg38.fa',
        ref_idx='/gpfs/project/projects/medbioinf/users/ebler/long-read-1kg/data/1KG_ONT_VIENNA_hg38.fa.fai'
    output:
        cfg='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/pav/config.json',
    run:
        import json
        import shutil
        pav_cfg = {
            'assembly_table': input.assm_tsv,
            'reference': output.ref
        }
        with open(output.cfg, 'w') as dump:
            _ = json.dump(pav_cfg, dump, ensure_ascii=True)

rule run_pav:
    input:
        cfg='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/pav/config.json'
    output:
        chk='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/pav/run.complete',
        summary='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/pav/pav.summary',
    params:
        dir='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/pav/'
    threads: 24
    log:
        pav='/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/pav/pav.log'
    benchmark:
        '/gpfs/project/projects/medbioinf/users/spani/results/1000GP/svarp-giggles/chm13-90c.r518/svarp/{sample}/pav/pav_run.rsrc'
    container:
        '/gpfs/project/projects/medbioinf/container/pav_v2.2.4.1.sif'
    shell:
        'snakemake --verbose --jobs {threads} -d {params.dir} -s {output.summary} --rerun-incomplete --keep-incomplete --restart-times 0 &> {log.pav} && touch {output.chk}'