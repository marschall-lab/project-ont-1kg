"""
rule prepare_svan_annot:
    input:
        insertions=config['path_to_svan_ins'],
        deletions=config['path_to_svan_del'],
        complex_insertions=config['path_to_svan_complex_ins'],
        complex_deletions=config['path_to_svan_complex_del'],
        gts=config['path_to_gts']
    output:
        ins_zip=temp('results/annotated-gts/svan-vcfs/ins.vcf.gz'),
        del_zip=temp('results/annotated-gts/svan-vcfs/del.vcf.gz'),
        ins_zip_tbi=temp('results/annotated-gts/svan-vcfs/ins.vcf.gz.tbi'),
        del_zip_tbi=temp('results/annotated-gts/svan-vcfs/del.vcf.gz.tbi'),
        complex_ins_zip=temp('results/annotated-gts/svan-vcfs/complex_ins.vcf.gz'),
        complex_del_zip=temp('results/annotated-gts/svan-vcfs/complex_del.vcf.gz'),
        complex_ins_zip_tbi=temp('results/annotated-gts/svan-vcfs/complex_ins.vcf.gz.tbi'),
        complex_del_zip_tbi=temp('results/annotated-gts/svan-vcfs/complex_del.vcf.gz.tbi'),
        svan_full=temp('results/annotated-gts/svan-vcfs/svan_full.vcf.gz'),
        svan_full_tbi=temp('results/annotated-gts/svan-vcfs/svan_full.vcf.gz.tbi'),
        gts_annotated='results/annotated-gts/annotated_gts.vcf'
    conda:
        '../envs/analysis.yml'
    shell:
        '''
        bcftools sort {input.insertions} | bgzip -c > {output.ins_zip}
        tabix -p vcf {output.ins_zip}
        bcftools sort {input.deletions} | bgzip -c > {output.del_zip}
        tabix -p vcf {output.del_zip}
        bcftools sort {input.complex_insertions} | bgzip -c > {output.complex_ins_zip}
        tabix -p vcf {output.complex_ins_zip}
        bcftools sort {input.complex_deletions} | bgzip -c > {output.complex_del_zip}
        tabix -p vcf {output.complex_del_zip}
        
        bcftools concat -a {output.complex_del_zip} {output.complex_ins_zip} {output.del_zip} {output.ins_zip} | bcftools sort | bgzip -c > {output.svan_full}
        tabix -p vcf {output.svan_full}

        bcftools annotate -c INFO -a {output.svan_full} {input.gts} > {output.gts_annotated}
        '''
"""

# extracting the VNTR sites from SVAN annotations and adding genotypes to them
rule prepare_svan_annot:
    input:
        svan=config['path_to_svan'],
        gts=config['path_to_gts']
    output:
        gts_annotated='results/annotated-gts/annotated_gts.vcf'
    conda:
        '../envs/analysis.yml'
    shell:
        'bcftools annotate -c INFO -a {input.svan} {input.gts} > {output.gts_annotated}'

# extract samples from the SVAN vntrs
rule extract_subsample_svan:
    input:
        'results/annotated-gts/annotated_gts.vcf'
    output:
        'results/annotated-gts/sample-wise-gts/{sample}.vcf'
    conda:
        '../envs/analysis.yml'
    shell:
        'bcftools view --samples {wildcards.sample} {input} | bcftools view --min-ac 1 > {output}'

# get stats for NUMT, PSD, etc each sample
rule get_sample_statistics:
    input:
        'results/annotated-gts/sample-wise-gts/{sample}.vcf'
    output:
        'results/annotated-gts/sample-wise-gts/{sample}.counts'
    shell:
        'cat {input} | grep -e PSD -e NUMT -e INV | wc -l > {output}'