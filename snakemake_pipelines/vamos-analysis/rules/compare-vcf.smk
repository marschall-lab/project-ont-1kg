# making the vcf compatible for the comparison pipeline
rule prepare_miller_vcf:
    input:
        vcf='results/miller-comparison/analysis/vamos-miller.sorted.vcf',
        ref='results/miller-comparison/reference-vntrs.bed'
    output:
        'results/miller-comparison/callset-comparison/miller-vcfs/{sample}.vcf'
    shell:
        'python scripts/subsample-vcf.py -vcf {input.vcf} -reference {input.ref} -sample {wildcards.sample} > {output}'

# adding SVLEN and SVTYPE tags to the miller vcf
rule add_tags_miller:
    input:
        'results/miller-comparison/callset-comparison/miller-vcfs/{sample}.vcf'
    output:
        'results/miller-comparison/callset-comparison/miller-vcfs/{sample}.tagged.vcf'
    log:
        'results/miller-comparison/callset-comparison/miller-vcfs/{sample}.tagged.log'
    conda:
        "../envs/comparison.yml"
    shell:
        "cat {input} | python scripts/add-svtags.py 2> {log} 1> {output}"

# extracting the VNTR sites from SVAN annotations and adding genotypes to them
rule prepare_svan_annot:
    input:
        insertions=config['path_to_svan_ins'],
        deletions=config['path_to_svan_del'],
        gts=config['path_to_gts']
    output:
        ins_zip=temp('results/miller-comparison/callset-comparison/svan-vcfs/ins.vcf.gz'),
        del_zip=temp('results/miller-comparison/callset-comparison/svan-vcfs/del.vcf.gz'),
        ins_zip_tbi=temp('results/miller-comparison/callset-comparison/svan-vcfs/ins.vcf.gz.tbi'),
        del_zip_tbi=temp('results/miller-comparison/callset-comparison/svan-vcfs/del.vcf.gz.tbi'),
        ins_gts=temp('results/miller-comparison/callset-comparison/svan-vcfs/ins-gts.vcf.gz'),
        del_gts=temp('results/miller-comparison/callset-comparison/svan-vcfs/del-gts.vcf.gz'),
        ins_gts_tbi=temp('results/miller-comparison/callset-comparison/svan-vcfs/ins-gts.vcf.gz.tbi'),
        del_gts_tbi=temp('results/miller-comparison/callset-comparison/svan-vcfs/del-gts.vcf.gz.tbi'),
        vntr_exp=temp('results/miller-comparison/callset-comparison/svan-vcfs/vntrs-exp.vcf.gz'),
        vntr_con=temp('results/miller-comparison/callset-comparison/svan-vcfs/vntrs-con.vcf.gz'),
        vntr_exp_tbi=temp('results/miller-comparison/callset-comparison/svan-vcfs/vntrs-exp.vcf.gz.tbi'),
        vntr_con_tbi=temp('results/miller-comparison/callset-comparison/svan-vcfs/vntrs-con.vcf.gz.tbi'),
        vntrs='results/miller-comparison/callset-comparison/svan-vcfs/vntrs.vcf'
    conda:
        '../envs/comparison.yml'
    shell:
        '''
        bcftools sort {input.insertions} | bgzip -c > {output.ins_zip}
        tabix -p vcf {output.ins_zip}
        bcftools sort {input.deletions} | bgzip -c > {output.del_zip}
        tabix -p vcf {output.del_zip}
        
        bcftools annotate -c INFO -a {output.ins_zip} {input.gts} | bgzip -c > {output.ins_gts}
        tabix -p vcf {output.ins_gts}
        bcftools annotate -c INFO -a {output.del_zip} {input.gts} | bgzip -c > {output.del_gts}
        tabix -p vcf {output.del_gts}
        
        bcftools view -i "ITYPE_N=\'VNTR\'" {output.ins_gts} | bcftools sort | bgzip -c > {output.vntr_exp}
        tabix -p vcf {output.vntr_exp}
        bcftools view -i "DTYPE_N=\'VNTR\'" {output.del_gts} | bcftools sort | bgzip -c > {output.vntr_con}
        tabix -p vcf {output.vntr_con}
        
        bcftools concat -a {output.vntr_exp} {output.vntr_con} > {output.vntrs}
        '''

# extract samples from the SVAN vntrs
rule extract_subsample_svan:
    input:
        'results/miller-comparison/callset-comparison/svan-vcfs/vntrs.vcf'
    output:
        'results/miller-comparison/callset-comparison/svan-vcfs/{sample}.vcf'
    conda:
        '../envs/comparison.yml'
    shell:
        'bcftools view --samples {wildcards.sample} {input} > {output}'

# adding SVLEN and SVTYPE tags to the SVAN vcf
rule add_tags_svan:
    input:
        'results/miller-comparison/callset-comparison/svan-vcfs/{sample}.vcf'
    output:
        'results/miller-comparison/callset-comparison/svan-vcfs/{sample}.tagged.vcf'
    log:
        'results/miller-comparison/callset-comparison/svan-vcfs/{sample}.tagged.log'
    conda:
        '../envs/comparison.yml'
    shell:
        'cat {input} | python scripts/add-svtags.py 2> {log} 1> {output}'

# intersecting the vcfs
# the order of the vcfs given to the intersection matters
rule intersect_vcfs:
    input:
        'results/miller-comparison/callset-comparison/svan-vcfs/{sample}.tagged.vcf',
        'results/miller-comparison/callset-comparison/miller-vcfs/{sample}.tagged.vcf'
    output:
        tsv="results/miller-comparison/callset-comparison/intersection/{sample}.tsv",
        vcf="results/miller-comparison/callset-comparison/intersection/{sample}.vcf",
        pdf="results/miller-comparison/callset-comparison/intersection/{sample}.pdf",
        plot="results/miller-comparison/callset-comparison/intersection/upsetplot_{sample}.pdf"
    conda:
        "../envs/comparison.yml"
    log:
        intersect="results/miller-comparison/callset-comparison/intersection/intersect-{sample}.log",
        plot="results/miller-comparison/callset-comparison/intersection/plot-{sample}.log"
    params:
        names = ['svan', 'miller'],
        columns = ['in_svan', 'in_miller']
    shell:
        '''
        python scripts/intersect_callsets.py intersect -c {input} -n {params.names} -t {output.tsv} -v {output.vcf} -p {output.pdf} --id-from-vcf &> {log.intersect}
        python scripts/plot-upset.py -t {output.tsv} -o {output.plot} -n {params.columns} &> {log.plot}
        '''