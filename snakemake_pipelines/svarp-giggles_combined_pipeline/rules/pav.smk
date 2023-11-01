###### PAV rules from Arda ######
rule all:
    'pav/run.complete'

rule prepare_pav_assembly_table:
    input:
        svtigs_h1=config['sample']+'_svtigs_H1.fa',
        svtigs_h2=config['sample']+'_svtigs_H2.fa',
    output:
        svtigs_h1='pav/svtigs_h1.asm.fa',
        svtigs_h2='pav/svtigs_h2.asm.fa',
        assm_tsv='pav/assemblies.tsv'
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
        assm_tsv='pav/assemblies.tsv',
        ref='/gpfs/project/projects/medbioinf/users/ebler/long-read-1kg/data/1KG_ONT_VIENNA_hg38.fa',
        ref_idx='/gpfs/project/projects/medbioinf/users/ebler/long-read-1kg/data/1KG_ONT_VIENNA_hg38.fa.fai'
    output:
        cfg='pav/config.json',
    run:
        import json
        import shutil
        pav_cfg = {
            'assembly_table': input.assm_tsv,
            'reference': input.ref
        }
        with open(output.cfg, 'w') as dump:
            _ = json.dump(pav_cfg, dump, ensure_ascii=True)

rule run_pav:
    input:
        cfg='pav/config.json'
    output:
        chk='pav/run.complete'
    threads: 8
    log:
        pav='pav/pav.log'
    benchmark:
        'pav/pav_run.rsrc'
    container:
        '/gpfs/project/projects/medbioinf/container/pav_v2.2.4.1.sif'
    shell:
        'snakemake --verbose --jobs {threads} -d pav/ -s /opt/pav/Snakefile --rerun-incomplete --keep-incomplete --restart-times 0 &> {log.pav} && touch {output.chk}'