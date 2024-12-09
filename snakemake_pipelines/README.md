# Notes of the Snakemake Pipelines

- The snakemake pipelines in this directory are prepared for the HILBERT cluster of HHU, Duesseldorf.
- The following files are specialized for HILBERT
    - <pipeline name>/hhu_hilbert.py
    - <pipeline name>/profile/hilbert_cluster.json
    - <pipeline name>/envs/*.yaml (channels set to local mirrors in HILBERT)
- Each pipeline is explained in <pipeline name>/README.md
- The pipelines given are interconnected. Some pipelines require input from other pipelines to work. Refer to the <pipeline name>/config.yaml to determine the inputs.
- Once the release files are available on IGSR, the pipeline config files will be updated to show which files can be downloaded and used.