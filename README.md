# 1000 Genomes ONT Sequencing
## Project overview
![overview_figure](https://github.com/marschall-lab/project-ont-1kg/blob/main/figures/ont-1kg-overview.png)

## Source of graph
```
$ wget https://zenodo.org/record/6499594/files/chm13-90c.r518.gfa.gz
$ gunzip chm13-90c.r518.gfa.gz
$ md5sum chm13-90c.r518.gfa
566c961c63e95b338b6c08b6737d2c93  chm13-90c.r518.gfa
$ awk 'BEGIN {OFS="\t"} /^S/ {$3="*"} {print}' chm13-90c.r518.gfa > chm13-90c.r518.noseq.gfa
$ md5sum chm13-90c.r518.noseq.gfa
7abf0a440392b4936bc4cffa890130a3  chm13-90c.r518.noseq.gfa
```
