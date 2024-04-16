# Description of the Pipeline

- The long read DNA sequences of the sample cohort can be tagged as "H1" or "H2" using `whatshap`.
- The haplotags were determined using the NYGC phased panels (PMID: 36055201) where 967 samples intersect between their study and our study.
- Due to WhatsHap's limitation at the pseudo-autosomal regions of the X and Y chromosomes, certain post processing steps were added to accurately tag the reads.