# Description of the Pipeline

- These are phasing experiments performed using WhatsHap to check the quality of the reads.
- The NYGC study (PMID: 36055201) produced raw genotypes along with the final phased panels. We attempted to phase the raw genotypes using various strategies and consequently compared it against the phased panels to investigate concordance.
- Three phasing strategies were applied:
    - Trio-based phasing: 6 trios are present in the sample set and the genotypes of the trios were used for phasing.
    - Longread-based phasing: For all the samples, phasing was done with only the longread data.
    - Trio Longread phasing: For the 6 trios, the trio data and the longread data was used together to phase.
- The three strategies are then compared to the statistical phasing.