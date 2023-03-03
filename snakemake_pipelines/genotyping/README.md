# WhatsHap genotyping

- Preprocessing the samples (In the current version, I am deleting the INFO field of the VCF. And then I am adding it back later.)
	- Need to fix this the next time I run for all the samples.
- Genotype all the samples.
- The genotyping is done for each sample separately. So all the VCFs are then merged for certain evaluations like HWE calculation, het-hom ratio, etc.
- Sample-wise genotype VCFs are used for genotype concordance calculation. It is checked against PanGenie-genotyped VCF of the 1000GP samples.
