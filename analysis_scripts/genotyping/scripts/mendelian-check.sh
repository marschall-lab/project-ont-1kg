trioPaternal=("NA19818" "HG01256" "HG00418" "NA19128" "NA12889" "NA12891")
trioMaternal=("NA19819" "HG01257" "HG00419" "NA19127" "NA12890" "NA12892")
trioChild=("NA19828" "HG01258" "HG00420" "NA19129" "NA12877" "NA12878")

# on the unfiltered genotypes

for i in 0 1 2 3 4 5; do
    python scripts/check-mendelian-consistency.py -map augmented-graph/map-bub_ids-to-allele_ids.tsv -vcf augmented-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz -child ${trioChild[$i]} -father ${trioPaternal[$i]} -mother ${trioMaternal[$i]} -out augmented-graph/mendelian-consistency/ > augmented-graph/mendelian-consistency/${trioChild[$i]}.log;
    python scripts/check-mendelian-consistency.py -map original-graph/map-bub_ids-to-allele_ids.tsv -vcf original-graph/multisample-genotypes/biallelic.ac0-filtered.vcf.gz -child ${trioChild[$i]} -father ${trioPaternal[$i]} -mother ${trioMaternal[$i]} -out original-graph/mendelian-consistency/ > original-graph/mendelian-consistency/${trioChild[$i]}.log;
done;

# on the filtered genotypes

for i in 0 1 2 3 4 5; do
    python scripts/check-mendelian-consistency.py -map augmented-graph/map-bub_ids-to-allele_ids.tsv -vcf genotype-filtering/augmented/biallelic.filtered.all-samples.vcf.gz -child ${trioChild[$i]} -father ${trioPaternal[$i]} -mother ${trioMaternal[$i]} -out genotype-filtering/augmented/mendelian-consistency/ > genotype-filtering/augmented/mendelian-consistency/${trioChild[$i]}.log;
    python scripts/check-mendelian-consistency.py -map original-graph/map-bub_ids-to-allele_ids.tsv -vcf genotype-filtering/original/biallelic.filtered.all-samples.vcf.gz -child ${trioChild[$i]} -father ${trioPaternal[$i]} -mother ${trioMaternal[$i]} -out genotype-filtering/original/mendelian-consistency/ > genotype-filtering/original/mendelian-consistency/${trioChild[$i]}.log;
done;