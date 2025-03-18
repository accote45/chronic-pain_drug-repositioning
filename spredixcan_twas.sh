ml python/2.7.16

for i in $(ls /data/elastic_net_models_gtexv8/* | xargs -n1 basename | sed 's/.db//g' | sed 's/.txt.gz//g' | sed 's/en_//g' | grep Brain | sort | uniq); do

./MetaXcan_update2/MetaXcan/software/MetaXcan.py  \
--weight_db_path /data/elastic_net_models_gtexv8/en_${i}.db \
--covariance /data/elastic_net_models_gtexv8/en_${i}.txt.gz \
--gwas_folder chronic_pain \
--gwas_file  chronic_pain-bgen.stats_noindels \
--snp_column SNP \
--effect_allele_column ALLELE1 \
--non_effect_allele_column ALLELE0 \
--se_column SE \
--beta_column BETA \
--pvalue_column P_BOLT_LMM_INF \
--output_file /chronic_pain/gtexv8_${i}.csv

rm -r intermediate
done
