#!/bin/bash
#run by ./hallifted_bed_to_count_STAR_new_only_uni_EG80.BZR1_paper_07_2022.sh genotype 
EXPECTED_ARGS=1
E_BADARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: ./hallifted_bed_to_count_STAR_new_only_uni_EG80.BZR1_paper_07_2022.sh genotype "
  exit $E_BADARGS
fi


g=$1



 cd  /netscratch/dep_psl/grp_frommer/Thomas/Results/HASCh_BZR1/custom/New_counts_paper

 
echo "bedtool mapping"


bedtools map -a B73.${g}.ID.hallifted.SNPs.bed -b /netscratch/dep_psl/grp_frommer/Thomas/Results/HASCh_BZR1/bamCoverage/EG_norm/BZR1_6.q255.RPGC.exScal.bin1.sm0.rm.bedgraph -g /netscratch/dep_psl/grp_frommer/Michael_Thomas/Genomes/Zea_mays/diploid/${g}/ref_B73${g}.fasta.size.new.txt -c 4 > B73.${g}.WW.counts.uni_EG80.BZR1.bed

bedtools map -a ${g}.B73.ID.hallifted.SNPs.bed -b /netscratch/dep_psl/grp_frommer/Thomas/Results/HASCh_BZR1/bamCoverage/EG_norm/BZR1_6.q255.RPGC.exScal.bin1.sm0.rm.bedgraph -g /netscratch/dep_psl/grp_frommer/Michael_Thomas/Genomes/Zea_mays/diploid/${g}/ref_B73${g}.fasta.size.new.txt -c 4 > ${g}.B73.WW.counts.uni_EG80.BZR1.bed



##This will give us : Chr start stop ID Value


