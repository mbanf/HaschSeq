#!/bin/bash
#run by ./deduplicate_and_translate_hallifted_SNPs_BZR1_paper.sh genotype 
EXPECTED_ARGS=1
E_BADARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: ./deduplicate_and_translate_hallifted_SNPs_BZR1_paper.sh genotype "
  exit $E_BADARGS
fi


g=$1

export PATH="/netscratch/dep_psl/grp_frommer/Thomas/bin/bedops/bin:$PATH"


 cd /netscratch/dep_psl/grp_frommer/Thomas/Results/HASCh_BZR1/custom/New_counts_paper


echo "get duplicates"


gawk -v OFS='\t' '{print $4}' ${g}.B73.hallifted.SNPs.bed | sort | uniq -d > ${g}.B73.hallifted.SNPs.dups_B73coord.txt

echo "remove duplicates"


gawk -v OFS='\t' 'NR==FNR{a[$1]=$1; next} !($4 in a){print $0}'  ${g}.B73.hallifted.SNPs.dups_B73coord.txt ${g}.B73.hallifted.SNPs.bed > ${g}.B73.hallifted.SNPs.clean2.bed



echo "renaming chromosomes"

gawk -v OFS='\t' -v ALTn=$g '{print ALTn"-"$1,$2,$3,$4}' ${g}.B73.hallifted.SNPs.clean2.bed | sortBed -g /netscratch/dep_psl/grp_frommer/Michael_Thomas/Genomes/Zea_mays/diploid/${g}/ref_B73${g}.fasta.size.new.txt > ${g}.B73.ID.hallifted.SNPs.bed

gawk -v OFS='\t'  '{{split($4,var,"."); print "B73-"var[1],var[2]-1,var[2],$1"."$3"."var[4]"."var[3]}}' ${g}.B73.hallifted.SNPs.clean2.bed  | sortBed -g /netscratch/dep_psl/grp_frommer/Michael_Thomas/Genomes/Zea_mays/diploid/${g}/ref_B73${g}.fasta.size.new.txt  > B73.${g}.ID.hallifted.SNPs.bed
