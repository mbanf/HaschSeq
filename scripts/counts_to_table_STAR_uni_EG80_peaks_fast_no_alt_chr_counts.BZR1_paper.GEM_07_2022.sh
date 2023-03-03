#!/bin/bash
#run by ./counts_to_table_STAR_uni_EG80_peaks_fast_no_alt_chr_counts.BZR1_paper.GEM_07_2022.sh genotype 
EXPECTED_ARGS=1
E_BADARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: ./counts_to_table_STAR_uni_EG80_peaks_fast_no_alt_chr_counts.BZR1_paper.GEM_07_2022.sh genotype "
  exit $E_BADARGS
fi


g=$1


 cd /netscratch/dep_psl/grp_frommer/Thomas/Results/HASCh_BZR1/custom/New_counts_paper

#Start with : Chr start stop ID Value

#To add peak information


echo "sort count files"

for tr in WW; do gawk -v OFS='\t' -v g=$g '{if($4!~"scaf-alt"){print $0} else{gsub("scaf-alt","scaf",$4);print $1,$2,$3,$4,"het"}}' B73.${g}.${tr}.counts.uni_EG80.BZR1.bed | sort -k1,1 -k2,2n  > B73.${g}.${tr}.counts.uni_EG80_sorted.BZR1.bed; done

for tr in WW;  do gawk -v OFS='\t' -v g=$g '{if($1!~"scaf-alt"){print $0} else{gsub("scaf-alt","scaf",$1);print $1,$2,$3,$4,"het"}}' ${g}.B73.${tr}.counts.uni_EG80.BZR1.bed | sort -k1,1 -k2,2n  > ${g}.B73.${tr}.counts.uni_EG80_sorted.BZR1.bed; done



#What we want: Chr(B73)  STop=POS(B73)   ID(NAMchr_NAMPOS)   REFal(B73)   ALTal(NAM)   Count(Either B73 or NAM)    and we want this B73 coordinate sorted
echo "intersect peaks and sort"

for tr in WW; do intersectBed -a B73.${g}.${tr}.counts.uni_EG80_sorted.BZR1.bed -b /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/custom/counts_STAR_80_EG/${g}/new_04_2021/BZR1_GEM/BZR1.peaks.names.sorted.bed -wa -wb -sorted -loj| gawk -v OFS='\t' '{split($4,var,"\."); print $1,$3,var[1]":"var[2],var[4], var[3], $5, $9}' >  B73.${g}.${tr}.ID.count.sort.uni_EG80.BZR1.csv; done

for tr in WW; do intersectBed -a ${g}.B73.${tr}.counts.uni_EG80_sorted.BZR1.bed -b /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/custom/counts_STAR_80_EG/${g}/new_04_2021/BZR1_GEM/BZR1.peaks.names.sorted.bed -wa -wb -sorted -loj|gawk -v OFS='\t' '{split($4,var,"\."); print var[1],var[2],$1":"$3,var[3], var[4], $5, $9}' | sort -k1,1 -k2,2n  > ${g}.B73.${tr}.ID.count.sort.uni_EG80.BZR1.csv; done


echo "preparing postfrequency files"

for tr in WW; do paste B73.${g}.${tr}.ID.count.sort.uni_EG80.BZR1.csv ${g}.B73.${tr}.ID.count.sort.uni_EG80.BZR1.csv | gawk -v OFS='\t' '{if($2==$9){if($6+$13==0){print "B73-"substr($1,5,length($1)-4),$2,$10,$4,$5,$6,$13,$7,$14,"n.r."} else {print "B73-"substr($1,5,length($1)-4),$2,$10,$4,$5,$6,$13,$7,$14,$6/($6+$13)}}}' > B73.${g}.${tr}.PF.uni_EG80.BZR1.csv; done


echo "adding genotypes and depth "

CW=`awk -v OFS="\t" 'BEGIN{n=1} {if($4>0 && $4<n){n=$4}} END { print n }' /netscratch/dep_psl/grp_frommer/Thomas/Results/HASCh_BZR1/bamCoverage/EG_norm/BZR1_6.q255.RPGC.exScal.bin1.sm0.rm.bedgraph`;

for tr in WW ; do gawk -v OFS='\t' -v g=$g -v CW=${CW} '{if($8~"Peak" || $9~"Peak"){print $1,$2,$4,$5,$3,"1/1",$6,$7,$8,$9,$10,int(($6/CW)+0.5),int(($7/CW)+0.5)} else {print $1,$2,$4,$5,$3,"1/1",$6,$7,$8,$9,"n.p."$10,int(($6/CW)+0.5),int(($7/CW)+0.5)}}' B73.${g}.${tr}.PF.uni_EG80.BZR1.csv  > B73.${g}.${tr}.PF.uni_EG80.GT.RN.BZR1.GEM.csv; done

