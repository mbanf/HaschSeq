#!/bin/bash
#run by ./halliftover_of_SNPs_Mo17_for_BZR1_paper.sh genotype 
EXPECTED_ARGS=1
E_BADARGS=1
if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: ./halliftover_of_SNPs_Mo17_for_BZR1_paper.sh genotype "
  exit $E_BADARGS
fi


g=$1

export PATH="/netscratch/dep_psl/grp_frommer/Thomas/bin/bedops/bin:$PATH"



 cd /netscratch/dep_psl/grp_frommer/Thomas/Results/HASCh_BZR1/custom/New_counts_paper
##tsv file contains 0-based positions

echo "halLiftover"
 
halLiftover  --hdf5InMemory /netscratch/dep_psl/grp_frommer/Thomas/Results/HybMoa_0819_WWvsDS/cactus/${g}/B73v5_${g}.hal B73 Mo17.SNPs.for_liftover.bed ${g} ${g}.B73.hallifted.SNPs.bed






