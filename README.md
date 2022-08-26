# HaschSeq


# Background

Variation in transcriptional regulation is a major cause of phenotypic diversity1
. Genomewide association studies (GWAS) have shown that most functional variants reside in
noncoding regions, where they potentially affect transcription factor (TF) binding and
chromatin accessibility to alter gene expression. Pinpointing such regulatory variations,
however, remains challenging2. Here, we developed a hybrid allele-specific chromatin
binding sequencing (HASCh-seq) approach and identified variations in target binding of
the brassinosteroid (BR) responsive transcription factor ZmBZR1 in maize. Chromatin
immunoprecipitation followed by sequencing (ChIP-seq) in B73xMo17 F1s identified
thousands of target genes of ZmBZR1. Allele-specific ZmBZR1 binding (ASB) was
observed for about 14.3% of target genes. It correlated with over 550 loci containing
sequence variation in BZR1-binding motifs and over 340 loci with haplotype-specific DNA
methylation, linking genetic and epigenetic variations to ZmBZR1 occupancy.
Comparison with GWAS data linked hundreds of ASB loci to important yield, growth and
disease-related traits. Our study provides a robust method for analyzing genome-wide
variations of transcription factor occupancy and identified genetic and epigenetic
variations of the BR response transcription network in maize.

###

Revision 

old counts, new peaks and input - 7817 ASBs





klaeren 

arabidopsis orthologs basierend auf Pythozome 13 oder Ensemble .. (3200 bzw. 1770 ..)


original gem daten fuer die BZR1 binding peaks 
im text basierend auf 5 out of 6
die gene partitioning distribution anders als im paper


### FIGURES ###

Figure 1
c) 


### OVERVIEW ### 

hachseq_functions:

addNearestGene () # add orthologs, add gene 


# gene partitioning 

Table 4 (orthologs of 17431 chip peaks... )


analyse_asb_methylation_vs_postfrequency()  # ... 



hachseq_expression:

bQTL_on_ArabidopsisHomolog_geneexpression_evaluation - figure 1 h, i



### datasets 

17463 high confidence ZmBZR1 binding peaks (Table S2)
- 5731 putative ZMBZR1 target genes (Table S3)
- with 580 and 469 (text and figure 1 d) BR repressed and activated genes 




#############
##########

Fig 2c - main prototype



HachSeq_functions.R

17463 high confidence ZmBZR1 binding peaks (Table S2)
6371 putative ZmBZR1 target genes (Table S3) 
including 580 and 469 BR repressed and activated genes (Fig 1d)


Overall, maize BZR1 target genes had 2060 Arabidopsis homologs, of which 1438 (34.1%) were previously identified as AtBZR1 targets 14,15 (Table S4).




analyse_asb_methylation_vs_postfrequency - HachSeq functions


TABLES SUPPLEMENT:

S1 - 2743 BR_up_and_Down - hachseq expression
S2 - 17463 high confidence ZmBZR1 binding peaks - main prototype
S3 - 6371 putative ZmBZR1 target 