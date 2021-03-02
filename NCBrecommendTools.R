

#P-Score (http://pound.med.utoronto.ca/~JFKlab/Software/psp.htm)

#step 1: software download:


#step 2 : data filter: proteins must be at least 140 residues:

library(tidyverse)
library(magrittr)
library(Biostrings)

setwd('~/Desktop/projects/ncba/litao/N_phase_separation/')

ncovPro <- readAAStringSet('~/Desktop/projects/phase_separation/data/sars-cov2/proteins_sars-cov2-29proteins.fa')

#only 16 proteins with length longer than 140:
nconPro140 <- ncovPro[width(ncovPro) >= 140]

writeXStringSet(nconPro140, '~/Desktop/projects/phase_separation/data/sars-cov2/proteins_sars-cov2-29proteinsLength140.fa')


#step 3:
#python2 elife_phase_separation_predictor.py  \
#    proteins_sars-cov2-29proteinsLength140.fa > P-Score-results_nCovLength140.txt

#python2 elife_phase_separation_predictor.py \
#	proteins_sars-cov2-29proteinsLength140.fa -residue_scores -score_components -output P-Score-ComResults_nCovLength140.txt



# R + Y -------------------------------------------------------------------

library(Biostrings)

N_seq <- readAAStringSet('~/Desktop/projects/phase_separation/results/sars-cov2/sep_fa/N.fa')
N_iupredScore <- readr::read_tsv('~/Desktop/projects/phase_separation/results/sars-cov2/iupred_cider/iupred2a_cider_N.txt')

N_iupredScore$meanScore <- (N_iupredScore$iupred2_score + N_iupredScore$achore2_score)/2

aa_freq <- dplyr::filter(N_iupredScore, meanScore >= 0.5) %$% table(amino_acid) %>%
	sort(decreasing = T) %>% tibble::as_tibble()

aa_freq$freq <- aa_freq$n / sum(aa_freq$n)

readr::write_tsv(aa_freq, '~/Desktop/projects/ncba/litao/N_phase_separation/results/N_IDR_AA_freq.txt')




#P-Score ------------------------------------------------------------------

#awk '{print $1 "\t" $2"\t" $3 "\t" $4"\t" $NF}'  P-Score-ComResults_nCovLength140.txt |
#	sed 's/>//' | grep -v '^P' > P-Score_nCovLength140_clean.txt
library(tidyverse)
pscore <- readr::read_tsv('/Users/xuejia/Desktop/projects/ncba/litao/N_phase_separation/results/llpsTools/P-Score/P-Score_nCovLength140_clean.txt', col_names = c('P-Score', 'locus', 'aa', 'pscore', 'protein'))

table(pscore$protein)


layout(matrix(1:18, nrow = 6, ncol = 3))


prots <- names(table(pscore$protein))

par(mar = c(4,4,2,1))

for (i in prots) {
	tmp <- filter(pscore, protein == i )
	plot(tmp$locus, tmp$pscore, type = 'l', xlab = i, ylab = "P-Score")
	abline(h =4, col = 'red', lty =2)
}







