
library(tidyverse)
library(magrittr)
library(Biostrings)
library(pheatmap)
library(RColorBrewer)

setwd('~/Desktop/projects/ncba/litao/N_phase_separation/')


#mutation frequency across nCoV strains (Fig 3D) ------------------------

conserved_nt <- readDNAStringSet('./results/exonerate/202009/N_nt_conserved.fa')
nt_freq_n <- read_csv('./results/exonerate/202009/Nt_freq_N.csv')

nt_freq_n_m <- as.matrix(nt_freq_n[2:1261])
rownames(nt_freq_n_m) <- nt_freq_n$nt

nt_freq_n_m[nt_freq_n_m > 50000] <- 0

pheatmap(nt_freq_n_m[1:4,], 
		 cluster_rows = F, cluster_cols = F, 
		 color = colorRampPalette((brewer.pal(n = 8, name = "Blues")))(100),
		 labels_col = conserved_nt, cex.axis= 0.8, las = 2,
		 show_colnames = F,
		 main = 'mutations of N genes of SARS-CoV2',
		 ylab ='conserved AA sequences')

pheatmap(nt_freq_n_m[1:4,], 
		 cluster_rows = F, cluster_cols = F, 
		 color = c("#F8F8F8",
		 		  colorRampPalette(c( gray(0.90), 
		 		  					rgb(0, 158, 143, maxColorValue = 255)))(n=599)),
		 labels_col = conserved_nt, cex.axis= 0.8, las = 2,
		 show_colnames = F,
		 main = 'mutations of N genes of SARS-CoV2',
		 ylab ='conserved AA sequences')



par(oma = c(4,4,4,4))
tes <- pheatmap(t(colSums(nt_freq_n_m) ), oma = c( 4,4,4,4), mar = c( 4,4,4,4),
				border_color = 'red',
				cluster_rows = F, cluster_cols = F, 
				color = c("#F8F8F8",
						  colorRampPalette(c( gray(0.90), 
						  					rgb(0, 158, 143, maxColorValue = 255)))(n=599)),
				cex.axis= 0.8, las = 2,
				show_colnames = F,
				main = 'mutations of N genes of SARS-CoV2',
				ylab ='conserved AA sequences')



#iupred2a score visualization (extended Fig 1b) ---------------------
gene_score_list <- list.files(
	path = '~/Desktop/projects/phase_separation/results/sars-cov2/iupred2a_results/',
	full.names = T)
names(gene_score_list) <- gsub('.iupred2a.txt', '', gene_score_list %>% basename)


source('~/Desktop/projects/phase_separation/scripts/functions_main.R')



library(RSkittleBrewer)
palette(RSkittleBrewer("tropical"))

n_ps <- read_tsv(gene_score_list[names(gene_score_list) == "N"], 
				 skip = 6, 
				 col_names = c('location', 'AA', 'iupred2a_score', 'anchor2_score'))

par(bg =(brewer.pal(n = 8, name = "Blues"))[2], oma= c(4,4,4,4))
plot(n_ps$location, n_ps$iupred2a_score, col =3, ann=F, axes =F,
	 type = 'l', lwd = 2, ylim = c(0, 1.1))
lines(n_ps$location, n_ps$anchor2_score,  lwd = 2, col = 4)
#abline(h =0.5, lty = 2)
axis(1, at = c(1, seq(50, 400, by = 50), 419))
axis(4, at = c(0, 0.2, 0.4, 0.6, 0.8, 1.0))

legend('topright',c('iupred2a','anchor2'), lwd =2, col = c(3,4), bty = 'n' )




#  trinucleotide mutations country distributions -------------------------------
#100849 sequences:
N_seqs <- readDNAStringSet('./results/exonerate/202009/N_nt_all.fa')
mutsite <- subseq(N_seqs, start = 608, width = 3)
SiteTibble <- tibble(trinucleotide=as.data.frame(mutsite)[,1], strain = names(mutsite))
table(SiteTibble$trinucleotide) %>% sort(decreasing = T) #GGG 63331; AAC 36941
#only keep GGG and AAC: 100272 cases
SiteTibble <- SiteTibble[SiteTibble$trinucleotide %in% c('GGG','AAC'), ]


nCoVmeta <- read_tsv('~/data/public/gisaid/20201002/metadata_2020-10-02_07-13.tsv')

nCoVmeta %>% group_by(region) %>% summarise(n())

# 313 strains not included in the meta data:
(SiteTibble$strain %in% nCoVmeta$strain) %>% table()


# region -------------------
MutByRegion <- left_join(SiteTibble, nCoVmeta) %>% select(trinucleotide,  region)
MutByRegion <- group_by(MutByRegion,trinucleotide, region) %>% count %>%
	filter( !is.na(region)) %>% spread( trinucleotide, n)

MutByRegionT <- t(as.matrix(MutByRegion[, 2:3]))
colnames(MutByRegionT) <- MutByRegion$region
barplot(MutByRegionT, legend.text = rownames(MutByRegionT))
mutate(MutByRegion, ratio = AAC/(AAC + GGG))

#country ---------------------

nCoVmeta <- read_tsv('~/data/public/gisaid/20201002/metadata_2020-10-02_07-13.tsv')

nCoVmeta <- select(nCoVmeta, region, country) %>% unique()


MutByRegion <- left_join(SiteTibble, nCoVmeta) %>% select(trinucleotide,  country)
MutByRegion <- group_by(MutByRegion,trinucleotide, country) %>% count %>%
	filter( !is.na(country)) %>% spread( trinucleotide, n) 

MutByRegionClean <- mutate(MutByRegion, sum = (AAC + GGG),ratio = AAC/(AAC + GGG) ) %>% arrange(desc(ratio))
MutByRegionClean1K <- filter(MutByRegionClean, sum > 1000)

#write_csv(MutByRegionClean, './results/MutByRegionClean.csv')
#write_csv(MutByRegionClean1K, './results/MutByRegionClean1K.csv')

MutByRegionClean <- read_csv('./results/MutByRegionClean.csv')
MutByRegionClean1K <- read_csv('./results/MutByRegionClean1K.csv')

#mutationDeathRatio <- read_csv('./NCB_revise/mutation_and_DeathRatio.csv')
mutationDeathRatio <- read_csv('./NCB_revise/MutByRegionClean all with mortality 20201014.csv')

mutationDeathRatio <- filter(mutationDeathRatio, !is.na(ratio)) %>% filter(PatientSum != '-')

mutationDeathRatio$deathratio <- as.numeric(mutationDeathRatio$PatientDeath) /as.numeric(mutationDeathRatio$PatientSum)

filter(mutationDeathRatio, sum > 100) %$% plot(ratio, deathratio)

mutationDeathRatio <- left_join(mutationDeathRatio, nCoVmeta) %>% arrange(region)
mutationDeathRatio$label <- rep(1:6, time = count(group_by(mutationDeathRatio, region))$n)
leginfo <- group_by(mutationDeathRatio, region, label) %>% count()


library(RSkittleBrewer)
palette(RSkittleBrewer("smarties"))
palette(c(brewer.pal(12, 'Set3'), brewer.pal(8, 'Pastel1')))

plot(mutationDeathRatio$ratio, mutationDeathRatio$deathratio, type = 'p',
	 col = mutationDeathRatio$label, pch = 16, cex = 2,
	 xlab = 'AAC frequency', ylab = 'death ratio', main = 'all country')

legend('topright', col = leginfo$label, pch =16, pt.cex =2,
	   legend = leginfo$region)
fm <- lm(mutationDeathRatio$deathratio ~mutationDeathRatio$ratio)
fm
summary(fm)  
abline(fm, lwd = 3)

#1k filteration:
mutationDeathRatio1K <- filter(mutationDeathRatio, sum > 100)


plot(mutationDeathRatio1K$ratio, mutationDeathRatio1K$deathratio, type = 'p',
	 col = mutationDeathRatio1K$label, pch = 16, cex = 2,
	 xlab = 'AAC frequency', ylab = 'death ratio', main = 'genomes > 100')

legend('topright', col = leginfo$label, pch =16, pt.cex =2,
	   legend = leginfo$region)
fm <- lm(mutationDeathRatio1K$deathratio ~mutationDeathRatio1K$ratio)
fm
summary(fm)  
abline(fm, lwd = 3)




# alignment between nCov and SARS -----------------------------------------
library(genbankr)
library(msa)
setwd('~/Desktop/projects/nCoV/')

sarsGb <- readGenBank('./data/genome/sars_gb/SARS_urbani.gb')
Nsars <- transcripts(sarsGb)[3]$translation


NnCoV <- readAAStringSet('~/Desktop/projects/ncba/litao/N_phase_separation/results/N_conserved.fa')

Ns <- append(NnCoV, Nsars)
names(Ns) <- c('N_nCoV', 'N_SARS')

msaPrettyPrint(
	msa(Ns, method = 'ClustalOmega',
		order = 'input'), 
	output="pdf", #paperWidth = 9,
	file = './N_sars_ncov.pdf',
	askForOverwrite=FALSE, verbose=FALSE)


# R+Y, DDX4 like ----------------------------------------------------------

library(RSkittleBrewer)
palette(RSkittleBrewer("smarties"))
palette(RSkittleBrewer("tropical"))

aa_cider_iupred2a_list <- readRDS('~/Desktop/projects/phase_separation/results/sars-cov2/rda/aa_cider_iupred2a_list.rda')
aa_cider_iupred2a_list$M

library(data.table)
source('~/Desktop/projects/phase_separation/scripts/functions_main.R')


plot_AAs_properties <- function(dt, gene_symbol){
	#Usage:
	#	plot the AA distributions, the charges and the phase separatioin scores
	#of given gene.
	#parameters:
	#	dt: data.table contain the amino acid properties of a protein.
	#		col names of dt:
	#			amino_acid	FCR	NCPR	iupred2_score	achore2_score.
	#	gene_symbol: the gene used for analysis.
	#AA charges:	
	#positive: R(Arg)	K(Lys)	H(His)
	#negative: D(Asp)	E(Glu)
	#polar: S(Ser)	T(Thr)	N(Asn)	Q(Gln)
	#special: C(Cys)	G(Gly)	P(Pro)
	#Hydrophobic: A(Ala)	V(Val)	I(Ile)	L(Leu)	M(Met)	F(Phe)	Y(Tyr)	W(Trp)
	
	layout(matrix(c(1,1,1,2,3,4)), 6,1)
	par(mar =c(0,4,1,1), oma = c(4, 3, 2, 1))
	
	plot(x =0, xlim=c(0,dt[,.N]),ylim=c(0,20), type = 'n', yaxt = 'n',  xaxt ='n',
		 ylab="amino acids", main = gene_symbol)
	#add y label:
	axis(2, at = seq(0.5, 20), tick = F, las =2, line = -0.7, 
		 labels = c('R','K','H', 
		 		   'D', 'E', 
		 		   'S','T', 'N','Q', 
		 		   'C', 'G', 'P',
		 		   'A', 'V', 'I', 'L','M','F', 'Y', 'W') %>% rev)
	
	#background color for the five classes of AAs.
	rect(0, 17, dt[, .N], 20, col = rgb(255,192,128, maxColorValue = 255),  border = NA)
	rect(0, 15, dt[, .N], 17, col = rgb(192,192,192, maxColorValue = 255),  border = NA)
	rect(0, 11, dt[, .N], 15, col = rgb(131,192,203, maxColorValue = 255),  border = NA)
	rect(0, 8, dt[, .N], 11, col = rgb(217,153,153, maxColorValue = 255),  border = NA)
	rect(0, 0, dt[, .N], 8, col = rgb(224,224,128, maxColorValue = 255),  border = NA)
	
	
	aa_y_location <- data.table(
		amino_acid = c('R','K','H',  'D', 'E', 'S','T', 'N','Q', 'C', 
					   'G','P','A', 'V', 'I', 'L','M','F', 'Y', 'W'),
		location = seq(20,1))
	
	for (i in seq(1,dt[, .N])) {
		AA <- dt[i, amino_acid]
		j <- aa_y_location[amino_acid == AA, location]
		rect(i, j-1, i+1, j, col = par("fg"))
	}
	
	#Fraction of charged residue(FCR):
	plot(1:dt[, .N], dt[, FCR], type = 'l', ylab = 'FCR', 
		 xaxt = 'n', xlab = '', lwd =2, col = 6)
	
	#Net charges per residue(NCPR):
	plot(dt[, NCPR], type = 'l', xaxt = 'n', ylab = 'NCPR',lwd =2, 
		 col = rgb(248,89,171, maxColorValue = 255))
	abline(h = 0, col = "lightgray")
	#Phase separation scores predicted by iupred2 software:
	plot(dt[, iupred2_score], type = 'l', ylab = 'phase separation',
		 xlab = '', col =  rgb(248,89,171, maxColorValue = 255), lwd =2, ylim = c(0,1))
	lines(dt[, achore2_score],type = 'l', 
		  col = rgb(0,158,143, maxColorValue = 255), lwd =2)
	abline(h = 0.5, col = "lightgray")
}

#par(bg = 'white')
plot_AAs_properties(aa_cider_iupred2a_list$N,'N')



Npro <- readAAStringSet('./results/N_conserved.fa')

N_AAfreq <- as.character(Npro[[1]]) %>% strsplit('abc', split = '') %>% 
	unlist %>% table %>% as.data.frame() %>% as_tibble()


names(N_AAfreq) <- c('AA', 'Num')

N_AAfreq$freq <- N_AAfreq$Num / sum(N_AAfreq$Num)

write_csv(N_AAfreq, './results/N_AA_freq.csv')


