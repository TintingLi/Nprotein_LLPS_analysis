
library(tidyverse)
library(magrittr)
library(Biostrings)

setwd('~/Desktop/projects/ncba/litao/N_phase_separation/')


#Update the N protein mutation analysis to 20200918 

#download updated data:-----------
#100850 start sequence
sars_cov2_genomes <- readDNAStringSet('./data/20200919/sequences_2020-09-18_07-12.fasta')
sars_cov2_genomes_abfreq <- alphabetFrequency(sars_cov2_genomes) %>% as_tibble()

table(sars_cov2_genomes_abfreq$`-` == 0)    #T 104289
table(sars_cov2_genomes_abfreq$`+` == 0)    #T 104289
table(sars_cov2_genomes_abfreq$`.` == 0)    #T 104289


colSums(sars_cov2_genomes_abfreq)
# 5 genome sequence contains "U": U is read as T:
#writeXStringSet(sars_cov2_genomes, './data/20200919/sequences_2020-09-18_07-12_removeUs.fasta')


sars_cov2_genomes <- readDNAStringSet('~/litt/projects/collection/ncba/litao/N_phase_separation/data/20200919/sequences_2020-09-18_07-12_removeUs.fasta')
sars_cov2_genomes_abfreq <- alphabetFrequency(sars_cov2_genomes) %>% as_tibble()


# map N protein to genome sequences by exonerate -------------------------------
#pro2genome.sh --> clean_align.txt

n2genome <- read_tsv('./results/exonerate/202009/clean_align.txt',
                     col_names = c('query_range', 'target_genome', 'targe_range'))
#only keep alignments containg complete N protein sequence: 101102 alignments
n2genome <- n2genome[n2genome$query_range == 'Query range: 0 -> 419',]
table(n2genome$query_range)
#write_tsv(n2genome, './results/exonerate/202009/clean_complete_align.txt')


#get N gene locus and subtract nt sequences   ---------------------------------
N_locus <- tibble( strain= gsub('Target: ', '', n2genome$target_genome), 
                   N_start = gsub('Target range: ', '', n2genome$targe_range) %>%
                       gsub(' ->.*$', '', x =.) %>% as.integer(),
                   N_end = gsub('Target range: ', '', n2genome$targe_range) %>%
                       gsub('^.*-> ', '', x =.) %>% as.integer())
#left locus: +1, right locus: +3
N_locus$N_start <- N_locus$N_start + 1
N_locus$N_end <- N_locus$N_end + 3
N_locus$width <- N_locus$N_end - N_locus$N_start + 1

table(N_locus$width) 

# 100898 seqs's length is 1260
#only keep gene_length == 1260:
N_locus_1260 <- filter(N_locus, width == 1260)

#subtract all N gene sequences: 
sars_cov2_genomes <- sars_cov2_genomes[names(sars_cov2_genomes) %in% N_locus_1260$strain]
sars_cov2_genomes <- sars_cov2_genomes[width(sars_cov2_genomes) > 29000] # 100850

(names(sars_cov2_genomes) %in% N_locus_1260$strain) %>% table()
(N_locus_1260$strain  %in% names(sars_cov2_genomes)) %>% table()
N_locus_1260 <- N_locus_1260[N_locus_1260$strain  %in% names(sars_cov2_genomes), ]


#discard"Netherlands/ZH-EMC-347/2020"
N_locus_1260 <- N_locus_1260[N_locus_1260$strain != "Netherlands/ZH-EMC-347/2020", ]



N_seqs <- DNAStringSet() # 100849 sequences
for (i in seq(1, dim(N_locus_1260)[1])) {
    print(i)
    gn <- N_locus_1260$strain[i]
    gs <- sars_cov2_genomes[names(sars_cov2_genomes) == gn]
    tmp <- subseq(gs, 
                  start = N_locus_1260$N_start[i],
                  end = N_locus_1260$N_end[i])
    N_seqs <- append(N_seqs, values = tmp)
}
#writeXStringSet(N_seqs, './results/exonerate/202009/N_nt_all.fa')




#calculate the freqs for each nt: -----------
#100849 sequences:
N_seqs <- readDNAStringSet('./results/exonerate/202009/N_nt_all.fa')
width(N_seqs) %>% table()


N_seqs_matrix <- as.matrix(N_seqs)
st_n <- tibble()
st_n <- N_seqs_matrix[, 1] %>% table() %>% as_tibble()
for (i in seq(2, width(N_seqs)[1])) {
    tmp <- N_seqs_matrix[, i] %>% table() %>% as_tibble()
    st_n <- full_join(st_n, tmp, by = '.')
}
names(st_n) <- c('nt', paste0('pos', 1:1260))

st_n_sort_m <- as.matrix(st_n[, 2:1261])
#turn NA into 0:
st_n_sort_m[is.na(st_n_sort_m)] <- 0
rownames(st_n_sort_m) <- st_n$nt
#write_csv(as_tibble(st_n_sort_m, rownames = 'nt'), "./results/exonerate/202009/Nt_freq_N.csv")
st_n_sort_m <- read_csv('./results/exonerate/202009/Nt_freq_N.csv')


#get conserved sequences:
conserved_nt <- apply(st_n_sort_m, 2, function(x){
    return(rownames(st_n_sort_m)[which.max(x)])
}) 

N_conserved <- paste(conserved_nt,collapse = '') %>% DNAStringSet()
#writeXStringSet(N_conserved, './results/exonerate/202009/N_nt_conserved.fa')
N_conserved <- readDNAStringSet('./results/exonerate/202009/N_nt_conserved.fa')

#get highest mutation sites:

mutation_freq_N <- apply(st_n_sort_m, 2, function(x){
    tmp <- x[order(x, decreasing = T)]
    names(tmp) <- rownames(st_n_sort_m)[order(x, decreasing = T)]
    num1 <- paste0(tmp[1], " (", names(tmp)[1], ')')
    num2 <- paste0(tmp[2], " (", names(tmp)[2], ')')
    num3 <- paste0(tmp[3], " (", names(tmp)[3], ')')
    return(c(num1, num2, num3))
}) %>% t() %>% as_tibble(rownames = 'position')


colnames(mutation_freq_N) <- c('postion', 'conserved', 'mut1', 'mut2')

mutation_freq_N

mutation_freq_N$mut1_num <- gsub(' (.*)$', '', mutation_freq_N$mut1) %>% as.numeric()
mutation_freq_N$mut1_AA <- gsub('^[0-9].* ', '', mutation_freq_N$mut1)
mutation_freq_N$mut2_num <- gsub(' (.*)$', '', mutation_freq_N$mut2) %>% as.numeric()
mutation_freq_N$mut2_AA <- gsub('^[0-9].* ', '', mutation_freq_N$mut2)
mutation_freq_N %<>% arrange(desc(mut1_num, mut2_num))

#write_csv(mutation_freq_N,  './results/exonerate/202009/mutation_freq_N_nt.csv' )

mutation_freq_N <- read_csv('./results/exonerate/202009//mutation_freq_N_nt.csv')





# calculate the trinucleotide mutations -----------------------------------
#100849 sequences:
N_seqs <- readDNAStringSet('./results/exonerate/202009/N_nt_all.fa')
mutsite <- subseq(N_seqs, start = 608, width = 3)
SiteTibble <- tibble(trinucleotide = as.data.frame(mutsite), strain = names(mutsite))
table(SiteTibble$trinucleotide) %>% sort(decreasing = T) #GGG 63331; AAC 36941




