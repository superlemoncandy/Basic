library("wiggleplotr")
library("dplyr")
library("GenomicRanges")
library("GenomicFeatures")
library("biomaRt")

TxDb = makeTxDbFromGFF(file = '/Users/shaoqianma/research_2022/research/expr_results/single_cell_AID/stringtie_gtf/annotations/gencode.vM22.annotation.gtf',
                       format = "gtf",
                       dataSource = "GENCODE",
                       organism = "Mus musculus")

transcripts <- transcriptsBy(TxDb, "gene")
transcript_lengths <- transcriptLengths(TxDb)
# 
#gene_ranges <- GenomicRanges::reduce(transcripts)
#gene_lengths <- width(gene_ranges)

# 
longest_transcript <-  lapply(transcripts, FUN = function(x)x[which.max(width(x))])
# #seq <- extractSeq("reference_genome.fa", longest_transcript)
#rd = splitAsList(transcripts, names(transcripts))
saveRDS(longest_transcript, file = './splicing/longest_transcript.rds')

# then collect the exons per gene id

exons <- exonicParts(TxDb, linked.to.single.gene.only=T)
#exons_re <- GenomicRanges::reduce(exons)

introns <- intronicParts(TxDb, linked.to.single.gene.only=T)
#introns_re <- GenomicRanges::reduce(introns)

IX <- data.frame(GeneID=introns$gene_id, Chr=as.character(seqnames(introns)), width= width(introns),
                 Start=start(introns), End=end(introns), Strand=strand(introns), stringsAsFactors=FALSE, intronic_part = introns$intronic_part)
EX <- data.frame(GeneID=exons$gene_id, Chr=as.character(seqnames(exons)), width = width(exons),
                 Start=start(exons), End=end(exons), Strand=strand(exons), stringsAsFactors=FALSE, exonic_part = exons$exonic_part)


# retrieve forward match and back ward match for each exon. ####
#only the internal exons
# 
# max_exon_rank<- as.data.frame(exons)%>%
#   group_by(gene_id)%>%
#   mutate(max_exon_rank = max(exonic_part))%>%ungroup()
# 
# internal_exon<- max_exon_rank[(max_exon_rank$exonic_rank != 1) & (max_exon_rank$exonic_rank != max_exon_rank$max_exon_rank), ]
# 
internal_exon_GRange<- GRanges(seqnames = internal_exon$seqnames,ranges = IRanges(internal_exon$start,internal_exon$end),strand = internal_exon$strand,
        tx_id = internal_exon$tx_id, tx_name =internal_exon$tx_name, gene_id = internal_exon$gene_id, exon_id = internal_exon$exon_id,
        exon_name = internal_exon$exon_name, exon_rank = internal_exon$exon_rank, exonic_part = internal_exon$exonic_part, max_exon_rank =internal_exon$max_exon_rank)

# EX_IX_match<- introns[precede(internal_exon_GRange, introns, select=c("first"), ignore.strand=FALSE)]
# IX_EX_match<- introns[follow(internal_exon_GRange, introns, select=c("last"), ignore.strand=FALSE)]

# only keep 2nd exon, query intron first ####
IX <- IX[IX$width > 10000,]
introns <- introns[width(introns) > 10000]

EX_exon_nofrist<- EX[EX$exonic_part != 1,]
exons_nofrist<- exons[exons$exonic_part !=1]
idx <- precede(exons_nofrist, introns, select=c("first"), ignore.strand=FALSE)
idx <- idx[!is.na(idx)]
IX_EX_match<- introns[idx]

idx <- follow(exons_nofrist, introns, select=c("last"), ignore.strand=FALSE)
idx <- idx[!is.na(idx)]
EX_IX_match<- introns[idx]

EX_IX_match<- as.data.frame(EX_IX_match)
IX_EX_match <- as.data.frame(IX_EX_match)

IX_EX_match<- unique(IX_EX_match)
EX_IX_match<- unique(EX_IX_match)

# then reverse query exon by EX_IX ####


EX_exon_nofrist_gr<- GRanges(seqnames = EX_exon_nofrist$Chr,ranges = IRanges(EX_exon_nofrist$Start,EX_exon_nofrist$End),strand = EX_exon_nofrist$Strand,
                             gene_id = EX_exon_nofrist$GeneID, exonic_part = EX_exon_nofrist$exonic_part)
EX_IX_match_gr<- GRanges(seqnames = EX_IX_match$seqnames,ranges = IRanges(EX_IX_match$start,EX_IX_match$end),strand = EX_IX_match$strand,
                         gene_id = EX_IX_match$gene_id, exonic_part = EX_IX_match$intronic_part)
IX_EX_match_gr<- GRanges(seqnames = IX_EX_match$seqnames,ranges = IRanges(IX_EX_match$start,IX_EX_match$end),strand = IX_EX_match$strand,
                         gene_id = IX_EX_match$gene_id, exonic_part = IX_EX_match$intronic_part)

idx <- precede(EX_IX_match_gr, EX_exon_nofrist_gr, select=c("first"), ignore.strand=FALSE)
#idx <- idx[!is.na(idx)]
EX_exon_1<- EX_exon_nofrist_gr[idx]
idx <- follow(IX_EX_match_gr, EX_exon_nofrist_gr, select=c("last"), ignore.strand=FALSE)
#idx <- idx[!is.na(idx)]
EX_exon_2<- EX_exon_nofrist_gr[idx]

# intersect: the exon must has both intron > 10kb ####
# 5254 exons 
co_splicing_EX<- intersect(EX_exon_1, EX_exon_2)
co_splicing_EX_df<- as.data.frame(co_splicing_EX)
co_splicing_EX_10kb_fwd <- co_splicing_EX_df[co_splicing_EX_df$strand == '+',] # 2688
co_splicing_EX_10kb_rev <- co_splicing_EX_df[co_splicing_EX_df$strand == '-',] # 2566
# output fwd and rev exon seperately ####
write.table(co_splicing_EX_10kb_fwd, file = './ref_bed/co_splicing_EX_10kb_fwd.bed', quote = F, sep = '\t', row.names = F, col.names = F)
write.table(co_splicing_EX_10kb_rev, file = './ref_bed/co_splicing_EX_10kb_rev.bed', quote = F, sep = '\t', row.names = F, col.names = F)

# for compute co-splicing score using deeptools
#computeMatrix scale-regions -p 10 --binSize 20 --regionBodyLength 20 -a 2000 -b 2000 -R co_splicing_EX_10kb_fwd.bed -S  --skipZeros -o RNA_co_splicing-fwd_scaleRegion-data.gz

library(ggplot2)
library(tidyverse)
library(reshape2)
library(ggsci)
library(Rmisc)
library(data.table)
pacman::p_load(BiocManager, rtracklayer, Rsamtools,
               GenomicRanges, IRanges,
               org.Mm.eg.db, biomaRt, AnnotationFilter,
               foreach, doParallel, 
               dplyr, tibble, matrixStats,
               ggplot2, grid, gridExtra, viridis, RColorBrewer)

colors_20 <- c('#bf4040', '#cc8800', '#808000', '#408000', '#006622', '#2d8659', '#5c8a8a',
               '#0073e6', '#4d4dff', '#5353c6', '#4d0099', '#660080', '#602060', '#c6538c',
               '#99003d', '#cc5200', '#b32400', '#663300', '#b3b300', '#4d9900')

setwd('/Users/shaoqianma/research_2022/research/expr_results/single_cell_AID/scFAST-IAA-data/TT-seq_mESC_pluripotency-master/fig3/')

# load data
scaler_fwd <- fread('../pol2data/RNA_co_splicing-fwd10kb_scaleRegion-data.gz',header = F)
scaler_rev <- fread('../pol2data/RNA_co_splicing-rev10kb_scaleRegion-data.gz',header = F)

#(2000 + 20 + 2000)/20*7

clean_scaler <- scaler_fwd %>% dplyr::select(c(-1:-3,-5,-6))
# merge replicate ###
for(i in 1:4){
  clean_scaler[,2:202] <- clean_scaler[,2:202] + clean_scaler[,(201*i+2):(201*(i+1)+1)]
}
clean_scaler[,(201*5+2):(201*6+1)] <- clean_scaler[,(201*5+2):(201*6+1)] + clean_scaler[,(201*6+2):(201*7+1)]

clean_scaler <- clean_scaler[,c(1:202,(201*5+2):(201*6+1))]
# clean <- scaler_rev %>% select(c(-1:-3,-5,-6)) 
clean_scaler[,2:202] <- clean_scaler[,2:202]/5
clean_scaler[,203:403] <- clean_scaler[,203:403]/2

# combine fwd and rev.  1:101, 103:202, 203:302, 304:403
# clean_scaler<- cbind(rbind(clean_scaler_fwd[,c(1:101, 203:302)], clean_scaler_rev[,c(1,103:202,304:403)], use.names =F), 
#                      rbind(clean_scaler_fwd[,c(103:202,304:403)], clean_scaler_rev[,c(2:101, 203:302)], use.names=F))


#clean_scaler<- cbind(rbind(clean_scaler_fwd[,1:201], clean_scaler_rev[,c(1,204:403)], use.names =F), rbind(clean_scaler_fwd[,204:403], clean_scaler_rev[,2:201], use.names=F))
long_scaler <- melt(clean_scaler_rev,id.vars = 'V4',value.name = 'signal') %>%
  dplyr::select(-variable)

# add sample name
long_scaler$sample <- rep(c("CTRL","IAA24h"),
                          each = nrow(clean_scaler_rev)*201) #200 for combined fwd and rev

# add x position
long_scaler$pos <- rep(c(1:201),each = nrow(clean_scaler_rev),times = 2)

long_scaler<- long_scaler[!is.infinite(long_scaler$signal),]

# calculate means
long_scaler <- long_scaler %>%
  # remove na
  tidyr::drop_na() %>%
  dplyr::group_by(sample,pos) %>%
  # mean and 95% interval confidence
  dplyr::summarise(mean_signal = mean(signal),
                   sd = sd(signal),
                   upper = CI(signal,ci = 0.95)[1],
                   lower = CI(signal,ci = 0.95)[3])

#long_scaler$roll_mean <- stats::filter(long_scaler$mean_signal, rep(1,10))
library(ggbreak) 
p_fwd <- ggplot(long_scaler,aes(x = pos,y = log1p(mean_signal))) +
  # geom_rect(aes(xmin = -Inf, xmax = 151, ymin = -Inf, ymax = Inf),
  #           fill = "grey", alpha = 0.01, linetype = 0) +
  # geom_rect(aes(xmin = 401, xmax = Inf, ymin = -Inf, ymax = Inf),
  #           fill = "grey", alpha = 0.01, linetype = 0) +
  geom_line(aes(color = sample),size = 1.5) +
  #geom_smooth(aes(color = sample),size = 1.5,method = 'loess') + 
  #geom_ribbon( aes(ymin = lower, ymax = upper), alpha =0.2 )+
  scale_color_manual(values = colors_20[c(13, 2)]) + 
  theme_classic(base_size = 14) +
  # x label
  #scale_y_break(breaks = c(2,3),space = 0.2,scales = 1.5,expand = c(0,0))+
  scale_x_continuous(breaks = c(0,100,200),
                     labels = c('-2 kb','Exon','+2 kb')) +
  xlab('') + ylab('log(Normalized signal)') + 
  theme(panel.grid = element_blank(),aspect.ratio = 0.8, axis.title.x = element_text(size = 18), 
        legend.title = element_text(size = 16),
        axis.title.y = element_text(size = 18), title = element_text(size = 16), legend.text = element_text(size = 14, lineheight = 2),
        plot.title = element_text(color="black", face="bold"),
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))+
  ggtitle('Total RNA-seq metaplot (N=2106)') 


library(cowplot)

plot_grid(p_fwd,p_rev)

library(export)
graph2ppt(file = "../pol2data/Figures/pol2_co_splicing.pptx")

library(export)
graph2ppt(file = "./reads_annotation.pptx")


if(FALSE){
  ##. #####
  
  
  # nrow(unique(EX_exon_nofrist))
  # nrow(unique(EX_IX_match))
  # nrow(unique(IX_EX_match))
  EX_exon_nofrist <- unique(EX_exon_nofrist)
  EX_IX_match <- unique(EX_IX_match)
  IX_EX_match <- unique(IX_EX_match)
  
  # 
  
  
  
  
  
  
  #extract the introns with exact match of both start and end positions
  idx <- EX_IX_match$end == (IX$Start -1) & IX_EX_match$start == (IX$End +1)
  IX<- IX[idx,]
  EX_IX_match<- EX_IX_match[idx,]
  IX_EX_match<- IX_EX_match[idx,]
  
  # Introns longer than 10 kb (n = 17,828) in the refseq gene list are grouped by their sizes, 
  #and the average RNA-seq read count per 100 bp bins are plotted by the distance from the 3′SS for each group.
  
  introns_re <- GenomicRanges::reduce(introns)
  table(width(introns_re) > 10000)
  
  # better filter out the last intron, extract intron with specific length #### 
  
  max_intron_rank<- as.data.frame(introns)%>%
    group_by(gene_id)%>%
    mutate(max_intron_rank = max(intronic_part))%>%ungroup()
  
  internal_intron<- max_intron_rank[(max_intron_rank$intronic_part != max_intron_rank$max_intron_rank), ]
  
  
  
  
  introns_re <- as.data.frame(introns_re[width(introns_re) > 10000]) #22018
  introns_re$id <- paste("gene","intron", 1:nrow(introns_re),sep = '_')
  introns_re_fwd <- introns_re[introns_re$strand == '+',] #11074 
  introns_re_rev <- introns_re[introns_re$strand == '-',] #10944
  
  
  write.table(introns_re_fwd[,c(1,2,3,4,6,5)],quote = F,sep = '\t',row.names = F,col.names = F,file = './pol2_ChIP/intron_1kb_fwd.bed')
  write.table(introns_re_rev[,c(1,2,3,4,6,5)],quote = F,sep = '\t',row.names = F,col.names = F,file = './pol2_ChIP/intron_1kb_rev.bed')
  
  introns_re <- introns_re[introns_re$width > 50000,]
  introns_re$id <- 1:nrow(introns_re)
  introns_re_fwd <- introns_re[introns_re$strand == '+',] # 2058
  introns_re_rev <- introns_re[introns_re$strand == '-',] #1971
  write.table(introns_re_fwd[,c(1,2,3,4,6,5)],quote = F,sep = '\t',row.names = F,col.names = F,file = './pol2_ChIP/intron_50kb_fwd.bed')
  write.table(introns_re_rev[,c(1,2,3,4,6,5)],quote = F,sep = '\t',row.names = F,col.names = F,file = './pol2_ChIP/intron_50kb_rev.bed')
  
  
  
  
  #retrieve ajacent introns
  # first filter out features with single intron
  # extract exon id and Txn id for exons, and then for each exon, finding precedes and next
  
  max_intron_rank<- IX%>%
    group_by(gene_id)%>%
    mutate(n = max(intronic_part))%>%ungroup()
  
  
  
  
  # only retain long introns for Read density in the windows is normalized to the level of reads at the 3′ SS, to compensate for expression differences between genes. The average and
  # the standard deviation of the slope are shown
  
  
  ### count intron number, length, gene length for each gene with the longest transcript, ####
  
  longest_transcript_id<- unlist(lapply(longest_transcript, function(x)x$tx_id))
  idx<- unlist(lapply(introns$tx_id, function(x)all(x%in% longest_transcript_id)))
  
  longest_transcript_introns<- introns[idx]
  
  
  
  saveRDS(longest_transcript_introns, file = './splicing/longest_transcript_introns.rds')
  
  load("~/research_2022/research/expr_results/single_cell_AID/published data analysis/xie_markers.rda")
  library(data.table)
  g2s <- fread('../IntergenicTranscription-master/IAA_gffcompare/GencodeReference/g2s_vm22_gencode.txt',header = F,data.table = F) #disable data.table mode
  colnames(g2s) <- c("geneid","symbol")
  ids <- data.frame(geneid=longest_transcript_introns$gene_id, symbol = 'NULL')
  table(ids$geneid %in% g2s$geneid) #
  ids <- ids[ids$geneid %in% g2s$geneid,] #
  ids$symbol <- g2s[match(ids$geneid,g2s$geneid),2]
  
  longest_transcript_introns$symbol <- ids[, ]$symbol
  
  longest_transcript_introns_df<- as.data.frame(longest_transcript_introns) %>%
    group_by(symbol)%>%
    count()
  
  longest_transcript_introns$length <- width(longest_transcript_introns)
  longest_transcript_introns_df<- as.data.frame(longest_transcript_introns) %>%
    group_by(symbol)%>%
    mutate(intron_length = sum(length))%>%ungroup()
  
  longest_transcript_introns_df <- longest_transcript_introns_df[!duplicated(longest_transcript_introns_df$gene_id),]
  
  # use of smart-seq markers
  # one_cell_minor_ZGA_intron<- longest_transcript_introns_df[longest_transcript_introns_df$symbol %in% Zygote_gene,]%>%mutate(gene_class = 'one_cell')
  # early_2cell_minor_ZGA_intron<- longest_transcript_introns_df[longest_transcript_introns_df$symbol %in% early_TwoCell_gene,]%>%mutate(gene_class = 'early_2cell')
  # late_2cell_ZGA_intron<- longest_transcript_introns_df[longest_transcript_introns_df$symbol %in% late_TwoCell_gene,]%>%mutate(gene_class = 'late_2cell')
  # #four_cell_ZGA_intron<- longest_transcript_introns_df[longest_transcript_introns_df$symbol %in% FourCell_gene,]%>%mutate(gene_class = '4cell')
  # #eight_cell_ZGA_intron<- longest_transcript_introns_df[longest_transcript_introns_df$symbol %in% eightCell_gene,]%>%mutate(gene_class = '8cell')
  # #sixteen_cell_ZGA_intron<- longest_transcript_introns_df[longest_transcript_introns_df$symbol %in% sixteenCell_gene,]%>%mutate(gene_class = '16cell')
  # pluri_minor_ZGA_intron<- longest_transcript_introns_df[longest_transcript_introns_df$symbol %in% ICM_gene,]%>%mutate(gene_class = 'ICM')
  
  # gencode<- gencode[!duplicated(gencode),]
  # 
  # one_cell_minor_ZGA_intron<- gencode[gencode$symbol %in% Zygote_gene,]%>%mutate(gene_class = 'one_cell')
  # early_2cell_minor_ZGA_intron<- gencode[gencode$symbol %in% early_TwoCell_gene,]%>%mutate(gene_class = 'early_2cell')
  # late_2cell_ZGA_intron<- gencode[gencode$symbol %in% late_TwoCell_gene,]%>%mutate(gene_class = 'late_2cell')
  # #four_cell_ZGA_intron<- longest_transcript_introns_df[longest_transcript_introns_df$symbol %in% FourCell_gene,]%>%mutate(gene_class = '4cell')
  # #eight_cell_ZGA_intron<- longest_transcript_introns_df[longest_transcript_introns_df$symbol %in% eightCell_gene,]%>%mutate(gene_class = '8cell')
  # #sixteen_cell_ZGA_intron<- longest_transcript_introns_df[longest_transcript_introns_df$symbol %in% sixteenCell_gene,]%>%mutate(gene_class = '16cell')
  # pluri_minor_ZGA_intron<- gencode[gencode$symbol %in% ICM_gene,]%>%mutate(gene_class = 'ICM')
  
  
  
  one_cell_minor_ZGA_intron<- longest_transcript_introns_df[longest_transcript_introns_df$symbol %in% one_cell_minor_ZGA,]%>%mutate(gene_class = 'one_cell')
  early_2cell_minor_ZGA_intron<- longest_transcript_introns_df[longest_transcript_introns_df$symbol %in% early_2cell_minor_ZGA,]%>%mutate(gene_class = 'early_2cell')
  ICM_df<- longest_transcript_introns_df[longest_transcript_introns_df$symbol %in% pluripotent_gene,]%>%mutate(gene_class = 'pluripotent')
  
  
  #选择xie one cell, xie early, smart late, smart pluri
  
  #db<- rbind(one_cell_minor_ZGA_intron,early_2cell_minor_ZGA_intron, pluri_minor_ZGA_intron)
  
  db<- rbind(one_cell_minor_ZGA_intron, early_2cell_minor_ZGA_intron, ICM_df)
  
  #combine xie class
  db[db$gene_class %in% c('one_cell','early_2cell'),]$gene_class <- 'totipotent'
  #db[db$gene_class %in% c('ICM'),]$gene_class <- 'pluripotent'
  db$gene_class <- factor(db$gene_class, levels = c('totipotent','pluripotent'))
  
  
  #db$gene_class <- factor(db$gene_class, levels = c('one_cell', 'early_2cell', 'late_2cell','ICM'))
  
  #db$gene_class <- factor(db$gene_class, levels = c('one_cell', 'early_2cell', 'pluripotent'))
  
  
  #table(db$gene_class %in% tmp2[tmp2$class == 'Up', ]$symbol)
  library(ggplot2)
  library(ggpubr)
  p1<- ggplot(db, aes(x = gene_class,y = length))+
    scale_y_log10()+
    geom_violin(position = 'dodge', draw_quantiles = c(0.5), aes(color = gene_class))+
    #scale_color_viridis_d()+
    geom_boxplot(cex=0.8,width=0.2,aes(color = gene_class))+
    theme_classic()+ 
    stat_compare_means(method = "wilcox.test",
                       label = "p.signif",##星号设置
                       #aes(label = "p.format"), #显示方式
                       label.x.npc ="left", label.y = 5, comparisons = list(c('totipotent','pluripotent')),
                       size = 5)+
    ylab('log10(intron length)')+
    xlab('gene class')+
    theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), legend.position="none",
          legend.title = element_text(size = 16),
          axis.title.y = element_text(size = 18), title = element_text(size = 16), legend.text = element_text(size = 14, lineheight = 2),
          plot.title = element_text(color="black", face="bold"),
          axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))
  
  p2 <- ggplot(db, aes(x = gene_class,y = n))+
    scale_y_log10()+
    geom_violin(position = 'dodge', draw_quantiles = c(0.5), aes(color = gene_class))+
    #scale_color_viridis_d()+
    geom_boxplot(cex=0.8,width=0.2,aes(color = gene_class))+
    theme_classic()+ 
    stat_compare_means(method = "wilcox.test",
                       label = "p.signif",##星号设置
                       #aes(label = "p.format"), #显示方式
                       label.x.npc ="left", label.y = 2, comparisons = list(c('totipotent','pluripotent')),
                       size = 5)+
    ylab('log10(intron number)')+
    xlab('gene class')+
    theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), legend.position="none",
          legend.title = element_text(size = 16),
          axis.title.y = element_text(size = 18), title = element_text(size = 16), legend.text = element_text(size = 14, lineheight = 2),
          plot.title = element_text(color="black", face="bold"),
          axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))
  
  library(cowplot)
  
  pdf(file = "../published data analysis/Figures/intron_num_len.pdf",   # The directory you want to save the file in
      width = 10, # The width of the plot in inches
      height = 8)
  
  plot_grid(p1,p2, nrow = 1)
  
  dev.off()
  library(eoffice)
  #topptx(p1 + p2,filename = "Figures/heatmap.pptx")
  library(export)
  graph2ppt(file = "../published data analysis/Figures/intron_num_len.pptx")
  
  
  # calculate gc content for gene ####
  
  library(BSgenome.Mmusculus.UCSC.mm10)
  names(longest_transcript)
  
  load("~/research_2022/research/expr_results/single_cell_AID/published data analysis/xie_markers.rda")
  
  totipotent <- c(one_cell_minor_ZGA,early_2cell_minor_ZGA)
  pluripotent <- pluripotent_gene
  
  
  totipotent<- longest_transcript[names(longest_transcript) %in% totipotent]
  pluripotent<- longest_transcript[names(longest_transcript) %in% pluripotent]
  
  library(Biostrings)
  
  totipotent.df<- do.call(rbind,(lapply(totipotent, FUN = function(x)data.frame(chr = seqnames(x), start = start(x), end = end(x), strand = strand(x)))))
  pluripotent.df<- do.call(rbind,(lapply(pluripotent, FUN = function(x)data.frame(chr = seqnames(x), start = start(x), end = end(x), strand = strand(x)))))
  
  totipotent.idx<- GRanges(seqnames = totipotent.df$chr,ranges = IRanges(start = totipotent.df$start,end = totipotent.df$end),strand = totipotent.df$strand)
  pluripotent.idx<- GRanges(seqnames = pluripotent.df$chr,ranges = IRanges(start = pluripotent.df$start,end = pluripotent.df$end),strand = pluripotent.df$strand)
  
  
  toti_seq<- getSeq(BSgenome.Mmusculus.UCSC.mm10,totipotent.idx)
  pluri_seq<- getSeq(BSgenome.Mmusculus.UCSC.mm10,pluripotent.idx)
  #rd = splitAsList(transcripts, names(transcripts))
  library(stringr)
  num_g <- str_count(toti_seq, pattern = "[GC]")
  num_g <- str_count(pluri_seq, pattern = "[GC]")
  
  toti_gc_content <- num_g / str_length(toti_seq) * 100
  pluri_gc_content <- num_g / str_length(pluri_seq) * 100
  
  hist(gc_content)
  
  gc_tab<- data.frame(gc_content = c(toti_gc_content, pluri_gc_content), 
                      gene_class = c(rep('totipotent',length(toti_gc_content)),rep('pluripotent',length(pluri_gc_content))))
  
  gc_tab$gene_class <- factor(gc_tab$gene_class, levels = c('totipotent', 'pluripotent'))
  
  pdf(file = "../published data analysis/Figures/gc_content.pdf",   # The directory you want to save the file in
      width = 10, # The width of the plot in inches
      height = 8)
  
  ggplot(gc_tab, aes(x = gc_content, fill = gene_class))+
    geom_density(linewidth = 1, alpha = 0.5)+
    theme_classic()+ 
    ylab('density')+
    xlab('gc content(%)')+guides(fill=guide_legend(title="gene class"))+
    theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), 
          legend.title = element_text(size = 16),
          axis.title.y = element_text(size = 18), title = element_text(size = 16), legend.text = element_text(size = 14, lineheight = 2),
          plot.title = element_text(color="black", face="bold"),
          axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))
  
  dev.off()
  library(eoffice)
  #topptx(p1 + p2,filename = "Figures/heatmap.pptx")
  library(export)
  graph2ppt(file = "../published data analysis/Figures/gc_content.pptx")
  
  
  
  
  
}


