library(PopGenome)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(rstatix)

# if you want to calculate the population genome statistics for all snps, please use core.vcf file output by snippy
# if you want do for non-recombination snps, please delete the recombination sites in core.vcf
# you need to name the gff file and the vcf file with the same prefix
# important note: if there is quote in gff file, you need to replace them first; and then, you need to delete the sequence content in the gff, this content is start with "##FASTA"
lineages = c('major.lineage', 'lineage2', 'lineage3', 'lineage5')
pics = list()
alldf = data.frame()

l = 'major.lineage'
wkdir = paste0('/home/qguo/projects/baoman/single/single_all/', l, '/nonrecomb_cgsnp/popgenome')
setwd(wkdir)
genenames = get_gff_info(gff.file = 'nonrecomb/gff/core.gff', extract.gene.names = T, chr = 'gnl|Prokka|EBJODJIK_1', feature = T)
df = data.frame(gene = genenames)

genome.class <- readData("nonrecomb/vcf/",format="VCF", gffpath = 'nonrecomb/gff/') 
genes <- splitting.data(genome.class, subsites="gene")
genes <- F_ST.stats(genes) 
genes <- neutrality.stats(genes) 
tajima_nonrecomb = as.vector(genes@Tajima.D)
df$gene_region = genes@region.names
df$nonrecomb_tajima = tajima_nonrecomb

rm('genome.class')
rm('genes')
genome.class <- readData("all/vcf/",format="VCF", gffpath = 'all/gff/') 
genes <- splitting.data(genome.class, subsites="gene")
genes <- F_ST.stats(genes) 
genes <- neutrality.stats(genes) 
tajima_all = as.vector(genes@Tajima.D)
df$all_tajima = tajima_all
df$CHROM = 'EBJODJIK_1'
write.table(df, 'popgenome_tajima_res.txt', sep = '\t', quote = F, row.names = F)

dat_for_plot = subset(df, select = c(nonrecomb_tajima, all_tajima, CHROM))
dat_melt = melt(dat_for_plot, id.vars = 'CHROM')
dat_melt = dat_melt[!is.na(dat_melt$value),]
dat_melt$variable = as.character(dat_melt$variable)
dat_melt$variable[dat_melt$variable == 'nonrecomb_tajima'] = 'Non-recombination SNVs'
dat_melt$variable[dat_melt$variable == 'all_tajima'] = 'All SNVs'
dat_melt$variable = factor(dat_melt$variable, levels = c('Non-recombination SNVs', 'All SNVs'))

p = ggplot(dat_melt, aes(x = variable, y = value, fill = variable)) + 
  geom_violin() + 
  #geom_jitter(width =0.2,shape = 21,size=2.5)+
  labs(x = '', y = "Tajima'D") + 
  #title() + 
  theme_classic()+theme_bw() +   
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(size = 14, colour = 'black'),
        axis.title.x = element_text(size = 14, colour = 'black'),
        axis.title.y = element_text(size = 14, colour = 'black'),
        legend.position = 'none',
        panel.grid=element_blank(),
        panel.border = element_rect(colour = "black", size=0.5))
ggsave(p, filename = 'all.and.nonrecomb.tajimaD.pdf', width = 4.5, height = 5, dpi = 1200)

p1 = ggplot(df[!is.na(df$nonrecomb_tajima),], aes(x = CHROM, y = nonrecomb_tajima)) + 
  geom_violin() + 
  #geom_jitter(width =0.2,shape = 21,size=2.5)+
  labs(x = 'Gene', y = "Tajima'D", title = 'Non-recombination SNVs') + 
  #title() + 
  theme_classic()+theme_bw() +   
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y = element_text(size = 12, colour = 'black'),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 14, colour = 'black'),
        axis.title.y = element_text(size = 14, colour = 'black'),
        legend.position = 'none',
        panel.grid=element_blank(),
        panel.border = element_rect(colour = "black", size=0.5))

p2 = ggplot(df[!is.na(df$all_tajima),], aes(x = CHROM, y = all_tajima)) + 
  geom_violin() + 
  labs(x = 'Gene', y = "Tajima'D", title = 'All SNVs') + 
  theme_classic()+theme_bw() +   
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y = element_text(size = 12, colour = 'black'),
        axis.text.x = element_blank(),
        axis.title.x = element_text(size = 14, colour = 'black'),
        axis.title.y = element_text(size = 14, colour = 'black'),
        legend.position = 'none',
        panel.grid=element_blank(),
        panel.border = element_rect(colour = "black", size=0.5))

p = plot_grid(p1, p2, nrow = 1)
ggsave(p, filename = 'tajimaD.pdf', width = 5, height = 6, dpi = 1200)  



pics[[i]] = p

df$lineage = l

alldf = rbind(alldf, df)
break


##################################### only used for compare
pout = plot_grid(pics[[1]], pics[[2]], pics[[3]], pics[[4]], nrow = 2, align = 'hv', axis = 'btlr')
ggsave(pout, filename = '/home/qguo/projects/baoman/single/single_all/all_lineages_tajimaD.pdf', width = 10, height = 10, dpi = 1200) 
p = ggplot(alldf[!is.na(alldf$nonrecomb_tajima),], aes(x = CHROM, y = nonrecomb_tajima, fill = lineage)) + 
  geom_violin() + 
  labs(x = 'SNP type', y = "Tajima'D") + 
  theme_classic()+theme_bw() +   
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(size = 12, colour = 'black'),
        axis.title.x = element_text(size = 14, colour = 'black'),
        axis.title.y = element_text(size = 14, colour = 'black'),
        legend.position = 'none',
        panel.grid=element_blank(),
        panel.border = element_rect(colour = "black", size=0.5))
ggsave(p, filename = '/home/qguo/projects/baoman/single/single_all/all_lineages_tajimaD_comp.pdf', width = 5, height = 6, dpi = 1200)  

tajima_test = alldf[!is.na(alldf$nonrecomb_tajima),] %>% wilcox_test(nonrecomb_tajima ~ lineage)
write.table(tajima_test, '/home/qguo/projects/baoman/single/single_all/all_lineages_tajimaD_wilcox.txt', sep = '\t', quote = F, row.names = F)

# the scripts guided for plotting coding :
# 1: https://www.jianshu.com/p/0ef757e2d19f
# 2: https://www.jianshu.com/p/a2c573002a25

outdf1 = alldf[!is.na(alldf$nonrecomb_tajima),] %>% group_by(lineage) %>% mutate(mean_nonrecomb = mean(nonrecomb_tajima), median_nonrecomb = median(nonrecomb_tajima))
outdf2 = subset(outdf1, select = c(lineage, mean_nonrecomb, median_nonrecomb))
outdf2 = outdf2[!duplicated(outdf2),]

#tajima_test = merge(tajima_test, outdf)

write.table(outdf2, '/home/qguo/projects/baoman/single/single_all/all_lineages_tajimaD_mean_median.txt', sep = '\t', quote = F, row.names = F)
  
###################################### not run 


p = ggboxplot(df[!is.na(df$nonrecomb_tajima),], x = "CHROM", y = "nonrecomb_tajima",add = "median_iqr") + 
  stat_pvalue_manual(chao1_stat_test, label = "p.adj.signif", hide.ns=T, tip.length = 0,  size = 10) + 
  theme_bw()+ labs(x ='', y = 'Chao1') + 
  theme(plot.title = element_text(hjust=0.5), axis.text.y = element_text(size = 12, colour = 'black'), 
        axis.title = element_text(size = 13), axis.text.x = element_text(size=13, vjust = 0.6, colour = 'black'),
        legend.position='none',legend.title = element_text(''), 
        axis.line = element_line(colour='black', size=0.2), 
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "white", size=0.2))

