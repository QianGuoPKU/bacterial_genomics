rm(list = ls())

library("genoPlotR")
library("seqinr")
library("geiger")
library("reshape2")
library("dplyr")
library(ggplot2)
library(DescTools)
library(scales)


######################### method 1: 

lineages = c('major.lineage', 'lineage2', 'lineage3', 'lineage5')
pics = list()


outgroup = read.csv('/home/qguo/projects/baoman/single/single_all/major.lineage/outgroup/pan_gene_pres_abs.txt',
                    sep = '\t', check.names = F, stringsAsFactors = F)
alldf = outgroup

for (i in c(1:length(lineages))){
    l = lineages[1]
    wkdir = paste0('/home/qguo/projects/baoman/single/single_all/', l, '/nonrecomb_cgsnp/selection/poisson')
    setwd(wkdir)

    indata = read.csv("all_snp_for_poison.txt", check.names = F, stringsAsFactors = F, sep = '\t')
    data = indata[indata$recomb != 'recombination',]
    data = data[data$FTYPE == 'CDS',]
    #data <- data[-which(data$LOCUS_TAG == ""), ]
    #data = data[data$EFFECT %in% c('synonymous_variant', 'missense_variant'),]
    data$nsub = 1
    
    data_dcast = dcast(data, LOCUS_TAG~nsub, sum)
    colnames(data_dcast)[2] = l
    nsub_table = table(data_dcast[,2])
    nsub_table = data.frame(nsub = c(0,as.numeric(names(nsub_table))), 
                            ngene = c(0, as.vector(nsub_table)))
    
    nsub_table$ngene[1] = 3771 - sum(nsub_table$ngene)
    #nsub_table = data.frame(nsub = as.numeric(names(nsub_table)), ngene = as.vector(nsub_table))
    
    
    fm1 = glm(ngene~nsub, data=nsub_table, family=poisson)
    pred <- predict(fm1, type="response", se.fit = TRUE)
    df2 = cbind(nsub_table, pred = pred$fit)
    df2 = cbind(df2, se = pred$se.fit) 
    #df2 = cbind(df2, ucl=exp(df2$pred + 1.96*df2$se))
    #df2 = cbind(df2, lcl=exp(df2$pred - 1.96*df2$se))
    df2 = cbind(df2, ucl=df2$pred + 1.96*df2$se)
    df2 = cbind(df2, lcl=df2$pred - 1.96*df2$se)
    df2$ngene_log = log10(df2$ngene + 0.05)
    df2$ucl_log = log10(df2$ucl)
    df2$lcl_log = log10(df2$lcl)

    write.table(df2, 'poisson_analysis_out.txt', sep = '\t', quote = F, row.names = F)

    p = ggplot() +
        geom_bar(data=df2, aes(x=nsub,y=ngene_log),
                 stat="identity",position = "dodge") +
    geom_point(data=df2,aes(x=nsub,y=ucl_log), colour = 'orange') +
        geom_line(data=df2,aes(x=nsub,y = ucl_log), colour = 'orange') + 
        geom_point(data=df2,aes(x=nsub,y=lcl_log), colour = 'black') +
        geom_line(data=df2,aes(x=nsub,y = lcl_log), colour = 'black') +
        labs(x = 'Number of mutations per gene', y = 'Number of genes') + 
        theme_classic()+theme_bw() +   
        theme(plot.title = element_text(hjust=0.5),
                                             axis.text.y = element_text(size = 12, colour = 'black'),
                                             axis.text.x = element_text(size = 12, colour = 'black'),
                                             axis.title.x = element_text(size = 14, colour = 'black'),
                                             axis.title.y = element_text(size = 14, colour = 'black'),
                                             legend.position = 'none',
                                             panel.grid=element_blank(),
                                             panel.border = element_rect(colour = "black", size=0.5)) + 
        scale_y_continuous(limits = c(0, max(df2$ngene_log)), labels = math_format(10^.x), expand = expansion(mult = c(0, 0.005))) 
        #coord_cartesian(ylim = c(0.1, max(df2$ngene_log))) #scale_y_continuous(limits = c(-0.1, max(df2$ngene_log)))
    # scale_y_log10()
    
    ggsave(p, filename = 'log_poisson.pdf', width = 10, height = 10, dpi = 1200)
    pics[[i]] = p
    
    alldf = merge(alldf, data_dcast, by.x = 'Gene', by.y = 'LOCUS_TAG', all.x = T)
    write.table(alldf, 'outgroup_nsub.txt', sep = '\t', quote = F, row.names = F)
    break
}

pout = plot_grid(pics[[1]], pics[[2]], pics[[3]], pics[[4]], nrow = 2, align = 'hv', axis = 'btlr')
ggsave(pout, filename = '/home/qguo/projects/baoman/single/single_all/all_lineages_poisson.pdf', width = 10, height = 10, dpi = 1200) 
p

write.table(alldf, '/home/qguo/projects/baoman/single/single_all/all_lineages_poisson.txt', sep = '\t', quote = F, row.names = F)

############################## not run #######################################################3
ggplot() +
    geom_bar(data=df2, aes(x=nsub,y=ngene),
             stat="identity",position = "dodge") +
    geom_point(data=df2,aes(x=nsub,y=ucl, color = 'orange')) +
    geom_line(data=df2,aes(x=nsub,y = ucl, color = 'orange')) + 
    geom_point(data=df2,aes(x=nsub,y=lcl, color = 'black')) +
    geom_line(data=df2,aes(x=nsub,y = lcl, color = 'black')) +     
    theme_classic()+theme_bw() + scale_y_log10() 
#+ scale_y_continuous( limits = c(-1, nsub_table$ngene[1]))



outgroup_nsub = merge(outgroup, data_dcast, by.x = 'Gene', by.y = 'LOCUS_TAG', all.x = T)
write.table(outgroup_nsub, 'outgroup_nsub.txt', sep = '\t', quote = F, row.names = F)


############ method 2:
subAnnotation <- indata[which(indata$FTYPE == "CDS"),]
#subAnnotation <- subAnnotation[-which(subAnnotation$GENE == ""), ]
subAnnotation$SNP = 1
# extra and reshape NS df
rawNS <- subAnnotation[,c('EFFECT','LOCUS_TAG', 'Gene_length', 'SNP')]
rawNS$EFFECT[which(rawNS$EFFECT!='synonymous_variant')] <-'nonsynonymous_variant'
gene_len = rawNS %>% group_by(LOCUS_TAG) %>% summarise(Gene_length = max(Gene_length))
rawNS$Gene_length = gene_len$Gene_length[match(rawNS$LOCUS_TAG, gene_len$LOCUS_TAG)]

rawNS <- dcast(rawNS,LOCUS_TAG +Gene_length~EFFECT,sum)
rownames(rawNS) = rawNS$LOCUS_TAG

# Total number of non-synonymous and synonymous substitutions for each gene
rawNS$nsub = rawNS$nonsynonymous_variant + rawNS$synonymous_variant

# Substitution rate per site
subRate = sum(rawNS$nsub)/sum(rawNS$Gene_length)
rawNS$expected_substitutions = subRate*rawNS$Gene_length
# Calculate the p-value for each gene
subRatePVal = sapply(1:length(rawNS$nsub),function(i){
    t = poisson.test(rawNS$nsub[i],r=(rawNS$expected_substitutions[i]))
    t$p.value
}) 
#wald_CI = sapply(1:length(nSubs),function(i){
#ci = PoissonCI(nSubs[i], n = 1, conf.level = 0.95, sides = "two.sided",method = "wald")})
rawNS$poisson_p.value = subRatePVal
rawNS$flag = 'less'
rawNS$flag[rawNS$nsub > rawNS$expected_substitutions] = 'more'
write.table(rawNS, 'method2_res.txt', sep = '\t', quote = F, row.names = F)
    
########################################################################3
############### not run

# Find all coding substitutions
subAnnotation <- variants[which(variants$FTYPE == "CDS"),]
subAnnotation$SNP = 1
# extra and reshape NS df
codingSubs <- subAnnotation[,c('EFFECT','LOCUS_TAG', 'Gene_length' , 'SNP')]
codingSubs$EFFECT[which(codingSubs$EFFECT!='synonymous_variant')] <-'nonsynonymous_variant'
codingSubs <- dcast(codingSubs,`LOCUS_TAG` + `Gene_length`~EFFECT,sum)
codingSubs$nSubs = codingSubs$nonsynonymous_variant + codingSubs$synonymous_variant

# Substitution rate per site
subRate = sum(codingSubs$nSubs)/sum(codingSubs$Gene_length)
codingSubs$predicted_nSubs = subRate * codingSubs$Gene_length

# hist
nSubs_table = table(codingSubs$nSubs)
nSubs_table = data.frame(nSubs = names(nSubs_table), nGenes = as.vector(nSubs_table))
predicted_nSub_nGene_table = aggregate(`predicted_nSubs` ~ nSubs, data = codingSubs, mean)
predicted_nSub_nGene_table = merge(predicted_nSub_nGene_table, nSubs_table)

# 两种计算均值的方法，一种是用poissonCI中exp的计算方法；一种是用rate * length的方法
ngenes = sum(predicted_nSub_nGene_table$nGenes)
nsubs = sum(predicted_nSub_nGene_table$nSubs * predicted_nSub_nGene_table$nGenes)
predicted_nSub_nGene_table$exp = dpois(as.numeric(predicted_nSub_nGene_table$nSubs), nsubs/ngenes) * ngenes

# exp & predicted_nSub_nGene_table
plot(predicted_nSub_nGene_table$nSubs, predicted_nSub_nGene_table$exp)
plot(predicted_nSub_nGene_table$nSubs, predicted_nSub_nGene_table$exp)


##################3 not run


p1<-ggplot(data = codingSubs,aes(x=nSubs))+
    geom_histogram(bins = 300)


# poisson distribution of 
# Calculate the p-value for each gene
subRatePVal = sapply(1:length(nSubs),function(i){
    t = poisson.test(nSubs[i],r=(subRate*geneLengths[i]))
    t$p.value
}) 
#wald_CI = sapply(1:length(nSubs),function(i){
    #ci = PoissonCI(nSubs[i], n = 1, conf.level = 0.95, sides = "two.sided",method = "wald")})
rawNS$poisson_p.value = subRatePVal
rawNS$substitutions = nSubs
expected_subs = sapply(1:length(nSubs),function(i){h=(subRate*geneLengths[i])})
rawNS$expected_substitutions = expected_subs
sigNS = rawNS[which(rawNS$poisson_p.value<0.01),]
sigNS = sigNS[which(sigNS$substitutions>0),]
plot = sigNS[order(sigNS$substitutions),]
write.table(file='all_snp_for_poison_rawNS.txt',rawNS,sep='\t',row.names=TRUE)
write.table(file='all_snp_for_poison_plot.txt',plot,sep='\t',row.names=TRUE)

#plot=cbind(GENE=row.names(plot), plot)
#row.names(plot)=NULL
#plot <- melt(plot,id.vars='GENE', measure.vars =c('substitutions',
            #'expected_substitutions'),variable.name = 'EFFECT',value.name='substitutions_number')


plot(plot$substitutions, type = "l",ylab="",col='red',xlab = "",xaxt = 'n',main='all_snp_for_poison')
lines(plot$expected_substitutions, col=rgb(89, 118, 186, 100, maxColorValue=186))
mtext(side = 2, text = "Number of substitutions", line = 2.2, cex = 1.1)
mtext(side = 1, text = "Gene", line = 0.5, cex = 1.1)
labs <- c("substitutions", "expected_substitutions")
legend("topright", legend = labs, cex = 0.8, lty = 1, lwd = 2, 
       col = c("red", rgb(89, 118, 186, 100, maxColorValue=186)), 
       inset = 0.01, horiz = TRUE, box.col = "white")

