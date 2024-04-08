library(circlize)
library(grid)

setwd('/Users/guoqian/projects/baoman_single/single_article/article/draft/code/card_vf_ld')
ref_card_vf <- read.table('card_vf_snp.txt', colClasses = c("character", "numeric", "numeric", "numeric","character","character"), 
                          sep = "\t",header=T)
seq_link <- read.table('card_vf_snp_linkage.txt', sep = "\t",header=F,check.names=FALSE, fileEncoding = 'utf-8' )

col <- c("#0066CC","#E41A1C","#0066CC","#E41A1C",rep("#0066CC",5),rep("#E41A1C",2),rep('#0066CC',10),rep("#E41A1C",2),rep("#0066CC",11))
anno_col <- read.table('anno_col.txt', sep = "\t",header=T,check.names=FALSE, fileEncoding = 'utf-8' )
anno_color = paste("#",anno_col$anno_color,sep='')

circos.par(gap.degree = 2, "track.height" = 0.1,start.degree = 30)
circos.genomicInitialize(ref_card_vf,track.height = 0.03, plotType = NULL)
circos.trackPlotRegion(ylim = c(0,1), bg.border = col, bg.col = col, panel.fun = function(x, y) {
    gene = get.cell.meta.data("sector.index")
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.text(mean(xlim), 1.3, labels = gene,adj = c(-0.2, 0), cex = 0.5,facing = "clockwise", niceFacing = TRUE)
    circos.axis(h = "bottom", labels = NULL, sector.index = gene, direction = "outside")
}, track.height = 0.05)


circos.trackPlotRegion(ylim = c(0,1), bg.border = anno_color, bg.col = anno_color, panel.fun = function(x, y) {
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
}, track.height = 0.05)

circos.genomicTrack(ref_card_vf,numeric.column = 4,
                    panel.fun = function(region, value, ...) {
                        circos.genomicPoints(
                            region, value, pch = 16, 
                            cex = 0.4, col = '#8A2BE2')
                    })


for ( i in 1:length(seq_link[,1])) {
    circos.link(
        seq_link[i,1], get.cell.meta.data("cell.xlim", sector.index = seq_link[i,1]),
        seq_link[i,2], get.cell.meta.data("cell.xlim", sector.index = seq_link[i,2]),
        col = paste("#",seq_link[i,3],"30",sep=''), border = NA)
}


                      