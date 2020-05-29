require(raster)
require(ggplot2)

setwd("/Users/niek/repos/ecgenetics/")
source("global.R")
source("helpers.R")
source("helpers.heatmap.R")
source("helpers.regionalplot.R")
source("get_nearest_gene.r")




   input= unique(dftable1$SNP)
   query <- process_user_input(input,mapping.proteincoding)
   tabix_query <- get_tabix_query(query,df.static.pos,df.static.rsid)
   data <- extract_multiple_variants(tabix_query,dir_data,
                                     f.data_p="unadjusted.logP.outfile.tsv.chr.gz.tophits.gz",
                                     f.data_beta="",
                                     f.data_se="",
                                     f.data.index="unadjusted.logP.outfile.index.tsv.gz.tophits.gz"
   )


   

d = as.matrix(data$df_snp_p)

rescale_amplitudes<-function(col){
   colrescale<-col
   colrescale[colrescale>0]=(col[col>0] - 0)/(max(col) - 0)
   colrescale[colrescale<0]=-1*(abs(col[col<0]) - 0)/(max(col) - 0)
   return(colrescale)
}
flipeffects <- function(row) {
   
   if (abs(min(row)) > max(row)){
      print("flip")
      return(-1*row)
   } else{
      print("don't flip")
      return(1*row)
   }
}
d <- t(apply(d, 1, flipeffects))
d <- t(apply(d, 1, rescale_amplitudes))
d <- d*25

library(gplots)
library(pheatmap)   
library(RColorBrewer)
colorscale= colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

pm <- pheatmap(d,color=colorscale, cluster_cols=FALSE,fontsize=5)
d <- d[pm$tree_row$order,]
d = data.frame(rasterToPoints(raster(d)))

d$elev = d$layer + d$y * 1200

p = ggplot() + 
   labs(x = '', y = NULL, title = '') + 
   theme_minimal() +
   theme(axis.ticks = element_blank(),
         panel.grid = element_blank(),
         axis.text = element_blank(),
         plot.background = element_rect(fill = "black"),
         plot.title = element_text(colour = 'white', size = 76, hjust = .1),
         axis.title.x = element_text(colour = 'white', size = 43))

for(g in unique(d$y)){
   dat = subset(d, y == g)
   p = p + geom_polygon(data=dat, aes(x, elev), fill='black') + 
      geom_line(data=dat, aes(x, elev), col='white', size=.5)
}
p
# save 10 x 15

#ggsave(filename = 'unknown_pleasures.png', plot = p, w=8, h = 8)
#browseURL('unknown_pleasures.png')

