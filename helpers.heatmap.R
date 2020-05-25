############################################################
###### HEATMAP
############################################################
###### SWITCH EFFECTS: 
################################################
## TODO? To the effect of trait; 
## TO MOST POSITIVE EFFECT

library(RColorBrewer)

#check.size(library(heatmaply))
flipeffects <- function(row) {
  
  if (abs(min(row)) > max(row)){
    #print("flip")
    return(-1*row)
  } else{
    #print("don't flip")
    return(1*row)
  }
}


rescale_amplitudes<-function(col){
  colrescale<-col
  colrescale[colrescale>0]=(col[col>0] - 0)/(max(col) - 0)
  colrescale[colrescale<0]=-1*(abs(col[col<0]) - 0)/(max(col) - 0)
  return(colrescale)
}


make_heatmap_plot <- function(df_snp_p = data$df_snp_p,vct_snp_info = data$df_snp_info,
                              df_ecg_unadjusted, pvaltreshold = 0.3) {
  print("heatmap")
  #print(pvaltreshold)
  colorscale = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  mx <- data.matrix(df_snp_p)
  rownames(mx)<-1:nrow(mx)
  pvaltreshold_selection = rowMaxs(abs(mx))> -log10(pvaltreshold)
  mx <- mx[pvaltreshold_selection,]
  
  ### TODO: HOW TO SHOW THIS AS MESSAGE?? 
  validate(need(nrow(mx)>1, 'Need more than 1 variant to make a heatplot, provide more variants or change P value treshold'))
  #validate(need(nrow(mx)<=100, 'Please have less than 100 variant, provide less SNPs, stricter P-value or smaller region; More filtering options will come, suggestions welcome .'))
  ####
  
  mx <- t(apply(mx, 1, flipeffects))
  mx <- t(apply(mx, 1, rescale_amplitudes))

  #rownames <- stringr::str_wrap(rownames,width=10, ",")
  hc <- hclust(dist(mx))
  # dend <- as.dendrogram(hc)
  # plot(dend, horiz=TRUE)
  mx <- mx[hc$order,]
  
  rownames <-paste0(vct_snp_info[pvaltreshold_selection,]$SNP," (",vct_snp_info[pvaltreshold_selection,]$Gene,")")
  rownames <- stringr::str_wrap(rownames,width=15)
  rownames <- rownames[hc$order]
  # hm <- heatmaply(mx,labRow=rownames,
  #                 showticklabels =c(FALSE,TRUE),color=colorscale, Colv=FALSE,
  #                 fontsize_row=7, 
  #                 margins = c(60,100,40,20))

  hm <- heatmaply(mx,labRow=rownames,
                  showticklabels =c(FALSE,TRUE),color=colorscale, Colv=FALSE,dendrogram="none",
                  fontsize_row=7, 
                  margins = c(60,100,40,20)) %>% layout(height=max(500,25*nrow(mx)))
  
  hm$order <- as.numeric( rownames(mx)) # hc$order
  # hm
  # hm <- heatmap(mx, colors = colorRamp(c("red", "green")))
  #out <- list()
  #out$hm <- hm
  #out$n <- nrow(mx)
  #return(hm)
  return(hm)
}


