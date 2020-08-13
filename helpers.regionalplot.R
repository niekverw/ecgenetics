#regional
library(scales)
loadrecombination<-function(chr){
  #recomb<-data.frame(fread(paste0("/data/p255753/References/ALL_1000G_phase1integrated_v3_impute/genetic_map_chr",chr,"_combined_b37.txt")))
  recomb<-data.frame(fread(paste0(datadir,"/recombination_rates/genetic_map_chr",chr,"_combined_b37.txt")))  
  return(recomb)  
}

get_ld <- function(data,LD=T){
  if(LD==F){
    data$df_snp_info$LD_BIN <- "NA"
    data$df_snp_info$LD_BIN = factor(data$df_snp_info$LD_BIN, levels = c("1.0-0.8", "0.8-0.6", "0.6-0.4", "0.4-0.2", "0.2-0.0", "NA"))
    return(data)
  }
  rsid = data$df_snp_info[ data$max.snp.index[,1][1],"SNP"]
  print(rsid)
  print(data$max.snp.index)
  cat("get Ld, sorry slow.. ")
  pops="EUR"
  ld_command = paste0("https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var=", rsid,"&pop=",pops,"&r2_d=r2&token=70722d8a3958")
  print(ld_command)
  ####### EDIT
  df.ld = as.data.frame(data.table::fread(ld_command))
  df.ld = df.ld[,c("RS_Number", "R2")]
  #df.ld = data$df_snp_info[,c("SNP","INFO")]
  ######
  names(df.ld) <- c("SNP","R2")
  data$df_snp_info <- plyr::join(data$df_snp_info, df.ld,by="SNP",match="first")
  
  data$df_snp_info$LD_BIN <- cut(data$df_snp_info$R2,
                           breaks=c(0,0.2,0.4,0.6,0.8,1.0),
                           labels=c("0.2-0.0","0.4-0.2","0.6-0.4","0.8-0.6","1.0-0.8"))
  data$df_snp_info$LD_BIN <- as.character(data$df_snp_info$LD_BIN)
  data$df_snp_info$LD_BIN[is.na(data$df_snp_info$LD_BIN)] <- "NA"
  data$df_snp_info$LD_BIN = as.factor(data$df_snp_info$LD_BIN)
  data$df_snp_info$LD_BIN = factor(data$df_snp_info$LD_BIN, levels = c("1.0-0.8", "0.8-0.6", "0.6-0.4", "0.4-0.2", "0.2-0.0", "NA"))

  return(data)
}

regional_plot <- function(data,for_plotly=F,LD=FALSE){
  data$df_snp_info$BP <- as.numeric(data$df_snp_info$BP)
  chr_in = unique(data$df_snp_info$CHR)
  
  
  gene_sub <- mapping.proteincoding[mapping.proteincoding$chromosome_name %in% chr_in,]
  colnames(gene_sub) = c("GENE_ID","GENE_NAME","BAND", "CHR", "TRX_START", "TRX_END", "TYPE",  "STRAND")
  gene_sub$LENGTH <- gene_sub$TRX_END -  gene_sub$TRX_START 
  
  start=data$query$startpos
  end=data$query$endpos
  # reading in gene data
  gene_sub = subset(gene_sub, gene_sub$TRX_START > (start-50000 ))
  gene_sub = subset(gene_sub, gene_sub$TRX_END < (end+50000))
  gene_sub = gene_sub[,c("GENE_NAME","TRX_START","TRX_END")]
  gene_sub = reshape2::melt(gene_sub,id.vars = "GENE_NAME")
  gene_sub$y_value = as.numeric(as.factor(gene_sub$GENE_NAME))
  plot_lab = subset(gene_sub, gene_sub$variable == "TRX_END")
  
  dfrecomb<-loadrecombination(chr_in)
  dfrecomb <- subset(dfrecomb, dfrecomb[,1] > start & dfrecomb[,1] < end)
  Pf=max(abs(dfrecomb[,2]))/max(abs(as.numeric(data$df_snp_info$maxLogP)),na.rm=TRUE)
  dfrecomb[,2] = dfrecomb[,2]/Pf
  
  trait_to_plot = as.numeric(data$max.snp.index[,2][1])
  rsid = data$df_snp_info[ data$max.snp.index[,1][1],"SNP"]
  
  data <- get_ld(data,LD=LD)
  
  df_assoc_plot <- data$df_snp_info[,c("i","SNP","BP","maxLogP","LD_BIN")] #
  df_assoc_plot$maxLogPTrait <- abs(data$df_snp_p[,trait_to_plot])
  df_assoc_plot = melt(df_assoc_plot, id.vars=c("i","SNP","BP","LD_BIN")) #,"LD_BIN"
  

  plt_genes = ggplot2::ggplot(gene_sub, ggplot2::aes_string(x = "value", y = "y_value")) +
    ggplot2::geom_line(ggplot2::aes_string(group = "GENE_NAME"), size = 2) + ggplot2::theme_bw() +
    ggplot2::geom_text(data = plot_lab, ggplot2::aes_string(x = "value", y = "y_value", label = "GENE_NAME"),
                       colour="red", fontface = "bold",
                       hjust = -0.1,vjust = 0.3, size = 2.5) +
    ggplot2::theme(axis.title.y = ggplot2::element_text(color = "white", size = 28),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank()) + 
    ggplot2::xlab(paste0("Position on chromosome ", chr_in, "")) +
    ggplot2::coord_cartesian(xlim = c(start,end), ylim = c(0,(max(gene_sub$y_value)+1)))
  
  plt_regional = ggplot2::ggplot(df_assoc_plot, aes(x = BP, y = value, colour=factor(LD_BIN) )) +
    ggplot2::geom_point(aes(shape=factor(variable)),size=2,alpha=0.75) + 
    ggplot2::theme_bw() + theme(legend.position = "none") +
    ggplot2::scale_colour_manual(name = "LD (r2)",
      values = c("1.0-0.8" = "#000004FF", "0.8-0.6" = "#51127CFF", "0.6-0.4" = "#B63679FF",
                 "0.4-0.2" = "#FB8861FF", "0.2-0.0" = "#5DC863FF", "NA" = "grey"), drop = FALSE) +
    #scale_color_manual(values=c("grey", "black")) +
    scale_shape_manual(values=c(1,19),guide=FALSE)+ #shapes
    ggplot2::xlab(paste0("Position on chromosome ", chr_in, "")) +
    ggplot2::ylab("-log10(p-value)") + 
    ggplot2::coord_cartesian(xlim = c(start, end), ylim = c(0,max(data$df_snp_info$maxLogP))) + 
    ggplot2::geom_line(data=dfrecomb,aes(x=position,y=COMBINED_rate.cM.Mb.),colour="blue",alpha=0.5) +
    geom_hline(yintercept = -log10(5*10^-8), linetype="dashed",colour="blue")+
    geom_hline(yintercept = -log10(1*10^-5), linetype="dashed",colour="grey")
    
  
  if(for_plotly==T){
    out = list()
    out$plt_regional <- plt_regional
    out$plt_genes <- plt_genes
    return(out)
  }
  
  ggpubr::ggarrange(plt_regional, plt_genes, heights = c(3,1), nrow = 2, ncol = 1)
  #common.legend = TRUE, legend = "right"
  }
#regional_plot(data,for_plotly = F)

regional_plotly <- function(data){
  plots = regional_plot(data,for_plotly = T)
  #return(ggplotly(plots$plt_regional,source="source_rp")
    p1 <- ggplotly(plots$plt_regional,source="source_rp")
    #attr(p1,"source") <- "asd"
    p2 <- ggplotly(plots$plt_genes,source="asd")
    subplot(p1, p2, nrows = 2, margin = 0.09, heights = c(0.8, 0.2))
    #p1
}

#regional_plotly(data)
