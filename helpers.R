

check.size <- function(x){
  format(object.size(x),units = "auto")
}


make_ecg_plot <- function(vct_snp_p,vct_snp_beta="",vct_snp_se="",vct_snp_info,df_ecg_stats=df_ecg_unadjusted,invert=FALSE){
  #print(nrow(vct_snp_p))
  dfecg <- cbind(as.data.frame(df_ecg_stats),
                 t(vct_snp_p)
                 )
  names(dfecg)[5:ncol(dfecg)] <- vct_snp_info$SNP
  dfecg = melt(dfecg, id.vars=names(dfecg)[1:4])
  names(dfecg)[5:6] <- c("SNP","logP")
  dfecg$logP <- as.numeric(dfecg$logP)
  if(invert==TRUE) dfecg$logP <- -1*dfecg$logP
  
  if(is.data.frame(vct_snp_beta) & is.data.frame(vct_snp_se) & length(vct_snp_p)==500 ) {
    dfecg$BETA <- as.numeric(vct_snp_beta)
    if(invert==TRUE) dfecg$BETA <- -1*dfecg$BETA
    dfecg$cil <-dfecg$BETA - 1.96*as.numeric(vct_snp_se)
    dfecg$ciu <-dfecg$BETA + 1.96*as.numeric(vct_snp_se)
    
    # 
    BETAf=max(abs(c(dfecg$cil,  dfecg$ciu)) )/max(abs(as.numeric(dfecg$logP)),na.rm=TRUE)
    dfecg$BETA  <- dfecg$BETA /BETAf
    dfecg$cil  <- dfecg$cil /BETAf
    dfecg$ciu  <- dfecg$ciu /BETAf
    
    # 
  }
  
  # dfecg <- cbind(as.data.frame(df_ecg_stats),
  #                vct_snp_p=as.numeric(as.character(unlist(vct_snp_p))))
  
  dfecg$time <- 2*(1:500-250)
  dfecg$SNP <- unname(unlist(as.data.frame(sapply( vct_snp_info$SNP, function(x) rep(x,500)))))
  
  Pf=max(abs(dfecg$mean0))/max(abs(as.numeric(dfecg$logP)),na.rm=TRUE)
  Pfdefault=max(abs(dfecg$mean0))/8
  Pf = min(Pfdefault,Pf)
  
  
  #print(Pf)
  #print(Pfdefault)
  dfecg$mean0  <- dfecg$mean0 /Pf
  #print(dfecg[order(dfecg$time),])
  #### FOR MULTIPLE LINES IN ONE PLOT, DONT KNOW YET HOW THIS IS DONE BEST WITH USER INPUT. 
  # ecgplot <- ggplot()+
  #   geom_line(data=dfecg, aes(x = time, y= logP ,group=SNP,colour = SNP,
  #                             text = paste0('Time(ms): ',time,
  #                                           '<br>P-value: ',formatC(10^(-abs(logP)),digits=3),
  #                                           '<br>-Log10(P-value): ',abs(logP) ) ))  +
  #   geom_line(data=dfecg,aes(x = time,y=mean0),linetype="dotted",colour="black") +
  #       scale_color_manual(name="Groups",values=c("red", "blue","#00AFBB", "#E7B800", "#FC4E07"))+
  #       labs(title=paste(vct_snp_info$SNP,collapse = ", " ), subtitle=vct_snp_info$Gene)
  
# 
  ecgplot <- ggplot(dfecg,aes(x=time))+
    geom_line(data=dfecg,colour="red", aes(x = time, y= logP ,
                                           text = paste0('Time(ms): ',time,
                                                         '<br>P-value: ',formatC(10^(-abs(logP)),digits=3),
                                                         '<br>-Log10(P-value): ',abs(logP) ),
                                           group=1  ))  +
    geom_line(data=dfecg,aes(x = time,y=mean0),linetype="dotted",colour="black") +
    scale_color_manual(name="Groups",values=c("red", "blue"))+
    ylab("Signed -log10(P-value)") + xlab("Time (ms)") +
    labs(title=vct_snp_info$SNP, subtitle=vct_snp_info$Gene)
  
  if(is.data.frame(vct_snp_beta) & is.data.frame(vct_snp_se) & length(vct_snp_p)==500 ) {
    ecgplot <- ecgplot +  geom_ribbon(aes(ymin=cil, ymax=ciu,x = time), linetype=2, alpha=0.1) + 
    geom_line(colour="grey",aes(y=BETA)) + 
    scale_y_continuous(
    
      sec.axis = sec_axis(~ . *BETAf, name = "Effect (per standard deviation) and CI")
    )
  }
   
  #ecgplot
  
  return(ecgplot) 
}


# extract_multiple_variants <- function(rsids,df.static_subset,dir_data) {
#   ptm <- proc.time()
#   
#   #df.static[ df.static$SNP %in% rsids,] <- doesnt work anymore if setkey()??
#   vct.missingsnps = df.static_subset[is.na(df.static_subset$uniqid)]$SNP
#   
#   df.static_subset = na.omit(df.static_subset)
#   pos = df.static_subset$i
#   chr = paste0("chr",df.static_subset$CHR)
#   tabix_queries = paste0(chr,":",pos,"-",pos)
#   tabix_f = df.static_subset[ df.static_subset$SNP %in% rsids]$files
# 
#   df_snp_info <- data.frame()
#   df_snp_p <- data.frame()
#   if(length(tabix_f)==0){
#     print("invalid snps, no matches found")
#     return(0)
#   }
#   
#   query_time_df <- proc.time() - ptm
#   ptm <- proc.time()
#   for (f in unique(tabix_f)){
#     tabix_file_logp = paste0(dir_data,"/",f,".signed_logP.tsv.gz.t.tsv.gz")
#     #print(paste0("tabix_file_logp: ",tabix_file_logp))
#     tabix_query = tabix_queries[tabix_f %in% f]
#     f_vct_snp_p = seqminer::tabix.read( tabix_file_logp, tabix_query)
#     f_vct_snp_p = strsplit(f_vct_snp_p, '\t')
#     f_df_snp_p = as.data.frame(do.call(rbind, f_vct_snp_p),stringsAsFactors=FALSE)
#     
#     x <- match(f_df_snp_p$V1,df.static_subset$i)
#     x <- x[!is.na(x)] 
#     df_snp_info <- rbind(df_snp_info,df.static_subset[x,])
#     df_snp_p <- rbind(df_snp_p,f_df_snp_p[,3:ncol(f_df_snp_p)])
#   }
#   if(nrow(df_snp_p)==1){
#     df_snp_info$minP <- formatC( min(10^-abs(as.numeric(df_snp_p))),digits = 3)
#   } else{
#     df_snp_info$minP <- formatC( rowMins(as.matrix(apply(df_snp_p,2, function (x) 10^-abs(as.numeric(x) ) )) ),digits=3)
#   }
# 
#   query_time_tabix <- proc.time() - ptm
#   
#   df_snp_info <- get_genes_df(df_snp_info,chr="CHR",position="BP",snp="SNP")
#   out=list()
#   out$df_snp_info = df_snp_info
#   out$df_snp_p = df_snp_p
#   out$num_snps_found = nrow(df_snp_p)
#   out$vct.missingsnps = vct.missingsnps
#   out$query_time_df = query_time_df
#   out$query_time_tabix = query_time_tabix
#   out$tabix_queries = tabix_queries
#   out$tabix_f = tabix_f
#   return(out)
# }




######################
#### VERSION 2.0
process_user_input <- function(entry,mapping.proteincoding,window=500000,subset="tophits",phenotype="unadjusted"){
  
  entry=unlist(strsplit( entry ,"\n|[, ]+"))
  entry=unlist(strsplit( entry ,"\n|[, ]+"))
  entry=unique(trimws(entry))
  
  if(length(entry )>1){
  #if multiple entries, keep only the ones that look like snps. 
  entry<-unlist(lapply(entry,function(x){
      if (length(unlist(strsplit( x ,"rs")))==2 | length(unlist(strsplit( x ,":")))==2) {x}
    }))
    
  }
  # raise error if entry ==0. 
  if (length(entry)==0){
    query$error <-TRUE
    return(query)
  }
  
  #initialize query object. 
  query=""
  if (length(entry)==1 ){
    #check if region;
    if( length(unlist(strsplit( entry ,"-|:")))==3 ){
      #region
      chr=unlist(strsplit( entry ,"-|:"))[1]
      startpos=as.numeric(unlist(strsplit( entry ,"-|:"))[2])-window
      endpos=as.numeric(unlist(strsplit( entry ,"-|:"))[3])+window
      query = list(entry="positional", snp="",chr=chr,startpos=startpos,endpos=endpos,gene="")
    } else if (length(unlist(strsplit( entry ,"rs")))==2 | length(unlist(strsplit( entry ,":")))==2)  {
      #1 snp;
      query = list(entry="snp",snp=entry,chr="",startpos="",endpos="",gene="")
    } else if (toupper(entry) %in% mapping.proteincoding$hgnc_symbol){
      # Gene. 
      entry=toupper(entry)
      query = list(entry="positional",snp="",chr = mapping.proteincoding[mapping.proteincoding$hgnc_symbol %in% entry,"chromosome_name"],
                   startpos = mapping.proteincoding[mapping.proteincoding$hgnc_symbol %in% entry,"start_position"]-window,
                   endpos = mapping.proteincoding[mapping.proteincoding$hgnc_symbol %in% entry,"end_position"]+window,
                   gene=entry)
      
    } 
  } else {
    query = list( entry="snp",snp=entry,chr="",startpos="",endpos="",gene="")
  }
  
  ### ERROR CHECK if position > 2.5 MB. 
  if( query$entry =="positional"  ){
    if((query$endpos-query$startpos)>2500000){
      showModal(modalDialog(
        title = "error",
        "Window too large, can't be larger than 2500000bp",
        easyClose = TRUE
      ))
      query$error <-TRUE
    return(query)
    }
  }
  
  ## ERROR IF Too many snps. 
  if( query$entry =="snp"  ){
    if(length(query$snp)>500){
      showModal(modalDialog(
        title = "error",
        "too many random snps, max is 500",
        easyClose = TRUE
      ))
      query$error <-TRUE
      return(query)
    }
  }
  
  query$error <-FALSE
  
  query$phenotype <- phenotype
  query$subset <- subset
  return(query)
}



get_tabix_query <- function(query,df.static.pos,df.static.rsid){
  #### query example for region:
  # query <-list(entry = "positional", snp = "", chr = "10", startpos = "1017954", 
  #              endpos = "1037954", gene = "")
  # 
  if(query$entry == "positional"){
    snps <- df.static.pos[.( as.numeric(query$chr) )][BP>as.numeric(query$startpos) & BP < as.numeric(query$endpos),c("CHR", "i")]
    tabix_query = paste0(snps$CHR[1],":",snps$i[1],"-",snps$i[nrow(snps)] )
  }
  
  ####### extract rsids or chr:pos
  #query <- list(entry = "snp", snp = c("rs575272151", "rs544419019", "rs540538026", "rs367896724",  "1:13110", "rs201106462", "20:2511230"), chr = "", startpos = "", endpos = "", gene = "")
  if(query$entry == "snp"){
    snps_rsids <- as.numeric(unlist(lapply(query$snp,function (x) if(startsWith(x,"rs" ) ){ sub("rs","",x)} )))
    snps_rsids <- df.static.rsid[.( snps_rsids )][,c("CHR", "i")]
    snps_positions <-strsplit(query$snp[grep(":",query$snp )],":")
    snps_positions <- t(sapply(snps_positions,function(x) as.numeric(gsub("_.*","",x) )) ) # remove _ if alleles are given and make numeric + matrix format
    if (length(snps_positions)>0){
      snps_positions <- df.static.pos[.(snps_positions[,1],snps_positions[,2] ),c("CHR","i")]
      snps <- unique(rbind(snps_rsids,snps_positions))
    } else { snps <-snps_rsids  }
    
    
    if(nrow(snps[!is.na(snps$i),])>0){
    tabix_query = paste0("chr", snps$CHR,":",snps$i,"-",snps$i)
    tabix_query = tabix_query[!grepl('NA-NA', tabix_query)]
    #cat(file=stderr(),"BEFORE BEDR: nrow(snps_rsids):",nrow(snps_rsids),"||| chromosome:",snps$CHR,"||| indices:",snps$i,"||| tabix_queries:",tabix_query,"\n")
    tabix_query = bedr.merge.region(tabix_query,verbose=FALSE)
    #cat(file=stderr(),"AFTER BEDR: nrow(snps_rsids):",nrow(snps_rsids),"||| chromosome:",snps$CHR,"||| indices:",snps$i,"||| tabix_queries:",tabix_query,"\n")
    tabix_query = as.character(sapply(tabix_query,function(x) gsub("chr","",x) ))
    } else{
      tabix_query = list(error=TRUE)
      return(tabix_query)
      }

  }
  ########
  
  tabix_query = list(query=query,
                     tabix_query=tabix_query,
                     indices = snps$i[!is.na(snps$i)],
                     error=FALSE)
  return(tabix_query)
}


extract_multiple_variants <- function(tabix_query,
                                      dir_data,
                                      f.data_p="stretch.logP.outfile.tsv.chr.gz.tophits.gz",
                                      f.data_beta="stretch.logP.outfile.tsv.chr.gz.tophits.gz",
                                      f.data_se="stretch.logP.outfile.tsv.chr.gz.tophits.gz",
                                      f.data.index="stretch.logP.outfile.index.tsv.gz.tophits.gz") {
  
  ptm <- proc.time()
  df_snp_info = seqminer::tabix.read( paste0(dir_data,"/",f.data.index), tabix_query$tabix_query)
  df_snp_info = strsplit(df_snp_info, '\t')
  df_snp_info = as.data.frame(do.call(rbind, df_snp_info),stringsAsFactors=FALSE)
  names(df_snp_info) <- c("i","uniqid","SNP","CHR","BP","EFAL","NEFAL","EAF","INFO")
  
  n.missingsnps  = sum(!df_snp_info$i %in% tabix_query$indices)
  
  df_snp_p = seqminer::tabix.read( tabixFile=paste0(dir_data,"/",f.data_p), tabixRange= tabix_query$tabix_query )
  df_snp_p = strsplit(df_snp_p, '\t')
  df_snp_p = as.data.frame(do.call(rbind, df_snp_p),stringsAsFactors=FALSE)
  df_snp_p = df_snp_p[3:502]
  df_snp_p = sapply(df_snp_p, as.numeric)
  if(length(df_snp_p)==500) df_snp_p = as.data.frame(t(as.data.frame(df_snp_p)))
  if (f.data_beta!="" & f.data_se!="" ){
    df_snp_beta = seqminer::tabix.read( tabixFile=paste0(dir_data,"/",f.data_beta), tabixRange= tabix_query$tabix_query )
    df_snp_beta = strsplit(df_snp_beta, '\t')
    df_snp_beta = as.data.frame(do.call(rbind, df_snp_beta),stringsAsFactors=FALSE)
    df_snp_beta = df_snp_beta[3:502]
    
    df_snp_se = seqminer::tabix.read( tabixFile=paste0(dir_data,"/",f.data_se), tabixRange= tabix_query$tabix_query )
    df_snp_se = strsplit(df_snp_se, '\t')
    df_snp_se = as.data.frame(do.call(rbind, df_snp_se),stringsAsFactors=FALSE)
    df_snp_se = df_snp_se[3:502]
  } else {
    df_snp_beta = ""
    df_snp_se = ""
  }
  
  if(length(df_snp_p)==500){
    df_snp_info$maxLogP <- max(abs(as.numeric(df_snp_p)))
    df_snp_info$minP <-  10^-df_snp_info$maxLogP
  } else{
    df_snp_info$maxLogP <- rowMaxs( as.matrix(apply(df_snp_p,2, function (x) abs(as.numeric(x) ) ) ))
    df_snp_info$minP <-   10^-df_snp_info$maxLogP
  }
  max.snp.index=which(abs(df_snp_p) == max(df_snp_info$maxLogP,na.rm=TRUE), arr.ind = TRUE)
  # lowestsnp=df_snp_info[max.snp.index[1,][1],]
  # lowesttrait=max.snp.index[,1][1]
  
  df_snp_info$minP <- as.numeric(formatC(as.numeric(df_snp_info$minP),digits=3))
  df_snp_info$INFO <- as.numeric(formatC(as.numeric(df_snp_info$INFO),digits=3))
  df_snp_info$EAF <- as.numeric(formatC(as.numeric(df_snp_info$EAF),digits=3))
  #df_snp_info$maxLogP <- as.numeric(formatC(as.numeric(df_snp_info$maxLogP),digits=3))
  
  query_time_tabix <- proc.time() - ptm
  
  df_snp_info <- get_genes_df(df_snp_info,chr="CHR",position="BP",snp="SNP")
  out=list()
  out$df_snp_info = df_snp_info
  out$df_snp_p =  as.data.frame(df_snp_p)
  out$df_snp_beta = df_snp_beta
  out$df_snp_se = df_snp_se
  out$num_snps_found = nrow(df_snp_p)
  out$n.missingsnps = n.missingsnps
  out$query_time_tabix = query_time_tabix
  out$query = tabix_query$query
  out$tabix_queries = tabix_query$tabix_query
  out$indices = tabix_query$indices
  out$max.snp.index = max.snp.index
  return(out)
}



