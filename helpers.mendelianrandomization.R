library(MendelianRandomization)


make_uniqID<-function(chr,pos,allele1,allele2){
   alleleA <- allele1
   alleleB <- allele2
   alleleA[allele1 > allele2] <- allele2[allele1 > allele2]
   alleleB[allele1 > allele2] <- allele1[allele1 > allele2]
   snpids=paste0(chr,":",pos,"_",alleleA,"_",alleleB)
   
   return(snpids)
}
harmonizedfs <- function(dfexposure,data){
   
   dfexposure
   # harmonize
   UNIQID1 = "uniqid"
   EFAL1 <- "EFAL"
   NEFAL1 <- "NEFAL"
   EAF1 <- "EAF"
   BETA1 <- "BETA"
   SE1 <- "SE"
   
   dfexposure <- dfexposure[,c(UNIQID1,EFAL1,NEFAL1,EAF1,BETA1,SE1)]
   colnames(dfexposure) <- c("uniqid","EFAL","NEFAL","EAF","BETA","SE")
   
   
   dfexposure <- dfexposure[match(data$df_snp_info$uniqid,dfexposure$uniqid),]
   
   available_in_ecgdata <- data$df_snp_info$uniqid %in% dfexposure$uniqid
   data$df_snp_info <- data$df_snp_info[available_in_ecgdata,]
   data$df_snp_p <- data$df_snp_p[available_in_ecgdata,]
   data$df_snp_beta <- data$df_snp_beta[available_in_ecgdata,]
   data$df_snp_se <- data$df_snp_se[available_in_ecgdata,]
   
   iswitch = dfexposure$EFAL != data$df_snp_info$EFAL & dfexposure$NEFAL == data$df_snp_info$EFAL
   dfexposure[iswitch,"BETA"] <- dfexposure[iswitch,"BETA"]*-1
   dfexposure[iswitch,"EAF"] <- 1-dfexposure[iswitch,"EAF"]
   dfexposure[iswitch,"tmp"] <- dfexposure[iswitch,"EFAL"]
   dfexposure[iswitch,"EFAL"] <- dfexposure[iswitch,"NEFAL"]
   dfexposure[iswitch,"NEFAL"] <- dfexposure[iswitch,"tmp"]
   dfexposure[,"tmp"] <- NULL
   #print(plot(dfexposure$EAF, data$df_snp_info$EAF))
   data$dfexposure<-dfexposure
   return(data)
}

ecg_wide_ivw <- function(data){
   pvalues=c()
   betas=c()
   ses=c()
   for (i in 1:500) {
      
      #data$dfexposure
      
      
      
      MRInputObject <- mr_input(snp=data$df_snp_info$SNP, bx = data$dfexposure$BETA, bxse = data$dfexposure$SE,
                                by = as.numeric(data$df_snp_beta[,i]),
                                byse = as.numeric(data$df_snp_se[,i]))
      ivw <- MendelianRandomization::mr_ivw(MRInputObject)
      
      pvalues=c(pvalues,-log10(ivw@Pvalue)*sign(ivw@Estimate))
      
      #beta=if(is.nan(-log10(ivw@Estimate)) ){0} else{(ivw@Estimate)}
      betas=c(betas,ivw@Estimate)
      ses=c(ses,ivw@StdError)
      
   }
   
   #plot(1:500,pvalues)
   
   d <- data.frame(p=pvalues,b=betas,s=ses)
   return(d)
}

