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
   # dfexposure <- dfmrexposures
   dfexposure <- dfexposure[!is.na(dfexposure$uniqid),]
   
   # harmonize
   UNIQID1 = "uniqid"
   EFAL1 <- "EFAL"
   NEFAL1 <- "NEFAL"
   EAF1 <- "EAF"
   BETA1 <- "BETA"
   SE1 <- "SE"

   dfexposure <- dfexposure[,c(UNIQID1,EFAL1,NEFAL1,EAF1,BETA1,SE1)]
   colnames(dfexposure) <- c("uniqid","EFAL","NEFAL","EAF","BETA","SE")
   
   available_in_ecgdata <- data$df_snp_info$uniqid %in% dfexposure$uniqid
   data$df_snp_info <- data$df_snp_info[available_in_ecgdata,]
   data$df_snp_p <- data$df_snp_p[available_in_ecgdata,]
   data$df_snp_beta <- data$df_snp_beta[available_in_ecgdata,]
   data$df_snp_se <- data$df_snp_se[available_in_ecgdata,]
   
   dfexposure <- dfexposure[match(data$df_snp_info$uniqid,dfexposure$uniqid),]
   

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

#TwoSampleMR
mr_ivw_fe <- function (b_exp, b_out, se_exp, se_out, parameters = default_parameters()) 
{
   if (sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & 
           !is.na(se_out)) < 2) 
      return(list(b = NA, se = NA, pval = NA, nsnp = NA))
   ivw.res <- summary(lm(b_out ~ -1 + b_exp, weights = 1/se_out^2))
   b <- ivw.res$coef["b_exp", "Estimate"]
   se <- ivw.res$coef["b_exp", "Std. Error"]/ivw.res$sigma
   pval <- 2 * pnorm(abs(b/se), lower.tail = FALSE)
   Q_df <- length(b_exp) - 1
   Q <- ivw.res$sigma^2 * Q_df
   Q_pval <- pchisq(Q, Q_df, lower.tail = FALSE)
   return(list(b = b, se = se, pval = pval, nsnp = length(b_exp), 
               Q = Q, Q_df = Q_df, Q_pval = Q_pval))
}

ecg_wide_ivw <- function(data){
   pvalues=c()
   betas=c()
   ses=c()
   for (i in 1:500) {
      
      #data$dfexposure
      
      t <- mr_ivw_fe(b_exp = data$dfexposure$BETA, 
                se_exp = data$dfexposure$SE,
                b_out = as.numeric(data$df_snp_beta[,i]),
                                           se_out = as.numeric(data$df_snp_se[,i]))
      if(is.infinite(abs(-log10(t$pval)) )){ t$pval=1e-307 }
      
      pvalues=c(pvalues,-log10(t$pval)*sign(t$b))
      betas=c(betas,t$b)
      ses=c(ses,t$se)
      
      # MRInputObject <- mr_input(snp=data$df_snp_info$SNP, bx = data$dfexposure$BETA, bxse = data$dfexposure$SE,
      #                           by = as.numeric(data$df_snp_beta[,i]),
      #                           byse = as.numeric(data$df_snp_se[,i]))
      # #ivw <- MendelianRandomization::mr_ivw (MRInputObject)
      # pvalues=c(pvalues,-log10(ivw@Pvalue)*sign(ivw@Estimate))
      # betas=c(betas,ivw@Estimate)
      # ses=c(ses,ivw@StdError)
      
   }
   
   
   #plot(1:500,pvalues)
   
   d <- data.frame(p=pvalues,b=betas,s=ses)
   return(d)
}

