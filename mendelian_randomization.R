require(raster)
require(ggplot2)

setwd("/Users/niek/repos/ecgenetics/")
source("global.R")
source("helpers.R")
source("helpers.heatmap.R")
source("helpers.regionalplot.R")
source("get_nearest_gene.r")

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
      
      data$dfexposure
      
      
      
      MRInputObject <- mr_input(snp=data$df_snp_info$SNP, bx = data$dfexposure$BETA, bxse = data$dfexposure$SE,
                                by = as.numeric(data$df_snp_beta[,i]),
                                byse = as.numeric(data$df_snp_se[,i]))
      ivw <- mr_ivw(MRInputObject)
      
      pvalues=c(pvalues,-log10(ivw@Pvalue)*sign(ivw@Estimate))
      
      #beta=if(is.nan(-log10(ivw@Estimate)) ){0} else{(ivw@Estimate)}
      betas=c(betas,ivw@Estimate)
      ses=c(ses,ivw@StdError)
      
   }
   
   #plot(1:500,pvalues)
   
   d <- data.frame(p=pvalues,b=betas,s=ses)
   return(d)
}
dftables13 <- data.frame(fread("tableS13.tsv"))
#dftables13 <- dftables13[dftables13$TRAIT %in% "PR-Interval",]
#dftables13 <- dftables13[dftables13$TRAIT %in% "QRS-Duration" & dftables13$Included.in.the.PRS %in% 1,]
dftables13 <- dftables13[dftables13$TRAIT %in% unique(dftables13$TRAIT)[6] & dftables13$Included.in.the.PRS %in% 1,]

dftables13$uniqid <- make_uniqID(dftables13$CHR,dftables13$position,dftables13$EFAL,dftables13$NEFAL)
   
input= unique(dftables13$uniqid)
query <- process_user_input(input,mapping.proteincoding)
tabix_query <- get_tabix_query(query,df.static.pos,df.static.rsid)

# unadjusted
data_unadjusted <- extract_multiple_variants(tabix_query,dir_data,
                                  f.data_p="unadjusted.logP.outfile.tsv.chr.gz.tophits.gz",
                                  f.data_beta="unadjusted.BETA.outfile.tsv.gz.tophits.gz",
                                  f.data_se="unadjusted.SE.outfile.tsv.gz.tophits.gz",
                                  f.data.index="unadjusted.logP.outfile.index.tsv.gz.tophits.gz"
)

data_unadjusted <- harmonizedfs(dftables13,data_unadjusted)
mr_unadjusted <- ecg_wide_ivw(data_unadjusted)
data$df_snp_info[1,]$SNP <- ""
data$df_snp_info[1,]$Gene <- "IVW-fixed effect unadjusted"
make_ecg_plot(vct_snp_p=t(data.frame(mr_unadjusted$p)),
              vct_snp_beta=t(data.frame(mr_unadjusted$b)), #data$df_snp_beta,
              vct_snp_se=t(data.frame(mr_unadjusted$s)), #data$df_snp_se,
              vct_snp_info=data$df_snp_info[1,],
              df_ecg_stats=df_ecg_unadjusted,
              invert=FALSE)




data_rradjusted <- extract_multiple_variants(tabix_query,dir_data,
                                             f.data_p="stretch.logP.outfile.tsv.chr.gz.tophits.gz",
                                             f.data_beta="stretch.BETA.outfile.tsv.gz.tophits.gz",
                                             f.data_se="stretch.SE.outfile.tsv.gz.tophits.gz",
                                             f.data.index="stretch.logP.outfile.index.tsv.gz.tophits.gz"
)

data_rradjusted <- harmonizedfs(dftables13,data_rradjusted)
mr_rradjusted <- ecg_wide_ivw(data_rradjusted)
data$df_snp_info[1,]$SNP <- ""
data$df_snp_info[1,]$Gene <- "IVW-fixed effect RR-adjusted"
make_ecg_plot(vct_snp_p=t(data.frame(mr_rradjusted$p)),
              vct_snp_beta=t(data.frame(mr_rradjusted$b)), #data$df_snp_beta,
              vct_snp_se=t(data.frame(mr_rradjusted$s)), #data$df_snp_se,
              vct_snp_info=data$df_snp_info[1,],
              df_ecg_stats=df_ecg_unadjusted,
              invert=FALSE)





