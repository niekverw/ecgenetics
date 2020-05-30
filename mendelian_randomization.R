require(raster)
require(ggplot2)
#install.packages("devtools")
#devtools::install_github("MRCIEU/TwoSampleMR")

setwd("/Users/niek/repos/ecgenetics/")
source("global.R")
source("helpers.R")
source("helpers.heatmap.R")
source("helpers.regionalplot.R")
source("helpers.mendelianrandomization.R")
source("get_nearest_gene.r")

## choose instruments; 
# classic ECG: 
dfmrexposures <- data.frame(fread("data/tableS13.tsv"))
dfmrexposures <- dfmrexposures[dfmrexposures$TRAIT %in% unique(dfmrexposures$TRAIT)[6] & dfmrexposures$Included.in.the.PRS %in% 1,]

# MR BASE: 
#devtools::install_github("mrcieu/ieugwasr")

library(ieugwasr)

#remotes::install_github("MarkEdmondson1234/googleAnalyticsR")
                        
#get_access_token()
ieugwases <- ieugwasr::gwasinfo()

gwasid = "ebi-a-GCST006414"
gwasid="ukb-d-12340_irnt"
gwasid = "ukb-b-19953"

dfmrexposures <- tophits(
  id=gwasid,
  pval = 5e-08,
  clump = 1,
  r2 = 0.001,
  kb = 10000,
  pop = "EUR",
  force_server = FALSE,
  access_token = check_access_token()
)
dfmrexposures<-dfmrexposures[,c("p", "se", "n", "beta", "position", "chr", "id", "rsid", "ea",  "nea", "eaf", "trait")]
names(dfmrexposures) <- c("PVAL", "SE", "N", "BETA", "BP", "CHR", "id", "SNP", "EFAL","NEFAL", "EAF", "TRAIT")

fwrite(x=dfmrexposures,file = paste0("/Users/niek/Downloads/",gwasid,".tsv"),quote = F,sep="\t" )

##############

dfmrexposures$uniqid <- make_uniqID(dfmrexposures$CHR,dfmrexposures$BP,dfmrexposures$EFAL,dfmrexposures$NEFAL)

input= unique(dfmrexposures$uniqid)
query <- process_user_input(input,mapping.proteincoding)
tabix_query <- get_tabix_query(query,df.static.pos,df.static.rsid)

data_unadjusted <- extract_multiple_variants(tabix_query,dir_data,
                                               f.data_p="unadjusted.logP.outfile.tsv.chr.gz.tophits.gz",
                                               f.data_beta="unadjusted.BETA.outfile.tsv.gz.tophits.gz",
                                               f.data_se="unadjusted.SE.outfile.tsv.gz.tophits.gz",
                                               f.data.index="unadjusted.logP.outfile.index.tsv.gz.tophits.gz"
)
  
  data_unadjusted <- harmonizedfs(dfmrexposures,data_unadjusted)
  mr_unadjusted <- ecg_wide_ivw(data_unadjusted)
  data_unadjusted$df_snp_info[1,]$SNP <- "IVW-fixed effect unadjusted"
  data_unadjusted$df_snp_info[1,]$Gene <- "IVW-fixed effect unadjusted"
  mrplot_unadjusted <- make_ecg_plot(vct_snp_p=t(data.frame(mr_unadjusted$p)),
                                     vct_snp_beta=t(data.frame(mr_unadjusted$b)), #data$df_snp_beta,
                                     vct_snp_se=t(data.frame(mr_unadjusted$s)), #data$df_snp_se,
                                     vct_snp_info=data_unadjusted$df_snp_info[1,],
                                     df_ecg_stats=df_ecg_unadjusted,
                                     invert=FALSE) #+ ggtitle(paste0("unadjusted"))
  
  mrplot_unadjusted
  
  
  data_rradjusted <- extract_multiple_variants(tabix_query,dir_data,
                                               f.data_p="stretch.logP.outfile.tsv.chr.gz.tophits.gz",
                                               f.data_beta="stretch.BETA.outfile.tsv.gz.tophits.gz",
                                               f.data_se="stretch.SE.outfile.tsv.gz.tophits.gz",
                                               f.data.index="stretch.logP.outfile.index.tsv.gz.tophits.gz"
  )
  
  data_rradjusted <- harmonizedfs(dfmrexposures,data_rradjusted)
  mr_rradjusted <- ecg_wide_ivw(data_rradjusted)
  data_rradjusted$df_snp_info[1,]$SNP <- "IVW-fixed effect RR-adjusted"
  data_rradjusted$df_snp_info[1,]$Gene <- "IVW-fixed effect RR-adjusted"
  
  mrplot_rradjusted <- make_ecg_plot(vct_snp_p=t(data.frame(mr_rradjusted$p)),
                                     vct_snp_beta=t(data.frame(mr_rradjusted$b)), #data$df_snp_beta,
                                     vct_snp_se=t(data.frame(mr_rradjusted$s)), #data$df_snp_se,
                                     vct_snp_info=data_rradjusted$df_snp_info[1,],
                                     df_ecg_stats=df_ecg_unadjusted,
                                     invert=FALSE) #+ ggtitle(paste0("RR-adjusted"))
  
  mrplot_rradjusted
  
  dfmrexposures_NA <- dfmrexposures[!dfmrexposures$uniqid %in% data_unadjusted$df_snp_info$uniqid,]
  dfmrexposures <- dfmrexposures[dfmrexposures$uniqid %in% data_unadjusted$df_snp_info$uniqid,]
  



